"""Regression checks for scientifically material verification edge cases."""
import importlib
import io
import pathlib
import sys
import tempfile
import unittest
from unittest.mock import patch
import numpy as np
import pandas as pd
from scipy.integrate import quad

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "scripts"))
obs = importlib.import_module("02_obs")
scores = importlib.import_module("07_scores")
stats = importlib.import_module("08_stats")
provenance = importlib.import_module("provenance")

class PrecipitationTests(unittest.TestCase):
    def resolve(self, aa="24,9999,9,9", rem="SYN000 333 70000=", report="FM-12"):
        row = dict(time=pd.Timestamp("2023-04-06"), REPORT_TYPE=report, AA1=aa, REM=rem)
        return obs.precipitation(pd.DataFrame([row]))

    def test_missing_amount_requires_raw_evidence(self):
        pr, _, audit = self.resolve(rem="SYN000 333 7////=")
        self.assertEqual(pr, {})
        self.assertEqual(audit["n_unresolved_missing_reports"], 1)
        pr, source, _ = self.resolve()
        self.assertEqual(pr[pd.Timestamp("2023-04-05")], 0.0)
        self.assertEqual(source[pd.Timestamp("2023-04-05")], "synop33370000")

    def test_weather_metar_and_later_sections_are_not_daily_totals(self):
        for rem, report in [("SYN000 70000 333 10177=", "FM-12"),
                            ("SYN000 333 10177 555 70000=", "FM-12"),
                            ("SYN000 333 70000 70001=", "FM-12"),
                            ("MET000 333 70000=", "FM-12"),
                            ("SYN000 333 70000=", "FM-15")]:
            with self.subTest(rem=rem, report=report):
                self.assertEqual(self.resolve(rem=rem, report=report)[0], {})

    def test_rejected_conditions_and_quality_cannot_be_recovered(self):
        for condition in "1345678IJ":
            self.assertEqual(self.resolve(aa=f"24,9999,{condition},9")[0], {})
            self.assertEqual(self.resolve(aa=f"24,0100,{condition},9")[0], {})
        for quality in "2367":
            self.assertEqual(self.resolve(aa=f"24,9999,9,{quality}")[0], {})

    def test_numeric_zero_trace_and_positive_amounts_are_valid(self):
        for field, expected in [("24,0000,9,9", 0), ("24,0000,2,1", 0),
                                ("24,0123,9,1", 12.3)]:
            pr, _, _ = self.resolve(aa=field, rem="")
            self.assertAlmostEqual(pr[pd.Timestamp("2023-04-05")], expected)

    def test_conflicting_duplicate_reports_are_missing(self):
        rows = [dict(time=pd.Timestamp("2023-04-06"), AA1=f"24,{amount},9,1")
                for amount in ("0100", "0200")]
        pr, _, audit = obs.precipitation(pd.DataFrame(rows))
        self.assertEqual(pr, {})
        self.assertEqual(audit["n_conflicting_days"], 1)

    def test_no_aa24_does_not_recover_an_unmatched_raw_group(self):
        self.assertEqual(self.resolve(aa="06,9999,9,9")[0], {})

    def test_flagged_numeric_amount_needs_matching_daily_source_group(self):
        pr, source, audit = self.resolve(aa="24,0011,3,1", rem="SYN000 333 60011 70011=")
        self.assertEqual(pr[pd.Timestamp("2023-04-05")], 1.1)
        self.assertEqual(source[pd.Timestamp("2023-04-05")], "synop333-confirmed-aa24")
        self.assertEqual(audit["n_flagged_amounts_confirmed_groups"], 1)
        self.assertEqual(self.resolve(aa="24,0011,3,1", rem="SYN000 333 70012=")[0], {})
        self.assertEqual(self.resolve(aa="24,9999,3,1", rem="SYN000 333 70011=")[0], {})

    def test_real_beijing_report_recovers_only_documented_zero(self):
        path = pathlib.Path(__file__).parent / "fixtures" / "isd_beijing_2023.csv"
        _, daily, audit = obs.parse(path.read_text(), "54511099999", 2023, return_audit=True)
        self.assertEqual(audit["n_zero_recovered_days"], 3)
        self.assertEqual(daily.pr24.tolist(), [0.0, 0.0, 0.0])

class StatisticsTests(unittest.TestCase):
    def test_valid_zero_reference_losses_are_included(self):
        a = np.ones(60); b = np.tile([1.0, 0.0], 30)
        skill, lo, hi = stats.block_boot(a, b)
        self.assertAlmostEqual(skill, -1.0)
        self.assertAlmostEqual(lo, -1.0)
        self.assertAlmostEqual(hi, -1.0)

    def test_missing_days_keep_calendar_spacing(self):
        a = np.ones(90); b = np.ones(90)
        available = np.ones(90, bool); available[30:45] = False
        class FixedStarts:
            def integers(self, low, high, size):
                self.high = high
                return np.zeros(size, dtype=int)
        fixed = FixedStarts()
        with patch.object(stats, "rng", fixed):
            self.assertEqual(stats.block_boot(a, b, available=available), (0, 0, 0))
        # The resampling range covers the full 90-day calendar, including the gap.
        self.assertEqual(fixed.high, 90 - 14 + 1)

    def test_zero_total_reference_is_undefined(self):
        self.assertTrue(all(np.isnan(v) for v in stats.block_boot(np.ones(60), np.zeros(60))))

    def test_calibration_has_no_future_observations(self):
        inits = pd.date_range("2023-12-01", "2024-01-10")
        for var, lead, last in [("t2", 240, "2023-12-21"), ("ws", 6, "2023-12-31"),
                                ("pr24", 10, "2023-12-21"), ("tmax", 9, "2023-12-22"),
                                ("tmin", 1, "2023-12-30")]:
            train, test = scores.oos_masks(inits, np.ones(len(inits), bool), "2023-12-31", var, lead)
            self.assertEqual(inits[train][-1], pd.Timestamp(last))
            self.assertEqual(inits[test][0], pd.Timestamp("2024-01-01"))
            self.assertTrue((scores.observation_end(var, lead, inits[train]) < inits[test][0]).all())

    def test_all_six_hourly_leads_and_precipitation_models(self):
        self.assertTrue(all(scores.save_lead("t2", lead) for lead in range(6, 241, 6)))
        self.assertEqual(scores.expected_models("pr24", "core"), ["gfs", "graphcast"])
        self.assertEqual(len(scores.expected_models("t2", "core")), 4)

    def test_step_crps_and_tail_against_independent_quadrature(self):
        z = np.array([-1.0, 0.5, 3.0]); F = np.array([[.2, .6, 1.0]])
        for y in (-3.0, -0.3, 4.0):
            for a, b in [(-5, 5), (1, 5), (-5, 0)]:
                def integrand(x):
                    k = np.searchsorted(z, x, side="right") - 1
                    cdf = 0.0 if k < 0 else F[0, k]
                    return (cdf - (y <= x)) ** 2
                points = [x for x in [*z, y] if a < x < b]
                expected = quad(integrand, a, b, points=points)[0]
                actual = scores.crps_step(z, F, np.array([y]), a, b)[0]
                self.assertAlmostEqual(actual, expected, places=9)
        self.assertEqual(scores.crps_step(np.array([0.]), np.ones((1, 1)), np.array([2.]))[0], 2.)

    def test_gaussian_threshold_weighted_crps_against_independent_quadrature(self):
        from scipy.stats import norm
        def reference(mu, sig, y, a, b):
            f = lambda x: (norm.cdf((x - mu) / sig) - (y <= x)) ** 2
            lo = max(a, min(mu - 40 * sig, y) - 1.0); hi = min(b, max(mu + 40 * sig, y) + 1.0)
            if hi <= lo:
                return 0.0
            pts = [p for p in (y, mu) if lo < p < hi]
            return quad(f, lo, hi, points=pts or None, limit=500, epsabs=1e-12, epsrel=1e-10)[0]
        cases = []
        for mu, sig in [(0.0, 1.0), (25.0, 0.4), (-3.0, 2.5)]:
            for k in (-3.0, 0.0, 1.5, 12.0):       # threshold position in predictive sigmas
                t = mu + k * sig
                for y in (t - 5.0, t - 0.1, t, t + 0.3, t + 6.0, mu):
                    cases += [(mu, sig, y, t, np.inf), (mu, sig, y, -np.inf, t)]
            cases += [(mu, sig, mu + 1.0, mu - sig, mu + 2 * sig)]
        for mu, sig, y, a, b in cases:
            with self.subTest(mu=mu, sig=sig, y=y, a=a, b=b):
                got = scores.gauss_twcrps(np.array([mu]), np.array([sig]), np.array([y]), a, b)[0]
                self.assertGreaterEqual(got, 0.0)
                self.assertAlmostEqual(got, reference(mu, sig, y, a, b), places=7)
        # threshold far above a sharp predictive and the observation: score ~0, never negative
        mu = np.array([20.0, 20.0]); sig = np.array([0.3, 0.3]); y = np.array([19.0, 31.0])
        s = scores.gauss_twcrps(mu, sig, y, 30.0, np.inf)
        self.assertTrue(np.all(s >= 0)); self.assertLess(s[0], 1e-12); self.assertAlmostEqual(s[1], 1.0, places=9)
        # unweighted limit equals the closed-form CRPS
        self.assertAlmostEqual(scores.gauss_twcrps(mu, sig, y, -np.inf, np.inf)[0], scores.gauss_crps(mu, sig, y)[0], places=12)

    def test_score_and_statistics_pipeline_on_paired_synthetic_cubes(self):
        with tempfile.TemporaryDirectory() as directory:
            root = pathlib.Path(directory); analysis = root / "data/analysis"; data = analysis / "t2"
            data.mkdir(parents=True); (root / "data/interim").mkdir()
            output = root / "scores"; output.mkdir(); tables = root / "tables"; tables.mkdir()
            inits = pd.date_range("2022-01-01", "2024-05-31")
            rng = np.random.default_rng(41)
            y = rng.normal(10, 5, (len(inits), 2, 1))
            y[20:25] = np.nan
            np.save(data / "obs.npy", y)
            np.save(data / "inits.npy", inits.values)
            np.save(data / "leads.npy", np.array([6, 24]))
            np.save(data / "stations.npy", np.array(["test-station"]))
            for i, model in enumerate(scores.expected_models("t2", "core")):
                np.save(data / f"{model}.npy", y + rng.normal(0, 2 + i, y.shape))
            pd.DataFrame([dict(station_id="test-station", zone="Temperate")]).to_csv(root / "data/interim/stations_analysis.csv", index=False)
            row = dict(station_id="test-station")
            for hour in (0, 6, 12, 18):
                row[f"t2_h{hour:02d}_q95"] = 18
                row[f"t2_h{hour:02d}_q05"] = 2
            pd.DataFrame([row]).to_csv(analysis / "climatology_2011_2020.csv", index=False)
            provenance.write_metadata(analysis / "provenance.json")
            with patch.object(scores, "AN", analysis), patch.object(scores, "ROOT", root), patch.object(scores, "OUT", output), \
                 patch.object(stats, "SC", output), patch.object(stats, "TB", tables):
                scores.main("t2", "core", workers=1)
                stats.main("t2", "core")
            station = pd.read_parquet(output / "core_t2_station.parquet")
            self.assertEqual(set(station.lead), {6, 24})
            self.assertEqual(station.loc[station.lead == 6, "n_train"].unique().tolist(), [725])
            self.assertEqual(station.loc[station.lead == 24, "n_train"].unique().tolist(), [724])
            bootstrap = pd.read_csv(tables / "core_t2_bootstrap.csv")
            self.assertEqual(set(bootstrap[bootstrap.score == "oos_tw2"].lead), {6, 24})
            provenance.require_current(tables / "core_t2_meta.json")

if __name__ == "__main__":
    unittest.main()
