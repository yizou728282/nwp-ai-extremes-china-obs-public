#!/bin/bash
# Full core design after the 2023 validation year: GFS-initialised AI + GFS for 2022, 2024, 2025, 2021,
# then IFS-initialised AI 2022-2025 (initial-condition decomposition). Every step is resumable.
cd "$(dirname "$0")"
for Y in 2022 2024 2025 2021; do
  echo "$(date) start $Y"; ./run_year.sh $Y; echo "$(date) done $Y"; df -h . | tail -1
done
cd scripts
echo "$(date) start IFS-init"
../.venv/bin/python -u 03_stream_mlwp.py graphcast_ifs pangu_ifs fourcastnetv2_ifs --years 2022-2025 --threads 16 >> ../logs/03_ai_ifs.log 2>&1
echo "$(date) all done"
