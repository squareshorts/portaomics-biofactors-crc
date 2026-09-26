from __future__ import annotations

import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent


def run(name: str):
    print(f"\n=== {name} ===")
    subprocess.run([sys.executable, str(HERE / name)], check=True)


def main():
    for s in ["00_download_geo.py", "01_build_manifest.py", "02_preprocess_expression.py", "03_stage0_report.py"]:
        run(s)
    print("\nStage 0 complete. No classifier or disease-association analysis was executed.")


if __name__ == "__main__":
    main()
