# Windows setup for Stage 0

Do not use the MSYS2/UCRT64 Python installation at `C:\msys64\...\python.exe` for this project. Its Python package tags differ from the standard Windows CPython ABI, so pip can fall back to source builds for NumPy and other scientific packages.

Use standard CPython 3.13.

From PowerShell in `C:\work\portaomics_biofactors_crc`:

```powershell
Remove-Item -Recurse -Force .venv -ErrorAction SilentlyContinue
winget install -e --id Python.Python.3.13 --scope user
& "$env:LOCALAPPDATA\Programs\Python\Python313\python.exe" --version
& "$env:LOCALAPPDATA\Programs\Python\Python313\python.exe" -m venv .venv
.\run_stage0.ps1
```

The launcher is fail-closed. It stops immediately if dependency installation or any Stage-0 script fails and prints the Stage-0 report path only after a successful run.
