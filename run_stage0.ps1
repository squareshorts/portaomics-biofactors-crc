$ErrorActionPreference = "Stop"
Set-StrictMode -Version Latest

function Invoke-Checked {
    param(
        [Parameter(Mandatory = $true)][string]$Exe,
        [Parameter(Mandatory = $true)][string[]]$CommandArgs
    )

    & $Exe @CommandArgs
    if ($LASTEXITCODE -ne 0) {
        throw "Command failed with exit code $LASTEXITCODE: $Exe $($CommandArgs -join ' ')"
    }
}

function Test-UsableCpython {
    param([string]$Exe)

    if (-not $Exe -or -not (Test-Path $Exe)) {
        return $false
    }

    $resolved = (Resolve-Path $Exe).Path
    if ($resolved -match '(?i)\\msys64\\') {
        return $false
    }

    & $resolved -c "import sys; raise SystemExit(0 if sys.implementation.name == 'cpython' and sys.version_info >= (3, 11) and sys.version_info < (3, 14) else 1)" | Out-Null
    return ($LASTEXITCODE -eq 0)
}

$VenvPython = $null

if (Test-Path ".venv\Scripts\python.exe") {
    $candidate = (Resolve-Path ".venv\Scripts\python.exe").Path
    if (Test-UsableCpython $candidate) {
        $VenvPython = $candidate
    }
    else {
        throw "The existing .venv was not created by a supported standard CPython installation. Delete .venv and recreate it with standard CPython 3.13. See WINDOWS_SETUP.md."
    }
}
elif (Test-Path ".venv\bin\python.exe") {
    $candidate = (Resolve-Path ".venv\bin\python.exe").Path
    throw "The existing .venv uses the MSYS2/POSIX-style bin layout ($candidate). PyPI scientific wheels are not available for this interpreter. Delete .venv and recreate it with standard CPython 3.13. See WINDOWS_SETUP.md."
}
else {
    $Candidates = @(
        "$env:LOCALAPPDATA\Programs\Python\Python313\python.exe",
        "$env:ProgramFiles\Python313\python.exe"
    )

    $BasePython = $null
    foreach ($candidate in $Candidates) {
        if (Test-UsableCpython $candidate) {
            $BasePython = $candidate
            break
        }
    }

    if (-not $BasePython) {
        $pythonCmd = Get-Command python -ErrorAction SilentlyContinue
        if ($pythonCmd -and (Test-UsableCpython $pythonCmd.Source)) {
            $BasePython = $pythonCmd.Source
        }
    }

    if (-not $BasePython) {
        throw "A supported standard CPython installation was not found. Install CPython 3.13 with: winget install -e --id Python.Python.3.13 --scope user ; then rerun this script. Do not use C:\msys64\...\python.exe."
    }

    Write-Host "Creating .venv with: $BasePython"
    Invoke-Checked -Exe $BasePython -CommandArgs @("-m", "venv", ".venv")
    $VenvPython = (Resolve-Path ".venv\Scripts\python.exe").Path
}

Write-Host "Using virtual-environment Python: $VenvPython"
Invoke-Checked -Exe $VenvPython -CommandArgs @("-m", "pip", "install", "--upgrade", "pip")
Invoke-Checked -Exe $VenvPython -CommandArgs @("-m", "pip", "install", "-r", "requirements-stage0.txt")
Invoke-Checked -Exe $VenvPython -CommandArgs @("scripts\run_stage0.py")

Write-Host ""
Write-Host "Stage 0 completed successfully."
Write-Host "Stage-0 report: reports\stage0_audit.md"
