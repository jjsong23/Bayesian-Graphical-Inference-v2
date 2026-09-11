param(
    [int]$Port = 8765,
    [switch]$Background
)

$workbenchDir = $PSScriptRoot
$bundledPython = Join-Path $env:USERPROFILE ".cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe"
$pythonExecutable = $null
if (Test-Path -LiteralPath $bundledPython) {
    $pythonExecutable = $bundledPython
} else {
    $pythonCommand = Get-Command python -ErrorAction SilentlyContinue
    if ($null -eq $pythonCommand) {
        $pythonCommand = Get-Command py -ErrorAction SilentlyContinue
    }
    if ($null -ne $pythonCommand) {
        $pythonExecutable = $pythonCommand.Source
    }
}
if ($null -eq $pythonExecutable) {
    throw "Python 3 was not found. Start the workbench from the Codex project environment."
}

if ($Background) {
    & $pythonExecutable "$workbenchDir\launch_workbench.py" --port $Port
} else {
    & $pythonExecutable "$workbenchDir\server.py" --port $Port
}
