# ============================================================
# MAFFT多序列比对脚本
# 用法: .\run_mafft.ps1
# ============================================================

Write-Host "============================================================" -ForegroundColor Cyan
Write-Host "MAFFT Multiple Sequence Alignment" -ForegroundColor Cyan
Write-Host "============================================================" -ForegroundColor Cyan

# 参数设置
$InputFasta = "Tree\phylogeny\all_sequences_for_tree.fasta"
$OutputFasta = "Tree\phylogeny\all_sequences_aligned.fasta"
$MAFFT_BAT = "Tree\mafft-win\mafft.bat"

# 检查输入文件
if (-Not (Test-Path $InputFasta)) {
    Write-Host "ERROR: Input file not found: $InputFasta" -ForegroundColor Red
    exit 1
}

Write-Host ""
Write-Host "Input file: $InputFasta" -ForegroundColor Yellow
Write-Host "Output file: $OutputFasta" -ForegroundColor Yellow
Write-Host ""

# 运行MAFFT
Write-Host "Running MAFFT..." -ForegroundColor Green
& $MAFFT_BAT --auto $InputFasta > $OutputFasta

# 检查输出
if (Test-Path $OutputFasta) {
    $seqCount = (Get-Content $OutputFasta | Where-Object { $_ -match "^>" }).Count
    Write-Host "SUCCESS: Alignment completed" -ForegroundColor Green
    Write-Host "Sequences aligned: $seqCount" -ForegroundColor Green
} else {
    Write-Host "ERROR: Alignment failed" -ForegroundColor Red
    exit 1
}

Write-Host ""
Write-Host "============================================================" -ForegroundColor Cyan
Write-Host "MAFFT completed!" -ForegroundColor Cyan
Write-Host "============================================================" -ForegroundColor Cyan