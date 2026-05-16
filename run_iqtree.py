"""
IQ-TREE最大似然建树脚本

功能说明：
  1. 使用trimAI修剪多序列比对结果
  2. 使用IQ-TREE构建最大似然树

注意：
  - MAFFT比对需要单独运行（参见 run_mafft.ps1 或 run_phylogeny_pipeline.ps1）
  - IQ-TREE参数：
    -m MFP: ModelFinder Plus自动选择最佳模型
    -B 1000: 1000次ultrafast bootstrap
    --bnni: Bootstrap + NNI优化提高树稳定性
    -T AUTO: 自动选择最佳线程数

使用方法：
  python run_iqtree.py

输出文件：
  Tree/phylogeny/full_length.treefile  - 最大似然树（Newick格式）
  Tree/phylogeny/full_length.contree   - Bootstrap共识树
  Tree/phylogeny/full_length.iqtree   - IQ-TREE详细报告
"""

import subprocess
import os

# ===================== 参数设置 =====================
TRIMAL_PATH = "Tree/trimAI/trimal.exe"
IQ_TREE_PATH = "Tree/iqtree-3.1.2-Windows/bin/iqtree3.exe"
ALIGNED_FASTA = "Tree/phylogeny/all_sequences_aligned.fasta"
OUTPUT_DIR = "Tree/phylogeny"

# ===================== 1. trimAI修剪 =====================
print("=== 1. Running trimAI ===")
trimmed_fasta = f"{OUTPUT_DIR}/all_sequences_aligned_trimmed.fasta"

if os.path.exists(trimmed_fasta):
    print(f"  Skipping trimAI (output exists): {trimmed_fasta}")
else:
    print(f"  Input: {ALIGNED_FASTA}")
    print(f"  Output: {trimmed_fasta}")
    print(f"  Parameters: -gt 0.75 -cons 80")
    subprocess.run(f'"{TRIMAL_PATH}" -in {ALIGNED_FASTA} -out {trimmed_fasta} -gt 0.75 -cons 80', shell=True)

# ===================== 2. IQ-TREE建树 =====================
print("\n=== 2. Running IQ-TREE ===")
print("  Parameters:")
print("    -m MFP      : ModelFinder Plus (automatic model selection)")
print("    -B 1000     : 1000 ultrafast bootstrap replicates")
print("    --bnni      : Bootstrap + NNI for better tree stability")
print("    -T AUTO     : Automatic thread selection")
print("  This will take several hours...")

tree_prefix = f"{OUTPUT_DIR}/full_length"

if os.path.exists(f"{tree_prefix}.treefile"):
    print(f"  Skipping IQ-TREE (output exists): {tree_prefix}.treefile")
else:
    cmd = f'"{IQ_TREE_PATH}" -s {trimmed_fasta} -m MFP -B 1000 --bnni -T AUTO --prefix {tree_prefix}'
    result = subprocess.run(cmd, shell=True)
    print(f"\n  IQ-TREE completed with return code: {result.returncode}")

# ===================== 完成 =====================
print("\n=== 完成 ===")
print(f"  Tree file: {tree_prefix}.treefile")
print(f"  Report: {tree_prefix}.iqtree")