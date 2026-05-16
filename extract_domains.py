import os
import subprocess
import pandas as pd

# ===================== 参数 =====================
HMMSEARCH_PATH = "Tree/hmmer3/hmmsearch.exe"
HMMPRESS_PATH = "Tree/hmmer3/hmmpress.exe"
SEQKIT_PATH = "Tree/seqkit.exe"
PF00069_HMM = "Tree/PF00069.hmm"
PF00433_HMM = "Tree/PF00433.hmm"
INPUT_FASTA = "Tree/phylogeny/all_sequences_for_tree.fasta"
OUTPUT_DIR = "Tree/phylogeny"

MAX_GAP = 50

# ===================== 1. 使用hmmpress预处理HMM文件 =====================
print("=== 1. 使用hmmpress预处理HMM文件 ===")

cmd_press1 = f'"{HMMPRESS_PATH}" {PF00069_HMM}'
print("  预处理PF00069.hmm...")
subprocess.run(cmd_press1, shell=True)

cmd_press2 = f'"{HMMPRESS_PATH}" {PF00433_HMM}'
print("  预处理PF00433.hmm...")
subprocess.run(cmd_press2, shell=True)

# ===================== 2. 使用hmmsearch扫描结构域 =====================
print("\n=== 2. 使用hmmsearch扫描结构域 ===")

pf00069_out = os.path.join(OUTPUT_DIR, "hmmsearch_PF00069.txt")
cmd1 = f'"{HMMSEARCH_PATH}" --domtblout {pf00069_out} {PF00069_HMM} {INPUT_FASTA}'
print("  扫描PF00069...")
subprocess.run(cmd1, shell=True)

pf00433_out = os.path.join(OUTPUT_DIR, "hmmsearch_PF00433.txt")
cmd2 = f'"{HMMSEARCH_PATH}" --domtblout {pf00433_out} {PF00433_HMM} {INPUT_FASTA}'
print("  扫描PF00433...")
subprocess.run(cmd2, shell=True)

# ===================== 3. 解析hmmsearch结果 =====================
print("\n=== 3. 解析hmmsearch结果 ===")

def parse_hmmsearch_output(filepath):
    domains = {}
    if not os.path.exists(filepath):
        return domains
    with open(filepath, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) >= 20:
                seq_id = parts[0]
                e_value = float(parts[6])
                ali_start = int(parts[19])
                ali_end = int(parts[20])
                
                if seq_id not in domains:
                    domains[seq_id] = []
                domains[seq_id].append({
                    'e_value': e_value,
                    'start': ali_start,
                    'end': ali_end
                })
    return domains

pf00069_domains = parse_hmmsearch_output(pf00069_out)
pf00433_domains = parse_hmmsearch_output(pf00433_out)

print(f"  PF00069找到 {len(pf00069_domains)} 个序列")
print(f"  PF00433找到 {len(pf00433_domains)} 个序列")

# ===================== 4. 合并结构域坐标 =====================
print("\n=== 4. 合并结构域坐标 ===")

all_seq_ids = set(pf00069_domains.keys()).union(set(pf00433_domains.keys()))
print(f"  总共有 {len(all_seq_ids)} 个序列")

core_bed = []
full_bed = []

for seq_id in all_seq_ids:
    has_pf00069 = seq_id in pf00069_domains
    has_pf00433 = seq_id in pf00433_domains
    
    if has_pf00069:
        best_pf00069 = min(pf00069_domains[seq_id], key=lambda x: x['e_value'])
        core_start = best_pf00069['start']
        core_end = best_pf00069['end']
        core_bed.append((seq_id, core_start, core_end))
        
        if has_pf00433:
            best_pf00433 = min(pf00433_domains[seq_id], key=lambda x: x['e_value'])
            pf00433_start = best_pf00433['start']
            pf00433_end = best_pf00433['end']
            
            gap = pf00433_start - core_end
            if gap <= MAX_GAP and gap > 0:
                full_start = core_start
                full_end = pf00433_end
                full_bed.append((seq_id, full_start, full_end))
            else:
                full_bed.append((seq_id, core_start, core_end))
        else:
            full_bed.append((seq_id, core_start, core_end))

core_bed_file = os.path.join(OUTPUT_DIR, "AGC_core_domain.bed")
with open(core_bed_file, 'w') as f:
    for seq_id, start, end in core_bed:
        f.write(f"{seq_id}\t{start}\t{end}\n")

full_bed_file = os.path.join(OUTPUT_DIR, "AGC_full_domain.bed")
with open(full_bed_file, 'w') as f:
    for seq_id, start, end in full_bed:
        f.write(f"{seq_id}\t{start}\t{end}\n")

print(f"  生成AGC_core_domain.bed: {len(core_bed)} 个序列")
print(f"  生成AGC_full_domain.bed: {len(full_bed)} 个序列")

# ===================== 5. 使用seqkit提取结构域序列 =====================
print("\n=== 5. 使用seqkit提取结构域序列 ===")

core_fasta = os.path.join(OUTPUT_DIR, "AGC_core_domain.fasta")
cmd3 = f'"{SEQKIT_PATH}" subseq --bed {core_bed_file} -o {core_fasta} {INPUT_FASTA}'
print("  提取核心结构域序列...")
subprocess.run(cmd3, shell=True)

full_fasta = os.path.join(OUTPUT_DIR, "AGC_full_domain.fasta")
cmd4 = f'"{SEQKIT_PATH}" subseq --bed {full_bed_file} -o {full_fasta} {INPUT_FASTA}'
print("  提取完整结构域序列...")
subprocess.run(cmd4, shell=True)

# ===================== 6. 统计输出 =====================
print("\n=== 6. 统计信息 ===")

def count_sequences(filepath):
    if not os.path.exists(filepath):
        return 0
    with open(filepath, 'r') as f:
        content = f.read()
        return content.count('>')

print(f"  AGC_core_domain.fasta: {count_sequences(core_fasta)} 条序列")
print(f"  AGC_full_domain.fasta: {count_sequences(full_fasta)} 条序列")

print("\n=== 结构域提取完成 ===")