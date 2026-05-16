import pandas as pd
import os
import re

# ===================== 参数 =====================
PROTEIN_DIR = "Tree/AGC_protein"
PAN_LIST_FILE = "AGC_pan_list_custom.csv"
OUTPUT_DIR = "Tree/phylogeny"
AT_SEQUENCE_FILE = "Tree/AtAGCsimpleID.fa"

# ===================== 品系名称映射 =====================
STRAIN_NAME_MAPPING = {
    'I114Hprotein.fa': 'Il14H',
    'KI11protein.fa': 'Ki11',
    'KI3protein.fa': 'Ki3',
    'KY21protein.fa': 'Ky21',
    'MS71protein.fa': 'Ms71',
    'Mo18Wprotein.fa': 'M018W',
    'TX303protein.fa': 'Tx303',
    'TZI8protein.fa': 'Tzi8',
}

def get_strain_name(filename):
    if filename in STRAIN_NAME_MAPPING:
        return STRAIN_NAME_MAPPING[filename]
    else:
        return filename.replace('protein.fa', '')

# ===================== 1. 读取65个exemplar的泛基因信息 =====================
print("=== 1. 读取65个exemplar信息 ===")
pan_df = pd.read_csv(PAN_LIST_FILE)
print("泛基因列表列名:", pan_df.columns.tolist())
print()

exemplars = pan_df['exemplar'].unique()
print("共有 %d 个独特的exemplar" % len(exemplars))

# 使用原始的panID作为ID
exemplar_to_panid = {}
for idx, row in pan_df.iterrows():
    exemplar = row['exemplar']
    pan_id = row['panID']
    if exemplar not in exemplar_to_panid:
        exemplar_to_panid[exemplar] = pan_id

print("Exemplar到Pan_ID映射示例:")
for ex, panid in list(exemplar_to_panid.items())[:5]:
    print("  %s -> %s" % (ex, panid))

# ===================== 2. 创建fam_gene.list文件 =====================
print("\n=== 2. 创建Fam_gene.list文件 ===")
os.makedirs(OUTPUT_DIR, exist_ok=True)

fam_gene_list = []
for exemplar in exemplars:
    pan_id = exemplar_to_panid[exemplar]
    trans_list = pan_df[pan_df['exemplar'] == exemplar]['trans_clean'].tolist()
    for trans in trans_list:
        fam_gene_list.append((pan_id, trans))

fam_df = pd.DataFrame(fam_gene_list, columns=['Pan_ID', 'Transcript'])
fam_df.to_csv(os.path.join(OUTPUT_DIR, "Fam_gene.list"), sep='\t', index=False)
print("已保存Fam_gene.list")

# ===================== 3. 从各品系蛋白文件中提取序列 =====================
print("\n=== 3. 从各品系蛋白文件中提取序列 ===")

def read_protein_file(filepath):
    sequences = {}
    current_id = None
    current_seq = []

    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_id is not None:
                    sequences[current_id] = ''.join(current_seq)
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        if current_id is not None:
            sequences[current_id] = ''.join(current_seq)

    return sequences

def extract_exemplar_sequences(protein_dir, pan_df, exemplar_to_panid):
    exemplar_sequences = {}
    
    # 预加载所有品系的蛋白序列
    all_strain_sequences = {}
    for filename in os.listdir(protein_dir):
        if not filename.endswith('.fa') and not filename.endswith('.fasta'):
            continue
        strain_name = get_strain_name(filename)
        filepath = os.path.join(protein_dir, filename)
        all_strain_sequences[strain_name] = read_protein_file(filepath)

    for exemplar in exemplar_to_panid.keys():
        pan_id = exemplar_to_panid[exemplar]
        trans_list = pan_df[pan_df['exemplar'] == exemplar]['trans_clean'].tolist()
        
        found = False
        for trans in trans_list:
            # 尝试不同的后缀
            for suffix in ['_P001', '_P002', '_P003', '_P004', '_P005', '']:
                protein_id = trans + suffix
                
                # 在所有品系中查找
                for strain_name, sequences in all_strain_sequences.items():
                    if protein_id in sequences:
                        exemplar_sequences[pan_id] = sequences[protein_id]
                        found = True
                        break
                if found:
                    break
            if found:
                break

    return exemplar_sequences

exemplar_sequences = extract_exemplar_sequences(PROTEIN_DIR, pan_df, exemplar_to_panid)
print("\n成功提取 %d 个exemplar的序列" % len(exemplar_sequences))

output_fasta = os.path.join(OUTPUT_DIR, "maize.pep_fulllength.fasta")
with open(output_fasta, 'w') as f:
    for pan_id in sorted(exemplar_sequences.keys()):
        f.write(">%s\n" % pan_id)
        f.write(exemplar_sequences[pan_id] + "\n")
print("已保存全长蛋白序列")

# ===================== 4. 准备拟南芥序列 =====================
print("\n=== 4. 准备拟南芥序列 ===")

at_sequences = {}
current_id = None
current_seq = []

with open(AT_SEQUENCE_FILE, 'r') as f:
    for line in f:
        line = line.strip()
        if line.startswith('>'):
            if current_id is not None:
                at_sequences[current_id] = ''.join(current_seq)
            current_id = line[1:].split()[0]
            current_seq = []
        else:
            current_seq.append(line)
    if current_id is not None:
        at_sequences[current_id] = ''.join(current_seq)

print("读取到 %d 个拟南芥AGC序列" % len(at_sequences))

at_output = os.path.join(OUTPUT_DIR, "At_AGC.pep.fasta")
with open(at_output, 'w') as f:
    for at_id in sorted(at_sequences.keys()):
        f.write(">AT_%s\n" % at_id)
        f.write(at_sequences[at_id] + "\n")
print("已保存拟南芥序列")

# ===================== 5. 合并所有序列用于建树 =====================
print("\n=== 5. 合并序列 ===")

all_sequences = {}
all_sequences.update(exemplar_sequences)
for at_id, seq in at_sequences.items():
    all_sequences["AT_" + at_id] = seq

combined_output = os.path.join(OUTPUT_DIR, "all_sequences_for_tree.fasta")
with open(combined_output, 'w') as f:
    for seq_id in sorted(all_sequences.keys(), key=lambda x: (x.startswith('AT_'), x)):
        f.write(">%s\n" % seq_id)
        f.write(all_sequences[seq_id] + "\n")

print("已保存合并序列")
print("共 %d 条序列（65个玉米exemplar + %d个拟南芥）" % (len(all_sequences), len(at_sequences)))

# ===================== 6. 输出统计信息 =====================
print("\n=== 6. 统计信息 ===")
print("玉米exemplar序列: %d" % len(exemplar_sequences))
print("拟南芥序列: %d" % len(at_sequences))
print("合并序列: %d" % len(all_sequences))

missing_exemplars = set(exemplar_to_panid.keys()) - set(exemplar_sequences.keys())
if missing_exemplars:
    print("\n警告: %d 个exemplar未能提取到序列:" % len(missing_exemplars))
    for ex in list(missing_exemplars)[:5]:
        print("  - %s" % ex)

print("\n=== 数据准备完成 ===")
print("输出目录: %s" % OUTPUT_DIR)