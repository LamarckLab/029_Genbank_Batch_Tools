import re
from pathlib import Path
import pandas as pd

# ===== 1) 输入 Excel 路径 & 列名（如与你文件不一致，请改这里） =====
excel_file = Path(r"C:\Users\Lamarck\Desktop\Accession_G_L.xlsx")
accession_col = "accession"   # 第一列：样本 accession（含版本号更好）
g_col = "G"                   # 第二列：G 蛋白氨基酸序列
l_col = "L"                   # 第三列：L 蛋白氨基酸序列
sheet_name = 0                # 读取第一个工作表，如需指定表名可填字符串

# ===== 2) 工具函数 =====
def wrap_fasta(seq: str, width: int = 60) -> str:
    """按指定宽度换行（FASTA 常用 60）"""
    return "\n".join(seq[i:i+width] for i in range(0, len(seq), width))

def clean_seq(seq: str) -> str:
    """清理序列：去空白字符，转大写。"""
    seq = re.sub(r"\s+", "", seq)  # 去掉所有空白（空格、换行、制表等）
    return seq.upper()

def _to_str(x):
    """把可能为 NaN 的对象安全转字符串；若为空返回空串。"""
    if pd.isna(x):
        return ""
    return str(x)

# ===== 3) 读取 Excel =====
df = pd.read_excel(excel_file, sheet_name=sheet_name)

# 检查列是否存在
missing_cols = [c for c in [accession_col, g_col, l_col] if c not in df.columns]
if missing_cols:
    raise ValueError(
        f"找不到这些列：{missing_cols}。当前表头为：{list(df.columns)}。\n"
        f"请修改脚本顶部的列名设置。"
    )

# ===== 4) 输出目录（与 Excel 同目录） =====
out_dir = excel_file.parent
g_fasta_path = out_dir / "G.fasta"
l_fasta_path = out_dir / "L.fasta"

# ===== 5) 写 G.fasta =====
with open(g_fasta_path, "w", encoding="utf-8") as g_out:
    for _, row in df.iterrows():
        acc = _to_str(row[accession_col]).strip()
        seq = clean_seq(_to_str(row[g_col]))
        if not acc or not seq:
            continue  # 跳过空 accession 或空序列
        g_out.write(f">{acc}\n{wrap_fasta(seq, 60)}\n")

# ===== 6) 写 L.fasta =====
with open(l_fasta_path, "w", encoding="utf-8") as l_out:
    for _, row in df.iterrows():
        acc = _to_str(row[accession_col]).strip()
        seq = clean_seq(_to_str(row[l_col]))
        if not acc or not seq:
            continue
        l_out.write(f">{acc}\n{wrap_fasta(seq, 60)}\n")

print(f"已生成：\n  {g_fasta_path}\n  {l_fasta_path}")
