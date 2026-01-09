#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import os
import math
import matplotlib.pyplot as plt
from scipy.stats import t as tdist
from matplotlib.backends.backend_pdf import PdfPages
from reportlab.platypus import SimpleDocTemplate, Table, TableStyle
from reportlab.lib import colors
from reportlab.lib.pagesizes import letter

# =========================
# 文件路径 
# =========================
expr_path = "/content/RNA_seq_composite_expression.csv"
drug_path = "/content/DTP_NCI60_ZSCORE.csv"

RadSig_genes = ["VEGFC", "SERPINE1", "CHEK1", "MEG3", "SLC16A1", "PRRX1"]
RadSig_coef = np.array([
    0.0004562371,0.0000420980,0.0003481120,
    0.0000614104,0.0001242109,0.0001237774
])

outdir = "/content"
os.makedirs(outdir, exist_ok=True)

# =========================
# 加载基因表达
# =========================
expr_full = pd.read_csv(expr_path, header=2, low_memory=False)
expr_cols = [c for c in expr_full.columns if ":" in str(c)]
expr_mat = (
    expr_full[["Gene name d"] + expr_cols]
    .rename(columns={"Gene name d": "gene"})
    .dropna(subset=["gene"])
    .set_index("gene")
)
expr_mat = expr_mat.apply(pd.to_numeric, errors="coerce")

# =========================
# 加载药物 Zscore
# =========================
drug_full = pd.read_csv(drug_path, header=1, low_memory=False)
drug_cols = [c for c in drug_full.columns if ":" in str(c)]
drug_mat = (
    drug_full[["Drug name"] + drug_cols]
    .rename(columns={"Drug name": "drug"})
    .dropna(subset=["drug"])
    .set_index("drug")
)
drug_mat = drug_mat.apply(pd.to_numeric, errors="coerce")

def get_drug_series(name):
    row = drug_mat.loc[name]
    return row if not isinstance(row, pd.DataFrame) else row.mean(axis=0)

# =========================
# 检测异常药物名
# =========================
def is_bad_drug_name(name):
    s = name.strip().lower()
    if s.replace(",", "").replace("h", "").replace(".", "").isdigit():
        return True
    if "h" in s and "," in s:
        return True
    if len(s) <= 2:
        return True
    return False

# =========================
# 自动药物简写
# =========================
def drug_shortname_v3(name):
    s = name.strip()
    if "(" in s and ")" in s:
        code = s[s.find("(")+1:s.find(")")].strip()
        if not is_bad_drug_name(code):
            return code

    parts = s.replace("-", " ").split()
    if len(parts) == 0:
        return s

    stop_words = {
        "hydrochloride","hydrate","solution","sodium","potassium",
        "monohydrate","dihydrate","trihydrate","acid","ester",
        "tablet","capsule","inhibitor"
    }
    parts_clean = [p for p in parts if p.lower() not in stop_words]
    if len(parts_clean) == 0:
        parts_clean = parts[:2]

    if len(parts_clean) >= 2:
        w1, w2 = parts_clean[:2]
    else:
        w1 = parts_clean[0]
        w2 = ""

    def shrink(w):
        if len(w) >= 12:
            return w[0] + str(len(w)-2) + w[-1]
        return w

    w1 = shrink(w1)
    w2 = shrink(w2) if w2 != "" else ""
    return w1 if w2 == "" else f"{w1} {w2}"

# =========================
# Top12 固定 drug 简写
# =========================
drug_shortname_fixed = {
    "zorbamycin, dihydrochloride": "Zorbamycin",
    "JNJ-3887618": "JNJ-3887618",
    "protein: pahiv1": "Protein-PAHIV1",
    "Diethyl 2-[[4-[3-phenyl-6-(trifluoromethyl)quinoxalin-2-yl]oxybenzoyl]amino]pentanedioate": "Diethyl-quinoxalinyl",
    "4-piperidinone, 3,5-dimethyl-1-nitroso-2,6-diphenyl-": "Nitroso-piperidinone",
    "LS-150387": "LS-150387",
    "1-[2-methoxyanilinomalonyl]-3,5-dimethyl-4-[2-fluorophenyl]": "Fluoro-anilide",
    "4(3h)-quinazolinone, 3-[[(2-nitrophenyl)methylene]amino]": "Nitrophenyl-quinazolinone",
    "4,8-ethenobenzo[1,2-c:4,5-c']dipyrrole-1,3,5,7(2h,6h)- t": "Etheno-dipyrrole",
    "STK673880": "STK673880",
    "3-benzyl-4-[(E)-(2,4-dimethoxyphenyl)methylideneamino]-1H-1,2,4-triazol-5-one": "Dimethoxy-triazolone"
}

# =========================
# 计算 RadSigScore
# =========================
expr_sub = expr_mat.loc[RadSig_genes]
scores = (expr_sub.T.values * RadSig_coef).sum(axis=1)
RadSigScore = pd.Series((scores - scores.mean())/scores.std(), index=expr_mat.columns)

# =========================
# 计算药物相关性
# =========================
cors=[]
ps=[]
drugs_used = []
for drug_name, row in drug_mat.iterrows():
    if is_bad_drug_name(drug_name):
        continue

    mask = row.notna() & RadSigScore.notna()
    if mask.sum() < 5:
        continue

    r = np.corrcoef(row[mask], RadSigScore[mask])[0,1]
    tval = r * math.sqrt((mask.sum()-2)/(1-r*r)) if abs(r) != 1 else 999
    p = 2*tdist.sf(abs(tval), df=mask.sum()-2) if abs(r) != 1 else 0

    drugs_used.append(drug_name)
    cors.append(r)
    ps.append(p)

cor_df = pd.DataFrame({"drug": drugs_used, "cor": cors, "p": ps}).dropna()
cor_df = cor_df.sort_values("p")
m = len(cor_df)
cor_df["FDR"] = (cor_df["p"] * m / (np.arange(1,m+1))).clip(upper=1)

sig_drugs = set(cor_df.loc[cor_df["FDR"] < 0.05, "drug"].tolist())

# =========================
# 选出每个基因 top2
# =========================
gene_top2 = {}
for gene in RadSig_genes:
    gexpr = expr_mat.loc[gene]
    pairs=[]
    for d in sig_drugs:
        y = get_drug_series(d)
        mask = y.notna() & gexpr.notna()
        if mask.sum()<5:
            continue

        xr = gexpr[mask]
        yr = y[mask]
        r = np.corrcoef(xr, yr)[0,1]
        tval = r * math.sqrt((mask.sum()-2)/(1-r*r)) if abs(r) != 1 else 999
        p = 2*tdist.sf(abs(tval), df=mask.sum()-2) if abs(r) != 1 else 0
        pairs.append((d,r,p))

    gdf = pd.DataFrame(pairs, columns=["drug","cor","p"])
    if len(gdf)>0:
        gdf["abs_cor"] = gdf["cor"].abs()
        gene_top2[gene] = gdf.sort_values("abs_cor", ascending=False).head(2)

# =========================
# 输出 Top12 CSV
# =========================
out_top12 = []
for gene in RadSig_genes:
    if gene not in gene_top2:
        continue
    df = gene_top2[gene][["drug","cor","p"]].copy()
    df.insert(0, "gene", gene)
    df["drug_short"] = df["drug"].apply(lambda x: drug_shortname_fixed.get(x, drug_shortname_v3(x)))
    out_top12.append(df)

top12_df = pd.concat(out_top12, axis=0)
csv_path = os.path.join(outdir, "Top12_Gene_Drug_Correlation.csv")
top12_df.to_csv(csv_path, index=False)

# =========================
# 输出 PDF 表格
# =========================
pdf_table_path = os.path.join(outdir, "Top12_Gene_Drug_Table.pdf")

doc = SimpleDocTemplate(pdf_table_path, pagesize=letter)
data = [["Gene","Drug","Short Name","Correlation","P-value"]]

for _, r in top12_df.iterrows():
    data.append([
        r["gene"],
        r["drug"],
        r["drug_short"],
        f"{r['cor']:.3f}",
        f"{r['p']:.3e}"
    ])

table = Table(data, repeatRows=1)
table.setStyle(TableStyle([
    ('BACKGROUND', (0,0), (-1,0), colors.lightgrey),
    ('ALIGN', (0,0), (-1,-1), 'CENTER'),
    ('FONTNAME', (0,0), (-1,0), 'Helvetica-Bold'),
    ('GRID', (0,0), (-1,-1), 0.5, colors.grey),
]))
doc.build([table])

# =========================
# 绘图：使用固定简写
# =========================
pdf_plot_path = os.path.join(outdir, "Top12_Gene_Drug_CorPlot.pdf")

with PdfPages(pdf_plot_path) as pdf:
    fig, axes = plt.subplots(3, 4, figsize=(12, 9))
    axes = axes.flatten()
    plt.subplots_adjust(hspace=0.9, wspace=0.25)

    plot_idx = 0

    for gene in RadSig_genes:
        gdf = gene_top2[gene]
        gexpr = expr_mat.loc[gene]

        for _, row in gdf.iterrows():
            drug = row["drug"]
            drugshort = drug_shortname_fixed.get(drug, drug_shortname_v3(drug))

            r_val = row["cor"]
            p_val = row["p"]

            y = get_drug_series(drug)
            mask = y.notna() & gexpr.notna()
            x = gexpr[mask]
            yv = y[mask]

            ax = axes[plot_idx]
            ax.scatter(x, yv, s=15, color="black", alpha=0.8)

            if len(x) > 2:
                mcoef, bcoef = np.polyfit(x, yv, 1)
                xx = np.linspace(x.min(), x.max(), 60)
                yy = mcoef*xx + bcoef
                ax.plot(xx, yy, color="#1f77b4", linewidth=1.5)
                ax.fill_between(xx, yy-0.2, yy+0.2, color="gray", alpha=0.15)

            ax.set_title(f"{gene}, {drugshort}\nCor={r_val:.3f}, p={p_val:.3f}", fontsize=11, pad=6)
            ax.set_xlabel("Expression", fontsize=10)
            ax.set_ylabel("IC50", fontsize=10)

            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)

            plot_idx += 1

    pdf.savefig(fig)
    plt.close(fig)

