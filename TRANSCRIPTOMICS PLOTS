"""
Transcriptomic Data - All 13 Publication-Ready Plots
======================================================
Usage: python transcriptomic_plots.py
Output: 13 PNG files saved in the current directory
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LogNorm
from collections import Counter
import warnings
warnings.filterwarnings('ignore')

# ── Load Data ──────────────────────────────────────────────────────────────────
df = pd.read_csv('UPDOWN.csv')

# Colour palette
UP_COLOR   = '#E8593C'   # coral-red  → upregulated
DOWN_COLOR = '#3B8BD4'   # steel-blue → downregulated
NS_COLOR   = '#AAAAAA'   # grey       → not significant / neutral

# Common style
plt.rcParams.update({
    'font.family'     : 'DejaVu Sans',
    'axes.spines.top' : False,
    'axes.spines.right': False,
    'axes.grid'       : False,
    'figure.dpi'      : 150,
})

print("Generating 13 plots from your transcriptomic data…\n")

# ==============================================================================
# PLOT 1 — Volcano Plot
# ==============================================================================
fig, ax = plt.subplots(figsize=(7, 6))

df['neglog10p'] = -np.log10(df['P Value'])
colors = df['Regulation'].map({'UP': UP_COLOR, 'DOWN': DOWN_COLOR}).fillna(NS_COLOR)

ax.scatter(df['log2FoldChange'], df['neglog10p'],
           c=colors, s=18, alpha=0.7, linewidths=0)

ax.axvline(x=0,  color='black', linestyle='--', linewidth=0.8, alpha=0.4)
ax.axhline(y=-np.log10(0.05), color='grey', linestyle=':', linewidth=0.8, alpha=0.6)

ax.set_xlabel('log₂ Fold Change', fontsize=12)
ax.set_ylabel('−log₁₀ (p-value)',  fontsize=12)
ax.set_title('Volcano Plot', fontsize=14, fontweight='bold', pad=10)

up_patch   = mpatches.Patch(color=UP_COLOR,   label=f'UP ({(df.Regulation=="UP").sum()})')
down_patch = mpatches.Patch(color=DOWN_COLOR, label=f'DOWN ({(df.Regulation=="DOWN").sum()})')
ax.legend(handles=[up_patch, down_patch], frameon=False, fontsize=10)

plt.tight_layout()
plt.savefig('Plot01_Volcano.png', dpi=150, bbox_inches='tight')
plt.close()
print("✓ Plot 01 — Volcano Plot")

# ==============================================================================
# PLOT 2 — MA Plot
# ==============================================================================
fig, ax = plt.subplots(figsize=(7, 6))

df['A'] = (np.log2(df['144 Expression'] + 1) + np.log2(df['187 Expression'] + 1)) / 2
df['M'] = df['log2FoldChange']

ax.scatter(df['A'], df['M'], c=colors, s=18, alpha=0.7, linewidths=0)
ax.axhline(y=0, color='black', linestyle='--', linewidth=0.8, alpha=0.5)

ax.set_xlabel('A  (average log₂ expression)', fontsize=12)
ax.set_ylabel('M  (log₂ fold change)',         fontsize=12)
ax.set_title('MA Plot', fontsize=14, fontweight='bold', pad=10)
ax.legend(handles=[up_patch, down_patch], frameon=False, fontsize=10)

plt.tight_layout()
plt.savefig('Plot02_MA_Plot.png', dpi=150, bbox_inches='tight')
plt.close()
print("✓ Plot 02 — MA Plot")

# ==============================================================================
# PLOT 3 — Expression Scatter (Sample 144 vs 187)
# ==============================================================================
fig, ax = plt.subplots(figsize=(6, 6))

ax.scatter(df['144 Expression'], df['187 Expression'],
           c=colors, s=18, alpha=0.65, linewidths=0)

lim = max(df['144 Expression'].max(), df['187 Expression'].max()) * 1.05
ax.plot([0, lim], [0, lim], 'k--', linewidth=0.8, alpha=0.4, label='y = x')
ax.set_xlim(0, lim); ax.set_ylim(0, lim)

ax.set_xlabel('Sample 144 Expression', fontsize=12)
ax.set_ylabel('Sample 187 Expression', fontsize=12)
ax.set_title('Expression Scatter\n(Sample 144 vs 187)', fontsize=14, fontweight='bold', pad=10)
ax.legend(handles=[up_patch, down_patch,
                   mpatches.Patch(color='none', label='')],
          frameon=False, fontsize=10)

plt.tight_layout()
plt.savefig('Plot03_Expression_Scatter.png', dpi=150, bbox_inches='tight')
plt.close()
print("✓ Plot 03 — Expression Scatter")

# ==============================================================================
# PLOT 4 — UP / DOWN Bar Chart
# ==============================================================================
fig, ax = plt.subplots(figsize=(5, 5))

counts = df['Regulation'].value_counts()
bars = ax.bar(['UP', 'DOWN'], [counts['UP'], counts['DOWN']],
              color=[UP_COLOR, DOWN_COLOR], width=0.5, edgecolor='white')

for bar, val in zip(bars, [counts['UP'], counts['DOWN']]):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 5,
            str(val), ha='center', va='bottom', fontsize=12, fontweight='bold')

ax.set_ylabel('Number of DEGs', fontsize=12)
ax.set_title('Upregulated vs Downregulated Genes', fontsize=13, fontweight='bold', pad=10)
ax.set_ylim(0, counts['UP'] * 1.15)

plt.tight_layout()
plt.savefig('Plot04_UP_DOWN_Bar.png', dpi=150, bbox_inches='tight')
plt.close()
print("✓ Plot 04 — UP/DOWN Bar Chart")

# ==============================================================================
# PLOT 5 — Heatmap (Top 50 DEGs by |log2FC|)
# ==============================================================================
fig, ax = plt.subplots(figsize=(5, 12))

top50 = df.nlargest(50, 'log2FoldChange').copy()
top50 = pd.concat([
    df.nlargest(25, 'log2FoldChange'),
    df.nsmallest(25, 'log2FoldChange')
]).drop_duplicates()

# Normalise each gene to [0,1] across the two samples
expr = top50[['144 Expression', '187 Expression']].copy()
expr_norm = expr.div(expr.max(axis=1), axis=0)

labels = top50['Gene Names'].fillna(top50['Unigene_ID']).str[:25].values
reg    = top50['Regulation'].values

im = ax.imshow(expr_norm.values, aspect='auto', cmap='RdYlBu_r',
               vmin=0, vmax=1, interpolation='nearest')

ax.set_xticks([0, 1])
ax.set_xticklabels(['Sample 144', 'Sample 187'], fontsize=11)
ax.set_yticks(range(len(labels)))
ax.set_yticklabels([f"{'↑' if r=='UP' else '↓'} {l}" for l, r in zip(labels, reg)],
                   fontsize=7)

plt.colorbar(im, ax=ax, label='Normalised Expression', shrink=0.4)
ax.set_title('Top 50 DEGs Heatmap\n(by |log₂FC|)', fontsize=13, fontweight='bold', pad=10)

plt.tight_layout()
plt.savefig('Plot05_Heatmap.png', dpi=150, bbox_inches='tight')
plt.close()
print("✓ Plot 05 — Heatmap")

# ==============================================================================
# PLOT 6 — Box / Violin Plot (expression per sample)
# ==============================================================================
fig, axes = plt.subplots(1, 2, figsize=(9, 6), sharey=False)

for ax, col, color, label in zip(
        axes,
        ['144 Expression', '187 Expression'],
        ['#5B9BD5', '#E67E22'],
        ['Sample 144', 'Sample 187']):

    data_log = np.log2(df[col] + 1)
    parts = ax.violinplot(data_log, positions=[0], widths=0.6,
                          showmedians=True, showextrema=True)
    for pc in parts['bodies']:
        pc.set_facecolor(color); pc.set_alpha(0.6)
    parts['cmedians'].set_color('black'); parts['cmedians'].set_linewidth(2)
    parts['cmaxes'].set_color('grey');  parts['cmins'].set_color('grey')
    parts['cbars'].set_color('grey')

    ax.set_xticks([0]); ax.set_xticklabels([label], fontsize=11)
    ax.set_ylabel('log₂ (Expression + 1)', fontsize=10)
    ax.set_title(label, fontsize=12, fontweight='bold')

fig.suptitle('Expression Distribution per Sample', fontsize=13, fontweight='bold', y=1.01)
plt.tight_layout()
plt.savefig('Plot06_Violin_Box.png', dpi=150, bbox_inches='tight')
plt.close()
print("✓ Plot 06 — Violin / Box Plot")

# ==============================================================================
# PLOT 7 — log2FC Distribution Histogram
# ==============================================================================
fig, ax = plt.subplots(figsize=(7, 5))

up_fc   = df.loc[df.Regulation=='UP',   'log2FoldChange']
down_fc = df.loc[df.Regulation=='DOWN', 'log2FoldChange']

bins = np.linspace(df['log2FoldChange'].min(), df['log2FoldChange'].max(), 35)
ax.hist(up_fc,   bins=bins, color=UP_COLOR,   alpha=0.75, label=f'UP ({len(up_fc)})',   edgecolor='white')
ax.hist(down_fc, bins=bins, color=DOWN_COLOR, alpha=0.75, label=f'DOWN ({len(down_fc)})', edgecolor='white')
ax.axvline(x=0, color='black', linestyle='--', linewidth=0.9, alpha=0.5)

ax.set_xlabel('log₂ Fold Change', fontsize=12)
ax.set_ylabel('Number of Genes',  fontsize=12)
ax.set_title('log₂FC Distribution', fontsize=14, fontweight='bold', pad=10)
ax.legend(frameon=False, fontsize=10)

plt.tight_layout()
plt.savefig('Plot07_log2FC_Histogram.png', dpi=150, bbox_inches='tight')
plt.close()
print("✓ Plot 07 — log2FC Distribution")

# ==============================================================================
# PLOT 8 — GO Biological Process Bar Chart (Top 15)
# ==============================================================================
def parse_go_terms(series):
    """Split semicolon-separated GO terms, return Counter of term names."""
    terms = []
    for val in series.dropna():
        for t in str(val).split(';'):
            t = t.strip()
            if '~' in t:
                terms.append(t.split('~', 1)[1].strip())
            elif t:
                terms.append(t)
    return Counter(terms)

fig, ax = plt.subplots(figsize=(9, 7))

go_up   = parse_go_terms(df.loc[df.Regulation=='UP',   'GO_BP'])
go_down = parse_go_terms(df.loc[df.Regulation=='DOWN', 'GO_BP'])
all_go  = Counter()
for k in set(list(go_up.keys()) + list(go_down.keys())):
    all_go[k] = go_up.get(k, 0) + go_down.get(k, 0)

top_terms = [t for t, _ in all_go.most_common(15)][::-1]
y = np.arange(len(top_terms))

ax.barh(y - 0.2, [go_up.get(t, 0)   for t in top_terms], 0.4,
        color=UP_COLOR,   label='UP',   alpha=0.85)
ax.barh(y + 0.2, [go_down.get(t, 0) for t in top_terms], 0.4,
        color=DOWN_COLOR, label='DOWN', alpha=0.85)

ax.set_yticks(y)
ax.set_yticklabels([t[:55] for t in top_terms], fontsize=8)
ax.set_xlabel('Gene Count', fontsize=11)
ax.set_title('Top 15 GO Biological Process Terms', fontsize=13, fontweight='bold', pad=10)
ax.legend(frameon=False, fontsize=10)

plt.tight_layout()
plt.savefig('Plot08_GO_BP_Bar.png', dpi=150, bbox_inches='tight')
plt.close()
print("✓ Plot 08 — GO-BP Bar Chart")

# ==============================================================================
# PLOT 9 — GO Molecular Function Dot Plot
# ==============================================================================
fig, ax = plt.subplots(figsize=(9, 7))

go_mf_all = parse_go_terms(df['GO_MF'])
top_mf    = [t for t, _ in go_mf_all.most_common(15)][::-1]

counts_mf = [go_mf_all[t] for t in top_mf]
ratios_mf = [c / len(df) * 100 for c in counts_mf]

sc = ax.scatter(ratios_mf,
                range(len(top_mf)),
                s=[c * 18 for c in counts_mf],
                c=ratios_mf, cmap='viridis', alpha=0.85, edgecolors='white', linewidths=0.5)

ax.set_yticks(range(len(top_mf)))
ax.set_yticklabels([t[:60] for t in top_mf], fontsize=8)
ax.set_xlabel('Gene Ratio (%)', fontsize=11)
ax.set_title('Top 15 GO Molecular Function Terms\n(dot size = gene count)', fontsize=13, fontweight='bold', pad=10)

plt.colorbar(sc, ax=ax, label='Gene Ratio (%)', shrink=0.5)
plt.tight_layout()
plt.savefig('Plot09_GO_MF_Dotplot.png', dpi=150, bbox_inches='tight')
plt.close()
print("✓ Plot 09 — GO-MF Dot Plot")

# ==============================================================================
# PLOT 10 — GO Annotation Coverage Pie Chart
# ==============================================================================
fig, ax = plt.subplots(figsize=(6, 6))

total    = len(df)
bp_count = df['GO_BP'].notna().sum()
mf_count = df['GO_MF'].notna().sum()
cc_count = df['GO_CC'].notna().sum()

# Genes with at least one GO annotation
annotated   = df[['GO_BP','GO_MF','GO_CC']].notna().any(axis=1).sum()
unannotated = total - annotated

sizes  = [bp_count, mf_count - bp_count, cc_count, unannotated]
labels = [f'GO-BP\n({bp_count})',
          f'GO-MF only\n({mf_count - bp_count})',
          f'GO-CC only\n({cc_count})',
          f'No annotation\n({unannotated})']
colors_pie = ['#E8593C', '#3B8BD4', '#9B59B6', '#AAAAAA']

wedges, texts, autotexts = ax.pie(
    sizes, labels=labels, colors=colors_pie,
    autopct='%1.1f%%', startangle=140,
    pctdistance=0.78, textprops={'fontsize': 9},
    wedgeprops={'edgecolor': 'white', 'linewidth': 1.5})

for at in autotexts:
    at.set_fontsize(8); at.set_color('white'); at.set_fontweight('bold')

ax.set_title('GO Annotation Coverage\n(489 DEGs)', fontsize=13, fontweight='bold', pad=10)

plt.tight_layout()
plt.savefig('Plot10_GO_Annotation_Coverage.png', dpi=150, bbox_inches='tight')
plt.close()
print("✓ Plot 10 — GO Annotation Coverage Pie")

# ==============================================================================
# PLOT 11 — % Identity Distribution (BLAST)
# ==============================================================================
fig, ax = plt.subplots(figsize=(7, 5))

id_up   = df.loc[df.Regulation=='UP',   '%Identity'].dropna()
id_down = df.loc[df.Regulation=='DOWN', '%Identity'].dropna()

bins_id = np.linspace(0, 100, 30)
ax.hist(id_up,   bins=bins_id, color=UP_COLOR,   alpha=0.75, label=f'UP   (n={len(id_up)})',   edgecolor='white')
ax.hist(id_down, bins=bins_id, color=DOWN_COLOR, alpha=0.75, label=f'DOWN (n={len(id_down)})', edgecolor='white')

ax.axvline(x=id_up.median(),   color=UP_COLOR,   linestyle='--', linewidth=1.2,
           label=f'Median UP={id_up.median():.1f}%')
ax.axvline(x=id_down.median(), color=DOWN_COLOR, linestyle='--', linewidth=1.2,
           label=f'Median DOWN={id_down.median():.1f}%')

ax.set_xlabel('BLAST % Identity', fontsize=12)
ax.set_ylabel('Number of Genes',  fontsize=12)
ax.set_title('BLAST % Identity Distribution', fontsize=14, fontweight='bold', pad=10)
ax.legend(frameon=False, fontsize=9)

plt.tight_layout()
plt.savefig('Plot11_BLAST_Identity.png', dpi=150, bbox_inches='tight')
plt.close()
print("✓ Plot 11 — BLAST % Identity Distribution")

# ==============================================================================
# PLOT 12 — E-value vs % Identity Scatter (BLAST quality)
# ==============================================================================
fig, ax = plt.subplots(figsize=(7, 6))

blast_df = df.dropna(subset=['E-value', '%Identity']).copy()
blast_df['neglog10_evalue'] = -np.log10(blast_df['E-value'] + 1e-300)
b_colors = blast_df['Regulation'].map({'UP': UP_COLOR, 'DOWN': DOWN_COLOR}).fillna(NS_COLOR)

ax.scatter(blast_df['%Identity'], blast_df['neglog10_evalue'],
           c=b_colors, s=20, alpha=0.65, linewidths=0)

ax.set_xlabel('BLAST % Identity', fontsize=12)
ax.set_ylabel('−log₁₀ (E-value)',  fontsize=12)
ax.set_title('E-value vs % Identity\n(BLAST Hit Quality)', fontsize=13, fontweight='bold', pad=10)

up_patch2   = mpatches.Patch(color=UP_COLOR,   label=f'UP   (n={(blast_df.Regulation=="UP").sum()})')
down_patch2 = mpatches.Patch(color=DOWN_COLOR, label=f'DOWN (n={(blast_df.Regulation=="DOWN").sum()})')
ax.legend(handles=[up_patch2, down_patch2], frameon=False, fontsize=10)

plt.tight_layout()
plt.savefig('Plot12_Evalue_vs_Identity.png', dpi=150, bbox_inches='tight')
plt.close()
print("✓ Plot 12 — E-value vs % Identity Scatter")

# ==============================================================================
# PLOT 13 — Species Hit Distribution (Top 10)
# ==============================================================================
fig, ax = plt.subplots(figsize=(9, 6))

species = df['Organism Name'].dropna()
# Shorten names: take everything before the first parenthesis
species_short = species.str.replace(r'\s*\(.*?\)', '', regex=True).str.strip()
top_species   = species_short.value_counts().head(10)

colors_sp = plt.cm.tab10(np.linspace(0, 1, 10))
bars = ax.barh(range(len(top_species)), top_species.values[::-1],
               color=colors_sp[::-1], edgecolor='white', alpha=0.88)

ax.set_yticks(range(len(top_species)))
ax.set_yticklabels([s[:50] for s in top_species.index[::-1]], fontsize=9)
ax.set_xlabel('Number of DEGs', fontsize=11)
ax.set_title('Top 10 Species in BLAST Hits', fontsize=13, fontweight='bold', pad=10)

for bar, val in zip(bars, top_species.values[::-1]):
    ax.text(val + 0.5, bar.get_y() + bar.get_height()/2,
            str(val), va='center', fontsize=9)

plt.tight_layout()
plt.savefig('Plot13_Species_Distribution.png', dpi=150, bbox_inches='tight')
plt.close()
print("✓ Plot 13 — Species Hit Distribution")

# ==============================================================================
# Summary
# ==============================================================================
print("\n" + "="*55)
print("ALL 13 PLOTS SAVED SUCCESSFULLY!")
print("="*55)
plots = [
    "Plot01_Volcano.png",
    "Plot02_MA_Plot.png",
    "Plot03_Expression_Scatter.png",
    "Plot04_UP_DOWN_Bar.png",
    "Plot05_Heatmap.png",
    "Plot06_Violin_Box.png",
    "Plot07_log2FC_Histogram.png",
    "Plot08_GO_BP_Bar.png",
    "Plot09_GO_MF_Dotplot.png",
    "Plot10_GO_Annotation_Coverage.png",
    "Plot11_BLAST_Identity.png",
    "Plot12_Evalue_vs_Identity.png",
    "Plot13_Species_Distribution.png",
]
for i, p in enumerate(plots, 1):
    print(f"  {i:>2}. {p}")
print("\nMake sure UPDOWN.csv is in the same folder as this script.")
