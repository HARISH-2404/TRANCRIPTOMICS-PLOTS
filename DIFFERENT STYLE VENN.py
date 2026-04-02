"""
Venn Diagrams for Transcriptomic DEG Data
==========================================
3 Venn diagrams drawn with pure matplotlib (no extra library needed):
  1. UP vs DOWN regulated DEGs
  2. GO annotation categories (BP / MF / CC) — 3-way
  3. UP/DOWN DEGs with GO annotation and BLAST hits — 3-way

Usage: python venn_diagrams.py
Output: Venn_DEGs.png  (all 3 panels in one figure)
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Circle
from matplotlib.collections import PatchCollection

# ── Load data ──────────────────────────────────────────────────────────────────
df = pd.read_csv('UPDOWN.csv')

up_set   = set(df.loc[df['Regulation'] == 'UP',   'Unigene_ID'])
down_set = set(df.loc[df['Regulation'] == 'DOWN',  'Unigene_ID'])
blast_set= set(df.loc[df['%Identity'].notna(),     'Unigene_ID'])
go_bp    = set(df.loc[df['GO_BP'].notna(),          'Unigene_ID'])
go_mf    = set(df.loc[df['GO_MF'].notna(),          'Unigene_ID'])
go_cc    = set(df.loc[df['GO_CC'].notna(),          'Unigene_ID'])
go_any   = go_bp | go_mf | go_cc

# ── Helpers ────────────────────────────────────────────────────────────────────
def draw_circle(ax, cx, cy, r, color, alpha=0.35, lw=2):
    c = Circle((cx, cy), r, color=color, alpha=alpha, linewidth=lw,
               edgecolor=color, zorder=2)
    ax.add_patch(c)

def label(ax, x, y, text, fontsize=12, color='#222222', bold=False):
    ax.text(x, y, text, ha='center', va='center', fontsize=fontsize,
            color=color, fontweight='bold' if bold else 'normal', zorder=5)

# ── Colour palette ─────────────────────────────────────────────────────────────
UP_C   = '#E8593C'
DOWN_C = '#3B8BD4'
BP_C   = '#9B59B6'
MF_C   = '#27AE60'
CC_C   = '#F39C12'
BLAST_C= '#2980B9'
GOAN_C = '#8E44AD'

# ══════════════════════════════════════════════════════════════════════════════
fig, axes = plt.subplots(1, 3, figsize=(18, 6))
fig.patch.set_facecolor('white')

# ─────────────────────────────────────────────────────────────────────────────
# PANEL 1 — UP vs DOWN (2-circle Venn)
# ─────────────────────────────────────────────────────────────────────────────
ax = axes[0]
ax.set_xlim(-2.5, 2.5); ax.set_ylim(-2, 2)
ax.set_aspect('equal'); ax.axis('off')

# These two groups are mutually exclusive (no overlap)
draw_circle(ax, -0.85, 0, 1.25, UP_C,   alpha=0.45)
draw_circle(ax,  0.85, 0, 1.25, DOWN_C, alpha=0.45)

# Set labels inside circles
label(ax, -1.3, 0.15, str(len(up_set)),   fontsize=22, color='white', bold=True)
label(ax, -1.3, -0.3, 'Unique',           fontsize=9,  color='white')
label(ax,  1.3, 0.15, str(len(down_set)), fontsize=22, color='white', bold=True)
label(ax,  1.3, -0.3, 'Unique',           fontsize=9,  color='white')

# Overlap label (0 — circles barely touch visually to show exclusivity)
label(ax, 0, 0, '0', fontsize=16, color='#555555')

# Circle titles above
label(ax, -1.3,  1.45, 'UP-regulated',   fontsize=12, color=UP_C,   bold=True)
label(ax,  1.3,  1.45, 'DOWN-regulated', fontsize=12, color=DOWN_C, bold=True)

ax.set_title('Venn 1: UP vs DOWN DEGs\n(489 total DEGs)',
             fontsize=13, fontweight='bold', pad=14)

# ─────────────────────────────────────────────────────────────────────────────
# PANEL 2 — GO annotation categories: BP / MF / CC (3-circle Venn)
# ─────────────────────────────────────────────────────────────────────────────
ax = axes[1]
ax.set_xlim(-2.5, 2.5); ax.set_ylim(-2.2, 2.4)
ax.set_aspect('equal'); ax.axis('off')

r = 1.25
# Equilateral triangle centres
cx_bp = -0.75; cy_bp =  0.55
cx_mf =  0.75; cy_mf =  0.55
cx_cc =  0.0;  cy_cc = -0.75

draw_circle(ax, cx_bp, cy_bp, r, BP_C, alpha=0.30)
draw_circle(ax, cx_mf, cy_mf, r, MF_C, alpha=0.30)
draw_circle(ax, cx_cc, cy_cc, r, CC_C, alpha=0.30)

# Compute set intersections
only_bp     = len(go_bp - go_mf - go_cc)
only_mf     = len(go_mf - go_bp - go_cc)
only_cc     = len(go_cc - go_bp - go_mf)
bp_mf_only  = len((go_bp & go_mf) - go_cc)
bp_cc_only  = len((go_bp & go_cc) - go_mf)
mf_cc_only  = len((go_mf & go_cc) - go_bp)
all_three   = len(go_bp & go_mf & go_cc)

# Region labels (approximate centres)
label(ax, -1.35,  1.15, str(only_bp),    fontsize=15, bold=True)
label(ax,  1.35,  1.15, str(only_mf),    fontsize=15, bold=True)
label(ax,  0.0,  -1.55, str(only_cc),    fontsize=15, bold=True)
label(ax,  0.0,   0.95, str(bp_mf_only), fontsize=13)
label(ax, -0.78, -0.38, str(bp_cc_only), fontsize=13)
label(ax,  0.78, -0.38, str(mf_cc_only), fontsize=13)
label(ax,  0.0,   0.18, str(all_three),  fontsize=14, bold=True)

# Circle titles
label(ax, -1.7,  1.90, 'GO-BP', fontsize=12, color=BP_C, bold=True)
label(ax,  1.7,  1.90, 'GO-MF', fontsize=12, color=MF_C, bold=True)
label(ax,  0.0, -2.05, 'GO-CC', fontsize=12, color=CC_C, bold=True)

ax.set_title('Venn 2: GO Annotation Categories\n(BP / MF / CC)',
             fontsize=13, fontweight='bold', pad=14)

# ─────────────────────────────────────────────────────────────────────────────
# PANEL 3 — UP-DEGs / DOWN-DEGs / BLAST-annotated (3-circle Venn)
# ─────────────────────────────────────────────────────────────────────────────
ax = axes[2]
ax.set_xlim(-2.5, 2.5); ax.set_ylim(-2.2, 2.4)
ax.set_aspect('equal'); ax.axis('off')

cx_up   = -0.75; cy_up   =  0.55
cx_dn   =  0.75; cy_dn   =  0.55
cx_bl   =  0.0;  cy_bl   = -0.75

draw_circle(ax, cx_up, cy_up, r, UP_C,    alpha=0.30)
draw_circle(ax, cx_dn, cy_dn, r, DOWN_C,  alpha=0.30)
draw_circle(ax, cx_bl, cy_bl, r, BLAST_C, alpha=0.30)

# Compute intersections
up_blast   = up_set & blast_set
dn_blast   = down_set & blast_set
up_dn      = up_set & down_set          # always 0 (mutually exclusive)
all_three2 = up_set & down_set & blast_set  # 0

only_up    = len(up_set - blast_set)
only_dn    = len(down_set - blast_set)
only_bl    = len(blast_set - up_set - down_set)
up_bl_only = len(up_blast - down_set)
dn_bl_only = len(dn_blast - up_set)
up_dn_only = len(up_dn - blast_set)    # 0

label(ax, -1.35,  1.15, str(only_up),    fontsize=15, bold=True)
label(ax,  1.35,  1.15, str(only_dn),    fontsize=15, bold=True)
label(ax,  0.0,  -1.55, str(only_bl),    fontsize=15, bold=True)
label(ax,  0.0,   0.95, str(up_dn_only), fontsize=13)
label(ax, -0.78, -0.38, str(up_bl_only), fontsize=13, bold=True)
label(ax,  0.78, -0.38, str(dn_bl_only), fontsize=13, bold=True)
label(ax,  0.0,   0.18, str(len(all_three2)), fontsize=14, bold=True)

label(ax, -1.80,  1.90, 'UP DEGs',     fontsize=12, color=UP_C,    bold=True)
label(ax,  1.80,  1.90, 'DOWN DEGs',   fontsize=12, color=DOWN_C,  bold=True)
label(ax,  0.0,  -2.05, 'BLAST Hits',  fontsize=12, color=BLAST_C, bold=True)

ax.set_title('Venn 3: UP / DOWN DEGs\nvs BLAST-annotated Genes',
             fontsize=13, fontweight='bold', pad=14)

# ── Final layout ───────────────────────────────────────────────────────────────
plt.suptitle('Venn Diagrams — Differentially Expressed Genes (489 DEGs)',
             fontsize=15, fontweight='bold', y=1.02)
plt.tight_layout(pad=2.5)
plt.savefig('Venn_DEGs.png', dpi=180, bbox_inches='tight', facecolor='white')
plt.close()

print("✓ Venn_DEGs.png saved successfully!")
print()
print("─── Set sizes ───────────────────────────────")
print(f"  UP-regulated         : {len(up_set)}")
print(f"  DOWN-regulated       : {len(down_set)}")
print(f"  UP ∩ DOWN            : {len(up_set & down_set)}  (mutually exclusive)")
print()
print(f"  GO-BP annotated      : {len(go_bp)}")
print(f"  GO-MF annotated      : {len(go_mf)}")
print(f"  GO-CC annotated      : {len(go_cc)}")
print(f"  BP ∩ MF ∩ CC         : {len(go_bp & go_mf & go_cc)}")
print()
print(f"  BLAST hits           : {len(blast_set)}")
print(f"  UP with BLAST        : {len(up_blast)}")
print(f"  DOWN with BLAST      : {len(dn_blast)}")
