"""
All Possible Venn Diagrams — Transcriptomic DEG Data
=====================================================
Plots 9 Venn diagrams across 3 figures:
  Figure 1 (2x3): Six 2-way Venn diagrams
  Figure 2 (2x3): Six 3-way Venn diagrams
  Figure 3 (1x1): One 4-way Venn diagram (ellipse style)

Usage : python all_venn_diagrams.py
Output: Venn_2way.png, Venn_3way.png, Venn_4way.png
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyBboxPatch, Ellipse
from matplotlib.colors import to_rgba
import warnings
warnings.filterwarnings('ignore')

# ── Load & build sets ──────────────────────────────────────────────────────────
df = pd.read_csv('UPDOWN.csv')

up    = set(df.loc[df['Regulation']=='UP',    'Unigene_ID'])
down  = set(df.loc[df['Regulation']=='DOWN',  'Unigene_ID'])
blast = set(df.loc[df['%Identity'].notna(),   'Unigene_ID'])
go_bp = set(df.loc[df['GO_BP'].notna(),        'Unigene_ID'])
go_mf = set(df.loc[df['GO_MF'].notna(),        'Unigene_ID'])
go_cc = set(df.loc[df['GO_CC'].notna(),        'Unigene_ID'])
go_any= go_bp | go_mf | go_cc

# ═══════════════════════════════════════════════════════════════════════════════
#  HELPERS
# ═══════════════════════════════════════════════════════════════════════════════
def add_circle(ax, cx, cy, r, color, alpha=0.40, lw=2.2):
    p = Circle((cx, cy), r, facecolor=to_rgba(color, alpha),
                edgecolor=color, linewidth=lw, zorder=2)
    ax.add_patch(p)

def add_rect(ax, x, y, w, h, color, alpha=0.40, lw=2.2):
    p = FancyBboxPatch((x, y), w, h,
                       boxstyle="round,pad=0.25",
                       facecolor=to_rgba(color, alpha),
                       edgecolor=color, linewidth=lw, zorder=2)
    ax.add_patch(p)

def add_ellipse(ax, cx, cy, rx, ry, angle, color, alpha=0.35, lw=2.2):
    p = Ellipse((cx, cy), 2*rx, 2*ry, angle=angle,
                facecolor=to_rgba(color, alpha),
                edgecolor=color, linewidth=lw, zorder=2)
    ax.add_patch(p)

def num(ax, x, y, val, fs=15, bold=True):
    ax.text(x, y, str(val), ha='center', va='center',
            fontsize=fs, fontweight='bold' if bold else 'normal',
            color='#111111', zorder=6)

def lbl(ax, x, y, text, color, fs=11):
    ax.text(x, y, text, ha='center', va='center',
            fontsize=fs, fontweight='bold', color=color, zorder=7)

def style_ax(ax, title):
    ax.set_aspect('equal'); ax.axis('off')
    ax.set_title(title, fontsize=12, fontweight='bold', pad=8)

# ═══════════════════════════════════════════════════════════════════════════════
#  FIGURE 1 — Six 2-way Venn Diagrams
# ═══════════════════════════════════════════════════════════════════════════════
PALETTES_2 = [
    ('#FF6B6B','#4ECDC4'),   # coral / teal
    ('#A855F7','#F97316'),   # purple / orange
    ('#22D3EE','#F43F5E'),   # cyan / rose
    ('#84CC16','#F59E0B'),   # lime / amber
    ('#6366F1','#EC4899'),   # indigo / pink
    ('#10B981','#EF4444'),   # emerald / red
]

two_way_data = [
    (up,    down,   'UP-DEGs',  'DOWN-DEGs', 'UP vs DOWN DEGs',            PALETTES_2[0]),
    (up,    blast,  'UP-DEGs',  'BLAST hits','UP-DEGs vs BLAST Annotated', PALETTES_2[1]),
    (down,  blast,  'DOWN-DEGs','BLAST hits','DOWN-DEGs vs BLAST Annotated',PALETTES_2[2]),
    (up,    go_any, 'UP-DEGs',  'GO Annot.', 'UP-DEGs vs GO Annotated',    PALETTES_2[3]),
    (down,  go_any, 'DOWN-DEGs','GO Annot.', 'DOWN-DEGs vs GO Annotated',  PALETTES_2[4]),
    (go_bp, go_mf,  'GO-BP',    'GO-MF',     'GO-BP vs GO-MF',             PALETTES_2[5]),
]

fig1, axes1 = plt.subplots(2, 3, figsize=(16, 11))
fig1.patch.set_facecolor('white')
fig1.suptitle('2-Way Venn Diagrams — DEG Transcriptomic Data',
              fontsize=16, fontweight='bold', y=1.01)

for ax, (A, B, nA, nB, title, (cA, cB)) in zip(axes1.flat, two_way_data):
    only_A = len(A - B)
    only_B = len(B - A)
    overlap= len(A & B)

    ax.set_xlim(0, 10); ax.set_ylim(0, 8)
    style_ax(ax, title)

    # Two overlapping circles
    add_circle(ax, 3.8, 4.0, 2.5, cA)
    add_circle(ax, 6.2, 4.0, 2.5, cB)

    # Numbers
    num(ax, 2.4, 4.0, only_A)
    num(ax, 7.6, 4.0, only_B)
    num(ax, 5.0, 4.0, overlap)

    # Labels
    lbl(ax, 2.4, 7.0, nA, cA, fs=11)
    lbl(ax, 7.6, 7.0, nB, cB, fs=11)

    # Total annotations
    ax.text(0.3, 0.4, f'n(A)={len(A)}', fontsize=8,
            color=cA, transform=ax.transAxes)
    ax.text(0.7, 0.4, f'n(B)={len(B)}', fontsize=8,
            color=cB, transform=ax.transAxes, ha='right')

plt.tight_layout(pad=2.5)
plt.savefig('Venn_2way.png', dpi=180, bbox_inches='tight', facecolor='white')
plt.close()
print("✓ Venn_2way.png  — 6 two-way Venn diagrams")

# ═══════════════════════════════════════════════════════════════════════════════
#  FIGURE 2 — Six 3-way Venn Diagrams (circle + 2 rects style like reference)
# ═══════════════════════════════════════════════════════════════════════════════
PALETTES_3 = [
    ('#FF6B6B','#4ECDC4','#FFE66D'),
    ('#A855F7','#F97316','#22D3EE'),
    ('#F43F5E','#10B981','#6366F1'),
    ('#EF4444','#84CC16','#F59E0B'),
    ('#EC4899','#06B6D4','#8B5CF6'),
    ('#FF6B35','#004E89','#1A936F'),
]

three_way_data = [
    (up,   down,  blast,  'UP-DEGs', 'DOWN-DEGs','BLAST hits', 'UP / DOWN / BLAST',          PALETTES_3[0]),
    (up,   down,  go_any, 'UP-DEGs', 'DOWN-DEGs','GO Annot.',  'UP / DOWN / GO Annotated',   PALETTES_3[1]),
    (up,   blast, go_any, 'UP-DEGs', 'BLAST hits','GO Annot.', 'UP / BLAST / GO Annotated',  PALETTES_3[2]),
    (down, blast, go_any, 'DOWN-DEGs','BLAST hits','GO Annot.','DOWN / BLAST / GO Annotated',PALETTES_3[3]),
    (go_bp,go_mf, go_cc,  'GO-BP',   'GO-MF',    'GO-CC',      'GO-BP / GO-MF / GO-CC',      PALETTES_3[4]),
    (up,   go_bp, go_mf,  'UP-DEGs', 'GO-BP',    'GO-MF',      'UP / GO-BP / GO-MF',         PALETTES_3[5]),
]

def draw_3way_reference(ax, A, B, C, nA, nB, nC, title, cA, cB, cC):
    """
    Reference style: circle (A) + top-right rounded rect (B) + bottom rounded rect (C)
    """
    only_A   = len(A - B - C)
    only_B   = len(B - A - C)
    only_C   = len(C - A - B)
    AB_only  = len((A & B) - C)
    AC_only  = len((A & C) - B)
    BC_only  = len((B & C) - A)
    all3     = len(A & B & C)

    ax.set_xlim(0, 10); ax.set_ylim(0, 10)
    style_ax(ax, title)

    # Rect B — top right
    add_rect(ax, 4.2, 4.0, 5.4, 5.6, cB, alpha=0.38)
    # Rect C — bottom spanning
    add_rect(ax, 0.4, 0.4, 9.2, 5.6, cC, alpha=0.38)
    # Circle A — centre
    add_circle(ax, 4.2, 5.5, 3.0, cA, alpha=0.42)

    # Region numbers
    num(ax, 2.2, 6.5, only_A)    # only circle
    num(ax, 8.4, 8.5, only_B)    # only rect-B
    num(ax, 1.5, 1.8, only_C)    # only rect-C
    num(ax, 6.2, 7.5, AB_only)   # circle ∩ rect-B
    num(ax, 3.5, 3.5, AC_only)   # circle ∩ rect-C
    num(ax, 7.8, 2.8, BC_only)   # rect-B ∩ rect-C
    num(ax, 5.5, 5.2, all3)      # all three

    # Labels
    lbl(ax, 2.0, 9.5, nA, cA, fs=10)
    lbl(ax, 8.5, 9.8, nB, cB, fs=10)
    lbl(ax, 1.8, 0.05,nC, cC, fs=10)

fig2, axes2 = plt.subplots(2, 3, figsize=(18, 13))
fig2.patch.set_facecolor('white')
fig2.suptitle('3-Way Venn Diagrams — DEG Transcriptomic Data',
              fontsize=16, fontweight='bold', y=1.01)

for ax, (A, B, C, nA, nB, nC, title, (cA, cB, cC)) in zip(axes2.flat, three_way_data):
    draw_3way_reference(ax, A, B, C, nA, nB, nC, title, cA, cB, cC)

plt.tight_layout(pad=3)
plt.savefig('Venn_3way.png', dpi=180, bbox_inches='tight', facecolor='white')
plt.close()
print("✓ Venn_3way.png  — 6 three-way Venn diagrams")

# ═══════════════════════════════════════════════════════════════════════════════
#  FIGURE 3 — 4-way Venn Diagram (ellipse style: UP / DOWN / BLAST / GO)
# ═══════════════════════════════════════════════════════════════════════════════
A, B, C, D = up, down, blast, go_any

only_A    = len(A - B - C - D)
only_B    = len(B - A - C - D)
only_C    = len(C - A - B - D)
only_D    = len(D - A - B - C)
AB        = len((A & B) - C - D)
AC        = len((A & C) - B - D)
AD        = len((A & D) - B - C)
BC        = len((B & C) - A - D)
BD        = len((B & D) - A - C)
CD        = len((C & D) - A - B)
ABC       = len((A & B & C) - D)
ABD       = len((A & B & D) - C)
ACD       = len((A & C & D) - B)
BCD       = len((B & C & D) - A)
ABCD      = len(A & B & C & D)

cA4 = '#FF6B6B'
cB4 = '#4ECDC4'
cC4 = '#FFE66D'
cD4 = '#A855F7'

fig3, ax = plt.subplots(1, 1, figsize=(11, 10))
fig3.patch.set_facecolor('white')
ax.set_xlim(0, 14); ax.set_ylim(0, 14)
ax.set_aspect('equal'); ax.axis('off')

# 4 overlapping ellipses at different angles
add_ellipse(ax, 5.5, 7.5, 3.8, 2.2, -40, cA4, alpha=0.28)
add_ellipse(ax, 8.5, 7.5, 3.8, 2.2,  40, cB4, alpha=0.28)
add_ellipse(ax, 5.5, 6.5, 3.8, 2.2,  40, cC4, alpha=0.28)
add_ellipse(ax, 8.5, 6.5, 3.8, 2.2, -40, cD4, alpha=0.28)

# Region labels (approximate visual positions for each of the 15 regions)
fs4 = 13
num(ax, 2.4, 10.5, only_A,  fs4)    # only UP
num(ax, 11.4,10.5, only_B,  fs4)    # only DOWN
num(ax, 2.4,  3.0, only_C,  fs4)    # only BLAST
num(ax, 11.4, 3.0, only_D,  fs4)    # only GO
num(ax, 5.2, 10.0, AB,      fs4)    # UP∩DOWN
num(ax, 3.5,  7.0, AC,      fs4)    # UP∩BLAST
num(ax, 4.5,  8.8, AD,      fs4)    # UP∩GO
num(ax, 10.5, 7.0, BC,      fs4)    # DOWN∩BLAST
num(ax, 9.5,  8.8, BD,      fs4)    # DOWN∩GO
num(ax, 7.0,  4.5, CD,      fs4)    # BLAST∩GO
num(ax, 5.8,  9.0, ABC,     fs4)    # UP∩DOWN∩BLAST
num(ax, 8.2,  9.0, ABD,     fs4)    # UP∩DOWN∩GO
num(ax, 5.2,  5.8, ACD,     fs4)    # UP∩BLAST∩GO
num(ax, 8.8,  5.8, BCD,     fs4)    # DOWN∩BLAST∩GO
num(ax, 7.0,  7.0, ABCD,   fs4+1)  # ALL 4

# Set labels
lbl(ax,  1.5, 12.0, 'UP-DEGs',   cA4, fs=12)
lbl(ax, 12.5, 12.0, 'DOWN-DEGs', cB4, fs=12)
lbl(ax,  1.5,  1.5, 'BLAST hits',cC4, fs=12)
lbl(ax, 12.5,  1.5, 'GO Annot.', cD4, fs=12)

ax.set_title('4-Way Venn Diagram\nUP-DEGs / DOWN-DEGs / BLAST Annotated / GO Annotated',
             fontsize=14, fontweight='bold', pad=14)
fig3.suptitle('489 Differentially Expressed Genes', fontsize=13,
              color='grey', y=0.02)

plt.tight_layout()
plt.savefig('Venn_4way.png', dpi=180, bbox_inches='tight', facecolor='white')
plt.close()
print("✓ Venn_4way.png  — 4-way Venn diagram (UP/DOWN/BLAST/GO)")

# ── Summary ────────────────────────────────────────────────────────────────────
print()
print("=" * 50)
print("ALL VENN DIAGRAMS GENERATED!")
print("=" * 50)
print("  Venn_2way.png  — 6 two-way  Venns")
print("  Venn_3way.png  — 6 three-way Venns")
print("  Venn_4way.png  — 1 four-way  Venn")
print()
print("Sets used:")
print(f"  UP-DEGs   : {len(up)}")
print(f"  DOWN-DEGs : {len(down)}")
print(f"  BLAST hits: {len(blast)}")
print(f"  GO-BP     : {len(go_bp)}")
print(f"  GO-MF     : {len(go_mf)}")
print(f"  GO-CC     : {len(go_cc)}")
print(f"  GO (any)  : {len(go_any)}")
