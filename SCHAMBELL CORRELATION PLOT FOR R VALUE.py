"""
Schemaball Correlation Plot — WHITE BACKGROUND, VIVID VISIBLE CHORDS
=====================================================================
- White background
- Bold, thick, clearly visible chord lines
- Strong saturated colors: deep red (positive) / deep blue (negative)
- Opacity & width both scale with |r|
- Colorbar on white background

Usage : python schemaball_v2.py
Output: Schemaball_White.png
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap, Normalize
import matplotlib.cm as cm
import warnings
warnings.filterwarnings('ignore')

# ── Load data & compute correlation ───────────────────────────────────────────
df = pd.read_csv('UPDOWN.csv')

cols = [
    '144 Expression', '187 Expression', 'FoldChange', 'log2FoldChange',
    'P Value', 'Adjusted P Value', 'Unigene_length(Query)', 'Protein_length',
    '%Identity', 'Alignment_length', 'Query_coverage', 'Subject_coverage',
    'No.Of_mismatches', 'No.Of_gaps', 'E-value'
]

short_names = [
    'Expr-144', 'Expr-187', 'FC', 'log2FC',
    'P-val', 'Adj.P', 'Unigene_len', 'Prot_len',
    '%Identity', 'Align_len', 'Q_cov', 'S_cov',
    'Mismatches', 'Gaps', 'E-value'
]

corr = df[cols].corr(method='pearson')
n = len(cols)

# ── Vivid segment colours — each variable gets a unique bold colour ────────────
seg_colors = [
    '#E53935',   # deep red       Expr-144
    '#8E24AA',   # deep purple    Expr-187
    '#1E88E5',   # strong blue    FC
    '#00897B',   # teal           log2FC
    '#F4511E',   # deep orange    P-val
    '#3949AB',   # indigo         Adj.P
    '#00ACC1',   # cyan           Unigene_len
    '#7CB342',   # green          Prot_len
    '#FB8C00',   # amber          %Identity
    '#D81B60',   # pink           Align_len
    '#6D4C41',   # brown          Q_cov
    '#039BE5',   # light blue     S_cov
    '#43A047',   # medium green   Mismatches
    '#F06292',   # light pink     Gaps
    '#FF7043',   # deep orange2   E-value
]

# ── Diverging colormap: dark blue → light blue → light red → dark red
# (avoids white in the middle so lines stay visible on white bg)
cmap = LinearSegmentedColormap.from_list(
    'schemaball_white',
    ['#1A237E',   # very dark blue   (r = -1)
     '#1565C0',   # dark blue
     '#1E88E5',   # medium blue
     '#90CAF9',   # light blue       (r ~ -0.2)
     '#FFCC80',   # light orange     (r ~ +0.2)
     '#EF5350',   # medium red
     '#B71C1C'],  # very dark red    (r = +1)
    N=512
)
norm = Normalize(vmin=-1, vmax=1)

# ── Figure ─────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(14, 14))
fig.patch.set_facecolor('white')
ax.set_facecolor('white')
ax.set_xlim(-1.65, 1.65)
ax.set_ylim(-1.65, 1.65)
ax.set_aspect('equal')
ax.axis('off')

# ── Parameters ─────────────────────────────────────────────────────────────────
R_outer  = 1.30
R_inner  = 1.08
gap_deg  = 2.5
total_arc= 360 - n * gap_deg
arc_per  = total_arc / n

def angle_for_seg(i):
    s = 90.0 + i * (arc_per + gap_deg)
    e = s + arc_per
    m = np.radians((s + e) / 2)
    return s, e, m

# ── Draw outer ring segments ───────────────────────────────────────────────────
for i in range(n):
    s_deg, e_deg, mid_rad = angle_for_seg(i)
    color = seg_colors[i]

    theta = np.linspace(np.radians(s_deg), np.radians(e_deg), 100)
    x_out = R_outer * np.cos(theta)
    y_out = R_outer * np.sin(theta)
    x_in  = R_inner * np.cos(theta[::-1])
    y_in  = R_inner * np.sin(theta[::-1])

    xs = np.concatenate([x_out, x_in])
    ys = np.concatenate([y_out, y_in])
    ax.fill(xs, ys, color=color, alpha=1.0, zorder=5)

    # Dark border for segment definition
    ax.plot(x_out, y_out, color='#333333', linewidth=0.6, alpha=0.7, zorder=6)
    ax.plot(x_in,  y_in,  color='#333333', linewidth=0.3, alpha=0.5, zorder=6)

# ── Draw chord curves (Quadratic Bezier) ──────────────────────────────────────
threshold = 0.04   # skip near-zero correlations

pairs = []
for i in range(n):
    for j in range(i+1, n):
        r = corr.iloc[i, j]
        if abs(r) >= threshold:
            pairs.append((i, j, r))

# Draw weakest first, strongest on top
pairs.sort(key=lambda x: abs(x[2]))

for (i, j, r) in pairs:
    _, _, mid_i = angle_for_seg(i)
    _, _, mid_j = angle_for_seg(j)

    x1 = R_inner * np.cos(mid_i)
    y1 = R_inner * np.sin(mid_i)
    x2 = R_inner * np.cos(mid_j)
    y2 = R_inner * np.sin(mid_j)

    # Control point pulled toward centre
    ctrl_scale = 0.30
    cx = ctrl_scale * (x1 + x2)
    cy = ctrl_scale * (y1 + y2)

    t  = np.linspace(0, 1, 200)
    bx = (1-t)**2 * x1 + 2*(1-t)*t * cx + t**2 * x2
    by = (1-t)**2 * y1 + 2*(1-t)*t * cy + t**2 * y2

    rgba = cmap(norm(r))

    # KEY CHANGE: much higher alpha floor + bigger linewidth so always visible
    alpha = 0.45 + 0.50 * abs(r)          # min 0.45, max 0.95
    lw    = 0.8  + 4.5  * abs(r) ** 1.4   # min 0.8, max ~5.3 for r=1

    ax.plot(bx, by, color=rgba, alpha=alpha, linewidth=lw,
            solid_capstyle='round', zorder=3)

# ── Variable labels ────────────────────────────────────────────────────────────
R_label = 1.43
for i in range(n):
    _, _, mid_rad = angle_for_seg(i)
    deg = np.degrees(mid_rad)

    x = R_label * np.cos(mid_rad)
    y = R_label * np.sin(mid_rad)

    if 90 < deg <= 270:
        rot = deg + 180
        ha  = 'right'
    else:
        rot = deg
        ha  = 'left'

    ax.text(x, y, short_names[i],
            ha=ha, va='center',
            rotation=rot, rotation_mode='anchor',
            fontsize=11, fontweight='bold',
            color=seg_colors[i], zorder=8,
            bbox=dict(boxstyle='round,pad=0.1', fc='white', ec='none', alpha=0.6))

# ── Colourbar ──────────────────────────────────────────────────────────────────
sm = cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
cbar_ax = fig.add_axes([0.83, 0.20, 0.022, 0.28])
cbar = fig.colorbar(sm, cax=cbar_ax)
cbar.set_label('Pearson  r', color='black', fontsize=11, labelpad=10)
cbar.ax.yaxis.set_tick_params(color='black', labelsize=9)
plt.setp(cbar.ax.yaxis.get_ticklabels(), color='black')
cbar.outline.set_edgecolor('#555555')
cbar.outline.set_linewidth(0.8)

# ── Title ──────────────────────────────────────────────────────────────────────
ax.text(0, 1.60, 'Schemaball Correlation Plot',
        ha='center', va='center', fontsize=18, fontweight='bold',
        color='#1A1A2E', zorder=10)
ax.text(0, 1.50, 'Transcriptomic DEG Parameters  (Pearson r)',
        ha='center', va='center', fontsize=12, color='#444444', zorder=10)

# ── Legend ────────────────────────────────────────────────────────────────────
pos_patch = mpatches.Patch(color='#B71C1C', label='Positive correlation')
neg_patch = mpatches.Patch(color='#1A237E', label='Negative correlation')
ax.legend(handles=[pos_patch, neg_patch],
          loc='lower center', bbox_to_anchor=(0.0, -0.06),
          frameon=True, facecolor='white', edgecolor='#CCCCCC',
          fontsize=11, ncol=2,
          labelcolor='black')

plt.tight_layout()
plt.savefig('Schemaball_White.png', dpi=200,
            bbox_inches='tight', facecolor='white')
plt.close()

print("✓ Schemaball_White.png saved!")
print(f"  Variables : {n}")
print(f"  Chords    : {len(pairs)}  (|r| ≥ {threshold})")
print(f"  Background: WHITE")
print(f"  Chord colors: dark blue (negative) → dark red (positive)")
