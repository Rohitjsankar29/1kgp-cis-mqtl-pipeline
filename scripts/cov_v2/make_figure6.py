#!/usr/bin/env python3
"""
Figure 6 | Multi-evidence convergence (4-panel composite)

Run on Gadi, where the data lives. Produces Figure6_composite.png at 300 dpi.

Panels
  (a) Convergence-score distribution           <- computed from the scored table (col 37)
  (b) UpSet of the five evidence axes          <- computed from ev flags (cols 32-36)
  (c) ENCODE cCRE class composition            <- needs a cCRE-class source (see CONFIG)
  (d) Imprinted loci as positive control       <- asm_frac_fisher (col 19) split by imprinted set

Only a and b are fully self-contained. c and d need the two CONFIG hooks below.

  micromamba activate /scratch/cy94/rs4477/micromamba/envs/mqtl_tensorqtl   # (has matplotlib)
  python3 make_figure6.py
"""

import gzip, math, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec

# ============================ CONFIG ============================
PRIOR_GZ = "/g/data/cy94/rs4477/downstream/genome_452/genome.prioritised_final_asm_rescored.txt.gz"

# 1-based column numbers in the header (confirmed layout)
C_NEAREST    = 10    # nearest_gene
C_CCRE       = 13    # ccre (binary)
C_ASM_FISHER = 19    # asm_frac_fisher   (y-value for panel d)
C_ASM_TESTED = 22    # asm_tested (binary)
C_SCORE      = 37    # score_convergence
EV_COLS      = (32, 33, 34, 35, 36)          # ev_pip, ev_asm, ev_ccre, ev_sv, ev_effect
EV_NAMES     = ["PIP", "ASM", "cCRE", "SV", "Effect"]

# ---- (c) cCRE class source -------------------------------------------------
# Option 1: a TSV with two columns  variant_id <TAB> class   (no header needed).
#           The script counts classes among cCRE-overlapping variants (ccre==1).
CCRE_CLASS_FILE = None    # e.g. "/g/data/cy94/rs4477/downstream/genome_452/ccre_class_by_variant.tsv"
# Option 2 (fallback, used only if CCRE_CLASS_FILE is None): paste counts here.
#   Get them, e.g.:  <your ccre-class annotation> | awk '...' | sort | uniq -c
CCRE_CLASS_COUNTS = {
    # "dELS": 632000, "pELS": 125000, "CA": 118000, "CA-CTCF": 49000,
    # "CA-H3K4me3": 46000, "TF": 42000, "PLS": 23000, "CA-TF": 11000,
}

# ---- (d) imprinted control -------------------------------------------------
# Loci are flagged imprinted when nearest_gene is in this set. Adjust to match
# your imprinted-region definition (your mock-up had n=525). Alternatively set
# IMPRINTED_CPG_FILE to a file of imprinted best_cpg ids (one per line) and it
# takes precedence.
IMPRINTED_CPG_FILE = None
IMPRINTED_GENES = {
    "PEG13","SNRPN","SNURF","UBE3A","MAGEL2","NDN","NPAP1","PWAR1",      # PWS/AS
    "KCNQ1","KCNQ1OT1","CDKN1C","PHLDA2","SLC22A18",                      # KCNQ1
    "IGF2","IGF2-AS","H19","INS",                                         # IGF2-H19
    "GNAS","GNAS-AS1","NESP","NESPAS",                                    # GNAS
    "ZNF331","PEG3","ZIM2","MEST","MESTIT1","PLAGL1","HYMAI",
    "GRB10","DLK1","MEG3","MEG8","RTL1","DIO3","NNAT","BLCAP",
    "L3MBTL1","NAP1L5","INPP5F","ZDBF2","GPR1","NAA60","ERLIN2",
}

OUT = "Figure6_composite.png"
MIN_COMBO = 1000     # UpSet shows combinations with >= this many variants
# ================================================================


def read_table():
    """Single pass over the scored table."""
    score_counts = {s: 0 for s in range(6)}
    combo_counts = {}
    asm_imp, asm_seq = [], []

    imp_cpgs = None
    if IMPRINTED_CPG_FILE:
        with open(IMPRINTED_CPG_FILE) as fh:
            imp_cpgs = {ln.strip() for ln in fh if ln.strip()}

    # optional cCRE class lookup
    cls_lookup = None
    if CCRE_CLASS_FILE:
        cls_lookup = {}
        op = gzip.open if CCRE_CLASS_FILE.endswith(".gz") else open
        with op(CCRE_CLASS_FILE, "rt") as fh:
            for ln in fh:
                p = ln.rstrip("\n").split("\t")
                if len(p) >= 2:
                    cls_lookup[p[0]] = p[1]
    ccre_class_counts = {}

    op = gzip.open if PRIOR_GZ.endswith(".gz") else open
    with op(PRIOR_GZ, "rt") as fh:
        fh.readline()  # header
        iNear, iCcre, iFish, iTest, iScore = (C_NEAREST-1, C_CCRE-1, C_ASM_FISHER-1,
                                              C_ASM_TESTED-1, C_SCORE-1)
        iev = [c-1 for c in EV_COLS]
        for line in fh:
            c = line.rstrip("\n").split("\t")
            vid, best_cpg = c[0], c[1]
            s = int(float(c[iScore]))
            if 0 <= s <= 5:
                score_counts[s] += 1
            ev = tuple(int(float(c[k])) for k in iev)
            combo_counts[ev] = combo_counts.get(ev, 0) + 1

            # (d) asm by imprinted/seq, among tested regions
            if c[iTest] not in ("", "0", "0.0") and c[iFish] not in ("", "NA"):
                try:
                    frac = float(c[iFish])
                except ValueError:
                    frac = None
                if frac is not None:
                    if imp_cpgs is not None:
                        is_imp = best_cpg in imp_cpgs
                    else:
                        is_imp = c[iNear] in IMPRINTED_GENES
                    (asm_imp if is_imp else asm_seq).append(frac)

            # (c) cCRE class among cCRE-overlapping variants
            if cls_lookup is not None and c[iCcre] in ("1", "1.0"):
                cl = cls_lookup.get(vid)
                if cl:
                    ccre_class_counts[cl] = ccre_class_counts.get(cl, 0) + 1

    if not cls_lookup:
        ccre_class_counts = dict(CCRE_CLASS_COUNTS)
    return score_counts, combo_counts, asm_imp, asm_seq, ccre_class_counts


def marginals(combo_counts):
    m = [0]*5
    for ev, n in combo_counts.items():
        for i in range(5):
            if ev[i]:
                m[i] += n
    return m


# --------------------------- panels ---------------------------
def panel_a(ax, score_counts):
    vals = [score_counts[s] for s in range(6)]
    cols = ["#9aa0a6", "#a9cde6", "#6ba7d0", "#357fb0", "#e8952e", "#c0392b"]
    x = np.arange(6)
    ax.bar(x, vals, width=0.72, color=cols, edgecolor="white", lw=0.8, zorder=3)
    ax.set_yscale("log"); ax.set_ylim(1e3, 3e6)
    for xi, v in zip(x, vals):
        ax.text(xi, v*1.25, f"{v:,}", ha="center", va="bottom", fontsize=7.5, color="#223")
    ax.set_xticks(x)
    ax.set_xticklabels(["0 axes","1 axis","2 axes","3 axes","4 axes","5 axes\n(all)"], fontsize=8.5)
    ax.set_xlabel("Number of evidence axes satisfied", fontsize=9.5)
    ax.set_ylabel("Variants (log scale)", fontsize=9.5)
    ax.set_title("a   Multi-evidence convergence", fontsize=11, fontweight="bold", loc="left")
    ax.text(0.98, 0.97, "Axes: PIP, ASM, cCRE, SV, Effect size", transform=ax.transAxes,
            ha="right", va="top", fontsize=7.5, color="#6b7177")
    for sp in ["top","right"]: ax.spines[sp].set_visible(False)
    ax.tick_params(labelsize=8); ax.grid(axis="y", lw=.4, color="#eceef0", zorder=0)


def panel_b(fig, subspec, combo_counts):
    m = marginals(combo_counts)
    order = sorted(range(5), key=lambda i: -m[i])          # rows by set size desc
    combos = sorted([kv for kv in combo_counts.items() if sum(kv[0])>0 and kv[1]>=MIN_COMBO],
                    key=lambda kv: -kv[1])
    sizes = [n for _, n in combos]; mem = [k for k, _ in combos]; n = len(combos)

    gs = GridSpecFromSubplotSpec(2, 2, subplot_spec=subspec,
                                 width_ratios=[1.1, 6.5], height_ratios=[3.0, 2.2],
                                 wspace=0.05, hspace=0.05)
    NAVY="#13537d"; OFF="#d9dde1"
    axbar = fig.add_subplot(gs[0, 1]); x = np.arange(n)
    axbar.bar(x, sizes, width=0.6, color=NAVY, zorder=3)
    axbar.set_yscale("log"); axbar.set_ylim(3e3, 2e6)
    for xi, s in zip(x, sizes):
        axbar.text(xi, s*1.15, f"{s:,}", ha="center", va="bottom", fontsize=5.6, rotation=90, color="#223")
    axbar.set_ylabel("Intersection", fontsize=8.5); axbar.set_xticks([]); axbar.set_xlim(-0.7, n-0.3)
    for sp in ["top","right"]: axbar.spines[sp].set_visible(False)
    axbar.tick_params(labelsize=7); axbar.grid(axis="y", lw=.35, color="#eceef0", zorder=0)
    axbar.set_title("b   Convergence across five axes", fontsize=11, fontweight="bold", loc="left")

    axm = fig.add_subplot(gs[1, 0]); yy = np.arange(5)[::-1]
    names_o = [EV_NAMES[i] for i in order]; vals_o = [m[i] for i in order]
    axm.barh(yy, vals_o, color="#6ba7d0", height=0.55, zorder=3)
    axm.set_xscale("log"); axm.set_xlim(3e6, 2e4); axm.invert_xaxis()
    for y, v in zip(yy, vals_o):
        axm.text(v*1.3, y, f"{v/1000:.0f}k", va="center", ha="right", fontsize=6.2, color="#334")
    axm.set_yticks(yy); axm.set_yticklabels(names_o, fontsize=8); axm.set_xticks([])
    for sp in ["top","right","bottom"]: axm.spines[sp].set_visible(False)
    axm.set_xlabel("Set size", fontsize=7.5)

    axd = fig.add_subplot(gs[1, 1]);
    for r, i in enumerate(order):
        yr = 4 - r
        if r % 2 == 0: axd.axhspan(yr-0.5, yr+0.5, color="#f5f6f7", zorder=0)
        for xi, k in zip(x, mem):
            axd.scatter(xi, yr, s=34, color=NAVY if k[i] else OFF, zorder=3)
    for xi, k in zip(x, mem):
        rows_on = [4-r for r, i in enumerate(order) if k[i]]
        if len(rows_on) > 1:
            axd.plot([xi, xi], [min(rows_on), max(rows_on)], color=NAVY, lw=1.4, zorder=2)
    axd.set_xlim(-0.7, n-0.3); axd.set_ylim(-0.6, 4.6); axd.set_xticks([]); axd.set_yticks([])
    for sp in axd.spines.values(): sp.set_visible(False)


def panel_c(ax, counts):
    if not counts:
        ax.text(0.5, 0.5, "cCRE class data not configured\n(see CONFIG: CCRE_CLASS_FILE / COUNTS)",
                ha="center", va="center", fontsize=8, color="#b0453a", transform=ax.transAxes)
        ax.set_title("c   Regulatory context: cCRE class", fontsize=11, fontweight="bold", loc="left")
        ax.axis("off"); return
    items = sorted(counts.items(), key=lambda kv: -kv[1])
    labels = [k for k, _ in items]; vals = [v for _, v in items]
    y = np.arange(len(labels))[::-1]
    ax.barh(y, vals, color="#3a9a5c", height=0.62, zorder=3)
    ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=8.5)
    ax.set_xlabel("Prioritised variants", fontsize=9.5)
    ax.set_title("c   Regulatory context: cCRE class", fontsize=11, fontweight="bold", loc="left")
    for sp in ["top","right"]: ax.spines[sp].set_visible(False)
    ax.tick_params(labelsize=8); ax.grid(axis="x", lw=.4, color="#eceef0", zorder=0)


def panel_d(ax, asm_seq, asm_imp):
    bins = np.linspace(0, 1, 41)
    if asm_seq:
        ax.hist(asm_seq, bins=bins, density=True, color="#6ba7d0", alpha=0.6, label="sequence-driven")
    if asm_imp:
        ax.hist(asm_imp, bins=bins, density=True, color="#c0556a", alpha=0.6,
                label=f"imprinted control (n={len(asm_imp)})")
    ax.set_xlabel("ASM fraction (Fisher)", fontsize=9.5); ax.set_ylabel("density", fontsize=9.5)
    ax.set_title("d   Imprinted loci as positive control", fontsize=11, fontweight="bold", loc="left")
    ax.legend(fontsize=8, frameon=False)
    for sp in ["top","right"]: ax.spines[sp].set_visible(False)
    ax.tick_params(labelsize=8)


def main():
    sc, cc, asm_imp, asm_seq, ccre_cls = read_table()
    tot = sum(sc.values())
    print("score distribution:", sc, " total:", f"{tot:,}")
    print("per-axis marginals:", dict(zip(EV_NAMES, marginals(cc))))
    print("ASM tested  seq:", len(asm_seq), " imprinted:", len(asm_imp))
    print("cCRE classes:", ccre_cls if ccre_cls else "(not configured)")

    fig = plt.figure(figsize=(13.5, 9.0), dpi=300)
    gs = GridSpec(2, 2, figure=fig, hspace=0.32, wspace=0.22)
    panel_a(fig.add_subplot(gs[0, 0]), sc)
    panel_b(fig, gs[0, 1], cc)
    panel_c(fig.add_subplot(gs[1, 0]), ccre_cls)
    panel_d(fig.add_subplot(gs[1, 1]), asm_seq, asm_imp)
    fig.savefig(OUT, bbox_inches="tight", facecolor="white")
    print("wrote", OUT)


if __name__ == "__main__":
    main()
