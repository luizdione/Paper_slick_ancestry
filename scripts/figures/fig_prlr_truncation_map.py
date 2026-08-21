#!/usr/bin/env python3
"""make_fig_prlr_truncation_map.py -- mapa funcional COMPARADO das tres truncagens do PRLR.

Leva os alelos do nivel do gene (Fig. dos alelos) para o nivel da PROTEINA e da
consequencia de sinalizacao. Nao ha DE de expressao aqui: e a arquitetura do receptor
e onde cada truncante corta a cauda citoplasmatica.

Ancoras (nao inventadas):
  - PRLR bovino, isoforma canonica UniProt **Q28172, 581 aa**:
      extracelular (ligacao a prolactina, dominios FN3)  25-237
      transmembrana                                      238-258
      cauda citoplasmatica                               259-581
  - Box1 (motivo rico em prolina que ancora JAK2): proximal a cauda, ~262-271
    (arquitetura conservada de receptor de citocina classe I; homologo a PRLR humano).
  - Tirosina distal de docking de STAT5: perto do C-terminal (homologo a Y580 do
    PRLR humano P16471); representada em ~560.
  - Posicoes das truncagens (config do locus + notacao ja usada no paper):
      SLICK1  c.1382del  frameshift  -> ~aa 461   (origem N'Dama/africana)
      SLICK3  c.1393C>A  p.Ser465*   -> aa 465    (origem a definir, n=2)
      SLICK2  c.1489C>T  p.Arg497*   -> aa 497    (origem taurina europeia)

Mensagem: as tres cortam DISTAL ao Box1 (JAK2 preservado) e removem o C-terminal
que recruta STAT5 -> convergencia funcional de alelos de origem geografica distinta;
o alelo europeu (SLICK2) conserva o coto citoplasmatico mais longo.
"""
from __future__ import annotations
import pathlib
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrowPatch

OUT = pathlib.Path(__file__).parent / "figuras"
OUT.mkdir(exist_ok=True)

L = 581                       # comprimento da isoforma canonica
EXTRA = (25, 237)             # extracelular (ligacao)
TM = (238, 258)              # transmembrana
CYTO = (259, 581)            # cauda citoplasmatica
BOX1 = (262, 271)            # ancoragem de JAK2
STAT5_Y = 560                # tirosina distal de docking de STAT5 (por homologia)

TRUNC = [  # nome, aa, c-notacao, p-notacao, origem, cor
    ("SLICK1", 461, "c.1382del", "frameshift", "N'Dama / African", "#1b7837"),
    ("SLICK3", 465, "c.1393C>A", "p.Ser465*", "origin undetermined (n=2)", "#8c6d31"),
    ("SLICK2", 497, "c.1489C>T", "p.Arg497*", "European taurine", "#2166ac"),
]


def draw_receptor(ax, y, height=0.5, label_domains=True, xmax=None):
    """Desenha o receptor linear (1..L) numa faixa y."""
    ax.add_patch(Rectangle((EXTRA[0], y), EXTRA[1] - EXTRA[0], height, facecolor="#dfe6ef",
                           edgecolor="black", lw=0.7, zorder=2))
    ax.add_patch(Rectangle((TM[0], y), TM[1] - TM[0], height, facecolor="#8a8f98",
                           edgecolor="black", lw=0.7, zorder=2))
    xend = xmax if xmax is not None else CYTO[1]
    ax.add_patch(Rectangle((CYTO[0], y), xend - CYTO[0], height, facecolor="#f4efe6",
                           edgecolor="black", lw=0.7, zorder=2))
    ax.add_patch(Rectangle((BOX1[0], y), BOX1[1] - BOX1[0], height, facecolor="#2ca25f",
                           edgecolor="black", lw=0.7, zorder=3))
    if label_domains:
        ax.text((EXTRA[0] + EXTRA[1]) / 2, y + height / 2, "extracellular\n(prolactin-binding)",
                ha="center", va="center", fontsize=8, zorder=4)
        ax.text((TM[0] + TM[1]) / 2, y + height + 0.14, "TM", ha="center", va="bottom", fontsize=7.5)
        ax.text(BOX1[1] + 4, y + height + 0.14, "Box1\n(JAK2)", ha="left", va="bottom",
                fontsize=7.5, color="#166b3a", fontweight="bold", linespacing=1.0)


def main():
    fig = plt.figure(figsize=(12.5, 7.6))
    gs = fig.add_gridspec(2, 1, height_ratios=[1.15, 1.0], hspace=0.42)
    ax = fig.add_subplot(gs[0])       # (a) receptor completo + cortes
    axb = fig.add_subplot(gs[1])      # (b) receptores retidos (convergencia)

    # ---------- (a) receptor completo com os 3 cortes ----------
    y0 = 0.0
    draw_receptor(ax, y0)
    # regiao removida (a partir do corte mais proximal) sombreada
    first_cut = min(t[1] for t in TRUNC)
    ax.add_patch(Rectangle((first_cut, y0 - 0.02), CYTO[1] - first_cut, 0.54, facecolor="#d62728",
                           alpha=0.10, edgecolor="none", zorder=1))
    ax.text((first_cut + CYTO[1]) / 2, y0 - 0.26, "C-terminus removed by all three alleles",
            ha="center", va="top", fontsize=8, color="#a11", style="italic")
    # tirosina distal de STAT5 (marcador + rotulo curto entre o bar e a fileira de rotulos)
    ax.plot(STAT5_Y, y0 + 0.5 + 0.03, marker="o", ms=8, mfc="#d62728", mec="black", zorder=5)
    ax.text(STAT5_Y, y0 + 0.70, "distal Tyr (STAT5)", ha="center", va="bottom", fontsize=7.5,
            color="#a11", fontweight="bold")
    # setas dos cortes (461 e 465 quase coincidem; 497 depois)
    for nome, aa, cnot, pnot, orig, cor in TRUNC:
        ax.add_patch(FancyArrowPatch((aa, y0 + 0.74), (aa, y0 + 0.52), arrowstyle="-|>",
                                     mutation_scale=13, color=cor, lw=1.8, zorder=6))
    # rotulos numa fileira acima, espalhados, com leader em L ate cada corte
    ORIG = {"SLICK1": "N'Dama / African", "SLICK3": "origin undet. (n=2)", "SLICK2": "European taurine"}
    anchors = {"SLICK1": 235, "SLICK3": 380, "SLICK2": 525}
    for nome, aa, cnot, pnot, orig, cor in TRUNC:
        ax_x = anchors[nome]
        ax.plot([ax_x, ax_x, aa], [2.02, 1.12, 0.80], color=cor, lw=0.9, zorder=5,
                solid_capstyle="round")
        ax.text(ax_x, 2.06, f"{nome}  {pnot}\n{cnot}\n{ORIG[nome]}", ha="center", va="bottom",
                fontsize=7.6, color=cor, fontweight="bold", linespacing=1.3,
                bbox=dict(boxstyle="round,pad=0.28", fc="white", ec=cor, lw=0.9))
    ax.set_xlim(0, L + 12)
    ax.set_ylim(-0.6, 3.05)
    ax.set_yticks([])
    ax.set_xlabel("PRLR residue (Bos taurus, 581 aa; UniProt Q28172)", fontsize=9)
    ax.set_title("(a) Where each slick allele truncates the prolactin receptor",
                 fontsize=11, fontweight="bold", loc="left")
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)

    # ---------- (b) receptores retidos por alelo (convergencia) ----------
    order = sorted(TRUNC, key=lambda t: t[1])
    for i, (nome, aa, cnot, pnot, orig, cor) in enumerate(order):
        y = (len(order) - 1 - i) * 1.0
        draw_receptor(axb, y, height=0.42, label_domains=False, xmax=aa)
        # marca o ponto de corte
        axb.plot(aa, y + 0.21, marker="|", ms=16, color=cor, mew=2.4, zorder=6)
        axb.text(aa + 6, y + 0.21, f"* {aa}", ha="left", va="center", fontsize=8,
                 color=cor, fontweight="bold")
        axb.text(-6, y + 0.21, f"{nome}", ha="right", va="center", fontsize=9,
                 color=cor, fontweight="bold")
        # residuos de cauda conservados
        axb.text(L + 20, y + 0.21, f"{aa - CYTO[0]} aa tail kept", ha="left", va="center",
                 fontsize=7.5, color="#555")
    # onde ficaria a tirosina de STAT5 (perdida por todos)
    axb.axvline(STAT5_Y, color="#d62728", ls="--", lw=1.0, zorder=1, ymin=0.05, ymax=0.95)
    axb.text(STAT5_Y, len(order) - 0.15, "STAT5 Tyr\n(lost in all)", ha="center", va="bottom",
             fontsize=7.5, color="#a11", linespacing=1.0)
    axb.set_xlim(-70, L + 130)
    axb.set_ylim(-0.35, len(order) + 0.15)
    axb.set_yticks([])
    axb.set_xlabel("PRLR residue", fontsize=9)
    axb.set_title("(b) Retained receptor per allele: Box1/JAK2 kept, STAT5 C-terminus lost by all "
                  "— European SLICK2 keeps the longest stub", fontsize=10.5, fontweight="bold", loc="left")
    for s in ("top", "right", "left"):
        axb.spines[s].set_visible(False)

    fig.suptitle("Three independent-origin truncations converge on the same receptor module",
                 fontsize=12.5, fontweight="bold", y=0.995)
    for ext in ("pdf", "png"):
        fig.savefig(OUT / f"Fig_PRLR_truncation_map.{ext}", dpi=300, bbox_inches="tight")
    print("->", OUT / "Fig_PRLR_truncation_map.pdf")


if __name__ == "__main__":
    main()
