#!/usr/bin/env python3
"""
fig_prlr_alleles.py — Figura 1: os alelos truncantes do PRLR.

Le genotipo real (via TSV exportado do VCF fasado incorporado) e a estrutura
de exons real (do GFF3), nao numero copiado a mao.

Tres paineis:
  (a) modelo do gene em BTA20 (exons do transcript ENSBTAT00000014437, GFF3
      Ensembl 113) com as posicoes de SLICK1/SLICK2/SLICK3_cand marcadas.
  (b) consequencia funcional dos sitios codificantes recuperados no
      faseamento (contagem fixa, ver Metodos: 24 sinonimos, 32 nao
      sinonimos, 1.360 intronicos, 38 fora das tres classes).
  (c) frequencia alelica de cada truncante por PAINEL geografico (nao por
      raca fina — o painel de 1.842 animais originais nao tem rotulo de raca
      carregado nesta sessao; ver nota no rodape da figura).
"""
from __future__ import annotations
import argparse
import os
import pathlib
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrowPatch

OUT = pathlib.Path(__file__).parent / "figuras"
OUT.mkdir(exist_ok=True)

DATA_DIR = pathlib.Path(os.environ.get("SLICK_DATA_DIR", "data"))

# --- (a) exons do transcript principal, direto do GFF3 (rank, inicio, fim) ---
EXONS = [
    (1, 39037656, 39037692), (2, 39077900, 39078016), (3, 39080773, 39080905),
    (4, 39082809, 39082978), (5, 39093142, 39093311), (6, 39094967, 39095108),
    (7, 39096738, 39096837), (8, 39097533, 39097602), (9, 39098688, 39100135),
]
SITIOS = [
    ("SLICK1", 39099213, "c.1382del", "frameshift", "#1b7837"),
    ("SLICK3_cand", 39099226, "c.1393C>A", "p.Ser465*", "#8c6d31"),
    ("SLICK2", 39099321, "c.1489C>T", "p.Arg497*", "#2166ac"),
]

# --- (b) consequencia funcional (Metodos: bcftools csq, Ensembl 113) ---
CONSEQUENCIAS = [("Synonymous", 24, "#7fbf7b"), ("Non-synonymous", 32, "#af8dc3"),
                  ("Intronic", 1360, "#999999"), ("UTR / splice / unclassified", 38, "#d9d9d9")]

PAINEL_COR = {"EUR_TAURINE": "#2166ac", "NDAMA": "#1b7837", "CRIOLLO": "#b2abd2",
              "ANKOLE": "#8c6d31", "BORAN": "#d6604d", "OGADEN": "#4d4d4d"}
PAINEL_ROTULO = {"EUR_TAURINE": "European\ntaurine", "NDAMA": "N'Dama", "CRIOLLO": "Criollo",
                  "ANKOLE": "Ankole", "BORAN": "Boran", "OGADEN": "Ogaden"}


def ler_frequencias(gt_dir: pathlib.Path, groups_tsv: pathlib.Path):
    """Le sample->painel e os 3 TSV de genotipo; devolve {alelo: {painel: (copias_alt, copias_totais)}}."""
    painel_de = {}
    with open(groups_tsv, encoding="utf-8") as f:
        for linha in f:
            partes = linha.rstrip("\n").split("\t")
            if len(partes) >= 2:
                painel_de[partes[0]] = partes[1]

    freqs = {}
    for alelo, arq in (("SLICK1", "gt_slick1.tsv"), ("SLICK2", "gt_slick2.tsv"),
                        ("SLICK3_cand", "gt_slick3.tsv")):
        contagem = {}
        with open(gt_dir / arq, encoding="utf-8") as f:
            for linha in f:
                amostra, gt = linha.rstrip("\n").split("\t")
                painel = painel_de.get(amostra)
                if painel is None or painel not in PAINEL_COR:
                    continue
                alelos = gt.replace("|", "/").split("/")
                alelos = [a for a in alelos if a in ("0", "1")]
                if not alelos:
                    continue
                alt, tot = contagem.get(painel, (0, 0))
                contagem[painel] = (alt + sum(1 for a in alelos if a == "1"), tot + len(alelos))
        freqs[alelo] = contagem
    return freqs


def montar_figura(freqs: dict, out_pdf: pathlib.Path):
    fig = plt.figure(figsize=(13, 9.5))
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 2.1], hspace=0.42, wspace=0.28)
    ax_gene = fig.add_subplot(gs[0, :])
    ax_cons = fig.add_subplot(gs[1, 0])
    ax_freq = fig.add_subplot(gs[1, 1])

    # ---------- (a) modelo do gene ----------
    g_ini, g_fim = EXONS[0][1], EXONS[-1][2]
    ax_gene.add_patch(FancyArrowPatch((g_ini, 0), (g_fim, 0), color="#888", lw=1.4,
                                       arrowstyle="-|>", mutation_scale=14, zorder=1))
    for rank, ini, fim in EXONS:
        largura = max(fim - ini, 400)
        ax_gene.add_patch(Rectangle((ini, -0.09), largura, 0.18, facecolor="#4d4d4d",
                                     edgecolor="black", linewidth=0.5, zorder=2))
    # os 3 sitios ficam a poucas dezenas de pb um do outro (mesmo exon final) --
    # rotulo direto colidiria. Espalha os rotulos num degrau horizontal com
    # linha de chamada ate a posicao real (leader line).
    vao = g_fim - g_ini
    deslocamentos_rotulo = [-0.11, 0.0, 0.11]
    for (nome, pos, cnot, pnot, cor), dx_frac in zip(SITIOS, deslocamentos_rotulo):
        ax_gene.axvline(pos, color=cor, lw=1.3, zorder=3, ymin=0.12, ymax=0.5)
        ax_gene.plot(pos, 0.09, marker="v", color=cor, ms=8, zorder=4)
        x_rotulo = pos + dx_frac * vao
        ax_gene.plot([pos, pos, x_rotulo], [0.5, 0.58, 0.58], color=cor, lw=0.9,
                      zorder=3, solid_capstyle="round")
        ax_gene.text(x_rotulo, 0.62, f"{nome}\n{cnot}\n{pnot}", ha="center", va="bottom",
                      fontsize=7.6, color=cor, fontweight="bold", linespacing=1.35)
    ax_gene.set_xlim(g_ini - 3000, g_fim + max(0.14 * vao, 3000))
    ax_gene.set_ylim(-0.5, 1.15)
    ax_gene.set_yticks([])
    ax_gene.set_xlabel("BTA20 (ARS-UCD1.2), posição em pb", fontsize=9)
    ax_gene.set_title("(a) Modelo do gene *PRLR* — 9 éxons, sítios focais no éxon final",
                       fontsize=10.5, fontweight="bold", loc="left")
    for side in ("top", "right", "left"):
        ax_gene.spines[side].set_visible(False)

    # ---------- (b) consequencia funcional ----------
    labels = [c[0] for c in CONSEQUENCIAS]
    valores = [c[1] for c in CONSEQUENCIAS]
    cores = [c[2] for c in CONSEQUENCIAS]
    y = range(len(labels))
    ax_cons.barh(list(y), valores, color=cores, edgecolor="black", linewidth=0.6)
    for yi, v in zip(y, valores):
        ax_cons.text(v * 1.02, yi, f"{v:,}", va="center", fontsize=9)
    ax_cons.set_yticks(list(y))
    ax_cons.set_yticklabels(labels, fontsize=9)
    ax_cons.set_xscale("log")
    ax_cons.set_xlabel("sítios (escala log)", fontsize=9)
    ax_cons.set_title("(b) Consequência funcional dos sítios do gene faseado",
                       fontsize=10.5, fontweight="bold", loc="left")
    ax_cons.spines[["top", "right"]].set_visible(False)

    # ---------- (c) frequencia por painel ----------
    alelos = list(freqs.keys())
    paineis = [p for p in PAINEL_COR if any(p in freqs[a] for a in alelos)]
    n_alelos = len(alelos)
    largura_barra = 0.8 / n_alelos
    x = range(len(paineis))
    for i, alelo in enumerate(alelos):
        cor_alelo = next(c for n, _, _, _, c in SITIOS if n == alelo)
        vals = []
        for p in paineis:
            alt, tot = freqs[alelo].get(p, (0, 0))
            vals.append(alt / tot if tot else 0.0)
        deslocamento = [xi + (i - n_alelos / 2 + 0.5) * largura_barra for xi in x]
        ax_freq.bar(deslocamento, vals, width=largura_barra * 0.92, color=cor_alelo,
                    edgecolor="black", linewidth=0.4, label=alelo)
    ax_freq.set_xticks(list(x))
    ax_freq.set_xticklabels([PAINEL_ROTULO[p] for p in paineis], fontsize=8.5)
    ax_freq.set_ylabel("frequência alélica no painel", fontsize=9)
    ax_freq.set_title("(c) Frequência alélica por painel geográfico\n"
                       "(não por raça fina — ver nota no rodapé)",
                       fontsize=10.5, fontweight="bold", loc="left")
    ax_freq.legend(fontsize=8, frameon=False, loc="upper right")
    ax_freq.spines[["top", "right"]].set_visible(False)
    ax_freq.grid(axis="y", alpha=0.25)

    fig.suptitle("Figura 1. Os alelos truncantes do receptor de prolactina",
                  fontsize=12.5, fontweight="bold", y=0.985)
    for ext in ("pdf", "png"):
        fig.savefig(out_pdf.with_suffix(f".{ext}"), dpi=300, bbox_inches="tight")
    print(f"-> {out_pdf.with_suffix('.pdf')}")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gt-dir", type=pathlib.Path,
                     default=DATA_DIR / "genotypes")
    ap.add_argument("--groups", type=pathlib.Path,
                     default=pathlib.Path(__file__).parents[1] /
                     "analises/03_slick/results/groups_afr_split.tsv")
    ap.add_argument("--out", type=pathlib.Path, default=OUT / "Fig5_PRLR_alleles.pdf")
    args = ap.parse_args()

    freqs = ler_frequencias(args.gt_dir, args.groups)
    montar_figura(freqs, args.out)


if __name__ == "__main__":
    main()
