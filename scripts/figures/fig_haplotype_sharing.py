#!/usr/bin/env python3
"""
fig_haplotype_sharing.py — figuras de compartilhamento de haplotipo, geradas
dos TSV de resultado.

Fig 1  o contraste que e o achado: comprimento do trato compartilhado pelos
       portadores do SLICK1 e do SLICK2, contra os grupos ancestrais.
Fig 2  a mesma medida convertida em tempo, com as ancoras historicas que
       importam (fundacao do Senepol, chegada do gado iberico).
"""
from __future__ import annotations
import pathlib
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

R = pathlib.Path(__file__).parent.parent / "results"
OUT = pathlib.Path(__file__).parent.parent / "PAPER" / "figuras"
OUT.mkdir(parents=True, exist_ok=True)

C1 = "#c0392b"   # SLICK1
C2 = "#2471a3"   # SLICK2
CG = "#95a5a6"   # grupos ancestrais

NOME = {"PORTADORES homozigotos": "carriers\n(homozygous)",
        "portadores heterozigotos": "carriers\n(heterozygous)",
        "AFR_TAURINE": "African\ntaurine", "AFR_SANGA": "African\nsanga",
        "CRIOLLO": "Criollo\n(non-carrier)", "EUR_TAURINE": "European\ntaurine",
        "ZEBU": "Zebu"}
ORD = ["PORTADORES homozigotos", "portadores heterozigotos", "AFR_TAURINE",
       "AFR_SANGA", "CRIOLLO", "EUR_TAURINE", "ZEBU"]

plt.rcParams.update({"font.size": 9, "axes.spines.top": False,
                     "axes.spines.right": False, "figure.dpi": 150})


def save(fig, stem):
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"{stem}.{ext}", bbox_inches="tight", dpi=300)
    plt.close(fig)
    print(f"  -> {stem}.png / .pdf")


def load(tag):
    d = pd.read_csv(R / f"sharing_{tag}.tsv", sep="\t").set_index("grupo")
    return d.reindex([g for g in ORD if g in d.index])


def fig1():
    d1, d2 = load("SLICK1"), load("SLICK2")
    labels = [NOME[g] for g in d1.index]
    y = np.arange(len(labels))
    fig, ax = plt.subplots(figsize=(7.4, 4.3))
    h = .38
    ax.barh(y + h/2, d1["kb_mediano"], h, color=C1, label="SLICK1  (Costeño con Cuernos + Limonero)")
    ax.barh(y - h/2, d2["kb_mediano"], h, color=C2, label="SLICK2  (Limonero only)")
    for i, (v1, v2) in enumerate(zip(d1["kb_mediano"], d2["kb_mediano"])):
        if v1 > 0: ax.text(v1*1.25, i + h/2, f"{v1:,.0f}", va="center", fontsize=7.5, color=C1)
        if v2 > 0: ax.text(v2*1.25, i - h/2, f"{v2:,.0f}", va="center", fontsize=7.5, color=C2)
    ax.set_xscale("symlog", linthresh=1)
    ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel("shared tract with the carrier haplotype (kb, log scale)")
    ax.legend(frameon=False, fontsize=8, loc="lower right")
    ax.set_title("Two truncating alleles, two haplotype histories",
                 loc="left", fontweight="bold", fontsize=10)
    save(fig, "Fig1_sharing")


def fig2():
    fig, ax = plt.subplots(figsize=(7.6, 3.4))
    dados = [("SLICK1", 5177.2, 19, C1), ("SLICK2", 155.7, 642, C2)]
    for nome, kb, gen, c in dados:
        ax.barh(nome, gen, color=c, height=.5)
        ax.text(gen*1.06, nome, f"~{gen} generations   ({kb:,.0f} kb tract)",
                va="center", fontsize=8.5, color=c, fontweight="bold")
    for x, lab, ls in [(22, "Senepol founded, 1918", "--"),
                       (105, "Iberian cattle reach the Americas, ~1500", ":")]:
        ax.axvline(x, color="0.35", ls=ls, lw=1.2)
        ax.text(x, 1.62, lab, rotation=90, fontsize=7.5, color="0.35",
                ha="right", va="top")
    ax.set_xscale("log")
    ax.set_xlim(8, 3000)
    ax.set_xlabel("generations since expansion (implied by tract length; 1 cM = 1 Mb)")
    ax.set_ylim(-.6, 1.9)
    ax.set_title("SLICK1 postdates the founding of the modern Criollo breeds; SLICK2 predates colonisation",
                 loc="left", fontweight="bold", fontsize=9.5)
    save(fig, "Fig2_timing")


if __name__ == "__main__":
    print("gerando figuras:")
    fig1(); fig2()
    print(f"\nfiguras em {OUT}")
