#!/usr/bin/env python3
"""
fig_genomic_ancestry.py — figuras da ancestralidade genomica, geradas dos
resultados reais.

Fig 1  calibracao do ensaio: dose conhecida x dose estimada, e o limite de
       deteccao em funcao da densidade de marcadores. Responde "o ensaio
       enxerga o que dizemos que enxerga?".
Fig 2  o achado no contexto: todas as populacoes ordenadas por ancestralidade
       africana aparente, coloridas pelo grupo. Mostra que o valor do Senepol
       so significa alguma coisa contra o controle certo.
Fig 3  os dois controles que sustentam a conclusao, lado a lado, e a
       concordancia entre os dois estimadores.

Idempotente. Le apenas TSV ja gravado pelas analises.
"""
from __future__ import annotations
import os
import pathlib
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

DATA_DIR = pathlib.Path(os.environ.get("SLICK_DATA_DIR", "data"))
FAIR = DATA_DIR / "senepol_widde/decker2014/plink/fair"
ZEBU = DATA_DIR / "senepol_widde/decker2014/plink/zebu"
OUT = pathlib.Path(__file__).parent / "figuras"
OUT.mkdir(exist_ok=True)

C_SEN = "#c0392b"      # Senepol
C_BRIT = "#2471a3"     # linhagem britanica (controle pareado)
C_CONT = "#7f8c8d"     # continentais e alpinas
C_COMP = "#f39c12"     # compostos com zebu, sem Africa
C_REF = "#27ae60"
C_CRIO = "#8e44ad"      # criollos de origem iberica

BRIT = {"REDP", "DEV", "SUSS", "LINC", "SDEV", "SH", "GALL", "DEX", "KERR"}
COMP = {"SGT", "BEFM"}
CRIOLLO = {"TXLH", "ROMO"}   # descendentes de gado iberico: Africa plausivel via Iberia
NOME = {"REDP": "Red Poll", "DEV": "Devon", "SUSS": "Sussex", "LINC": "Lincoln Red",
        "SDEV": "South Devon", "SH": "Shorthorn", "GALL": "Galloway", "DEX": "Dexter",
        "KERR": "Kerry", "SGT": "Santa Gertrudis", "BEFM": "Beefmaster",
        "SENP": "SENEPOL", "JER": "Jersey", "PIED": "Piedmontese", "GAS": "Gascon",
        "BSW": "Brown Swiss", "GNS": "Guernsey", "TAR": "Tarentaise", "VOS": "Vosgienne",
        "SIM": "Simmental", "AUB": "Aubrac", "NORM": "Normande", "GEL": "Gelbvieh",
        "CHA": "Charolais", "MAAN": "Maine Anjou", "MRI": "Meuse-Rhine-Issel",
        "BBLU": "Belgian Blue", "AYR": "Ayrshire", "ANR": "Angus Red",
        "TXLH": "Texas Longhorn", "ROMO": "Romosinuano"}

plt.rcParams.update({"font.size": 9, "axes.spines.top": False,
                     "axes.spines.right": False, "figure.dpi": 150})


def save(fig, stem):
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"{stem}.{ext}", bbox_inches="tight", dpi=300)
    plt.close(fig)
    print(f"  -> {stem}.png / .pdf")


# ---------------------------------------------------------------- Figura 1
def fig1():
    t = pd.read_csv(FAIR / "titration.tsv", sep="\t")
    fig, (a, b) = plt.subplots(1, 2, figsize=(9.5, 3.9),
                               gridspec_kw={"width_ratios": [1.35, 1]})

    x, y = t["dose_nominal_pct"], t["aft_estimado_pct"]
    a.fill_between(x, t["ic_low"], t["ic_high"], color=C_REF, alpha=.25, lw=0)
    a.plot(x, y, "o-", color=C_REF, ms=5, lw=1.8, label="estimated (95% CI)")
    lim = float(t.iloc[0]["ic_high"])
    a.axhline(float(t.iloc[0]["aft_estimado_pct"]), color="0.4", ls=":", lw=1.2)
    a.text(26, float(t.iloc[0]["aft_estimado_pct"]) + .8,
           f"zero-dose baseline {t.iloc[0]['aft_estimado_pct']:.2f}%", color="0.35", fontsize=8)
    a.plot([0, 50], [t.iloc[0]["aft_estimado_pct"], t.iloc[0]["aft_estimado_pct"] + 50],
           ls="--", lw=1, color="0.6", label="baseline + nominal dose")
    a.scatter([0.5], [t.loc[t.dose_nominal_pct == 0.5, "aft_estimado_pct"]],
              s=110, facecolor="none", edgecolor=C_SEN, lw=1.8, zorder=5)
    a.annotate("detection limit 0.5%", xy=(0.5, 7.97), xytext=(7, 4.2),
               arrowprops=dict(arrowstyle="->", color=C_SEN, lw=1.2),
               color=C_SEN, fontsize=8.5)
    a.set_xlabel("N'Dama introgressed (nominal %)")
    a.set_ylabel("African ancestry estimated (%)")
    a.set_title("a  Assay calibration", loc="left", fontweight="bold")
    a.legend(frameon=False, fontsize=8, loc="upper left")

    dens, lims = [2000, 5000, 10000, 29723], [5.0, 1.0, 1.0, 0.5]
    b.plot(dens, lims, "s-", color=C_BRIT, ms=7, lw=1.8)
    for d, l in zip(dens, lims):
        b.annotate(f"{l}%", (d, l), textcoords="offset points", xytext=(0, 9),
                   ha="center", fontsize=8.5, color=C_BRIT)
    b.axhline(0.6, color=C_SEN, ls="--", lw=1.3)
    b.text(2400, 0.66, "0.6% — published Senepol estimate", color=C_SEN, fontsize=8)
    b.set_xscale("log")
    b.set_xticks(dens)
    b.set_xticklabels(["2k", "5k", "10k", "30k"])
    b.set_ylim(0, 5.9)
    b.set_xlabel("markers")
    b.set_ylabel("detection limit (%)")
    b.set_title("b  Resolution vs marker density", loc="left", fontweight="bold")
    save(fig, "Fig1_calibration")


# ---------------------------------------------------------------- Figura 2
def fig2():
    d = pd.read_csv(ZEBU / "zebu_control.tsv", sep="\t").sort_values("AFT")

    def grupo(p):
        return ("SENP" if p == "SENP" else "COMP" if p in COMP
                else "CRIO" if p in CRIOLLO
                else "BRIT" if p in BRIT else "CONT")
    d["g"] = d["pop"].map(grupo)
    cor = {"SENP": C_SEN, "BRIT": C_BRIT, "CONT": C_CONT, "COMP": C_COMP, "CRIO": C_CRIO}

    fig, ax = plt.subplots(figsize=(6.2, 8.2))
    y = np.arange(len(d))
    ax.barh(y, d["AFT"], color=[cor[g] for g in d["g"]],
            edgecolor="white", lw=.6, height=.78)
    ax.set_yticks(y)
    ax.set_yticklabels([NOME.get(p, p) for p in d["pop"]], fontsize=8)
    for i, (v, g) in enumerate(zip(d["AFT"], d["g"])):
        ax.text(v + .12, i, f"{v:.2f}", va="center", fontsize=7.5,
                color=cor[g], fontweight="bold" if g == "SENP" else "normal")

    brit_m = d.loc[d.g == "BRIT", "AFT"].mean()
    ax.axvline(brit_m, color=C_BRIT, ls="--", lw=1.4)
    ax.text(brit_m + .1, -1.6, f"British mean {brit_m:.2f}%", color=C_BRIT, fontsize=8)
    ax.set_xlabel("apparent African taurine ancestry (%)")
    ax.set_ylim(-2.2, len(d) - .3)
    ax.set_title("Every population estimated against the same fixed panel",
                 loc="left", fontweight="bold", fontsize=10)
    ax.legend(handles=[
        Line2D([], [], color=C_SEN, lw=7, label="Senepol (N'Dama × Red Poll)"),
        Line2D([], [], color=C_BRIT, lw=7, label="British lineage — matched control"),
        Line2D([], [], color=C_CONT, lw=7, label="continental / alpine European"),
        Line2D([], [], color=C_COMP, lw=7, label="British × Brahman, no African founder"),
        Line2D([], [], color=C_CRIO, lw=7, label="Iberian-derived Criollo")],
        frameon=False, fontsize=8, loc="lower right")
    save(fig, "Fig2_ranking")


# ---------------------------------------------------------------- Figura 3
def fig3():
    d = pd.read_csv(ZEBU / "zebu_control.tsv", sep="\t").set_index("pop")
    s = pd.read_csv(ZEBU / "snmf_labeled.tsv", sep="\t")
    fig, (a, b) = plt.subplots(1, 2, figsize=(9.6, 4.0),
                               gridspec_kw={"width_ratios": [1.25, 1]})

    grupos = ["British\nlineage", "Santa\nGertrudis", "Beefmaster", "SENEPOL"]
    brit = d.loc[[p for p in BRIT if p in d.index]]
    vals = np.array([
        [brit["EUT"].mean(), brit["AFT"].mean(), brit["ZEB"].mean()],
        [d.loc["SGT", "EUT"], d.loc["SGT", "AFT"], d.loc["SGT", "ZEB"]],
        [d.loc["BEFM", "EUT"], d.loc["BEFM", "AFT"], d.loc["BEFM", "ZEB"]],
        [d.loc["SENP", "EUT"], d.loc["SENP", "AFT"], d.loc["SENP", "ZEB"]]])
    x = np.arange(4)
    a.bar(x, vals[:, 0], .62, color="#d5d8dc", edgecolor="white", label="European")
    a.bar(x, vals[:, 2], .62, bottom=vals[:, 0], color=C_COMP, edgecolor="white",
          alpha=.85, label="Zebu")
    a.bar(x, vals[:, 1], .62, bottom=vals[:, 0] + vals[:, 2], color=C_SEN,
          edgecolor="white", label="African")
    for i in range(4):
        a.text(i, vals[i, 0] + vals[i, 2] + vals[i, 1] + 1.6, f"AFT\n{vals[i,1]:.2f}%",
               ha="center", fontsize=8.5, color=C_SEN, fontweight="bold")
    a.set_xticks(x)
    a.set_xticklabels(grupos, fontsize=8.5)
    a.set_ylabel("ancestry (%)")
    a.set_ylim(0, 118)
    a.legend(frameon=False, fontsize=8, ncol=3, loc="lower center")
    a.set_title("a  Zebu does not drag the African component",
                loc="left", fontweight="bold")

    def snmf(pops):
        return s.loc[s.FID.isin(pops), "AFT"].mean() * 100
    pares = [("British\nlineage", brit["AFT"].mean(), snmf(list(BRIT))),
             ("Santa\nGertrudis", d.loc["SGT", "AFT"], snmf(["SGT"])),
             ("Beefmaster", d.loc["BEFM", "AFT"], snmf(["BEFM"])),
             ("SENEPOL", d.loc["SENP", "AFT"], snmf(["SENP"]))]
    xs = np.arange(len(pares))
    w = .36
    b.bar(xs - w / 2, [p[1] for p in pares], w, color=C_BRIT, label="supervised")
    b.bar(xs + w / 2, [p[2] for p in pares], w, color=C_REF, label="sNMF (unsupervised)")
    for i, p in enumerate(pares):
        b.text(i - w / 2, p[1] + .18, f"{p[1]:.2f}", ha="center", fontsize=7.5)
        b.text(i + w / 2, p[2] + .18, f"{p[2]:.2f}", ha="center", fontsize=7.5)
    b.set_xticks(xs)
    b.set_xticklabels([p[0] for p in pares], fontsize=8.5)
    b.set_ylabel("African ancestry (%)")
    b.set_ylim(0, 11.4)
    b.legend(frameon=False, fontsize=8, loc="upper left")
    b.set_title("b  Two estimators agree", loc="left", fontweight="bold")
    save(fig, "Fig3_controls")



# ---------------------------------------------------------------- Figura 4
def fig4():
    """Sensibilidade ao numero de eixos ancestrais: 3 vs 4 componentes."""
    pops = ["Beefmaster", "Santa\nGertrudis", "British\nlineage", "SENEPOL",
            "Romosinuano", "Texas\nLonghorn"]
    tres = [0.00, 0.99, 2.46, 8.16, 12.42, 13.78]      # 3 componentes
    quatro = [0.00, 0.07, 0.46, 3.85, 5.42, 6.45]      # 4 componentes (eixo NAF)
    cor = [C_COMP, C_COMP, C_BRIT, C_SEN, C_CRIO, C_CRIO]
    fig, (a, b) = plt.subplots(1, 2, figsize=(9.6, 4.0),
                               gridspec_kw={"width_ratios": [1.3, 1]})
    x = np.arange(len(pops)); w = .38
    a.bar(x - w/2, tres, w, color=cor, alpha=.55, edgecolor="white",
          label="3 ancestries (West African panel only)")
    a.bar(x + w/2, quatro, w, color=cor, edgecolor="white",
          label="4 ancestries (+ North African axis)")
    for i, (v1, v2) in enumerate(zip(tres, quatro)):
        a.text(i - w/2, v1 + .25, f"{v1:.2f}", ha="center", fontsize=7)
        a.text(i + w/2, v2 + .25, f"{v2:.2f}", ha="center", fontsize=7, fontweight="bold")
    a.set_xticks(x); a.set_xticklabels(pops, fontsize=8)
    a.set_ylabel("African taurine ancestry (%)")
    a.legend(frameon=False, fontsize=8, loc="upper left")
    a.set_title("a  The absolute value depends on the axes offered",
                loc="left", fontweight="bold")
    exc3 = [t - 2.46 for t in tres]; exc4 = [q - 0.46 for q in quatro]
    b.axhline(0, color="0.6", lw=1)
    b.plot(x, exc3, "o--", color="0.5", ms=6, label="3 ancestries")
    b.plot(x, exc4, "s-", color=C_SEN, ms=6, label="4 ancestries")
    b.annotate(f"+{exc3[3]:.2f}", (3, exc3[3]), textcoords="offset points",
               xytext=(6, 6), fontsize=8, color="0.4")
    b.annotate(f"+{exc4[3]:.2f}", (3, exc4[3]), textcoords="offset points",
               xytext=(6, -14), fontsize=8, color=C_SEN, fontweight="bold")
    b.set_xticks(x); b.set_xticklabels(pops, fontsize=8)
    b.set_ylabel("excess over matched control (pp)")
    b.legend(frameon=False, fontsize=8, loc="upper left")
    b.set_title("b  The excess over control is what holds",
                loc="left", fontweight="bold")
    save(fig, "Fig4_sensitivity")

if __name__ == "__main__":
    print("gerando figuras a partir dos resultados reais:")
    fig1(); fig2(); fig3(); fig4()
    print(f"\nfiguras em {OUT}")
