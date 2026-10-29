# -*- coding: utf-8 -*-
"""Figuras (2026-10-26). Substitui as cinco figuras da versao Animal Genetics.

Diagnostico que motivou a reconstrucao (auditoria do autor + inventario de 2026-10-20): das
cinco figuras antigas, so a rede e o ranking de ancestralidade tinham dado denso. A Fig. 2 era a
Tabela 2 desenhada; a Fig. 1b eram 4 numeros; a 1c tinha 18 barras, 15 delas zero; a Fig. 5
eram 3 numeros ja escritos no texto. Aqui cada figura carrega dado que uma tabela nao mostra.

  Figura 1  O locus e os seis alelos publicados (funde as antigas 1a e 5; acrescenta SLICK4-6).
  Figura 2  Rede de haplotipos codificantes, fundo UNICO, quatro alelos e o painel iberico.
  Figura 3  Robustez da afinidade: haplotipo a haplotipo, rarefacao dos paineis, lote.
  Figura 4  Codificante x flanco: paisagem de informatividade e atribuicao por distancia.

Todo numero vem de `tables/derived/resultado_final.json`, nunca digitado aqui, exceto as
coordenadas e os dominios da Figura 1, que sao anotacao e estao citados na legenda.

Uso:  python make_figuras_jabg.py [1 2 3 4]   (sem argumento, gera todas)
"""
from __future__ import annotations

import collections
import json
import pathlib
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.patheffects as pe  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch, Rectangle  # noqa: E402

AQUI = pathlib.Path(__file__).resolve().parents[2]   # raiz do repositorio
RAIZ = AQUI
OUT = AQUI / "figures"
RES = AQUI / "tables" / "derived" / "resultado_final.json"
sys.path.insert(0, str(RAIZ / "scripts" / "lib"))

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8, "axes.linewidth": 0.6,
                     "pdf.fonttype": 42, "svg.fonttype": "none"})

# Cores por painel: as mesmas da rede antiga, mais o iberico.
COR = {"EUR_TAURINE": "#2166ac", "IBERIAN": "#e08214", "NDAMA": "#1b7837",
       "CRIOLLO": "#b2abd2", "ANKOLE": "#8c6d31", "BORAN": "#d6604d", "OGADEN": "#4d4d4d"}
ROT = {"EUR_TAURINE": "Continental European taurine", "IBERIAN": "Iberian taurine",
       "NDAMA": "N'Dama", "CRIOLLO": "Criollo", "ANKOLE": "Ankole", "BORAN": "Boran",
       "OGADEN": "Ogaden"}
# Cor de cada alelo (anel na rede, marcador na Fig. 1)
ALELO_COR = {"SLICK1": "#1b7837", "SLICK2": "#2166ac", "SLICK3": "#e08214",
             "SLICK6": "#a6611a", "SLICK4": "#999999", "SLICK5": "#999999"}


def salvar(fig, nome):
    OUT.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(OUT / f"{nome}.{ext}", dpi=600, bbox_inches="tight")
    plt.close(fig)
    print(f"  -> {OUT / nome}.pdf / .png")


def letra(ax, txt, x=-0.08, y=1.04):
    """Letra de painel em negrito, sem titulo: o titulo pertence a legenda."""
    ax.text(x, y, txt, transform=ax.transAxes, fontsize=10, fontweight="bold",
            va="bottom", ha="left")


# =====================================================================================
# FIGURA 1 - O locus e os seis alelos
# =====================================================================================
# Exons do transcrito ENSBTAT00000014437 (GFF3 Ensembl, ARS-UCD1.2): (rank, inicio, fim).
EXONS = [(1, 39037656, 39037692), (2, 39077900, 39078016), (3, 39080773, 39080905),
         (4, 39082809, 39082978), (5, 39093142, 39093311), (6, 39094967, 39095108),
         (7, 39096738, 39096837), (8, 39097533, 39097602), (9, 39098688, 39100135)]
# Os seis alelos publicados (OMIA:001372-9913; coordenadas ARS-UCD1.2, HGVS em NM_001039726.2).
# ultimo residuo nativo retido -> comprimento do coto citoplasmatico (cauda = 259-581).
ALELOS = [  # nome, pos, c., p., ultimo_aa_retido, status
    ("SLICK4", 39099113, "c.1281C>G", "p.Tyr427Ter", 426, "monomorphic"),
    ("SLICK1", 39099213, "c.1382del", "p.Ala461ValfsTer2", 461, "carriers"),
    ("SLICK3", 39099226, "c.1394C>A", "p.Ser465Ter", 464, "carriers"),
    ("SLICK5", 39099228, "c.1396A>T", "p.Lys466Ter", 465, "candidate"),
    ("SLICK6", 39099267, "c.1435C>T", "p.Gln479Ter", 478, "carriers"),
    ("SLICK2", 39099321, "c.1489C>T", "p.Arg497Ter", 496, "carriers"),
]
L_PROT, EXTRA, TM, CYTO = 581, (25, 237), (238, 258), (259, 581)
BOX1 = (262, 271)      # por homologia com o PRLR humano; nao e anotacao UniProt
STAT5_Y = 560          # tirosina distal de docking de STAT5, por homologia


def fig1(res):
    n_hap = {a: res["alelos"][a]["n_haplotipos"] for a in res["alelos"]}
    fig = plt.figure(figsize=(7.2, 6.4))
    gs = fig.add_gridspec(3, 1, height_ratios=[0.55, 1.15, 1.35], hspace=0.62)

    # ---- (a) modelo do gene
    ax = fig.add_subplot(gs[0])
    g0, g1 = EXONS[0][1], EXONS[-1][2]
    ax.plot([g0, g1], [0, 0], color="#555555", lw=1.0, zorder=1)
    for rank, ini, fim in EXONS:
        final = rank == EXONS[-1][0]
        ax.add_patch(Rectangle((ini, -0.32), max(fim - ini, 260), 0.64,
                               facecolor="#c51b7d" if final else "#bdbdbd",
                               edgecolor="black", lw=0.5, zorder=2))
        ax.text((ini + fim) / 2, 0.52, str(rank), ha="center", va="bottom", fontsize=6.5)
    Z0, Z1 = 39099090, 39099345            # janela do painel b
    # 255 pb somem na escala do gene: marca-se a posicao da janela do painel b com uma seta
    ax.annotate("b", xy=((Z0 + Z1) / 2, -0.36), xytext=((Z0 + Z1) / 2 - 2500, -0.9),
                fontsize=7, fontweight="bold", ha="right", va="center",
                arrowprops=dict(arrowstyle="-|>", color="black", lw=0.7))
    ax.set_xlim(g0 - 1500, g1 + 1500)
    ax.set_ylim(-1.0, 1.1)
    ax.set_yticks([])
    for s in ("left", "right", "top"):
        ax.spines[s].set_visible(False)
    ax.set_xticks(np.arange(39040000, 39100001, 20000))
    ax.set_xticklabels([f"{x/1e6:.2f}" for x in np.arange(39040000, 39100001, 20000)])
    ax.set_xlabel("BTA20 position (Mb, ARS-UCD1.2); exons of ENSBTAT00000014437", fontsize=7.5)
    letra(ax, "a", y=1.0)

    # ---- (b) zoom no exon terminal: os seis alelos, rotulos em faixas equidistantes
    ax = fig.add_subplot(gs[1])
    ax.add_patch(Rectangle((Z0, -0.12), Z1 - Z0, 0.24, facecolor="#f3c6dd",
                           edgecolor="#c51b7d", lw=0.6, zorder=1))
    slots = np.linspace(Z0 + 22, Z1 - 22, len(ALELOS))
    for (nome, p, c, prot, _, st), xs in zip(ALELOS, slots):
        cheio = st == "carriers"
        ax.plot([p, p], [0.12, 0.42], color="#333333", lw=0.7, zorder=2)
        ax.plot([p, xs], [0.42, 0.95], color="#333333", lw=0.7, zorder=2)
        ax.plot(xs, 0.95, marker="o", ms=8 if cheio else 6.5, zorder=3,
                mfc=ALELO_COR[nome] if cheio else "white",
                mec="black" if cheio else "#777777", mew=0.8 if cheio else 1.0)
        if st == "carriers":
            nota = f"{n_hap.get(nome, 0)} carrier hap."
        elif st == "monomorphic":
            nota = "absent (0 copies)"
        else:
            nota = "not assigned"
        ax.text(xs, 1.12, "\n".join((nome, c, prot, nota)), ha="center", va="bottom",
                fontsize=6.2, linespacing=1.2, fontweight="normal",
                color="black" if cheio else "#666666")
    ax.set_xlim(Z0, Z1)
    ax.set_ylim(-0.3, 2.35)
    ax.set_yticks([])
    for s in ("left", "right", "top"):
        ax.spines[s].set_visible(False)
    tk = np.arange(39099100, 39099341, 40)
    ax.set_xticks(tk)
    ax.set_xticklabels([f"{x:,}" for x in tk], fontsize=6.5)
    ax.set_xlabel("Final exon, BTA20 bp: the six published truncating alleles lie within 208 bp",
                  fontsize=7.5)
    letra(ax, "b", y=1.0)

    # ---- (c) receptor e coto retido por alelo
    ax = fig.add_subplot(gs[2])
    y = len(ALELOS) + 0.6
    ax.add_patch(Rectangle((1, y), EXTRA[1], 0.5, facecolor="#d9d9d9", edgecolor="black", lw=0.5))
    ax.add_patch(Rectangle((TM[0], y), TM[1] - TM[0], 0.5, facecolor="#8a8f98", edgecolor="black", lw=0.5))
    ax.add_patch(Rectangle((CYTO[0], y), CYTO[1] - CYTO[0], 0.5, facecolor="#f7f7f7", edgecolor="black", lw=0.5))
    ax.add_patch(Rectangle((BOX1[0], y), BOX1[1] - BOX1[0], 0.5, facecolor="#4dac26", edgecolor="black", lw=0.4))
    ax.plot(STAT5_Y, y + 0.25, marker="o", ms=6, mfc="#d62728", mec="black", zorder=5)
    ax.text(120, y + 0.25, "extracellular", ha="center", va="center", fontsize=6.5)
    ax.text(TM[0] + 10, y + 0.72, "TM", ha="center", fontsize=6.5)
    ax.text(BOX1[0] + 5, y - 0.22, "Box1", ha="center", va="top", fontsize=6.2, color="#2b7a0b")
    ax.text(STAT5_Y, y + 0.72, "STAT5 Tyr", ha="center", fontsize=6.2, color="#d62728")
    ax.text(420, y + 0.72, "cytoplasmic tail", ha="center", fontsize=6.5)
    for k, (nome, _, _, prot, ultimo, st) in enumerate(ALELOS):
        yy = len(ALELOS) - 1 - k
        cheio = st == "carriers"
        ax.add_patch(Rectangle((CYTO[0], yy), ultimo - CYTO[0], 0.55,
                               facecolor=ALELO_COR[nome] if cheio else "white",
                               edgecolor="black" if cheio else "#888888",
                               hatch=None if cheio else "////", lw=0.6))
        ax.text(CYTO[0] - 6, yy + 0.27, nome, ha="right", va="center", fontsize=7,
                fontweight="bold" if cheio else "normal", color="black" if cheio else "#666666")
        ax.text(ultimo + 4, yy + 0.27, f"{ultimo} aa (−{L_PROT - ultimo})", ha="left",
                va="center", fontsize=6.3, zorder=4,
                bbox=dict(facecolor="white", edgecolor="none", pad=0.6))
    ax.axvline(STAT5_Y, color="#d62728", ls="--", lw=0.7, ymin=0.02, ymax=0.84, zorder=0)
    ax.axvspan(BOX1[0], BOX1[1], ymin=0.02, ymax=0.84, color="#4dac26", alpha=0.18, zorder=0)
    ax.set_xlim(-10, L_PROT + 40)
    ax.set_ylim(-0.4, len(ALELOS) + 1.6)
    ax.set_yticks([])
    for s in ("left", "right", "top"):
        ax.spines[s].set_visible(False)
    ax.set_xlabel("PRLR residue (wild type 581 aa, UniProt Q28172); bars end at the last residue "
                  "before the premature stop",
                  fontsize=7.5)
    letra(ax, "c", y=1.0)

    leg = [Line2D([], [], marker="o", ls="", ms=7, mfc="#666666", mec="black",
                  label="carrier haplotypes assigned in this study"),
           Line2D([], [], marker="o", ls="", ms=6, mfc="white", mec="#777777",
                  label="not assigned (see text)")]
    fig.legend(handles=leg, loc="lower center", ncol=2, frameon=False, fontsize=7,
               bbox_to_anchor=(0.5, -0.01))
    salvar(fig, "Figure_1")


# =====================================================================================
# Dados das analises novas (uma fonte por figura; nada digitado aqui)
# =====================================================================================
LAI = AQUI / "tables" / "derived" / "local_ancestry"
RFM = LAI
COR.update({"INDICINE": "#b2182b", "ASIAN_TAURINE": "#762a83",
            "TAURINO_EUROPEU": "#2166ac", "TAURINO_AFRICANO": "#1b7837"})
ROT.update({"INDICINE": "Indicine", "ASIAN_TAURINE": "North/East Asian taurine",
            "TAURINO_EUROPEU": "European taurine", "TAURINO_AFRICANO": "African taurine"})
CURTO = {"EUR_TAURINE": "Cont. Eur.", "IBERIAN": "Iberian", "NDAMA": "N'Dama",
         "INDICINE": "Indicine", "ASIAN_TAURINE": "Asian taur.", "CRIOLLO": "Criollo",
         "TAURINO_EUROPEU": "Eur. taurine", "TAURINO_AFRICANO": "Afr. taurine",
         "ANKOLE": "Ankole", "BORAN": "Boran", "OGADEN": "Ogaden"}
ALELOS_SEG = ["SLICK1", "SLICK3", "SLICK6", "SLICK2"]      # por coordenada


def _json(p):
    return json.load(open(p, encoding="utf-8"))


def le_msp(run):
    """(haplotipos, janelas) do .msp.tsv do RFMix."""
    linhas = (RFM / f"{run}.msp.tsv").read_text().splitlines()
    cod = {int(v): k for k, v in (x.split("=") for x in linhas[0].split(":", 1)[1].split())}
    haps = linhas[1].lstrip("#").split("\t")[6:]
    jan = []
    for l in linhas[2:]:
        c = l.split("\t")
        jan.append({"spos": int(c[1]), "epos": int(c[2]), "call": [cod[int(x)] for x in c[6:]]})
    return haps, jan


# =====================================================================================
# FIGURA 2 - validacao da inferencia de ancestralidade
# =====================================================================================
def fig4(_res):
    val = _json(RFM / "validacao.json")["runs"]
    fig = plt.figure(figsize=(7.2, 5.6))
    gs = fig.add_gridspec(2, 3, height_ratios=[1.0, 0.82], width_ratios=[1, 1.25, 0.95],
                          hspace=0.70, wspace=0.78)

    # ---- (a) e (b): matrizes de confusao da simulacao interna
    for col, run in enumerate(("m3", "m4")):
        ax = fig.add_subplot(gs[0, col])
        sim = val[run]["simulacao_interna"]
        pops, M = sim["populacoes"], np.array(sim["matriz"])
        im = ax.imshow(M, cmap="Blues", vmin=0, vmax=100, aspect="auto")
        for i in range(M.shape[0]):
            for j in range(M.shape[1]):
                if M[i, j] >= 1:
                    ax.text(j, i, f"{M[i, j]:.0f}", ha="center", va="center", fontsize=6.3,
                            color="white" if M[i, j] > 55 else "black",
                            fontweight="bold" if i == j else "normal")
        ax.set_xticks(range(len(pops)))
        ax.set_xticklabels([CURTO.get(p, p) for p in pops], rotation=45, ha="right", fontsize=6.5)
        ax.set_yticks(range(len(pops)))
        ax.set_yticklabels([CURTO.get(p, p) for p in pops], fontsize=6.5)
        ax.set_xlabel("simulated ancestry", fontsize=7)
        if col == 0:
            ax.set_ylabel("inferred ancestry", fontsize=7)
        pior = val[run]["pior_painel"]
        ok = val[run]["passa_todos_os_paineis"]
        ax.set_title(f"{len(pops)}-way reference\n"
                     + ("all panels $\\geq$70%" if ok else
                        f"{CURTO.get(pior['painel'], pior['painel'])} only {pior['acerto']:.0f}%"),
                     fontsize=7.5, color="#1b7837" if ok else "#b2182b", pad=4)
        for s in ax.spines.values():
            s.set_visible(False)
        letra(ax, "ab"[col], x=-0.30 if col == 0 else -0.22, y=1.22)
    # sem barra de cores: cada celula ja traz o numero, e a barra roubava o painel c

    # ---- (c) controle por omissao
    ax = fig.add_subplot(gs[0, 2])
    itens = []
    for run, rot in (("m4", "4-way"), ("c4", "4-way,\nreduced")):
        for pan, v in (val[run].get("controle_por_omissao") or {}).items():
            itens.append((f"{CURTO.get(pan, pan)}\n({rot})", v["certo_no_prlr"], v["n_hap"], pan))
    y = np.arange(len(itens))
    frac = [c / n for _, c, n, _ in itens]
    ax.barh(y, frac, color=[COR[p] for *_, p in itens], edgecolor="black", lw=0.5, height=0.65)
    for k, (rot, c, n, _) in enumerate(itens):
        ax.text(min(frac[k] + 0.03, 0.97), k, f"{c}/{n}", va="center", fontsize=6.3)
    ax.axvline(1.0, color="#333333", lw=0.7, ls=":")
    ax.set_yticks(y)
    ax.set_yticklabels([x[0] for x in itens], fontsize=6.2)
    ax.set_xlim(0, 1.22)
    ax.set_xticks([0, 0.5, 1.0])
    ax.set_xlabel("held-out reference haplotypes\nassigned to their own panel", fontsize=7)
    ax.invert_yaxis()
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    letra(ax, "c", x=-0.42, y=1.22)

    # ---- (d) controle positivo: africanos de mistura conhecida
    ax = fig.add_subplot(gs[1, :2])
    mist = val["m3"]["controle_positivo_mistura"]
    pops = ["TAURINO_EUROPEU", "TAURINO_AFRICANO", "INDICINE"]
    nomes = list(mist)
    x = np.arange(len(nomes))
    base = np.zeros(len(nomes))
    for p in pops:
        v = np.array([mist[n].get(p, 0) for n in nomes])
        ax.bar(x, v, 0.6, bottom=base, color=COR[p], edgecolor="black", lw=0.4,
               label=ROT.get(p, p))
        for k, vv in enumerate(v):
            if vv > 0.07:
                ax.text(k, base[k] + vv / 2, f"{100*vv:.0f}", ha="center", va="center",
                        fontsize=6.3, color="white")
        base += v
    ax.set_xticks(x)
    cons = collections.Counter(l.split("\t")[1] for l in
                               (LAI / "consulta.tsv").read_text().split("\n") if l.strip())
    ax.set_xticklabels([f"{CURTO.get(n, n)}\n({2 * cons[n]} hap.)" for n in nomes], fontsize=6.8)
    ax.set_ylim(0, 1)
    ax.set_yticks([0, 0.5, 1])
    ax.set_yticklabels(["0", "50", "100"])
    ax.set_ylabel("% of the 10 Mb window", fontsize=7)
    ax.set_xlabel("African breeds of known taurine $\\times$ indicine admixture, run as queries",
                  fontsize=7)
    ax.legend(fontsize=6.3, frameon=False, ncol=3, loc="upper center",
              bbox_to_anchor=(0.5, 1.22), handlelength=1.1)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    letra(ax, "d", x=-0.11, y=1.24)

    # ---- (e) resolucao da janela e peso do CRF, por configuracao
    ax = fig.add_subplot(gs[1, 2])
    runs = ["m3", "m4", "m5"]
    nj = [val[r].get("n_janelas") or len(le_msp(r)[1]) for r in runs]
    kb = [10_000 / n for n in nj]
    ax.bar(range(3), kb, 0.6, color="#999999", edgecolor="black", lw=0.5)
    for k, (n, v) in enumerate(zip(nj, kb)):
        ax.text(k, v * 1.08, f"{v:.0f} kb\n{n} win.", ha="center", fontsize=6.2)
    ax.set_yscale("log")
    ax.set_xticks(range(3))
    ax.set_xticklabels(["3-way", "4-way", "5-way"], fontsize=6.8)
    ax.set_ylabel("mean window size (kb, log)", fontsize=7)
    ax.set_ylim(1, 400)
    ax.set_xlabel("reference configuration", fontsize=7)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    letra(ax, "e", x=-0.42, y=1.24)
    salvar(fig, "Figure_4")


# =====================================================================================
# FIGURA 3 - o classificador codificante repousa numa variante ligada do exon 2
# =====================================================================================
def fig3(res):
    d = _json(AQUI / "tables" / "derived" / "exon2_dados.json")
    est = _json(AQUI / "tables" / "derived" / "estabilidade_hierarquica.json")
    p1 = str(d["posicoes"][0])
    fig = plt.figure(figsize=(7.2, 6.1))
    gs = fig.add_gridspec(3, 2, height_ratios=[0.95, 1.0, 1.0], width_ratios=[1.1, 1],
                          hspace=0.95, wspace=0.36)

    # ---- (a) frequencia do alelo do exon 2 por painel, e o que os portadores levam
    ax = fig.add_subplot(gs[0, :])
    ordem = ["OGADEN", "ANKOLE", "BORAN", "NDAMA", "CRIOLLO", "EUR_TAURINE", "IBERIAN"]
    f = [d["frequencia"][p1][p]["freq_alt"] for p in ordem]
    n = [d["frequencia"][p1][p]["n_hap"] for p in ordem]
    x = np.arange(len(ordem))
    ax.bar(x, f, 0.62, color=[COR[p] for p in ordem], edgecolor="black", lw=0.5)
    for k, v in enumerate(f):
        ax.text(k, v + 0.04, f"{v:.2f}", ha="center", fontsize=6.6)
    ax.set_xticks(x)
    ax.set_xticklabels([CURTO.get(p, p) + chr(10) + f"{nn} hap." for p, nn in zip(ordem, n)],
                       fontsize=6.6)
    ax.set_ylim(0, 1.2)
    ax.set_yticks([0, 0.5, 1.0])
    ax.set_ylabel("frequency of the exon 2" + chr(10) + "alternative allele", fontsize=7)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    txt = "   ".join(f"{a} {d['alelos'][a]['alelo_exon2']['ALT']}/{d['alelos'][a]['n_hap']}"
                     for a in ALELOS_SEG)
    ax.text(0.5, 1.05, "carrier haplotypes carrying that allele:   " + txt,
            transform=ax.transAxes, ha="center", va="bottom", fontsize=6.7)
    letra(ax, "a", x=-0.075, y=1.20)

    # ---- (b) quantos sitios sustentam a atribuicao: contribuicao acumulada
    ax = fig.add_subplot(gs[1, 0])
    for a in ALELOS_SEG:
        s2 = d["alelos"][a]["configuracao_dominante"]["sitios"]
        v = np.array([abs(t["log10RV"]) for t in s2])
        cum = np.cumsum(np.sort(v)[::-1]) / v.sum()
        ax.step(np.arange(1, len(cum) + 1), cum, where="post", color=ALELO_COR[a], lw=1.4, label=a)
        sh = sum(abs(t["log10RV"]) for t in s2 if t["exon2"]) / v.sum()
        ax.plot(2, sh, "o", ms=4.5, mfc=ALELO_COR[a], mec="black", mew=0.5, zorder=5)
    ax.axvline(2, color="#444444", lw=0.6, ls=":")
    ax.annotate("the two exon 2 sites carry" + chr(10) + "59-85% of the evidence",
                xy=(2.25, 0.60), xytext=(4.2, 0.17), fontsize=6.3,
                arrowprops=dict(arrowstyle="->", lw=0.6, color="#444444"))
    ax.set_xscale("log")
    ax.set_xlim(1, 56)
    ax.set_xticks([1, 2, 5, 10, 20, 56])
    ax.set_xticklabels(["1", "2", "5", "10", "20", "56"], fontsize=6.5)
    ax.set_ylim(0, 1.04)
    ax.set_yticks([0, 0.5, 1.0])
    ax.set_yticklabels(["0", "50", "100"])
    ax.set_xlabel("coding sites, ranked by contribution (log scale)", fontsize=7)
    ax.set_ylabel("cumulative % of the total" + chr(10) + "likelihood-ratio evidence", fontsize=7)
    ax.legend(fontsize=6.0, frameon=False, ncol=4, loc="upper center",
              bbox_to_anchor=(0.5, 1.21), handlelength=1.0, columnspacing=0.8)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    letra(ax, "b", x=-0.28, y=1.21)

    # ---- (c) controle por remocao: a quem cada alelo e atribuido com e sem os dois sitios
    ax = fig.add_subplot(gs[1, 1])
    larg, x = 0.36, np.arange(len(ALELOS_SEG))
    vistos = set()
    for j2, k in enumerate(("com", "sem")):
        base = np.zeros(len(ALELOS_SEG))
        for pan in ("EUR_TAURINE", "IBERIAN", "NDAMA", "ANKOLE", "BORAN", "OGADEN", "CRIOLLO"):
            v = np.array([d["alelos"][a]["classes"]["so nao sinonimas"][k].get(pan, 0)
                          / d["alelos"][a]["n_hap"] for a in ALELOS_SEG])
            if v.sum() == 0:
                continue
            ax.bar(x + (j2 - 0.5) * larg, v, larg * 0.9, bottom=base, color=COR[pan],
                   edgecolor="black", lw=0.35, label=None if pan in vistos else ROT.get(pan, pan))
            vistos.add(pan)
            base += v
    for k2 in x:
        ax.text(k2 - 0.5 * larg, -0.04, "all", ha="center", va="top", fontsize=5.6, color="#555555")
        ax.text(k2 + 0.5 * larg, -0.04, "minus" + chr(10) + "exon 2", ha="center", va="top",
                fontsize=5.6, color="#555555")
    ax.set_xticks(x)
    ax.set_xticklabels(ALELOS_SEG, fontsize=6.8, fontweight="bold")
    ax.tick_params(axis="x", pad=18, length=0)
    ax.set_ylim(0, 1)
    ax.set_yticks([0, 0.5, 1])
    ax.set_yticklabels(["0", "50", "100"])
    ax.set_ylabel("% of carrier haplotypes", fontsize=7)
    ax.set_xlabel("non-synonymous sites used", fontsize=7, labelpad=2)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.legend(fontsize=5.8, frameon=False, ncol=3, loc="upper center",
              bbox_to_anchor=(0.45, 1.33), handlelength=0.9, columnspacing=0.7)
    letra(ax, "c", x=-0.19, y=1.33)

    # ---- (d) estabilidade por reamostragem de sitios, no nivel do continente
    ax = fig.add_subplot(gs[2, :])
    for k, a in enumerate(ALELOS_SEG):
        for j2, c in enumerate(("so sinonimas", "so nao sinonimas")):
            e = est["alelos"][a][c]
            obs = e["grupo_observado"]
            v = e["consenso_por_grupo"].get(obs, 0)
            ax.barh(k + (0.5 - j2) * 0.32, v, 0.28,
                    color=COR["EUR_TAURINE"] if obs == "europeu" else COR["NDAMA"],
                    alpha=1.0 if j2 else 0.5, edgecolor="black", lw=0.4)
            ax.text(v + 0.008, k + (0.5 - j2) * 0.32, f"{100 * v:.0f}%", va="center", fontsize=6.1)
    ax.axvline(0.95, color="#b2182b", lw=0.8, ls="--")
    ax.text(0.95, -0.85, "95%", color="#b2182b", fontsize=6.3, ha="center")
    ax.set_yticks(range(len(ALELOS_SEG)))
    ax.set_yticklabels(ALELOS_SEG, fontsize=6.9, fontweight="bold")
    ax.set_xlim(0, 1.06)
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
    ax.set_xticklabels(["0", "25", "50", "75", "100"])
    ax.set_xlabel("% of 2,000 site-resampling replicates keeping the observed continental group"
                  + chr(10) + "(pale bar: synonymous sites; solid bar: non-synonymous sites)",
                  fontsize=7)
    ax.invert_yaxis()
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    letra(ax, "d", x=-0.075, y=1.10)
    salvar(fig, "Figure_3")


# =====================================================================================
# FIGURA 4 - ancestralidade local e distancia de haplotipo sem modelo
# =====================================================================================
FOCO = {"SLICK1": 39099213, "SLICK3": 39099226, "SLICK6": 39099267, "SLICK2": 39099321}
SEM_SUP = {"SLICK2": {"ERR12035125", "ERR12035133"}, "SLICK3": {"ERR12035159"},
           "SLICK6": {"ERR12035135"}}


def _portadores():
    lin = (RFM / "focais_faseados.tsv").read_text().splitlines()
    am = lin[0].split("\t")[1:]
    gt = {int(l.split("\t")[0]): {s: g.split("|") for s, g in zip(am, l.split("\t")[1:])}
          for l in lin[1:]}
    cons = dict(l.split("\t")[:2] for l in (LAI / "consulta.tsv").read_text().split("\n") if l.strip())
    out = []
    for a in ALELOS_SEG:
        for s in am:
            for i in (0, 1):
                if gt[FOCO[a]][s][i] == "1" and cons.get(s) == "CRIOLLO" and s not in SEM_SUP.get(a, set()):
                    out.append((a, s, i))
    return out, cons


def fig5(res):
    port, cons = _portadores()
    viz = _json(RFM / "vizinho.json")
    fig = plt.figure(figsize=(7.2, 7.0))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.55, 1.0], width_ratios=[1, 0.95],
                          hspace=0.46, wspace=0.34)

    # ---- (a) cariograma na configuracao VALIDADA (3 vias), portadores e fundo nao portador
    ax = fig.add_subplot(gs[0, :])
    haps, jan = le_msp("m3")
    idx = {}
    for k, h in enumerate(haps):
        s, i = h.rsplit(".", 1)
        idx[(s, int(i))] = k
    nao_port = [(s, i) for (s, i) in idx
                if cons.get(s) == "CRIOLLO" and all(not (a == aa and s == ss and i == ii)
                                                    for (aa, ss, ii) in port for a in [aa])]
    usados = {(s, i) for (_, s, i) in port}
    fundo = sorted({(s, i) for (s, i) in idx if cons.get(s) == "CRIOLLO"} - usados)[:12]

    def faixa(k, y):
        ini, atual = jan[0]["spos"], jan[0]["call"][k]
        for j in jan[1:] + [None]:
            if j is None or j["call"][k] != atual:
                fim = jan[-1]["epos"] if j is None else j["spos"]
                ax.add_patch(Rectangle((ini / 1e6, y - 0.38), (fim - ini) / 1e6, 0.76,
                                       facecolor=COR.get(atual, "#cccccc"), lw=0))
                if j is not None:
                    ini, atual = j["spos"], j["call"][k]

    y, marcas = 0, []
    for a in ALELOS_SEG:
        sub = [(s, i) for (aa, s, i) in port if aa == a and (s, i) in idx]
        y0 = y
        for (s, i) in sub:
            faixa(idx[(s, i)], y)
            y += 1
        marcas.append((a, y0, y - 1, ALELO_COR[a]))
        y += 0.7
    y += 0.5
    y0 = y
    for (s, i) in fundo:
        faixa(idx[(s, i)], y)
        y += 1
    marcas.append(("Criollo" + chr(10) + "non-carriers", y0, y - 1, "#555555"))

    for nome, a0, a1, c in marcas:
        ax.text(33.95, (a0 + a1) / 2, nome, ha="right", va="center", fontsize=6.8,
                fontweight="bold", color=c, linespacing=1.1)
    ax.axvline(39.0993, color="black", lw=0.9, zorder=5)
    ax.annotate("slick alleles", xy=(39.0993, y - 0.5), xytext=(40.8, y + 0.6), fontsize=6.6,
                ha="left", va="center",
                arrowprops=dict(arrowstyle="->", lw=0.6, color="black"))
    ax.set_yticks([])
    ax.set_xlim(34.1, 44.1)
    ax.set_ylim(y + 1.2, -0.9)
    ax.set_xlabel("BTA20 position (Mb): local ancestry over the 10 Mb window, 3-way reference",
                  fontsize=7.5)
    for sp in ("left", "right", "top"):
        ax.spines[sp].set_visible(False)
    ax.legend(handles=[Patch(facecolor=COR[p], edgecolor="none", label=ROT.get(p, p))
                       for p in ("TAURINO_EUROPEU", "TAURINO_AFRICANO", "INDICINE")],
              fontsize=6.5, frameon=False, ncol=3, loc="upper center",
              bbox_to_anchor=(0.5, 1.07), handlelength=1.1)
    letra(ax, "a", x=-0.085, y=1.06)

    # ---- (b) distancia ao vizinho mais proximo de cada painel (+-50 kb)
    ax = fig.add_subplot(gs[1, 0])
    pans = ["CRIOLLO", "EUR_TAURINE", "IBERIAN", "NDAMA", "INDICINE"]
    w = "50000"
    for k, a in enumerate(ALELOS_SEG):
        for j, p in enumerate(pans):
            v = [h["por_painel"][p]["dist"] for h in viz[w][a] if p in h["por_painel"]]
            if not v:
                continue
            xx = k + (j - 2) * 0.165
            ax.plot([xx] * len(v), v, "o", ms=2.6, mfc=COR[p], mec="black", mew=0.25, alpha=0.85)
            ax.plot([xx - 0.06, xx + 0.06], [np.median(v)] * 2, color=COR[p], lw=1.8,
                    solid_capstyle="butt", zorder=3)
    ax.set_xticks(range(len(ALELOS_SEG)))
    ax.set_xticklabels(ALELOS_SEG, fontsize=6.8, fontweight="bold")
    ax.set_ylabel("differences to the closest haplotype of" + chr(10)
                  + "each panel (" + chr(177) + "50 kb, 766 SNPs)", fontsize=7)
    ax.set_ylim(0, None)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.legend(handles=[Line2D([], [], marker="o", ls="", ms=4, mfc=COR[p], mec="black",
                              mew=0.3, label=CURTO.get(p, p)) for p in pans],
              fontsize=5.9, frameon=False, ncol=5, loc="upper center",
              bbox_to_anchor=(0.5, 1.17), handlelength=0.6, columnspacing=0.6,
              handletextpad=0.35)
    letra(ax, "b", x=-0.22, y=1.17)

    # ---- (c) chamada no locus por configuracao de referencia
    ax = fig.add_subplot(gs[1, 1])
    runs = [("m3", "3-way"), ("m4", "4-way"), ("m5", "5-way")]
    larg = 0.26
    for j, (run, rot) in enumerate(runs):
        hs, jn = le_msp(run)
        i0 = next(i for i, x in enumerate(jn) if x["spos"] <= 39099267 <= x["epos"])
        ix = {}
        for k, h in enumerate(hs):
            s, i = h.rsplit(".", 1)
            ix[(s, int(i))] = k
        cont = []
        for a in ALELOS_SEG:
            cont.append(collections.Counter(jn[i0]["call"][ix[(s, i)]] for (aa, s, i) in port
                                            if aa == a and (s, i) in ix))
        base = np.zeros(len(ALELOS_SEG))
        for pan in ("TAURINO_EUROPEU", "EUR_TAURINE", "IBERIAN", "ASIAN_TAURINE",
                    "TAURINO_AFRICANO", "NDAMA", "INDICINE"):
            v = np.array([c.get(pan, 0) / max(sum(c.values()), 1) for c in cont])
            if v.sum() == 0:
                continue
            ax.bar(np.arange(len(ALELOS_SEG)) + (j - 1) * larg, v, larg * 0.9, bottom=base,
                   color=COR[pan], edgecolor="black", lw=0.3)
            base += v
        for k in range(len(ALELOS_SEG)):
            ax.text(k + (j - 1) * larg, -0.03, rot, ha="center", va="top", fontsize=5.1,
                    rotation=90, color="#555555")
    ax.set_xticks(range(len(ALELOS_SEG)))
    ax.set_xticklabels(ALELOS_SEG, fontsize=6.8, fontweight="bold")
    ax.tick_params(axis="x", pad=24, length=0)
    ax.set_ylim(0, 1)
    ax.set_yticks([0, 0.5, 1])
    ax.set_yticklabels(["0", "50", "100"])
    ax.set_ylabel("% of carrier haplotypes", fontsize=7)
    ax.set_xlabel("call at the locus, by reference configuration", fontsize=7, labelpad=2)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.legend(handles=[Patch(facecolor=COR[p], edgecolor="none", label=ROT.get(p, p))
                       for p in ("TAURINO_EUROPEU", "IBERIAN", "ASIAN_TAURINE", "NDAMA",
                                 "INDICINE")],
              fontsize=5.6, frameon=False, ncol=2, loc="upper center",
              bbox_to_anchor=(0.5, 1.21), handlelength=0.9, columnspacing=0.6)
    letra(ax, "c", x=-0.20, y=1.21)
    salvar(fig, "Figure_5")


# ============================================================================================
# Figura 2 e Figuras S1-S3 (2026-10-29). Dados: tables/derived/figuras_novas.json (dados_figuras_novas.py)
# e tables/derived/endogamia_SLICK{1,2}.json (inbreeding_control.py --json). Nada digitado aqui.
# ============================================================================================
ALELO_MARC = {"SLICK1": "^", "SLICK3": "o", "SLICK6": "D", "SLICK2": "s"}
NOVAS = AQUI / "tables" / "derived" / "figuras_novas.json"


def fig2(_res):
    """Evidencia nas leituras: suporte de cada chamada, fase lida no portador composto e
    portadores por raca."""
    d = _json(NOVAS)
    fig = plt.figure(figsize=(7.2, 3.15))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 1.45, 1.05], wspace=0.62, left=0.075,
                          right=0.985, bottom=0.30, top=0.90)
    # (a) suporte de leitura das 18 chamadas conferidas
    ax = fig.add_subplot(gs[0])
    prof_max = max(int(r["leituras_ref"]) + int(r["leituras_alt"]) for r in d["leituras"])
    xs = np.linspace(3, prof_max + 3, 200)
    ax.plot(xs, 3 / xs, ls=":", lw=0.8, color="#777777")
    ax.axhline(0.25, ls="--", lw=0.8, color="#777777")
    sobrepos = collections.Counter()
    for r in d["leituras"]:
        a = r["alelo"]
        ref, alt = int(r["leituras_ref"]), int(r["leituras_alt"])
        prof = ref + alt
        frac = alt / prof if prof else 0.0
        usado = r["usado_na_atribuicao"].strip() == "sim"
        ax.scatter(prof, frac, marker=ALELO_MARC[a], s=30, zorder=3, linewidths=1.0,
                   facecolors=ALELO_COR[a] if usado else "white", edgecolors=ALELO_COR[a])
        sobrepos[(prof, round(frac, 4))] += 1
    for (x_, y_), n_ in sobrepos.items():   # chamadas com o mesmo ponto: nenhuma escondida
        if n_ > 1:
            ax.text(x_ + 0.7, y_ - 0.035, f"\u00d7{n_}", fontsize=6.5, va="top", color="#333333")
    ax.set_xlim(-1, prof_max + 3)
    ax.set_ylim(-0.04, 1.06)
    ax.set_xlabel("Reads at the focal site")
    ax.set_ylabel("Fraction with the alternative base")
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    leg = [Line2D([], [], marker=ALELO_MARC[a], ls="", ms=5, mfc=ALELO_COR[a], mec=ALELO_COR[a],
                  label=a) for a in ("SLICK2", "SLICK3", "SLICK6")]
    leg += [Line2D([], [], marker="o", ls="", ms=5, mfc="white", mec="#555555", label="rejected"),
            Line2D([], [], ls="--", lw=0.8, color="#777777", label="fraction 0.25"),
            Line2D([], [], ls=":", lw=0.8, color="#777777", label="3 alt. reads")]
    ax.legend(handles=leg, fontsize=6.3, frameon=False, loc="upper center",
              bbox_to_anchor=(0.45, -0.22), ncol=2, handletextpad=0.2, columnspacing=0.8)
    letra(ax, "a", x=-0.30, y=1.04)

    # (b) leituras do portador composto ERR12035149 sobre SLICK6 e SLICK2
    ax = fig.add_subplot(gs[1])
    p6, p2 = 39099267, 39099321
    reads = sorted(d["fase"], key=lambda r: (int(r["inicio"]), int(r["fim"])))
    for k, r in enumerate(reads):
        ini, fim = int(r["inicio"]), int(r["fim"])
        b6, b2 = r["base_SLICK6"], r["base_SLICK2"]
        ambas = bool(b6) and bool(b2)
        ax.plot([ini, fim], [k, k], lw=3.4, solid_capstyle="butt",
                color="#a9a9a9" if ambas else "#e3e3e3", zorder=1)
        for pos, b, al in ((p6, b6, "SLICK6"), (p2, b2, "SLICK2")):
            if not b:
                continue
            if b == "T":
                ax.plot(pos, k, marker=ALELO_MARC[al], ms=6.5, mfc=ALELO_COR[al], mec="black",
                        mew=0.5, zorder=3)
            else:
                ax.plot(pos, k, marker="|", ms=7, mew=1.2, color="#4d4d4d", zorder=3)
    for pos, al in ((p6, "SLICK6"), (p2, "SLICK2")):
        ax.axvline(pos, ls=":", lw=0.7, color=ALELO_COR[al], zorder=0)
        ax.text(pos + (-3 if al == "SLICK6" else 3), len(reads) - 0.3, al,
                ha="right" if al == "SLICK6" else "left", va="bottom", fontsize=7,
                color=ALELO_COR[al], fontweight="bold")
    yb = -1.1
    ax.annotate("", xy=(p6, yb), xytext=(p2, yb), arrowprops=dict(arrowstyle="<->", lw=0.7))
    ax.text((p6 + p2) / 2, yb - 0.3, "54 bp", ha="center", va="top", fontsize=6.5)
    ax.set_ylim(yb - 1.1, len(reads) + 0.7)
    ax.set_yticks([])
    lo = min(int(r["inicio"]) for r in reads) - 5
    hi = max(int(r["fim"]) for r in reads) + 5
    ax.set_xlim(lo, hi)
    tk = list(range((lo // 100 + 1) * 100, hi, 100))
    ax.set_xticks(tk)
    ax.set_xticklabels([f"{x:,}" for x in tk], fontsize=6.3)
    ax.set_xlabel("Position on BTA20 (bp)")
    ax.set_ylabel("Reads, ERR12035149")
    for sp in ("top", "right", "left"):
        ax.spines[sp].set_visible(False)
    leg = [Line2D([], [], marker="|", ls="", ms=7, mew=1.2, color="#4d4d4d", label="reference base"),
           Line2D([], [], color="#a9a9a9", lw=3.4, label="spans both sites"),
           Line2D([], [], color="#e3e3e3", lw=3.4, label="spans one site")]
    ax.legend(handles=leg, fontsize=6.3, frameon=False, loc="upper center",
              bbox_to_anchor=(0.5, -0.22), ncol=2, handletextpad=0.3, columnspacing=0.8)
    letra(ax, "b", x=-0.08, y=1.04)

    # (c) portadores por raca (animais) e animais triados
    ax = fig.add_subplot(gs[2])
    tri = d["racas"]["triados"]
    por = d["racas"]["portadores"]
    racas = sorted(tri, key=lambda r: (-sum(por[a].get(r, 0) for a in por), -tri[r]))
    ordem = ["SLICK2", "SLICK1", "SLICK6", "SLICK3"]
    for i, rc in enumerate(racas):
        esq = 0
        for a in ordem:
            n = por.get(a, {}).get(rc, 0)
            if n:
                ax.barh(i, n, left=esq, color=ALELO_COR[a], edgecolor="white", lw=0.5, height=0.62)
                esq += n
        ax.text(esq + 0.4, i, f"of {tri[rc]}", va="center", fontsize=6.3, color="#444444")
    ax.set_yticks(range(len(racas)))
    ax.set_yticklabels(racas, fontsize=6.8)
    ax.invert_yaxis()
    ax.set_xlabel("Carrier animals, by allele")
    tot = max(sum(por[a].get(r, 0) for a in por) for r in racas)
    ax.set_xlim(0, tot + 7)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    from matplotlib.patches import Patch
    leg = [Patch(facecolor=ALELO_COR[a], label=a) for a in ordem]
    ax.legend(handles=leg, fontsize=6.3, frameon=False, loc="upper center",
              bbox_to_anchor=(0.4, -0.22), ncol=2, handletextpad=0.3, columnspacing=0.8)
    letra(ax, "c", x=-0.62, y=1.04)
    salvar(fig, "Figure_2")


def figS1(_res):
    """Distribuicao da troca do SLICK5 (39.099.228) por populacao."""
    d = _json(NOVAS)["slick5"]
    an = d["animais"]
    paineis = ("European taurine", "Criollo", "N'Dama", "Ogaden")
    itens = sorted(an.items(), key=lambda kv: -kv[1])
    fica = [(k, v) for k, v in itens if v >= 2 or k in paineis]
    resto = [(k, v) for k, v in itens if not (v >= 2 or k in paineis)]
    if resto:
        fica.append((f"other declared breeds ({len(resto)})", sum(v for _, v in resto)))
    cor_de = {"European taurine": COR["EUR_TAURINE"], "Criollo": COR["CRIOLLO"],
              "N'Dama": COR["NDAMA"], "Ogaden": COR["OGADEN"]}
    fig, ax = plt.subplots(figsize=(5.2, 3.6))
    for i, (k, v) in enumerate(fica):
        c = cor_de.get(k, "#d9d9d9" if k == "no declared breed" else "#9ecae1")
        ax.barh(i, v, color=c, edgecolor="black", lw=0.4, height=0.7)
        ax.text(v + 0.4, i, str(v), va="center", fontsize=7)
    ax.set_yticks(range(len(fica)))
    ax.set_yticklabels([k for k, _ in fica], fontsize=7)
    ax.invert_yaxis()
    ax.set_xlabel(f"Animals carrying the SLICK5 change ({d['n_animais']} animals, "
                  f"{d['n_haplotipos']} haplotypes)")
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    leg = [Line2D([], [], color=c, lw=6, label=l) for l, c in
           (("reference panel (panel colour)", COR["EUR_TAURINE"]), ("declared breed, no panel", "#9ecae1"),
            ("no declared breed", "#d9d9d9"))]
    ax.legend(handles=leg, fontsize=6.5, frameon=False, loc="lower right")
    fig.tight_layout()
    salvar(fig, "Figure_S1")


def figS2(_res):
    """Atribuicao intronica por faixa de distancia ao sitio causal (mesmo criterio em todas)."""
    d = _json(NOVAS)["introns"]
    alelos = ["SLICK1", "SLICK3", "SLICK6", "SLICK2"]
    faixas = [b["faixa"] for b in d["SLICK1"]]
    curto = {"EUR_TAURINE": "Cont. European", "IBERIAN": "Iberian", "NDAMA": "N'Dama",
             "CRIOLLO": "Criollo", "ANKOLE": "Ankole", "BORAN": "Boran", "OGADEN": "Ogaden"}
    fig, ax = plt.subplots(figsize=(6.6, 2.9))
    for i, a in enumerate(alelos):
        for j, b in enumerate(d[a]):
            if "contagem" not in b:
                ax.add_patch(plt.Rectangle((j, i), 1, 1, facecolor="white", edgecolor="#bbbbbb"))
                ax.text(j + 0.5, i + 0.5, "too few\nsites", ha="center", va="center", fontsize=6)
                continue
            mx = max(b["contagem"].values())
            tops = [k for k, v in b["contagem"].items() if v == mx]
            tot = sum(b["contagem"].values())
            adj = b["posterior_mediana"] >= 0.50
            cor = COR.get(tops[0], "#999999") if len(tops) == 1 else "#bdbdbd"
            ax.add_patch(plt.Rectangle((j, i), 1, 1, facecolor=cor, alpha=0.95 if adj else 0.25,
                                       edgecolor="white", lw=1.2, hatch=None if adj else "////"))
            rot = " / ".join(curto.get(k, k) for k in tops)
            claro = adj and cor not in (COR["CRIOLLO"], "#bdbdbd")
            ax.text(j + 0.5, i + 0.5, f"{rot}\n{mx}/{tot}  ({b['posterior_mediana']:.2f})",
                    ha="center", va="center", fontsize=5.8, color="white" if claro else "black")
    ax.set_xlim(0, len(faixas))
    ax.set_ylim(len(alelos), 0)
    ax.set_xticks([j + 0.5 for j in range(len(faixas))])
    ax.set_xticklabels(faixas, fontsize=7)
    ax.set_yticks([i + 0.5 for i in range(len(alelos))])
    ax.set_yticklabels(alelos, fontsize=7.5)
    ax.set_xlabel("Distance of the intronic sites from the causal site")
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.tick_params(length=0)
    fig.tight_layout()
    salvar(fig, "Figure_S2")


def figS3(_res):
    """Controle de endogamia: trato no locus contra posicoes aleatorias da mesma janela."""
    fig, axs = plt.subplots(1, 2, figsize=(6.6, 2.8))
    for ax, a, le in zip(axs, ("SLICK1", "SLICK2"), "ab"):
        e = _json(AQUI / "tables" / "derived" / f"endogamia_{a}.json")
        base = np.array(e["linha_de_base_kb"])
        foco = e["trato_foco_mediana_kb"]
        lo_, hi_ = min(base.min(), foco), max(base.max(), foco)
        log_ = hi_ / max(lo_, 0.5) > 5            # faixa larga: escala log; estreita: linear
        if log_:
            bins = np.logspace(np.log10(max(lo_, 0.5)), np.log10(hi_ * 1.15), 30)
        else:
            m = (hi_ - lo_) * 0.06
            bins = np.linspace(lo_ - m, hi_ + m, 30)
        cont, _, _ = ax.hist(base, bins=bins, color="#bdbdbd", edgecolor="white", lw=0.4)
        ax.axvline(e["linha_de_base_mediana_kb"], ls="--", lw=0.8, color="#555555")
        ax.axvline(foco, lw=1.6, color=ALELO_COR[a])
        if log_:
            ax.set_xscale("log")
        ax.set_ylim(0, cont.max() * 1.6)
        ax.set_xlabel("Shared tract around the position (kb)")
        ax.set_ylabel(f"Random positions (of {e['n_random']})")
        nh = len(e["homozigotos"])
        ax.text(0.98, 0.97, f"{a}, {nh} homozygote{'s' if nh > 1 else ''}\n"
                f"at the locus {foco:,.1f} kb\nbaseline {e['linha_de_base_mediana_kb']:,.1f} kb\n"
                f"P = {e['p_empirico']:.3f}", transform=ax.transAxes, ha="right", va="top",
                fontsize=6.4, bbox=dict(facecolor="white", edgecolor="none", alpha=0.85, pad=1.5))
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        letra(ax, le, x=-0.2, y=1.04)
    fig.tight_layout()
    salvar(fig, "Figure_S3")


def main(quais):
    res = json.load(open(RES, encoding="utf-8"))
    for n in quais:                      # 1-5 (texto principal) e S1-S3 (suplementar)
        globals()[f"fig{n}"](res)


if __name__ == "__main__":
    q = sys.argv[1:] or ["1", "2", "3", "4", "5", "S1", "S2", "S3"]
    main(q)
