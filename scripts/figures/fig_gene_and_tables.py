#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Figuras e tabelas do eixo de descoberta do locus (PRLR / locus SLICK) a
partir dos resultados REAIS.

Idempotente e INCREMENTAL: cada figura verifica se sua entrada existe; se os
dados ainda nao foram produzidos, ela e pulada com aviso e o build segue.
Nada aqui inventa numero: tudo vem de results/*.

Uso:  python make_figures.py [--proj RAIZ]
"""
import argparse, json, sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd

PRLR = dict(chrom=20, start=39037656, end=39117082, cds_end=39099575,
            stop_codon=39099576)
# Dois alelos truncantes, nao um. O SLICK2 e um stop_gained (C>T) e por isso
# ficou invisivel enquanto a classificacao era por diferenca de comprimento.
TRUNCATING = [
    dict(pos=39099213, nome="SLICK1", hgvs="c.1382del", cor="#C1442E"),
    dict(pos=39099321, nome="SLICK2", hgvs="c.1489C>T", cor="#8C3D8C"),
]
FRAMESHIFT_POS = 39099213

GROUP_ORDER = ["TROP_TAURINE", "AFR_TAURINE", "AFR_SANGA_ZEBU",
               "IBERIAN", "TEMP_TAURINE", "INDICINE"]
GROUP_LABEL = {
    "TROP_TAURINE": "Tropical taurine\n(Criollo)",
    "AFR_TAURINE": "African taurine\n(N'Dama)",
    "AFR_SANGA_ZEBU": "African sanga/zebu",
    "IBERIAN": "Iberian",
    "TEMP_TAURINE": "Temperate taurine",
    "INDICINE": "Indicine",
}
GROUP_COLOR = {
    "TROP_TAURINE": "#C1442E", "AFR_TAURINE": "#E08A1E",
    "AFR_SANGA_ZEBU": "#B8912F", "IBERIAN": "#5B8C5A",
    "TEMP_TAURINE": "#3D6B99", "INDICINE": "#6A579B",
}

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 8,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.linewidth": 0.7, "xtick.major.width": 0.7,
    "ytick.major.width": 0.7, "savefig.dpi": 400,
    "figure.dpi": 150,
})

LOG = []
def ok(m):   LOG.append(("OK", m));   print(f"  [OK]   {m}")
def skip(m): LOG.append(("SKIP", m)); print(f"  [SKIP] {m}")


def read_tsv(p, **kw):
    p = Path(p)
    if not p.exists() or p.stat().st_size == 0:
        return None
    try:
        return pd.read_csv(p, sep="\t", **kw)
    except Exception as e:
        print(f"  [WARN] falha lendo {p.name}: {e}")
        return None


# ---------------------------------------------------------------- Figura 1
def fig1_prlr(proj, figs):
    """Modelo do PRLR + frameshift + frequencia por grupo ancestral."""
    cls = read_tsv(proj / "results/axis1/prlr_cds_variants_classified.tsv")
    afb = read_tsv(proj / "results/axis1/truncating_af_by_breed.tsv")
    grp = read_tsv(proj / "data/cohort/cohort_groups.tsv",
                   names=["sample", "breed", "group"])
    if cls is None or afb is None:
        return skip("Fig1: faltam resultados do Eixo 1")

    cds = read_tsv(proj / "results/axis1/prlr_cds.bed",
                   names=["c", "s", "e"])

    fig = plt.figure(figsize=(7.2, 5.0))
    gs = fig.add_gridspec(3, 1, height_ratios=[1.0, 1.0, 1.9], hspace=0.55)

    # (a) modelo do gene
    ax = fig.add_subplot(gs[0])
    ax.plot([PRLR["start"], PRLR["end"]], [0, 0], color="#999", lw=1.0, zorder=1)
    if cds is not None:
        for _, r in cds.iterrows():
            ax.add_patch(Rectangle((r.s, -0.22), r.e - r.s, 0.44,
                                   facecolor="#3D6B99", edgecolor="none", zorder=2))
    for i, t in enumerate(TRUNCATING):
        ax.axvline(t["pos"], color=t["cor"], lw=1.1, zorder=3)
    # Os dois distam 108 pb: na escala do gene caem no mesmo ponto, entao uma
    # chamada unica com as duas linhas evita falsa resolucao.
    ax.annotate("SLICK1  c.1382del  (frameshift)\nSLICK2  c.1489C>T  (stop gained)",
                xy=(TRUNCATING[0]["pos"], 0.30),
                xytext=(TRUNCATING[0]["pos"] - 30000, 0.74),
                fontsize=6.8, color="#333", ha="center", linespacing=1.5,
                arrowprops=dict(arrowstyle="->", color="#666", lw=0.8))
    ax.annotate("stop", xy=(PRLR["stop_codon"], -0.30), fontsize=6.5,
                color="#444", ha="center", va="top")
    ax.set_xlim(PRLR["start"] - 3000, PRLR["end"] + 3000)
    ax.set_ylim(-0.85, 1.05); ax.set_yticks([])
    ax.set_xlabel(f"BTA{PRLR['chrom']} position (bp)", fontsize=7.5)
    ax.set_title("a   PRLR gene model and the two truncating variants",
                 loc="left", fontsize=9, fontweight="bold")
    ax.ticklabel_format(axis="x", style="plain")
    ax.xaxis.set_major_formatter(
        matplotlib.ticker.FuncFormatter(lambda v, _: f"{int(v):,}"))
    ax.tick_params(axis="x", labelsize=6.5)

    # (b) consequencia funcional no CDS, pela anotacao do bcftools csq.
    # NAO usar a classificacao por diferenca de comprimento: ela nao enxerga
    # stop_gained (uma troca C>T nao muda o comprimento) e foi por isso que o
    # SLICK2 passou despercebido no primeiro rastreio.
    ax = fig.add_subplot(gs[1])
    csq = read_tsv(proj / "results/axis1/csq/csq_all.tsv",
                   names=["POS", "REF", "ALT", "AC", "AN", "FILTER", "MQ", "BCSQ"],
                   header=None)
    ORD = [("synonymous", "Synonymous", "#B9C6D4"),
           ("missense", "Missense", "#8FA9C4"),
           ("intron", "Intron", "#D3D9E2"),
           ("splice_region", "Splice region", "#7E93AC"),
           ("inframe", "In-frame indel", "#5B7fa5"),
           ("frameshift", "Frameshift", "#C1442E"),
           ("stop_gained", "Stop gained", "#8C3D8C")]
    if csq is not None:
        prim = (csq["BCSQ"].astype(str).str.split("|").str[0]
                .str.replace(r"^@\d+$", "continuation", regex=True))
        vals, labs, cols = [], [], []
        for key, lab, c in ORD:
            n = int(prim.str.contains(key, na=False).sum())
            if n:
                vals.append(n); labs.append(lab); cols.append(c)
        bars = ax.barh(labs, vals, color=cols, height=0.6)
        for b, v in zip(bars, vals):
            ax.text(b.get_width() + max(vals) * 0.015,
                    b.get_y() + b.get_height() / 2, str(v),
                    va="center", fontsize=7.5)
        ax.set_xlim(0, max(vals) * 1.18)
    ax.set_xlabel("Variants in PRLR coding sequence", fontsize=7.5)
    ax.set_title("b   Functional consequence of all coding variants (1,842 animals)",
                 loc="left", fontsize=9, fontweight="bold")

    # (c) frequencia dos DOIS alelos truncantes, por raca
    ax = fig.add_subplot(gs[2])
    if grp is not None and not afb.empty:
        b2g = dict(zip(grp["breed"], grp["group"]))
        afb = afb.copy()
        afb["AF"] = pd.to_numeric(afb["AF"], errors="coerce").fillna(0.0)
        afb["group"] = afb["BREED"].map(b2g)
        afb = afb[afb["group"].notna()]
        ordem = {g: i for i, g in enumerate(GROUP_ORDER)}
        # racas com n=1..3 nao sustentam frequencia e poluem o eixo; ficam de fora,
        # exceto se carregarem algum dos alelos (portador nunca e escondido)
        port = set(afb.loc[afb["AF"] > 0, "BREED"])
        breeds = (afb[["BREED", "N", "group"]].drop_duplicates()
                  .assign(_o=lambda d: d["group"].map(ordem))
                  .sort_values(["_o", "BREED"]))
        breeds = breeds[(breeds["N"] >= 4) | (breeds["BREED"].isin(port))]
        xs, labs, ticks = [], [], []
        x = 0.0; ultimo = None
        W = 0.38
        for _, b in breeds.iterrows():
            if ultimo is not None and b["group"] != ultimo:
                x += 0.9
            for k, t in enumerate(TRUNCATING):
                sel = afb[(afb["BREED"] == b["BREED"]) & (afb["POS"] == t["pos"])]
                af = float(sel["AF"].iloc[0]) if len(sel) else 0.0
                ax.bar(x + (k - 0.5) * W, af, width=W, color=t["cor"],
                       edgecolor="none")
            ticks.append(x)
            labs.append(f"{b['BREED'].replace('_',' ')} (n={int(b['N'])})")
            ultimo = b["group"]; x += 1.0
        ax.set_xticks(ticks)
        ax.set_xticklabels(labs, rotation=45, ha="right", fontsize=7.4)
        ax.set_ylabel("Allele frequency", fontsize=7.5)
        ax.set_ylim(0, 0.85)
        ax.set_xlim(-1.2, x)
        h = [plt.Rectangle((0, 0), 1, 1, color=t["cor"]) for t in TRUNCATING]
        ax.legend(h, [f"{t['nome']}  {t['hgvs']}" for t in TRUNCATING],
                  fontsize=6.8, frameon=False, ncol=2, loc="upper right")
    ax.set_title("c   Frequency of both truncating alleles, by breed",
                 loc="left", fontsize=9, fontweight="bold")

    out = figs / "Fig1_PRLR_alleles.png"
    fig.savefig(out, bbox_inches="tight"); plt.close(fig)
    fig_pdf = figs / "Fig1_PRLR_alleles.pdf"
    fig.savefig(fig_pdf, bbox_inches="tight")
    ok(f"Fig1 -> {out.name}")
    return True


# ---------------------------------------------------------------- Figura 2
def fig2_bta4(proj, figs):
    """Regiao BTA4 (janela de introgressao Zayas et al. 2025): dAF tropical vs
    temperado + trilha de genes."""
    # NOTA: no axis4 o `sort` roda DEPOIS do awk que imprime o cabecalho, entao a
    # linha de cabecalho foi ordenada para dentro do corpo. Ler sem header e
    # descartar linhas nao numericas resolve sem precisar reprocessar o VCF.
    daf = read_tsv(proj / "results/axis4/bta4_dAF_trop_vs_temp.tsv",
                   names=["SITE", "AF_TROP", "AF_TEMP", "dAF"], header=None)
    gen = read_tsv(proj / "results/axis4/bta4_genes_all.tsv",
                   names=["chrom", "start", "end", "strand", "name",
                          "biotype", "gene_id"])
    if daf is None or gen is None:
        return skip("Fig2: faltam resultados do Eixo 4")

    daf = daf.copy()
    daf["dAF"] = pd.to_numeric(daf["dAF"], errors="coerce")
    daf["POS"] = pd.to_numeric(
        daf["SITE"].astype(str).str.split(":").str[-1], errors="coerce")
    daf = daf.dropna(subset=["dAF", "POS"])
    daf["POS"] = daf["POS"].astype(int)
    if daf.empty:
        return skip("Fig2: tabela de dAF vazia")

    pc = gen[gen["biotype"] == "protein_coding"].copy()
    pc = pc[pc["name"].notna() & (pc["name"] != ".")]

    fig, (ax, axg) = plt.subplots(
        2, 1, figsize=(7.2, 3.9), sharex=True,
        gridspec_kw=dict(height_ratios=[3.1, 1.0], hspace=0.12))

    ax.scatter(daf["POS"] / 1e6, daf["dAF"], s=3.5, c="#8FA9C4",
               edgecolors="none", rasterized=True)
    thr = daf["dAF"].quantile(0.999)
    top = daf[daf["dAF"] >= thr]
    ax.scatter(top["POS"] / 1e6, top["dAF"], s=11, c="#C1442E",
               edgecolors="none", zorder=3, label=f"top 0.1% (n={len(top)})")
    ax.set_ylabel(r"$|\Delta$AF$|$  tropical vs temperate taurine", fontsize=7.5)
    ax.legend(fontsize=6.5, frameon=False, loc="upper left")
    ax.set_title(
        f"Zayas et al. (2025) BTA4 introgression window "
        f"({len(daf):,} SNPs, {len(pc)} named protein-coding genes)",
        loc="left", fontsize=9, fontweight="bold")

    for i, (_, g) in enumerate(pc.iterrows()):
        y = 0.62 - (i % 3) * 0.30
        axg.add_patch(Rectangle((g.start / 1e6, y - 0.055),
                                max((g.end - g.start) / 1e6, 0.006), 0.11,
                                facecolor="#3D6B99", edgecolor="none"))
        axg.text((g.start + g.end) / 2e6, y + 0.10, g["name"], fontsize=5.8,
                 ha="center", va="bottom", style="italic")
    axg.set_ylim(-0.12, 1.05); axg.set_yticks([])
    axg.set_xlabel("BTA4 position (Mb)", fontsize=7.5)
    axg.spines["left"].set_visible(False)

    out = figs / "Fig2_BTA4_Zayas_region.png"
    fig.savefig(out, bbox_inches="tight"); plt.close(fig)
    ok(f"Fig2 -> {out.name}")
    return True


# ---------------------------------------------------------------- Tabelas
def tables(proj, tabs):
    n = 0
    cls = read_tsv(proj / "results/axis1/prlr_cds_variants_classified.tsv")
    if cls is not None:
        sub = cls[cls["CLASSE"] != "substituicao"].copy()
        sub.to_csv(tabs / "Table1_PRLR_indels.tsv", sep="\t", index=False)
        ok(f"Table1 ({len(sub)} indels no CDS)"); n += 1
    car = read_tsv(proj / "results/axis1/indel_carriers.tsv")
    if car is not None:
        car.to_csv(tabs / "Table2_carriers.tsv", sep="\t", index=False)
        ok(f"Table2 ({len(car)} portadores)"); n += 1
    gen = read_tsv(proj / "results/axis4/bta4_genes_all.tsv",
                   names=["chrom", "start", "end", "strand", "name",
                          "biotype", "gene_id"])
    if gen is not None:
        pc = gen[gen["biotype"] == "protein_coding"]
        pc.to_csv(tabs / "Table3_BTA4_genes.tsv", sep="\t", index=False)
        ok(f"Table3 ({len(pc)} genes codificantes BTA4)"); n += 1
    grp = read_tsv(proj / "data/cohort/cohort_groups.tsv",
                   names=["sample", "breed", "group"])
    if grp is not None:
        comp = (grp.groupby(["group", "breed"]).size()
                .reset_index(name="n").sort_values(["group", "n"],
                                                   ascending=[True, False]))
        comp.to_csv(tabs / "TableS1_cohort.tsv", sep="\t", index=False)
        ok(f"TableS1 (coorte: {len(grp)} animais, {grp['breed'].nunique()} racas)")
        n += 1
    return n


# ---------------------------------------------------------------- pendentes
def pending(proj):
    """Relata o que ainda falta, para o texto final saber o que nao afirmar."""
    checks = {
        "Eixo 2/3 (coorte genome-wide)": proj / "data/cohort/chr1.cohort.vcf.gz",
        "Senepol (RNA-seq PRLR)": proj / "results/senepol/senepol_prlr.vcf",
        "WGS criollo (fatias)": proj / "data/wgs_slices/criollo",
        "WGS taurino africano": proj / "data/wgs_slices/afr_taurine",
        "Metilacao GSE101796": proj / "data/GSE101796/GSE101796_RAW.tar",
    }
    st = {}
    for k, p in checks.items():
        st[k] = p.exists() and (p.is_dir() and any(p.iterdir()) if p.is_dir()
                                else p.stat().st_size > 0)
    return st


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--proj", default=str(Path(__file__).resolve().parents[1]))
    a = ap.parse_args()
    proj = Path(a.proj)
    figs = proj / "PAPER/figuras"; figs.mkdir(parents=True, exist_ok=True)
    tabs = proj / "PAPER/tabelas"; tabs.mkdir(parents=True, exist_ok=True)

    print(f"== make_figures :: {proj}")
    fig1_prlr(proj, figs)
    fig2_bta4(proj, figs)
    ntab = tables(proj, tabs)

    st = pending(proj)
    print("\n== estado das entradas ==")
    for k, v in st.items():
        print(f"  {'PRONTO ' if v else 'PENDENTE'}  {k}")

    (proj / "PAPER/_build_status.json").write_text(json.dumps({
        "figuras": sorted(p.name for p in figs.glob("*.png")),
        "tabelas": sorted(p.name for p in tabs.glob("*.tsv")),
        "entradas": st,
        "log": LOG,
    }, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\n  {len(list(figs.glob('*.png')))} figuras, {ntab} tabelas")
    return 0


if __name__ == "__main__":
    sys.exit(main())
