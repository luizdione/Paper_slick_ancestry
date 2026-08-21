#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""panel_sanity_run9.py -- sanidade de painel (regra #1) para os rotulos Run9 v2, na JANELA do PRLR.

IMPORTANTE (regra #4): run9_chr20_prlr.vcf.gz e uma janela UNICA de ~0,4 Mb no locus PRLR, sob o
sweep do SLICK. Portanto isto NAO e recuperacao de raca genome-wide -- e um teste de COERENCIA:
os grupos coarse (african/european/indicine/east_asian/iberian) se separam mesmo aqui? indicino x
taurino deve separar forte; africano x europeu e mais sutil e separacao imperfeita NAO reprova.
Serve para pegar ERRO GROSSO de rotulo (ex.: 'Holstein' que cai no meio dos N'Dama).

Faz:
  1. matriz de dosagem 0/1/2 (1038 x nvar) do VCF da janela
  2. filtro MAF, imputacao por media, padroniza, PCA (top 10 PCs) -- so amostras REFERENCIA com grupo resolvido
  3. PC1xPC2 colorido por grupo coarse -> PNG
  4. leave-one-out nearest-centroid em espaco de PC -> matriz de confusao + % que recupera o proprio rotulo

Saidas: panel_sanity_run9_pca.png, panel_sanity_run9_report.txt
"""
import gzip, sys, time
from pathlib import Path
import numpy as np

BASE = Path(__file__).resolve().parent
VCF  = BASE / "run9_chr20_prlr.vcf.gz"
LAB  = BASE / "run9_labels_v2.tsv"
PNG  = BASE / "panel_sanity_run9_pca.png"
REP  = BASE / "panel_sanity_run9_report.txt"
MAF_MIN = 0.02
TOPPC = 10

def log(m, fh):
    print(m, flush=True); fh.write(m + "\n")

def gt_dosage(field):
    """'0|1:...' -> 1 ; './.' -> nan. Conta o alelo ALT."""
    g = field.split(":", 1)[0]
    if not g or g[0] == ".":
        return np.nan
    a = g.replace("|", "/").split("/")
    try:
        return float(sum(1 for x in a if x not in (".", "0")) if len(a) > 1 else (a[0] not in (".", "0")))
    except Exception:
        return np.nan

def main():
    fh = open(REP, "w", encoding="utf-8")
    log(f"[{time.strftime('%F %T')}] sanidade de painel Run9 (janela PRLR)", fh)

    # labels: sample -> coarse (so referencia candidata resolvida)
    coarse = {}
    for i, ln in enumerate(open(LAB, encoding="utf-8")):
        if i == 0: continue
        p = ln.rstrip("\n").split("\t")
        sid, cg, src = p[0], p[3], p[5]
        if cg in ("unresolved", "exclude"): continue
        coarse[sid] = cg

    # parse VCF -> dosage matrix
    samples = None
    rows = []; keep_idx = None
    with gzip.open(VCF, "rt") as f:
        for ln in f:
            if ln.startswith("##"): continue
            if ln.startswith("#CHROM"):
                cols = ln.rstrip("\n").split("\t")
                samples = cols[9:]
                keep_idx = [i for i, s in enumerate(samples) if s in coarse]
                kept_samples = [samples[i] for i in keep_idx]
                log(f"amostras no VCF: {len(samples)}; com grupo resolvido (referencia): {len(kept_samples)}", fh)
                continue
            p = ln.rstrip("\n").split("\t")
            gts = p[9:]
            row = np.array([gt_dosage(gts[i]) for i in keep_idx], dtype=float)
            # MAF
            m = np.nanmean(row) / 2.0
            if np.isnan(m) or m < MAF_MIN or m > 1 - MAF_MIN:
                continue
            rows.append(row)
    X = np.array(rows).T  # (n_samples, n_var)
    log(f"variantes apos MAF>= {MAF_MIN}: {X.shape[1]}; matriz {X.shape}", fh)

    # imputa por media da coluna, padroniza
    col_mean = np.nanmean(X, axis=0)
    inds = np.where(np.isnan(X))
    X[inds] = np.take(col_mean, inds[1])
    X = X - X.mean(0)
    sd = X.std(0); sd[sd == 0] = 1.0
    X = X / sd

    # PCA por SVD
    U, S, Vt = np.linalg.svd(X, full_matrices=False)
    PC = U[:, :TOPPC] * S[:TOPPC]
    ev = (S**2) / (S**2).sum()
    log(f"variancia explicada PC1..PC5: {', '.join(f'{100*e:.1f}%' for e in ev[:5])}", fh)

    y = np.array([coarse[s] for s in kept_samples])
    groups = sorted(set(y))

    # leave-one-out nearest-centroid em espaco de PC (regra #1: recupera o proprio rotulo?)
    log("\n== LEAVE-ONE-OUT nearest-centroid (regra #1) ==", fh)
    from collections import defaultdict
    conf = defaultdict(lambda: defaultdict(int))
    correct = defaultdict(int); total = defaultdict(int)
    P = PC.copy()
    for i in range(len(y)):
        best = None; bestd = None
        for g in groups:
            mask = (y == g); mask[i] = False
            if mask.sum() == 0: continue
            c = P[mask].mean(0)
            d = ((P[i] - c) ** 2).sum()
            if bestd is None or d < bestd:
                bestd = d; best = g
        conf[y[i]][best] += 1
        total[y[i]] += 1
        if best == y[i]: correct[y[i]] += 1
    log(f"{'grupo':14} {'n':>5} {'recupera_proprio':>17}", fh)
    overall_c = overall_n = 0
    for g in groups:
        overall_c += correct[g]; overall_n += total[g]
        log(f"{g:14} {total[g]:5d} {100*correct[g]/max(1,total[g]):16.1f}%", fh)
    log(f"{'TOTAL':14} {overall_n:5d} {100*overall_c/max(1,overall_n):16.1f}%", fh)
    log("\nmatriz de confusao (linha=rotulo verdadeiro, coluna=atribuido):", fh)
    log("row\\col       " + " ".join(f"{g[:8]:>9}" for g in groups), fh)
    for g in groups:
        log(f"{g:14}" + " ".join(f"{conf[g][gg]:9d}" for gg in groups), fh)

    # figura
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        cmap = {"african":"#c0392b","european":"#2980b9","indicine":"#27ae60",
                "east_asian":"#8e44ad","iberian":"#e67e22"}
        fig, ax = plt.subplots(figsize=(8, 6.5))
        for g in groups:
            mask = (y == g)
            ax.scatter(PC[mask,0], PC[mask,1], s=14, alpha=0.7,
                       c=cmap.get(g,"#555"), label=f"{g} (n={mask.sum()})", edgecolors="none")
        ax.set_xlabel(f"PC1 ({100*ev[0]:.1f}%)"); ax.set_ylabel(f"PC2 ({100*ev[1]:.1f}%)")
        ax.set_title("Run9 PRLR-window PCA por grupo coarse (sanidade regra #1)\njanela unica sob sweep — coerencia, nao ancestralidade genome-wide")
        ax.legend(fontsize=8, framealpha=0.9)
        fig.tight_layout(); fig.savefig(PNG, dpi=130)
        log(f"\n-> figura {PNG.name}", fh)
    except Exception as e:
        log(f"figura falhou: {e}", fh)
    log(f"-> relatorio {REP.name}", fh)
    fh.close()

if __name__ == "__main__":
    main()
