#!/usr/bin/env python3
"""
loo_control.py — controle negativo por leave-one-population-out.

Pergunta: quando uma raca taurina EUROPEIA PURA e estimada contra um painel
que NAO a contem, quanta ancestralidade africana/zebuina ESPURIA ela recebe?

Se uma raca sabidamente sem ancestralidade africana recebe X% de AFT so por
nao ter painel proximal, entao o AFT observado no Senepol tem de ser lido
contra esse X — nao contra zero. E o controle que decide se os 11,3% do
Senepol sao sinal ou artefato.
"""
from __future__ import annotations
import os, sys, pathlib
import numpy as np
import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).parent))
from supervised_ancestry import load_ref_freqs, load_dosages, fit_one

DATA_DIR = pathlib.Path(os.environ.get("SLICK_DATA_DIR", "data"))
BASE = DATA_DIR / "senepol_widde" / "admixture"
GROUPS = ["EUT", "AFT", "ZEB"]
POPS = ["ANG", "HOL", "HFD", "LMS", "MON", "SAL"]

print("[loo] carregando dosagens de todos os individuos...")
meta, dos, counted = load_dosages(str(BASE / "all_dosage.raw"))

rows = []
for pop in POPS:
    ref, _ = load_ref_freqs(str(BASE / "loo" / f"freq_no_{pop}.frq.strat"), GROUPS)
    shared = [s for s in dos.columns if s in ref.index]
    r = ref.loc[shared]
    d = dos[shared]
    flip = (counted.loc[shared].astype(str).values != r["A1"].astype(str).values)
    if flip.any():
        d = d.copy()
        d.iloc[:, flip] = 2 - d.iloc[:, flip]

    idx = meta.index[meta["FID"] == pop]
    F = r[GROUPS].to_numpy(dtype=float)
    X = d.iloc[idx].to_numpy(dtype=float)
    Q = np.vstack([fit_one(X[i], F, len(GROUPS)) for i in range(X.shape[0])])
    m = Q.mean(axis=0) * 100
    s = Q.std(axis=0) * 100
    rows.append({"pop": pop, "n": len(idx), "marcadores": len(shared),
                 "EUT": m[0], "AFT": m[1], "ZEB": m[2],
                 "AFT_dp": s[1], "AFT_max": Q[:, 1].max() * 100})
    print(f"  {pop:4s} n={len(idx):3d}  EUT {m[0]:5.1f}%  AFT {m[1]:5.2f}%  ZEB {m[2]:5.2f}%")

df = pd.DataFrame(rows)
df.to_csv(BASE / "loo_control.tsv", sep="\t", index=False, float_format="%.3f")

print("\n=== REFERENCIA ===")
print(f"  AFT espurio medio em taurinos europeus puros : {df['AFT'].mean():.2f}%")
print(f"  AFT espurio maximo                           : {df['AFT'].max():.2f}%  ({df.loc[df['AFT'].idxmax(),'pop']})")
print(f"  AFT observado no Senepol                     : 11.31%")
print(f"  AFT reportado por Flori et al. 2012          :  0.60%")
