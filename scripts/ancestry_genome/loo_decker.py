#!/usr/bin/env python3
"""
loo_decker.py — controle negativo (leave-one-population-out) no dataset Decker et al. 2014.

Mede o piso de deteccao: quanta ancestralidade africana ESPURIA recebe uma raca
taurina europeia pura quando estimada contra um painel que nao a contem.

Inclui o proprio Red Poll no LOO — e o controle mais informativo, porque o
Red Poll e o fundador europeu do Senepol: mesma linhagem britanica, mesma
situacao de "raca sem painel proximal".

O AFT do Senepol so e interpretavel como herança se estiver ACIMA desse piso.
"""
from __future__ import annotations
import os, sys, pathlib
import numpy as np
import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).parent))
from supervised_ancestry import load_ref_freqs, load_dosages, fit_one

DATA_DIR = pathlib.Path(os.environ.get("SLICK_DATA_DIR", "data"))
BASE = DATA_DIR / "senepol_widde" / "decker2014" / "plink" / "anc"
GROUPS = ["EUT", "AFT", "ZEB"]
POPS = ["AN", "HO", "HFD", "SAL", "LM", "MONT", "REDP"]
NOME = {"AN": "Angus", "HO": "Holstein", "HFD": "Hereford", "SAL": "Salers",
        "LM": "Limousin", "MONT": "Montbeliarde", "REDP": "Red Poll"}

print("[loo] carregando dosagens...")
meta, dos, counted = load_dosages(str(BASE / "all_nd.raw"))

rows = []
for pop in POPS:
    ref, _ = load_ref_freqs(str(BASE / "loo" / f"f_no_{pop}.frq.strat"), GROUPS)
    shared = [s for s in dos.columns if s in ref.index]
    r, d = ref.loc[shared], dos[shared]
    flip = (counted.loc[shared].astype(str).values != r["A1"].astype(str).values)
    if flip.any():
        d = d.copy()
        d.iloc[:, flip] = 2 - d.iloc[:, flip]

    idx = meta.index[meta["FID"] == pop]
    F = r[GROUPS].to_numpy(dtype=float)
    X = d.iloc[idx].to_numpy(dtype=float)
    Q = np.vstack([fit_one(X[i], F, len(GROUPS)) for i in range(X.shape[0])])
    m, s = Q.mean(axis=0) * 100, Q.std(axis=0) * 100
    rows.append({"pop": pop, "nome": NOME[pop], "n": len(idx),
                 "EUT": m[0], "AFT": m[1], "ZEB": m[2],
                 "AFT_dp": s[1], "AFT_min": Q[:, 1].min() * 100,
                 "AFT_max": Q[:, 1].max() * 100})
    print(f"  {NOME[pop]:13s} n={len(idx):3d}   EUT {m[0]:5.1f}%   "
          f"AFT {m[1]:5.2f}% (max {Q[:,1].max()*100:4.2f})   ZEB {m[2]:5.2f}%")

df = pd.DataFrame(rows)
df.to_csv(BASE / "loo_decker.tsv", sep="\t", index=False, float_format="%.4f")

sen = pd.read_csv(BASE / "SENP_comRedPoll.tsv", sep="\t")
piso_med, piso_max = df["AFT"].mean(), df["AFT"].max()
sen_aft = sen["AFT"].mean() * 100

print("\n" + "=" * 62)
print(f"  piso de artefato (AFT espurio, media)   : {piso_med:5.2f}%")
print(f"  piso de artefato (maximo)               : {piso_max:5.2f}%  ({df.loc[df['AFT'].idxmax(),'nome']})")
print(f"  AFT do Red Poll (fundador, sem si mesmo): {df.loc[df['pop']=='REDP','AFT'].iloc[0]:5.2f}%")
print(f"  AFT do SENEPOL (painel com Red Poll)    : {sen_aft:5.2f}%")
print(f"  razao Senepol / piso medio              : {sen_aft/piso_med:5.2f}x")
print("=" * 62)
veredito = ("SINAL acima do piso" if sen_aft > piso_max else
            "INDISTINGUIVEL do artefato")
print(f"  VEREDITO: {veredito}")
