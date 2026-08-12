#!/usr/bin/env python3
"""
fair_control.py — controle NEGATIVO justo, com painel de referencia fixo.

Problema do LOO anterior: ao remover uma raca do painel, o painel EUT encolhia
de 6 para 5 racas — logo a raca testada era avaliada contra um painel MENOR
que o usado para o Senepol. Comparacao enviesada.

Aqui o painel de referencia e FIXO (EUT = AN, HO, HFD, SAL, LM, MONT; AFT;
ZEB) e testamos racas europeias que simplesmente NAO fazem parte dele —
exatamente a mesma condicao do Senepol. Nada e removido, nada encolhe.

Se racas europeias puras (Jersey, Devon, Sussex, Red Poll, Galloway...) fora
do painel recebem X% de ancestralidade africana espuria, entao o AFT do
Senepol so conta como herança se ultrapassar X.
"""
from __future__ import annotations
import os, sys, pathlib
import numpy as np
import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).parent))
from supervised_ancestry import load_ref_freqs, load_dosages, fit_one

DATA_DIR = pathlib.Path(os.environ.get("SLICK_DATA_DIR", "data"))
BASE = DATA_DIR / "senepol_widde" / "decker2014" / "plink" / "fair"
GROUPS = ["EUT", "AFT", "ZEB"]

grp = pd.read_csv(BASE / "groups.txt", sep="\t", header=None,
                  names=["FID", "IID", "grupo"])
ref, _ = load_ref_freqs(str(BASE / "freq.frq.strat"), GROUPS)
meta, dos, counted = load_dosages(str(BASE / "all.raw"))

shared = [s for s in dos.columns if s in ref.index]
ref, dos = ref.loc[shared], dos[shared]
flip = (counted.loc[shared].astype(str).values != ref["A1"].astype(str).values)
if flip.any():
    dos = dos.copy()
    dos.iloc[:, flip] = 2 - dos.iloc[:, flip]
print(f"[fair] {len(shared):,} marcadores | painel EUT fixo (6 racas)\n")

F = ref[GROUPS].to_numpy(dtype=float)
X = dos.to_numpy(dtype=float)
alvo = grp[grp["grupo"].isin(["TEST", "SENP"])]

rows = []
for pop, sub in alvo.groupby("FID"):
    idx = sub.index.to_numpy()
    Q = np.vstack([fit_one(X[i], F, len(GROUPS)) for i in idx])
    m = Q.mean(axis=0) * 100
    rows.append({"pop": pop, "n": len(idx), "EUT": m[0], "AFT": m[1], "ZEB": m[2],
                 "AFT_max": Q[:, 1].max() * 100,
                 "eh_senepol": pop == "SENP"})

df = pd.DataFrame(rows).sort_values("AFT")
df.to_csv(BASE / "fair_control.tsv", sep="\t", index=False, float_format="%.4f")

print(f"{'raca':8s} {'n':>4s}  {'EUT':>7s} {'AFT':>7s} {'ZEB':>7s}")
print("-" * 40)
for _, r in df.iterrows():
    marca = "  <<< SENEPOL" if r["eh_senepol"] else ""
    print(f"{r['pop']:8s} {int(r['n']):4d}  {r['EUT']:6.1f}% {r['AFT']:6.2f}% {r['ZEB']:6.2f}%{marca}")

eur = df[~df["eh_senepol"]]
sen = df[df["eh_senepol"]].iloc[0]
print("\n" + "=" * 60)
print(f"  racas europeias fora do painel : n={len(eur)}")
print(f"    AFT espurio  mediana {eur['AFT'].median():5.2f}%   "
      f"faixa {eur['AFT'].min():.2f}–{eur['AFT'].max():.2f}%")
print(f"  SENEPOL                        : AFT {sen['AFT']:5.2f}%")
pct = (eur["AFT"] < sen["AFT"]).mean() * 100
print(f"  racas europeias com AFT MENOR que o Senepol: {pct:.0f}%")
print("=" * 60)
