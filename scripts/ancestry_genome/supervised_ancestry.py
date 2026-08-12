#!/usr/bin/env python3
"""
supervised_ancestry.py — ancestralidade global supervisionada por maxima verossimilhanca.

Substitui o ADMIXTURE supervisionado, que nao executa neste ambiente
(binario 1.3.0 de 2015 segfalha sob Docker/WSL2, com .bed comprovadamente
integro e byte-identico entre PLINK 1.9 e PLINK 2).

Modelo (identico ao do ADMIXTURE supervisionado):
  para o individuo i, marcador j, com genotipo g_ij em {0,1,2} contando A1,
      p_ij(q) = sum_k q_k * f_jk
      L(q)    = sum_j log[ Binom(g_ij ; 2, p_ij(q)) ]
  maximizado sobre o simplex { q_k >= 0, sum_k q_k = 1 }.

f_jk = frequencia de A1 no painel de referencia k, estimada dos individuos
rotulados (nunca do proprio alvo — evita circularidade).

Uso:
  python supervised_ancestry.py --freq freq_by_group.frq.strat \
                                --raw sen_dosage.raw --out ancestry_SEN.tsv
"""
from __future__ import annotations
import argparse, sys
import numpy as np
import pandas as pd
from scipy.optimize import minimize

EPS = 1e-6  # evita log(0) em marcadores fixados num painel


def load_ref_freqs(path: str, groups: list[str]) -> tuple[pd.DataFrame, list[str]]:
    """Le .frq.strat do PLINK -> matriz F (marcadores x grupos) de freq. de A1."""
    df = pd.read_csv(path, sep=r"\s+")
    df = df[df["CLST"].isin(groups)]
    # MAF na saida --freq --within e a frequencia de A1 dentro do cluster
    wide = df.pivot(index="SNP", columns="CLST", values="MAF")
    wide = wide.dropna(how="any")          # exige o marcador nos 3 paineis
    a1 = df.drop_duplicates("SNP").set_index("SNP")["A1"]
    wide["A1"] = a1.reindex(wide.index)
    return wide, list(groups)


def load_dosages(path: str) -> tuple[pd.DataFrame, pd.Series]:
    """Le .raw do PLINK (--recode A). Colunas de SNP tem sufixo _<A1>."""
    raw = pd.read_csv(path, sep=r"\s+")
    meta = raw[["FID", "IID"]].copy()
    dos = raw.drop(columns=["FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE"])
    # separa "SNPID_A1" -> SNPID e o alelo contado
    snp_ids, counted = [], []
    for c in dos.columns:
        sid, _, al = c.rpartition("_")
        snp_ids.append(sid)
        counted.append(al)
    dos.columns = snp_ids
    return meta, dos, pd.Series(counted, index=snp_ids)


def neg_loglik(q: np.ndarray, G: np.ndarray, F: np.ndarray) -> float:
    """-log L para um individuo. G: (m,) genotipos; F: (m,k) freq de referencia."""
    p = np.clip(F @ q, EPS, 1 - EPS)
    return -np.sum(G * np.log(p) + (2 - G) * np.log(1 - p))


def fit_one(G: np.ndarray, F: np.ndarray, k: int) -> np.ndarray:
    ok = ~np.isnan(G)
    Gv, Fv = G[ok], F[ok]
    cons = ({"type": "eq", "fun": lambda q: q.sum() - 1.0},)
    res = minimize(neg_loglik, np.full(k, 1.0 / k), args=(Gv, Fv),
                   method="SLSQP", bounds=[(0.0, 1.0)] * k, constraints=cons,
                   options={"maxiter": 500, "ftol": 1e-10})
    q = np.clip(res.x, 0, 1)
    return q / q.sum()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--freq", required=True)
    ap.add_argument("--raw", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--groups", default="EUT,AFT,ZEB")
    a = ap.parse_args()

    groups = a.groups.split(",")
    ref, groups = load_ref_freqs(a.freq, groups)
    meta, dos, counted = load_dosages(a.raw)

    shared = [s for s in dos.columns if s in ref.index]
    print(f"[ancestry] marcadores no painel de referencia : {len(ref):,}")
    print(f"[ancestry] marcadores no alvo                 : {dos.shape[1]:,}")
    print(f"[ancestry] compartilhados                     : {len(shared):,}")

    ref = ref.loc[shared]
    dos = dos[shared]
    counted = counted.loc[shared]

    # coerencia de alelo contado: se o .raw conta o outro alelo, inverte a dosagem.
    # AMBOS os lados como str: com alelos numericos (codificacao Illumina A/B = 1/2)
    # o pandas le o A1 do .frq como int e o sufixo do .raw como str, e '1' != 1
    # seria verdadeiro em todos os marcadores — invertendo o dataset inteiro.
    flip = (counted.astype(str).values != ref["A1"].astype(str).values)
    n_flip = int(flip.sum())
    if n_flip:
        print(f"[ancestry] alelo contado invertido em {n_flip:,} marcadores — corrigindo")
        dos.iloc[:, flip] = 2 - dos.iloc[:, flip]

    F = ref[groups].to_numpy(dtype=float)
    X = dos.to_numpy(dtype=float)

    Q = np.vstack([fit_one(X[i], F, len(groups)) for i in range(X.shape[0])])

    out = meta.copy()
    for j, g in enumerate(groups):
        out[g] = Q[:, j]
    out.to_csv(a.out, sep="\t", index=False, float_format="%.5f")

    print(f"\n[ancestry] n = {len(out)} individuos, {len(shared):,} marcadores\n")
    print("=== proporcao ancestral media ===")
    for g in groups:
        v = out[g]
        print(f"  {g}: {100*v.mean():6.2f}%   (dp {100*v.std():.2f}  "
              f"min {100*v.min():.2f}  max {100*v.max():.2f})")


if __name__ == "__main__":
    main()
