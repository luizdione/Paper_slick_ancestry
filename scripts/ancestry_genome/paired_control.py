#!/usr/bin/env python3
"""
paired_control.py — o controle pareado que faltava.

Comparar o Senepol com "as racas europeias" e um controle ruim: elas vao de
0,01% a 9,63% de africano aparente, e essa faixa nao e ruido sem estrutura.
O vies parece acompanhar a DISTANCIA GENETICA ao painel de referencia. Uma raca
britanica proxima do Angus e do Hereford (que estao no painel) tem pouco a
"sobrar"; uma raca alpina ou italiana, distante de todas as seis, tem muito.

O Senepol foi fundado sobre o Red Poll — taurino britanico. Seu controle pareado
nao e "a Europa", e sim a linhagem britanica vermelha da qual ele descende, que
compartilha a base genetica e nao tem ancestralidade africana documentada.

Este script:
  1. mede a distancia genetica (Fst de Hudson) de cada raca ao painel EUT
  2. testa se o africano espurio cresce com essa distancia
  3. compara o Senepol contra o grupo britanico pareado, com bootstrap
"""
from __future__ import annotations
import os, sys, pathlib
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(pathlib.Path(__file__).parent))
from supervised_ancestry import load_ref_freqs, load_dosages, fit_one

DATA_DIR = pathlib.Path(os.environ.get("SLICK_DATA_DIR", "data"))
BASE = DATA_DIR / "senepol_widde" / "decker2014" / "plink" / "fair"
GROUPS = ["EUT", "AFT", "ZEB"]
SEED = 20260803

# linhagem britanica: mesma base genetica do fundador europeu do Senepol
BRITANICAS = {"REDP": "Red Poll", "DEV": "Devon", "SUSS": "Sussex",
              "LINC": "Lincoln Red", "SDEV": "South Devon", "SH": "Shorthorn",
              "GALL": "Galloway", "DEX": "Dexter", "KERR": "Kerry"}


def hudson_fst(p1: np.ndarray, p2: np.ndarray, n1: int, n2: int) -> float:
    """Fst de Hudson na forma razao-de-medias (Bhatia et al. 2013).
    Insensivel a n desigual, ao contrario de Weir-Cockerham por SNP promediado."""
    num = (p1 - p2) ** 2 - p1 * (1 - p1) / (n1 - 1) - p2 * (1 - p2) / (n2 - 1)
    den = p1 * (1 - p2) + p2 * (1 - p1)
    ok = np.isfinite(num) & np.isfinite(den) & (den > 0)
    return float(num[ok].sum() / den[ok].sum())


def main() -> None:
    rng = np.random.default_rng(SEED)
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
    F = ref[GROUPS].to_numpy(dtype=float)
    X = dos.to_numpy(dtype=float)
    print(f"[paired] {len(shared):,} marcadores\n")

    p_eut = F[:, 0]
    n_eut = int((grp["grupo"] == "EUT").sum())

    alvo = grp[grp["grupo"].isin(["TEST", "SENP"])]
    rows = []
    for pop, sub in alvo.groupby("FID"):
        idx = sub.index.to_numpy()
        Q = np.vstack([fit_one(X[i], F, len(GROUPS)) for i in idx])
        p_pop = np.nanmean(X[idx], axis=0) / 2.0
        rows.append({"pop": pop, "n": len(idx),
                     "AFT": Q[:, 1].mean() * 100, "AFT_sd": Q[:, 1].std() * 100,
                     "ZEB": Q[:, 2].mean() * 100,
                     "fst_ao_painel_EUT": hudson_fst(p_pop, p_eut, len(idx), n_eut),
                     "britanica": pop in BRITANICAS,
                     "senepol": pop == "SENP"})
    df = pd.DataFrame(rows).sort_values("fst_ao_painel_EUT")
    df.to_csv(BASE / "paired_control.tsv", sep="\t", index=False, float_format="%.5f")

    print("=== o vies acompanha a distancia genetica ao painel? ===")
    eur = df[~df["senepol"]]
    r, p = stats.pearsonr(eur["fst_ao_painel_EUT"], eur["AFT"])
    rho, prho = stats.spearmanr(eur["fst_ao_painel_EUT"], eur["AFT"])
    print(f"  Pearson  r = {r:+.3f}  (P = {p:.2e})")
    print(f"  Spearman rho = {rho:+.3f}  (P = {prho:.2e})")
    print(f"  n = {len(eur)} racas europeias\n")

    brit = df[df["britanica"] & ~df["senepol"]]
    cont = df[~df["britanica"] & ~df["senepol"]]
    sen = df[df["senepol"]].iloc[0]

    print("=== grupo britanico (controle pareado do Senepol) ===")
    for _, x in brit.sort_values("AFT").iterrows():
        print(f"  {BRITANICAS[x['pop']]:13s} n={int(x['n']):3d}  "
              f"Fst {x['fst_ao_painel_EUT']:.4f}  AFT {x['AFT']:5.2f}%")
    print(f"\n  britanicas  : media {brit['AFT'].mean():5.2f}%  "
          f"faixa {brit['AFT'].min():.2f}–{brit['AFT'].max():.2f}%  "
          f"Fst medio {brit['fst_ao_painel_EUT'].mean():.4f}")
    print(f"  continentais: media {cont['AFT'].mean():5.2f}%  "
          f"faixa {cont['AFT'].min():.2f}–{cont['AFT'].max():.2f}%  "
          f"Fst medio {cont['fst_ao_painel_EUT'].mean():.4f}")
    print(f"  SENEPOL     : {sen['AFT']:5.2f}%  "
          f"Fst {sen['fst_ao_painel_EUT']:.4f}")

    print("\n=== Senepol contra o controle pareado ===")
    t = stats.ttest_1samp(brit["AFT"], sen["AFT"])
    z = (sen["AFT"] - brit["AFT"].mean()) / brit["AFT"].std(ddof=1)
    bs = np.array([rng.choice(brit["AFT"].to_numpy(), len(brit), replace=True).mean()
                   for _ in range(10000)])
    print(f"  media britanica  : {brit['AFT'].mean():.2f}%  "
          f"(IC95 bootstrap {np.percentile(bs,2.5):.2f}–{np.percentile(bs,97.5):.2f})")
    print(f"  Senepol          : {sen['AFT']:.2f}%")
    print(f"  excesso          : {sen['AFT']-brit['AFT'].mean():+.2f} pontos percentuais")
    print(f"  z                : {z:+.2f} desvios-padrao do grupo britanico")
    print(f"  t (1-amostra)    : t = {t.statistic:+.2f}, P = {t.pvalue:.2e}")
    print(f"  britanicas acima do Senepol: {int((brit['AFT']>sen['AFT']).sum())}/{len(brit)}")

    print("\n=== ressalva: o Senepol tambem carrega zebu ===")
    print(f"  Senepol ZEB {sen['ZEB']:.2f}%  |  britanicas ZEB medio {brit['ZEB'].mean():.2f}%")
    print("  taurinos africanos tem afinidade parcial com zebu, entao parte do")
    print("  excesso africano pode ser arrasto do componente zebuino.")


if __name__ == "__main__":
    main()
