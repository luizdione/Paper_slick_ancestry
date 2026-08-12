#!/usr/bin/env python3
"""
haplotype_ancestry_window.py — de que ancestral vem o haplotipo que carrega o slick?

Medir COMPRIMENTO de trato falha aqui por dois motivos. Um seculo e pouco tempo
para o relogio da recombinacao, e o comprimento confunde idade do alelo com
consanguinidade do animal que o carrega — foi o que aconteceu com o unico
homozigoto do SLICK1.

A pergunta pode ser feita de outro jeito. Em vez de perguntar HA QUANTO TEMPO,
perguntar DE QUEM: os SNPs que viajam com o alelo, numa janela ao redor dele,
parecem europeus, africanos, zebuinos ou criollos? Isso nao depende de tempo
decorrido, so de com quem o haplotipo se parece.

Para cada haplotipo portador e cada janela, calculamos a log-verossimilhanca de
observar aqueles alelos sob as frequencias de cada painel ancestral, e
normalizamos em probabilidade posterior. Varias janelas (10 kb a 2 Mb) mostram se
a atribuicao e estavel ou se muda de escala — uma janela pequena reflete o
haplotipo fundador, uma janela grande dilui no fundo da raca.
"""
from __future__ import annotations
import argparse, gzip
import numpy as np
import pandas as pd

EPS = 1e-4


def load_phased(path: str, centre: int, flank: int):
    pos, rows, samples = [], [], []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            f = line.rstrip("\n").split("\t")
            if line.startswith("#CHROM"):
                samples = f[9:]
                continue
            p = int(f[1])
            if abs(p - centre) > flank:
                continue
            r = np.empty(2 * len(samples), dtype=np.int8)
            for i, g in enumerate(f[9:]):
                a = g.split(":")[0]
                s = "|" if "|" in a else "/"
                x = a.split(s)
                r[2 * i] = int(x[0]) if x[0].isdigit() else -1
                r[2 * i + 1] = int(x[1]) if len(x) > 1 and x[1].isdigit() else -1
            pos.append(p)
            rows.append(r)
    return np.array(pos), np.vstack(rows), samples


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--groups", required=True)
    ap.add_argument("--carriers", required=True, help="portadores (homo + het)")
    ap.add_argument("--pos", type=int, required=True)
    ap.add_argument("--label", default="SLICK")
    ap.add_argument("--windows", default="10000,50000,200000,1000000,2000000")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    wins = [int(x) for x in a.windows.split(",")]
    pos, H, samples = load_phased(a.vcf, a.pos, max(wins))
    idx = {s: i for i, s in enumerate(samples)}
    foc = int(np.argmin(np.abs(pos - a.pos)))

    grp = pd.read_csv(a.groups, sep="\t", header=None, names=["sample", "grupo"])
    g = dict(zip(grp["sample"], grp["grupo"]))
    hap_grp = np.array([g.get(s, "NA") for s in samples for _ in (0, 1)])
    GRUPOS = [x for x in ["EUR_TAURINE", "AFR_TAURINE", "AFR_SANGA", "ZEBU", "CRIOLLO"]
              if (hap_grp == x).sum() >= 8]

    carriers = [s for s in (l.strip() for l in open(a.carriers)) if s in idx]
    car_haps = {s: (2 * idx[s], 2 * idx[s] + 1) for s in carriers}
    is_car = np.zeros(H.shape[1], dtype=bool)
    for h in car_haps.values():
        is_car[list(h)] = True

    print(f"[{a.label}] {len(pos):,} sitios na janela maxima | portadores: {len(carriers)}")
    print(f"[{a.label}] paineis: " + ", ".join(f"{x} ({(hap_grp==x).sum()} hap)" for x in GRUPOS) + "\n")

    # frequencia do alelo 1 por painel, EXCLUINDO os proprios portadores
    F = {}
    for x in GRUPOS:
        m = (hap_grp == x) & ~is_car
        sub = H[:, m].astype(float)
        sub[sub < 0] = np.nan
        F[x] = np.clip(np.nanmean(sub, axis=1), EPS, 1 - EPS)

    # o haplotipo portador de cada animal: o que mais concorda com o modal dos
    # homozigotos e definido no script de compartilhamento; aqui usamos ambos os
    # haplotipos e reportamos o de maior verossimilhanca criolla como candidato
    rows = []
    for w in wins:
        sel = np.where(np.abs(pos - a.pos) <= w)[0]
        if len(sel) < 5:
            continue
        for s in carriers:
            for k, h in enumerate(car_haps[s]):
                obs = H[sel, h]
                ok = obs >= 0
                if ok.sum() < 5:
                    continue
                ll = {}
                for x in GRUPOS:
                    p = F[x][sel][ok]
                    o = obs[ok]
                    ll[x] = float(np.sum(np.where(o == 1, np.log(p), np.log(1 - p))))
                m = max(ll.values())
                post = {x: np.exp(v - m) for x, v in ll.items()}
                tot = sum(post.values())
                post = {x: v / tot for x, v in post.items()}
                best = max(post, key=post.get)
                rows.append({"janela_kb": w // 1000, "amostra": s, "hap": k,
                             "n_snp": int(ok.sum()), "melhor": best,
                             **{f"p_{x}": post[x] for x in GRUPOS}})

    d = pd.DataFrame(rows)
    d.to_csv(a.out, sep="\t", index=False, float_format="%.4f")

    print("=== atribuicao por janela (probabilidade posterior do melhor painel) ===")
    for w in sorted(d["janela_kb"].unique()):
        sub = d[d["janela_kb"] == w]
        cnt = sub["melhor"].value_counts()
        pmed = sub.groupby("melhor").apply(
            lambda x: x[f"p_{x.name}"].median(), include_groups=False)
        n_snp = int(sub["n_snp"].median())
        parts = [f"{k} {v}hap (p~{pmed.get(k, float('nan')):.2f})" for k, v in cnt.items()]
        print(f"  {w:5d} kb  ({n_snp:4d} SNPs)  " + " | ".join(parts))

    print("\n=== por animal, na janela de 200 kb ===")
    sub = d[d["janela_kb"] == 200]
    for s in carriers:
        ss = sub[sub["amostra"] == s]
        if len(ss) == 0:
            continue
        txt = "  ".join(f"hap{int(r.hap)}: {r.melhor} ({max(r[c] for c in ss.columns if c.startswith('p_')):.2f})"
                        for _, r in ss.iterrows())
        print(f"  {s:16s} {g.get(s,'NA'):12s} {txt}")

    print(f"\n  -> {a.out}")


if __name__ == "__main__":
    main()
