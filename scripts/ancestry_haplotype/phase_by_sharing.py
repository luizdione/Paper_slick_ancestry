#!/usr/bin/env python3
"""
phase_by_sharing.py — recuperar a datacao do SLICK1 usando portadores menos consanguineos.

O problema: o unico homozigoto do SLICK1 e um Costeno con Cuernos cujo genoma
inteiro esta em blocos de ~5 Mb. O trato no locus nao excede a linha de base dele,
entao nao informa idade.

A saida: os HETEROZIGOTOS. Eles nao fixam a fase por si, mas o alelo, sendo raro e
presumivelmente identico por descendencia, deve estar no mesmo haplotipo de fundo
do homozigoto. Entao:

  1. o homozigoto define o haplotipo portador de referencia H*;
  2. em cada heterozigoto, dos dois haplotipos, o que acompanha H* na vizinhanca
     imediata do sitio e o portador — a assimetria entre os dois e a propria
     evidencia de que a atribuicao funcionou;
  3. mede-se o trato desse haplotipo contra a linha de base DAQUELE animal.

Se algum heterozigoto for pouco consanguineo, o trato dele informa idade onde o
homozigoto nao informava.
"""
from __future__ import annotations
import argparse, gzip
from collections import Counter
import numpy as np
import pandas as pd

SEED = 20260803


def load_phased(path: str):
    pos, rows, samples = [], [], []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            f = line.rstrip("\n").split("\t")
            if line.startswith("#CHROM"):
                samples = f[9:]
                continue
            r = np.empty(2 * len(samples), dtype=np.int8)
            for i, g in enumerate(f[9:]):
                a = g.split(":")[0]
                s = "|" if "|" in a else "/"
                x = a.split(s)
                r[2 * i] = int(x[0]) if x[0].isdigit() else -1
                r[2 * i + 1] = int(x[1]) if len(x) > 1 and x[1].isdigit() else -1
            pos.append(int(f[1]))
            rows.append(r)
    return np.array(pos), np.vstack(rows), samples


def tract_vs(H, pos, foc, ref, h):
    """bp de concordancia continua entre o haplotipo h e o vetor de referencia ref."""
    total = 0
    for step in (1, -1):
        i = foc
        while 0 <= i < len(pos):
            if ref[i] < 0 or H[i, h] < 0 or H[i, h] != ref[i]:
                break
            i += step
        i = min(max(i, 0), len(pos) - 1)
        total += abs(pos[i] - pos[foc])
    return total


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--homozygotes", required=True)
    ap.add_argument("--heterozygotes", required=True)
    ap.add_argument("--pos", type=int, required=True)
    ap.add_argument("--groups", required=True)
    ap.add_argument("--n-random", type=int, default=150)
    ap.add_argument("--label", default="SLICK1")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    rng = np.random.default_rng(SEED)
    pos, H, samples = load_phased(a.vcf)
    idx = {s: i for i, s in enumerate(samples)}
    foc = int(np.argmin(np.abs(pos - a.pos)))
    homo = [s for s in (l.strip() for l in open(a.homozygotes)) if s in idx]
    het = [s for s in (l.strip() for l in open(a.heterozygotes)) if s in idx]
    grp = pd.read_csv(a.groups, sep="\t", header=None, names=["sample", "grupo"])
    g = dict(zip(grp["sample"], grp["grupo"]))

    print(f"[{a.label}] sitio focal {pos[foc]:,} | homozigotos {len(homo)} | heterozigotos {len(het)}\n")

    # haplotipo portador de referencia, a partir do(s) homozigoto(s)
    ref_h = [2 * idx[s] + k for s in homo for k in (0, 1)]
    ref = np.array([Counter(H[i, ref_h][H[i, ref_h] >= 0].tolist()).most_common(1)[0][0]
                    if (H[i, ref_h] >= 0).any() else -1 for i in range(len(pos))], dtype=np.int8)

    def baseline(hs):
        ok = np.where(np.abs(pos - a.pos) > 1_000_000)[0]
        sel = rng.choice(ok, size=min(a.n_random, len(ok)), replace=False)
        vals = []
        for f2 in sel:
            r2 = np.array([Counter(H[i, ref_h][H[i, ref_h] >= 0].tolist()).most_common(1)[0][0]
                           if (H[i, ref_h] >= 0).any() else -1 for i in (f2,)], dtype=np.int8)
            # linha de base = homozigosidade do proprio animal ao redor de f2
            vals.append(np.median([tract_self(H, pos, f2, h) for h in hs]))
        return np.array(vals)

    def tract_self(H, pos, f2, h):
        """trato de homozigosidade do proprio animal (os dois haplotipos iguais)."""
        par = h - (h % 2)
        h1, h2 = par, par + 1
        total = 0
        for step in (1, -1):
            i = f2
            while 0 <= i < len(pos):
                if H[i, h1] < 0 or H[i, h2] < 0 or H[i, h1] != H[i, h2]:
                    break
                i += step
            i = min(max(i, 0), len(pos) - 1)
            total += abs(pos[i] - pos[f2])
        return total

    rows = []
    todos = [(s, "homozigoto") for s in homo] + [(s, "heterozigoto") for s in het]
    for s, tipo in todos:
        h0, h1 = 2 * idx[s], 2 * idx[s] + 1
        t0, t1 = tract_vs(H, pos, foc, ref, h0), tract_vs(H, pos, foc, ref, h1)
        port, outro = (h0, h1) if t0 >= t1 else (h1, h0)
        t_port, t_outro = max(t0, t1), min(t0, t1)
        # consanguinidade do animal: homozigosidade mediana em 150 posicoes aleatorias
        ok = np.where(np.abs(pos - a.pos) > 1_000_000)[0]
        sel = rng.choice(ok, size=min(a.n_random, len(ok)), replace=False)
        base = np.median([tract_self(H, pos, f2, h0) for f2 in sel])
        rows.append({"amostra": s, "tipo": tipo, "grupo": g.get(s, "NA"),
                     "trato_portador_kb": t_port / 1e3,
                     "trato_outro_kb": t_outro / 1e3,
                     "assimetria": (t_port + 1) / (t_outro + 1),
                     "consanguinidade_kb": base / 1e3,
                     "razao_locus_base": (t_port + 1) / (base + 1)})

    d = pd.DataFrame(rows).sort_values("consanguinidade_kb")
    print("=== fase resolvida por compartilhamento, e consanguinidade de cada portador ===")
    print(f"{'amostra':16s} {'tipo':13s} {'portador':>10s} {'outro':>9s} {'assim':>7s} {'ROH base':>10s} {'razao':>7s}")
    print("-" * 78)
    for _, r in d.iterrows():
        print(f"{r.amostra:16s} {r.tipo:13s} {r.trato_portador_kb:9.1f}k {r.trato_outro_kb:8.1f}k "
              f"{r.assimetria:6.1f}x {r.consanguinidade_kb:9.1f}k {r.razao_locus_base:6.2f}x")

    d.to_csv(a.out, sep="\t", index=False, float_format="%.2f")
    menos = d.iloc[0]
    print(f"\n=== portador menos consanguineo: {menos.amostra} ({menos.grupo}) ===")
    print(f"  homozigosidade de fundo : {menos.consanguinidade_kb:,.0f} kb")
    print(f"  trato portador no locus : {menos.trato_portador_kb:,.0f} kb")
    print(f"  razao locus / fundo     : {menos.razao_locus_base:.2f}x")
    if menos.razao_locus_base > 2:
        cm = menos.trato_portador_kb / 1000.0
        print(f"  -> o trato EXCEDE o fundo; idade implicada ~{100/cm:.0f} geracoes "
              f"(~{100/cm*5:.0f} anos)")
    else:
        print(f"  -> ainda dentro do fundo; a datacao segue impossivel neste painel")


if __name__ == "__main__":
    main()
