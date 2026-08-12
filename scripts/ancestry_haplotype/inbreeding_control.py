#!/usr/bin/env python3
"""
inbreeding_control.py — o trato longo do SLICK1 é idade do alelo ou endogamia do animal?

O homozigoto que ancora o SLICK1 é um Costeno con Cuernos, e essa raça tem duas
amostras no painel. Um animal de populacao pequena e consanguinea carrega tratos
longos de homozigosidade em TODO o genoma, nao so no locus de interesse. Se for
esse o caso, os 5.177 kb medem o pedigree do animal, nao a idade da mutacao.

O controle: medir, nos MESMOS animais, o comprimento de trato compartilhado em
posicoes focais SORTEADAS ao longo da mesma janela, longe do PRLR. Isso da a
linha de base de cada animal. O trato do sitio focal so significa idade do alelo
se exceder a propria linha de base do portador.
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


def tract(H, pos, foc, ref_haps, h):
    """bp de concordancia continua com o haplotipo modal do grupo de referencia."""
    total = 0
    for step in (1, -1):
        i = foc
        while 0 <= i < len(pos):
            col = H[i, ref_haps]
            col = col[col >= 0]
            if len(col) == 0 or H[i, h] < 0:
                break
            if H[i, h] != Counter(col.tolist()).most_common(1)[0][0]:
                break
            i += step
        i = min(max(i, 0), len(pos) - 1)
        total += abs(pos[i] - pos[foc])
    return total


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--homozygotes", required=True)
    ap.add_argument("--pos", type=int, required=True)
    ap.add_argument("--label", default="SLICK")
    ap.add_argument("--n-random", type=int, default=200)
    ap.add_argument("--min-dist", type=int, default=1_000_000,
                    help="distancia minima do sitio focal para os controles")
    a = ap.parse_args()

    rng = np.random.default_rng(SEED)
    pos, H, samples = load_phased(a.vcf)
    idx_of = {s: i for i, s in enumerate(samples)}
    homo = [s for s in (l.strip() for l in open(a.homozygotes)) if s in idx_of]
    ref_haps = [2 * idx_of[s] + h for s in homo for h in (0, 1)]
    foc = int(np.argmin(np.abs(pos - a.pos)))

    print(f"[{a.label}] {len(pos):,} sitios | homozigotos: {len(homo)} ({', '.join(homo)})")

    obs = np.array([tract(H, pos, foc, ref_haps, h) for h in ref_haps])
    print(f"\n  trato NO SITIO FOCAL      : mediana {np.median(obs)/1e3:9.1f} kb")

    # controles: posicoes sorteadas na mesma janela, longe do foco
    ok = np.where(np.abs(pos - a.pos) > a.min_dist)[0]
    sel = rng.choice(ok, size=min(a.n_random, len(ok)), replace=False)
    base = []
    for f2 in sel:
        base.append(np.median([tract(H, pos, f2, ref_haps, h) for h in ref_haps]))
    base = np.array(base)
    print(f"  linha de base do(s) animal(is): mediana {np.median(base)/1e3:9.1f} kb "
          f"(p90 {np.percentile(base,90)/1e3:.1f} kb, n={len(base)} posicoes)")

    p_emp = float((base >= np.median(obs)).mean())
    razao = np.median(obs) / max(np.median(base), 1)
    print(f"\n  razao foco / linha de base : {razao:8.1f}x")
    print(f"  P empirico (fracao de posicoes aleatorias com trato >= o do foco): {p_emp:.4f}")
    if p_emp < 0.05:
        print(f"  -> o trato do sitio focal EXCEDE a linha de base do proprio portador")
    else:
        print(f"  -> o trato do sitio focal NAO se distingue da linha de base:")
        print(f"     o comprimento reflete o pedigree do animal, nao a idade do alelo")


if __name__ == "__main__":
    main()
