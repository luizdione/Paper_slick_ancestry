#!/usr/bin/env python3
"""
haplotype_sharing.py — quem compartilha o haplotipo que carrega o alelo *slick*,
e por quanto?

Os sitios truncantes do PRLR nao sobrevivem ao faseamento: SLICK1 e um indel raro
(AC baixo) e SLICK2 e raro, e ambos caem nos filtros que tornam o faseamento
tratavel. Isso nao impede a analise — impede apenas a versao ingenua dela.

A estrategia aqui:

  1. Os **homozigotos** para o alelo derivado tem, por definicao, os DOIS
     haplotipos portadores. Eles definem o haplotipo de referencia sem ambiguidade
     de fase.
  2. Para cada haplotipo do painel, medimos ate onde ele acompanha esse haplotipo
     de referencia a partir do sitio focal, em cada direcao.
  3. Agregamos por grupo ancestral. Se o haplotipo portador for jovem e local, o
     compartilhamento longo fica restrito aos criollos. Se for antigo e de entrada
     africana, aparecem tratos compartilhados — mais curtos — em taurinos africanos.

Haplotipo longo = pouco tempo para recombinar = alelo recente.
"""
from __future__ import annotations
import argparse, gzip
from collections import Counter
import numpy as np
import pandas as pd

CM_PER_MB = 1.0     # Beagle rodou sem mapa; aproximacao declarada


def load_phased(path: str, focal: int, flank: int):
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
            if abs(p - focal) > flank:
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
    ap.add_argument("--homozygotes", required=True,
                    help="amostras homozigotas para o alelo derivado, uma por linha")
    ap.add_argument("--heterozygotes", default="", help="opcional, uma por linha")
    ap.add_argument("--pos", type=int, required=True)
    ap.add_argument("--flank", type=int, default=5_000_000)
    ap.add_argument("--label", default="SLICK")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    pos, H, samples = load_phased(a.vcf, a.pos, a.flank)
    foc = int(np.argmin(np.abs(pos - a.pos)))
    print(f"[{a.label}] {len(pos):,} sitios | {len(samples)} amostras")
    print(f"[{a.label}] sitio focal do faseamento: {pos[foc]:,} "
          f"(alvo {a.pos:,}; desvio {abs(pos[foc]-a.pos):,} pb)\n")

    idx_of = {s: i for i, s in enumerate(samples)}
    homo = [l.strip() for l in open(a.homozygotes) if l.strip()]
    homo = [s for s in homo if s in idx_of]
    if not homo:
        print("  [!] nenhum homozigoto presente no VCF faseado — sem referencia de fase")
        return
    het = []
    if a.heterozygotes:
        het = [l.strip() for l in open(a.heterozygotes) if l.strip() and l.strip() in idx_of]
    print(f"  homozigotos usados como referencia : {len(homo)}  ({', '.join(homo)})")
    print(f"  heterozigotos (fase ambigua)       : {len(het)}\n")

    ref_haps = [2 * idx_of[s] + h for s in homo for h in (0, 1)]

    grp = pd.read_csv(a.groups, sep="\t", header=None, names=["sample", "grupo"])
    g = dict(zip(grp["sample"], grp["grupo"]))
    hap_grp = np.array([g.get(s, "NA") for s in samples for _ in (0, 1)])
    hap_smp = np.array([s for s in samples for _ in (0, 1)])

    # haplotipo de referencia = alelo modal dos homozigotos em cada sitio
    ref = np.empty(len(pos), dtype=np.int8)
    for i in range(len(pos)):
        col = H[i, ref_haps]
        col = col[col >= 0]
        ref[i] = Counter(col.tolist()).most_common(1)[0][0] if len(col) else -1

    def shared_len(h: int) -> float:
        """bp de concordancia continua com o haplotipo de referencia, ambas as direcoes."""
        total = 0
        for step in (1, -1):
            i = foc
            while 0 <= i < len(pos):
                if H[i, h] < 0 or ref[i] < 0 or H[i, h] != ref[i]:
                    break
                i += step
            i = min(max(i, 0), len(pos) - 1)
            total += abs(pos[i] - pos[foc])
        return total

    L = np.array([shared_len(h) for h in range(H.shape[1])])
    is_ref = np.zeros(H.shape[1], dtype=bool)
    is_ref[ref_haps] = True
    is_het = np.isin(hap_smp, het)

    print("=== comprimento do trato compartilhado com o haplotipo portador ===")
    rows = []
    for nome, mask in [("PORTADORES homozigotos", is_ref),
                       ("portadores heterozigotos", is_het & ~is_ref)]:
        if mask.sum():
            v = L[mask] / 1e3
            rows.append({"grupo": nome, "n_hap": int(mask.sum()),
                         "kb_mediano": np.median(v), "kb_p90": np.percentile(v, 90)})
    for gg in sorted(set(hap_grp)):
        if gg == "NA":
            continue
        m = (hap_grp == gg) & ~is_ref & ~is_het
        if m.sum() < 4:
            continue
        v = L[m] / 1e3
        rows.append({"grupo": gg, "n_hap": int(m.sum()),
                     "kb_mediano": np.median(v), "kb_p90": np.percentile(v, 90)})
    df = pd.DataFrame(rows)
    for _, r in df.iterrows():
        print(f"  {r['grupo']:26s} n={int(r['n_hap']):5d}   "
              f"mediana {r['kb_mediano']:8.1f} kb   p90 {r['kb_p90']:9.1f} kb")

    # quem mais compartilha: os 20 haplotipos nao-portadores mais proximos
    print("\n=== os 20 haplotipos nao-portadores com maior compartilhamento ===")
    cand = np.where(~is_ref & ~is_het & (hap_grp != "NA"))[0]
    top = cand[np.argsort(-L[cand])][:20]
    print(Counter(hap_grp[top]).most_common())

    med = df.iloc[0]["kb_mediano"] if len(df) else 0
    if med > 0:
        cm = (med / 1000.0) * CM_PER_MB
        gen = 100.0 / cm if cm > 0 else float("inf")
        print(f"\n=== idade implicada pelo trato dos homozigotos ===")
        print(f"  mediana : {med:,.0f} kb (~{cm:.2f} cM)")
        print(f"  geracoes: ~{gen:.0f}   |   anos (5 a/ger): ~{gen*5:.0f}")
        print(f"  referencia: fundacao do Senepol em 1918 = ~{(2026-1918)/5:.0f} geracoes;")
        print(f"              chegada de gado iberico as Americas ~1500 = ~{(2026-1500)/5:.0f} geracoes")
        print("  (recombinacao uniforme de 1 cM/Mb — aproximacao, Beagle rodou sem mapa)")

    df.to_csv(a.out, sep="\t", index=False, float_format="%.2f")
    print(f"\n  -> {a.out}")


if __name__ == "__main__":
    main()
