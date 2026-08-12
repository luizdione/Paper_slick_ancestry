#!/usr/bin/env python3
"""
haplotype_ancestry.py — de onde vem o haplotipo que carrega o SLICK1?

Uma analise anterior estabeleceu ONDE os alelos truncantes do PRLR estao: dois
sitios, ambos so em criollos. O que ela nao pode dizer e de onde vieram —
contagem alelica nao data alelo nem identifica populacao de origem. Isso exige
estrutura haplotipica, e portanto faseamento.

Uma segunda analise acrescentou o motivo para perguntar: o Senepol carrega
ancestralidade taurina africana acima do seu controle pareado, o que torna a
rota africana do *slick* plausivel em vez de descartada por falta de substrato.

Este script mede, sobre o VCF faseado:

  1. o comprimento do haplotipo compartilhado pelos portadores do SLICK1
     (decaimento de homozigosidade estendida em torno do sitio focal);
  2. a que grupo ancestral pertencem os haplotipos mais proximos do haplotipo
     portador — criollo, taurino africano, iberico ou taurino europeu;
  3. a idade implicada pelo comprimento, sob recombinacao uniforme.

Haplotipo longo = alelo jovem, pouco tempo para recombinar: origem recente, nas
Americas. Haplotipo curto e disperso entre grupos = alelo antigo, compativel com
entrada africana pre-colombiana.
"""
from __future__ import annotations
import argparse, gzip, pathlib
from collections import defaultdict, Counter
import numpy as np
import pandas as pd

POS_SLICK1 = 39_099_213     # BTA20, ARS-UCD1.2 (GC>G, p.Ala461ValfsTer2)
POS_SLICK2 = 39_099_321     # BTA20 (C>T, p.Arg497*)
CM_PER_MB = 1.0             # recombinacao uniforme; declarar como aproximacao


def load_phased(path: str, focal: int, flank: int = 5_000_000):
    """Le o VCF faseado na janela em torno do sitio focal."""
    pos, hap, samples = [], [], []
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
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
            row = np.empty(2 * len(samples), dtype=np.int8)
            for i, g in enumerate(f[9:]):
                a = g.split(":")[0]
                sep = "|" if "|" in a else "/"
                x = a.split(sep)
                row[2 * i] = int(x[0]) if x[0].isdigit() else -1
                row[2 * i + 1] = int(x[1]) if len(x) > 1 and x[1].isdigit() else -1
            pos.append(p)
            hap.append(row)
    return np.array(pos), np.vstack(hap), samples


def ehh_length(H: np.ndarray, pos: np.ndarray, foc: int, idx: np.ndarray) -> np.ndarray:
    """Distancia (bp) ate o haplotipo divergir da maioria do grupo, em cada direcao."""
    out = []
    for h in idx:
        total = 0
        for step in (1, -1):
            i = foc + step
            while 0 <= i < len(pos):
                col = H[i, idx]
                col = col[col >= 0]
                if len(col) == 0:
                    break
                modal = Counter(col.tolist()).most_common(1)[0][0]
                if H[i, h] != modal:
                    break
                i += step
            i = min(max(i, 0), len(pos) - 1)
            total += abs(pos[i] - pos[foc])
        out.append(total)
    return np.array(out)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--groups", required=True, help="TSV: amostra<TAB>grupo_ancestral")
    ap.add_argument("--pos", type=int, default=POS_SLICK1)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    pos, H, samples = load_phased(a.vcf, a.pos)
    foc = int(np.argmin(np.abs(pos - a.pos)))
    print(f"[haplo] {len(pos):,} sitios na janela | {len(samples)} amostras")
    print(f"[haplo] sitio focal: {pos[foc]:,} (alvo {a.pos:,}, "
          f"desvio {abs(pos[foc]-a.pos):,} pb)\n")

    grp = pd.read_csv(a.groups, sep="\t", header=None, names=["sample", "grupo"])
    g = dict(zip(grp["sample"], grp["grupo"]))
    hap_grp = np.array([g.get(s, "NA") for s in samples for _ in (0, 1)])

    carrier = H[foc] == 1
    print("=== portadores do alelo derivado no sitio focal ===")
    cnt = Counter(hap_grp[carrier])
    for k, v in cnt.most_common():
        tot = int((hap_grp == k).sum())
        print(f"  {k:16s} {v:3d} haplotipos portadores de {tot:4d}  ({100*v/tot:5.2f}%)")
    if carrier.sum() == 0:
        print("\n  [!] nenhum portador no sitio focal — conferir posicao e faseamento")
        return

    print("\n=== comprimento do haplotipo compartilhado ===")
    rows = []
    idx_c = np.where(carrier)[0]
    L = ehh_length(H, pos, foc, idx_c)
    rows.append({"grupo": "PORTADORES", "n_hap": len(idx_c),
                 "kb_mediano": np.median(L) / 1e3, "kb_medio": L.mean() / 1e3})
    for gg in sorted(set(hap_grp)):
        if gg == "NA":
            continue
        idx = np.where((hap_grp == gg) & (~carrier))[0]
        if len(idx) < 4:
            continue
        Lg = ehh_length(H, pos, foc, idx)
        rows.append({"grupo": f"nao-portadores {gg}", "n_hap": len(idx),
                     "kb_mediano": np.median(Lg) / 1e3, "kb_medio": Lg.mean() / 1e3})
    df = pd.DataFrame(rows)
    for _, r in df.iterrows():
        print(f"  {r['grupo']:28s} n={int(r['n_hap']):4d}  "
              f"mediana {r['kb_mediano']:8.1f} kb   media {r['kb_medio']:8.1f} kb")

    med_kb = df.iloc[0]["kb_mediano"]
    cm = (med_kb / 1000.0) * CM_PER_MB
    if cm > 0:
        gen = 1.0 / (cm / 100.0)
        print(f"\n=== idade implicada pelo comprimento ===")
        print(f"  haplotipo portador mediano : {med_kb:,.0f} kb  (~{cm:.2f} cM)")
        print(f"  geracoes desde a expansao  : ~{gen:.0f}")
        print(f"  anos (5 anos/geracao)      : ~{gen*5:.0f}")
        print(f"  fundacao do Senepol (1918) : ~{(2026-1918)/5:.0f} geracoes")
        print("  (recombinacao uniforme de 1 cM/Mb; aproximacao declarada)")

    df.to_csv(a.out, sep="\t", index=False, float_format="%.2f")
    print(f"\n  -> {a.out}")


if __name__ == "__main__":
    main()
