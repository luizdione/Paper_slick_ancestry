# -*- coding: utf-8 -*-
"""Fase por LEITURA dos portadores compostos de alelos slick (sem BEAGLE, sem inferencia).

Os sitios slick ficam a 13-108 pb uns dos outros; uma leitura de 150 pb cobre dois deles. Cada
leitura que cobre os dois sitios mostra, diretamente, quais alelos estao no MESMO fragmento de
DNA. Isso confirma ou refuta a fase estatistica do BEAGLE nos animais compostos.

Roda dentro do container (samtools na PATH). Le `samtools view` e percorre a CIGAR de cada
leitura para achar a base nas duas posicoes de referencia.

ATENCAO: os BAMs do lote novo usam o contig `20:34100000-44100000` com coordenada
LOCAL. global = 34.100.000 + local - 1, logo local = global - 34.100.000 + 1.

uso: python3 fase_por_leitura.py <bam> <nomeA:posGlobalA:REF:ALT> <nomeB:posGlobalB:REF:ALT>
"""
import collections
import subprocess
import sys

CONTIG = "20:34100000-44100000"
DESLOC = 34100000 - 1


def base_em(pos_ref, inicio, cigar, seq):
    """Base da leitura na posicao de referencia `pos_ref` (1-based), ou None."""
    rpos, qpos, num = inicio, 0, ""
    for ch in cigar:
        if ch.isdigit():
            num += ch
            continue
        n = int(num)
        num = ""
        if ch in "M=X":
            if rpos <= pos_ref < rpos + n:
                return seq[qpos + (pos_ref - rpos)]
            rpos += n
            qpos += n
        elif ch in "IS":
            qpos += n
        elif ch in "DN":
            if rpos <= pos_ref < rpos + n:
                return "-"
            rpos += n
        # H, P: nao consomem
    return None


def main():
    bam = sys.argv[1]
    sitios = []
    for arg in sys.argv[2:4]:
        nome, pg, ref, alt = arg.split(":")
        sitios.append((nome, int(pg) - DESLOC, ref, alt))
    (na, pa, ra, aa), (nb, pb, rb, ab) = sitios
    lo, hi = min(pa, pb), max(pa, pb)
    cmd = ["samtools", "view", "-F", "0xF04", "-q", "20", bam, f"{CONTIG}:{lo}-{hi}"]
    out = subprocess.run(cmd, capture_output=True, text=True).stdout
    combos = collections.Counter()
    so_um = 0
    for linha in out.splitlines():
        c = linha.split("\t")
        inicio, cigar, seq = int(c[3]), c[5], c[9]
        ba, bb = base_em(pa, inicio, cigar, seq), base_em(pb, inicio, cigar, seq)
        if ba is None or bb is None:
            so_um += 1
            continue

        def rot(b, ref, alt, nome):
            if b == alt or (alt == "del" and b == "-"):
                return f"{nome}=ALT"
            if b == ref:
                return f"{nome}=ref"
            return f"{nome}=?{b}"
        combos[(rot(ba, ra, aa, na), rot(bb, rb, ab, nb))] += 1
    print(f"{bam.rsplit('/', 1)[-1]}  {na}@{pa + DESLOC} x {nb}@{pb + DESLOC}")
    print(f"  leituras cobrindo os dois sitios: {sum(combos.values())}  (cobrindo so um: {so_um})")
    for (x, y), n in combos.most_common():
        print(f"    {x:<14} {y:<14} {n:>4}")
    cis = combos.get((f"{na}=ALT", f"{nb}=ALT"), 0)
    trans = combos.get((f"{na}=ALT", f"{nb}=ref"), 0) + combos.get((f"{na}=ref", f"{nb}=ALT"), 0)
    veredicto = "CIS" if cis > 2 * trans else ("TRANS" if trans > 2 * cis else "AMBIGUO")
    print(f"  leituras ALT+ALT={cis}  ALT+ref ou ref+ALT={trans}  -> {veredicto}")


if __name__ == "__main__":
    main()
