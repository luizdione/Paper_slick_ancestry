#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""02_genotipa.py -- genotipa 1 BAM da isca de 10 Mb nos sitios do Run 8, alelos FORCADOS.

Mesmo principio do `criollo/genotipa_forcado.py` (o `call -C alleles` segfalta no bcftools
1.19), mas com o genotipo escolhido por verossimilhanca em vez de "qualquer leitura alternativa
vira heterozigoto": com erro de base e = 0,01,
    L(0/0) = (1-e)^r e^a,   L(0/1) = 0,5^(r+a),   L(1/1) = e^r (1-e)^a,
e o genotipo so e emitido com profundidade >= MIN_DP e vantagem de log10 >= MIN_LOD sobre o
segundo melhor. O resto fica ausente: o BEAGLE preenche no faseamento, e e melhor um ausente
que um heterozigoto falso, que o faseamento propaga.

Uso (no container de ferramentas (ver README)):
  02_genotipa.py <bait.fa> <sitios_local.tsv.gz> <BAM> <OUT.vcf.gz> <AMOSTRA>
Saida em coordenada GLOBAL (contig 20), na ordem e no esquema REF/ALT do Run 8.
"""
import gzip
import math
import os
import subprocess
import sys

FA, SITIOS, BAM, OUT, AMOSTRA = sys.argv[1:6]
OFF = 34099999
E, MIN_DP, MIN_LOD = 0.01, 4, 1.0
L1E, LE, LH = math.log10(1 - E), math.log10(E), math.log10(0.5)

alelos = []
with gzip.open(SITIOS, "rt") as f:
    for ln in f:
        c = ln.rstrip("\n").split("\t")
        alelos.append((int(c[1]), c[2], c[3]))

mp = subprocess.Popen(
    f"bcftools mpileup -f {FA} -T {SITIOS} -I -a AD -q 20 -Q 20 -Ou {BAM} 2>/dev/null"
    " | bcftools query -f '%POS\\t%REF\\t%ALT\\t[%AD]\\n'",
    shell=True, stdout=subprocess.PIPE, text=True)
ad = {}
for ln in mp.stdout:
    p, ref, alt, a = ln.rstrip("\n").split("\t")
    ad[int(p)] = (ref, alt.split(","), a)
mp.wait()

n = {"0/0": 0, "0/1": 0, "1/1": 0, "./.": 0}
tmp = OUT + ".txt"
with open(tmp, "w") as o:
    o.write("##fileformat=VCFv4.2\n##contig=<ID=20>\n"
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
            f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{AMOSTRA}\n")
    for p, pref, palt in alelos:
        gt = "./."
        if p in ad:
            ref, alts, a = ad[p]
            try:
                v = [int(x) for x in a.split(",")]
            except ValueError:
                v = []
            if v and ref == pref:
                r = v[0]
                x = v[alts.index(palt) + 1] if palt in alts and alts.index(palt) + 1 < len(v) else 0
                if r + x >= MIN_DP:
                    ll = {"0/0": r * L1E + x * LE, "0/1": (r + x) * LH, "1/1": r * LE + x * L1E}
                    o2 = sorted(ll.values(), reverse=True)
                    if o2[0] - o2[1] >= MIN_LOD:
                        gt = max(ll, key=ll.get)
        n[gt] += 1
        o.write(f"20\t{p + OFF}\t.\t{pref}\t{palt}\t.\t.\t.\tGT\t{gt}\n")
subprocess.run(f"bgzip -c {tmp} > {OUT} && bcftools index -f {OUT}", shell=True, check=True)
os.remove(tmp)
tot = sum(n.values())
print(f"{AMOSTRA}\t{tot}\t{tot - n['./.']}\t{100 * (tot - n['./.']) / tot:.1f}\t"
      f"{n['0/1']}\t{n['1/1']}")
