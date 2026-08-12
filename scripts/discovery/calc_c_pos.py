#!/usr/bin/env python3
"""Deriva a posicao c. de BTA20:39.099.213 no transcrito canonico do PRLR
(ENSBTAT00000014437, PRLR-201, ARS-UCD1.3 / Ensembl 113). Tudo verificavel."""
import subprocess, sys

FA = "/data/genome_reference/ARS-UCD1.3/Bos_taurus.ARS-UCD1.3.dna.toplevel.fa"

# CDS do transcrito canonico, 1-based inclusivo, fita +, na ordem do GTF
CDS = [(39077947, 39078016),
       (39080773, 39080905),
       (39082809, 39082978),
       (39093142, 39093311),
       (39094967, 39095108),
       (39096738, 39096837),
       (39097533, 39097602),
       (39098688, 39099575)]
STOP = (39099576, 39099578)
VAR_ANCHOR = 39099213   # posicao do G no REF=GC / ALT=G (VCF)
VAR_DEL    = 39099214   # base efetivamente deletada (o C)


def faidx(chrom, s, e):
    out = subprocess.run(["samtools", "faidx", FA, f"{chrom}:{s}-{e}"],
                         capture_output=True, text=True, check=True).stdout
    return "".join(l.strip() for l in out.splitlines() if not l.startswith(">")).upper()


def g2c(pos):
    """posicao genomica -> posicao c. (CDS), fita +. None se intronico/fora."""
    cum = 0
    for i, (s, e) in enumerate(CDS, 1):
        if s <= pos <= e:
            return cum + (pos - s + 1), i, cum
        cum += e - s + 1
    return None, None, None


def c2g(cpos):
    """posicao c. -> genomica."""
    cum = 0
    for i, (s, e) in enumerate(CDS, 1):
        L = e - s + 1
        if cum < cpos <= cum + L:
            return s + (cpos - cum - 1), i
        cum += L
    return None, None


CODON = {}
bases = "TCAG"
aas = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
for i, b1 in enumerate(bases):
    for j, b2 in enumerate(bases):
        for k, b3 in enumerate(bases):
            CODON[b1 + b2 + b3] = aas[i * 16 + j * 4 + k]


def translate(seq):
    prot = []
    for i in range(0, len(seq) - 2, 3):
        prot.append(CODON.get(seq[i:i + 3], "X"))
    return "".join(prot)


print("=" * 78)
print("1) TRANSCRITO CANONICO: ENSBTAT00000014437 (PRLR-201), tag Ensembl_canonical")
print("=" * 78)
cum = 0
print(f"{'bloco':>5} {'inicio':>10} {'fim':>10} {'len':>5} {'c.inicio':>9} {'c.fim':>9}")
for i, (s, e) in enumerate(CDS, 1):
    L = e - s + 1
    print(f"{i:>5} {s:>10} {e:>10} {L:>5} {cum+1:>9} {cum+L:>9}")
    cum += L
TOT = cum
print(f"CDS total = {TOT} pb  ({TOT/3:.2f} codons)  + stop {STOP[0]}-{STOP[1]}")

print()
print("=" * 78)
print("2) CONTA DA POSICAO c. DE 39.099.213")
print("=" * 78)
cpos, blk, cumbefore = g2c(VAR_ANCHOR)
print(f"39.099.213 cai no bloco CDS #{blk} ({CDS[blk-1][0]}-{CDS[blk-1][1]})")
print(f"soma dos CDS anteriores = {cumbefore} pb")
print(f"offset dentro do bloco  = 39099213 - {CDS[blk-1][0]} + 1 = {VAR_ANCHOR-CDS[blk-1][0]+1}")
print(f"=> c.{cumbefore} + {VAR_ANCHOR-CDS[blk-1][0]+1} = c.{cpos}   (base ancora G)")
cpos_del, blk2, _ = g2c(VAR_DEL)
print(f"base DELETADA (o C, 39.099.214) => c.{cpos_del}")

print()
print("=" * 78)
print("3) VERIFICACAO DE SEQUENCIA")
print("=" * 78)
cds_seq = "".join(faidx(20, s, e) for s, e in CDS)
stop_seq = faidx(20, *STOP)
print(f"len(CDS montado) = {len(cds_seq)}  (esperado {TOT})")
print(f"primeiros 12 nt : {cds_seq[:12]}   -> inicia com ATG? {cds_seq[:3]=='ATG'}")
print(f"ultimos 12 nt   : {cds_seq[-12:]}")
print(f"stop_codon anotado ({STOP[0]}-{STOP[1]}) = {stop_seq} -> e stop? {stop_seq in ('TAA','TAG','TGA')}")
prot = translate(cds_seq)
print(f"proteina = {len(prot)} aa; sem stop interno? {'*' not in prot}")
print(f"prot[:20]  = {prot[:20]}")
print(f"prot[-20:] = {prot[-20:]}")
print()
print(f"base em 39.099.213 = {faidx(20, VAR_ANCHOR, VAR_ANCHOR)}  (esperado G)")
print(f"base em 39.099.214 = {faidx(20, VAR_DEL, VAR_DEL)}  (esperado C)")
print(f"cds_seq[c.{cpos}] = {cds_seq[cpos-1]} ; cds_seq[c.{cpos_del}] = {cds_seq[cpos_del-1]}")
win = faidx(20, VAR_ANCHOR - 20, VAR_ANCHOR + 20)
print()
print("41 pb genomicos 39.099.193 - 39.099.233 (ancora ^ em 39.099.213, deletada v):")
print(win)
print(" " * 20 + "^v")
print()
print("contexto c. (c.1362-c.1402 do CDS montado):")
print(cds_seq[1361:1402])
print(" " * 20 + "^v")

print()
print("=" * 78)
print("4) POSICAO GENOMICA DE c.1382 (SLICK1) E DISTANCIA")
print("=" * 78)
g1382, b1382 = c2g(1382)
print(f"c.1382 -> genomico {g1382:,} (bloco CDS #{b1382}); base = {faidx(20,g1382,g1382)}")
print(f"c.1381 -> genomico {c2g(1381)[0]:,} ; base = {faidx(20,*(c2g(1381)[0],)*2)}")
print(f"nossa base deletada = {VAR_DEL:,} (c.{cpos_del})")
print(f"distancia |c.1382 - base deletada| = {abs(g1382 - VAR_DEL)} pb")

print()
print("=" * 78)
print("5) EFEITO NO QUADRO DE LEITURA (delecao de 1 pb em c.1382)")
print("=" * 78)
mut = cds_seq[:cpos_del - 1] + cds_seq[cpos_del:]  # remove a base c.1382
# a delecao e no ultimo exon codificante -> continuar lendo pelo genoma alem do stop
downstream = faidx(20, CDS[-1][1] + 1, CDS[-1][1] + 3000)  # 3' apos ultimo CDS (mesmo exon)
mut_ext = mut + downstream
prot_mut = translate(mut_ext)
stop_idx = prot_mut.find("*")
print(f"CDS mutante ate o fim do CDS-WT: {len(mut)} nt (perdeu 1)")
codon_afetado = (cpos_del - 1) // 3 + 1
print(f"c.1382 esta no codon {codon_afetado} (aa {codon_afetado}); aa WT = {prot[codon_afetado-1]}")
print(f"proteina mutante ate o 1o stop: {stop_idx} aa; stop no codon {stop_idx+1}")
print(f"prot WT   {codon_afetado-6}-{codon_afetado+14}: {prot[codon_afetado-7:codon_afetado+14]}")
print(f"prot MUT  {codon_afetado-6}-{codon_afetado+14}: {prot_mut[codon_afetado-7:codon_afetado+14]}")
print(f"prot MUT trecho ate stop (a partir do codon {codon_afetado}): {prot_mut[codon_afetado-1:stop_idx+1]}")
# posicao c. e genomica do novo stop
new_stop_c = stop_idx * 3 + 1  # 1a base do codon stop, coordenada no CDS MUTANTE
print(f"1a base do novo stop = posicao {new_stop_c} do CDS mutante")
# converter para coordenada genomica: no CDS mutante, apos a delecao, indice i (1-based)
# corresponde a c.(i+1) do WT para i >= 1382
wt_equiv = new_stop_c + 1 if new_stop_c >= cpos_del else new_stop_c
if wt_equiv <= TOT:
    g_new_stop, _ = c2g(wt_equiv)
else:
    g_new_stop = CDS[-1][1] + (wt_equiv - TOT)
print(f"=> equivale a c.{wt_equiv} do WT -> genomico {g_new_stop:,}")
print(f"novo stop esta {'DENTRO' if wt_equiv<=TOT else 'ALEM'} do CDS canonico (CDS termina em c.{TOT} = {CDS[-1][1]:,})")
print(f"aa perdidos vs WT: {len(prot)} - {stop_idx} = {len(prot)-stop_idx}")
