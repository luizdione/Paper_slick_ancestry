#!/bin/bash
# ============================================================================
# ITEM G - teste de haplotipo em torno de BTA20:39.099.214 (SLICK1)
#
# Objetivo: medir o comprimento do haplotipo compartilhado pelos portadores e
# compara-lo com haplotipos africanos e ibericos. Origem pos-1500 nas Americas
# e origem africana antiga produzem comprimentos distinguiveis - e so esse teste
# autoriza a palavra "origem" na descricao do resultado.
#
# Retomavel: cada etapa pula se a saida ja existe.
# ============================================================================
set -u
REF=/data/genome_reference/ARS-UCD1.3/Bos_taurus.ARS-UCD1.3.dna.toplevel.fa
VCF=/data/1000bg_vcf/Chr20-Run8-TAUIND-public.vcf.gz
R=/proj/results/haplo
COH=/proj/data/cohort
THREADS=${THREADS:-4}
WIN="20:34100000-44100000"     # +- 5 Mb em torno de 39.099.214
POS=39099214
mkdir -p "$R"

echo "########## ITEM G - haplotipo/idade do alelo SLICK1 ##########"
date

# ---------- 1. Janela do painel Run 8 ----------
if [ ! -s "$R/bta20_win10Mb.vcf.gz" ]; then
  echo ">> [1] Extraindo janela $WIN (medido: ~6 min, ~44 MB)"
  bcftools view --threads "$THREADS" -r "$WIN" -m2 -M2 -v snps -Ou "$VCF" \
    | bcftools +fill-tags -Ou -- -t AF,AC,AN,F_MISSING,MAF \
    | bcftools view -i 'MAF>=0.01 && F_MISSING<0.1' -Ou \
    | bcftools annotate --threads "$THREADS" -x '^FORMAT/GT' -Oz -o "$R/bta20_win10Mb.vcf.gz"
  bcftools index -f "$R/bta20_win10Mb.vcf.gz"
else
  echo ">> [1] janela ja extraida"
fi
echo "   SNPs: $(bcftools view -H "$R/bta20_win10Mb.vcf.gz" | wc -l)"

# ---------- 2. Incorporar as fatias WGS novas (criollo + oeste-africano) ----------
# As fatias guardam BTA20:34,1-44,1 Mb, exatamente a mesma janela.
if [ ! -s "$R/merged.vcf.gz" ]; then
  echo ">> [2] Chamando variantes nas fatias WGS e mesclando com o painel"
  ls /proj/data/wgs_slices/*/*.target.bam > "$R/bams.txt" 2>/dev/null || true
  NB=$(wc -l < "$R/bams.txt" 2>/dev/null || echo 0)
  echo "   BAMs de fatia: $NB"
  if [ "$NB" -gt 0 ]; then
    bcftools mpileup -f "$REF" -r "$WIN" -b "$R/bams.txt" \
        -a AD,DP -q 20 -Q 20 -L 1000000 -d 2000 --threads "$THREADS" -Ou \
      | bcftools call -m -Oz -o "$R/wgs_slices.vcf.gz" --threads "$THREADS"
    bcftools index -f "$R/wgs_slices.vcf.gz"
    bcftools merge --threads "$THREADS" -0 -Oz -o "$R/merged.vcf.gz" \
      "$R/bta20_win10Mb.vcf.gz" "$R/wgs_slices.vcf.gz"
  else
    echo "   [!] sem fatias WGS; seguindo so com o painel Run 8"
    cp "$R/bta20_win10Mb.vcf.gz" "$R/merged.vcf.gz"
  fi
  bcftools index -f "$R/merged.vcf.gz"
else
  echo ">> [2] merge ja feito"
fi
echo "   amostras no merge: $(bcftools query -l "$R/merged.vcf.gz" | wc -l)"

# ---------- 3. Faseamento (Beagle) ----------
BEAGLE=/opt/beagle.jar
if [ ! -s "$BEAGLE" ]; then
  echo ">> [3a] Baixando Beagle"
  for u in \
    "https://faculty.washington.edu/browning/beagle/beagle.27Feb25.75f.jar" \
    "https://faculty.washington.edu/browning/beagle/beagle.22Jul22.46e.jar" ; do
    wget -q -O "$BEAGLE" "$u" && break
  done
fi
if [ ! -s "$BEAGLE" ]; then
  echo "   [BLOQUEADO] nao foi possivel obter o Beagle. Instale manualmente em /opt/beagle.jar"
  echo "   (o estagio sera retomado na proxima passagem do supervisor)"
  exit 1
fi

if [ ! -s "$R/phased.vcf.gz" ]; then
  echo ">> [3b] Faseando (Beagle) - etapa mais longa, ~1-2 h"
  java -Xmx12g -jar "$BEAGLE" \
    gt="$R/merged.vcf.gz" out="$R/phased" nthreads="$THREADS" impute=false \
    > "$R/beagle.log" 2>&1 || { echo "   ERRO no Beagle"; tail -20 "$R/beagle.log"; exit 1; }
  bcftools index -f "$R/phased.vcf.gz"
else
  echo ">> [3] ja faseado"
fi

# ---------- 4. EHH / comprimento de haplotipo ----------
echo ">> [4] Comprimento do haplotipo compartilhado ao redor de $POS"
python3 - <<'PY' > "$R/ehh_summary.tsv"
import gzip, sys
from collections import defaultdict
V="/proj/results/haplo/phased.vcf.gz"; POS=39099214
samples=[]; rows=[]
with gzip.open(V,'rt') as f:
    for line in f:
        if line.startswith('##'): continue
        p=line.rstrip('\n').split('\t')
        if line.startswith('#CHROM'):
            samples=p[9:]; continue
        rows.append((int(p[1]), p[9:]))
if not rows: sys.exit("sem variantes")
# indice do sitio focal (mais proximo)
foc=min(range(len(rows)), key=lambda i: abs(rows[i][0]-POS))
def hap(gt,h):
    g=gt.split(':')[0]
    if '|' in g: a=g.split('|')
    else: a=g.split('/')
    return a[h] if len(a)>h else '.'
# haplotipos portadores vs nao portadores no sitio focal
carr=[]; ctrl=[]
for si,s in enumerate(samples):
    for h in (0,1):
        a=hap(rows[foc][1][si],h)
        (carr if a=='1' else ctrl).append((si,h))
print("grupo\tn_haplotipos\tcomprimento_mediano_kb\tcomprimento_medio_kb")
def decay(hs):
    # distancia ate o primeiro sitio em que o haplotipo difere da maioria
    lens=[]
    for si,h in hs:
        L=0
        for d in (1,-1):
            i=foc+d
            while 0<=i<len(rows):
                a=hap(rows[i][1][si],h)
                # compara com o alelo modal do grupo naquele sitio
                cnt=defaultdict(int)
                for sj,hj in hs: cnt[hap(rows[i][1][sj],hj)]+=1
                modal=max(cnt,key=cnt.get)
                if a!=modal: break
                i+=d
            L+=abs(rows[min(max(i,0),len(rows)-1)][0]-POS)
        lens.append(L/1000.0)
    lens.sort()
    if not lens: return (0,0)
    return (lens[len(lens)//2], sum(lens)/len(lens))
for nome,hs in (("portadores",carr),("controles",ctrl)):
    if hs:
        m,a=decay(hs)
        print(f"{nome}\t{len(hs)}\t{m:.1f}\t{a:.1f}")
PY
cat "$R/ehh_summary.tsv"

echo ""
echo "########## ITEM G CONCLUIDO ##########"
date
