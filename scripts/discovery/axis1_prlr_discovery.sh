#!/bin/bash
# ============================================================================
# EIXO 1 — Descoberta nao-dirigida de variantes de truncamento do PRLR
# Origem dos alelos SLICK / rota taurina tropical
# Entrada : coorte de descoberta (montada em /data) | Saida: resultados do eixo (montado em /proj)
# ============================================================================
set -u
GTF=/data/genome_reference/ARS-UCD1.3/Bos_taurus.ARS-UCD1.3.113.gtf
VCF=/data/1000bg_vcf/Chr20-Run8-TAUIND-public.vcf.gz
META=/data/1000bg_metadata
R=/proj/results/axis1
mkdir -p "$R"

GENE_CHR=20; GENE_START=39037656; GENE_END=39117082
REG="${GENE_CHR}:${GENE_START}-${GENE_END}"

echo "############ EIXO 1 — PRLR truncation discovery ############"
date

# ---------- 1. Intervalos CDS do PRLR (todos os transcritos, unificados) ----------
echo ">> [1] Extraindo CDS do PRLR do GTF"
awk -F'\t' -v OFS='\t' '$3=="CDS" && $9 ~ /gene_name "PRLR"/ {print $1,$4-1,$5}' "$GTF" \
  | sort -k1,1 -k2,2n | bedtools merge -i - 2>/dev/null > "$R/prlr_cds.bed" \
  || awk -F'\t' -v OFS='\t' '$3=="CDS" && $9 ~ /gene_name "PRLR"/ {print $1,$4-1,$5}' "$GTF" \
     | sort -k1,1 -k2,2n -u > "$R/prlr_cds.bed"
echo "   intervalos CDS: $(wc -l < "$R/prlr_cds.bed")"
awk '{s+=$3-$2} END{print "   bases codificantes: "s}' "$R/prlr_cds.bed"

# bed 1-based para bcftools -T
awk -v OFS='\t' '{print $1,$2+1,$3}' "$R/prlr_cds.bed" > "$R/prlr_cds.tsv"

# ---------- 2. Extrai o gene inteiro uma vez ----------
echo ">> [2] Extraindo regiao do gene PRLR ($REG)"
if [ ! -s "$R/prlr_gene.vcf.gz" ]; then
  bcftools view -r "$REG" -Oz -o "$R/prlr_gene.vcf.gz" "$VCF" 2>/dev/null
  bcftools index -f "$R/prlr_gene.vcf.gz" 2>/dev/null
fi
echo "   variantes no gene: $(bcftools view -H "$R/prlr_gene.vcf.gz" 2>/dev/null | wc -l)"

# ---------- 3. Restringe ao CDS e classifica ----------
echo ">> [3] Variantes no CDS + classificacao de frameshift"
bcftools view -T "$R/prlr_cds.tsv" -Oz -o "$R/prlr_cds.vcf.gz" "$R/prlr_gene.vcf.gz" 2>/dev/null
bcftools index -f "$R/prlr_cds.vcf.gz" 2>/dev/null
echo "   variantes no CDS: $(bcftools view -H "$R/prlr_cds.vcf.gz" 2>/dev/null | wc -l)"

# TIPO + delta de comprimento por alelo ALT -> frameshift se (len(ALT)-len(REF)) %% 3 != 0
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%FILTER\t%AC\t%AN\n' "$R/prlr_cds.vcf.gz" 2>/dev/null \
| awk -F'\t' -v OFS='\t' '
BEGIN{print "CHROM","POS","REF","ALT","QUAL","FILTER","AC","AN","DELTA","CLASSE"}
{
  n=split($4,alts,",");
  for(i=1;i<=n;i++){
    a=alts[i]; if(a=="*"||a=="<NON_REF>") continue;
    d=length(a)-length($3);
    cls=(d==0)?"substituicao":((d%3!=0)?"FRAMESHIFT":"indel_in_frame");
    print $1,$2,$3,a,$5,$6,$7,$8,d,cls
  }
}' > "$R/prlr_cds_variants_classified.tsv"

echo "   -- distribuicao por classe --"
awk -F'\t' 'NR>1{c[$10]++} END{for(k in c) printf "      %-18s %d\n",k,c[k]}' "$R/prlr_cds_variants_classified.tsv"

echo "   -- FRAMESHIFT candidatos (alelos SLICK potenciais) --"
awk -F'\t' 'NR==1||$10=="FRAMESHIFT"' "$R/prlr_cds_variants_classified.tsv" \
  | tee "$R/prlr_FRAMESHIFT_candidates.tsv" | column -t -s$'\t'

# ---------- 4. Anotacao funcional (SnpEff se disponivel) ----------
echo ">> [4] Anotacao funcional"
if command -v snpEff >/dev/null 2>&1; then
  snpEff -noStats ARS-UCD1.2.99 "$R/prlr_cds.vcf.gz" > "$R/prlr_cds.snpeff.vcf" 2>"$R/snpeff.log" \
    && echo "   SnpEff OK" || echo "   SnpEff falhou (ver snpeff.log) - classificacao por delta permanece valida"
else
  echo "   SnpEff ausente no container; usando classificacao por delta de comprimento (valida p/ frameshift)"
fi

# ---------- 5. Frequencia por raca de TODOS os frameshift ----------
echo ">> [5] Frequencia por raca (todas as racas com n>=1)"
POSLIST=$(awk -F'\t' 'NR>1{print $2}' "$R/prlr_FRAMESHIFT_candidates.tsv" | sort -un | tr '\n' ',' | sed 's/,$//')
echo "   posicoes frameshift: ${POSLIST:-nenhuma}"

OUTLONG="$R/frameshift_af_by_breed.tsv"
printf "BREED\tN_SAMPLES\tCHROM\tPOS\tREF\tALT\tAC\tAN\tAF\n" > "$OUTLONG"

if [ -n "${POSLIST:-}" ]; then
  bcftools view -i "POS=${POSLIST//,/ || POS=}" -Oz -o "$R/prlr_fs.vcf.gz" "$R/prlr_cds.vcf.gz" 2>/dev/null
  bcftools index -f "$R/prlr_fs.vcf.gz" 2>/dev/null

  for f in "$META"/samples_*.txt; do
    b=$(basename "$f"); b=${b#samples_}; b=${b%.txt}
    n=$(grep -c . "$f" 2>/dev/null)
    [ "${n:-0}" -ge 1 ] || continue
    bcftools view -S "$f" --force-samples "$R/prlr_fs.vcf.gz" 2>/dev/null \
      | bcftools +fill-tags -- -t AF,AC,AN 2>/dev/null \
      | bcftools query -f "${b}\t${n}\t%CHROM\t%POS\t%REF\t%ALT\t%AC\t%AN\t%AF\n" 2>/dev/null \
      >> "$OUTLONG"
  done
  echo "   linhas geradas: $(( $(wc -l < "$OUTLONG") - 1 ))"
  echo "   -- racas PORTADORAS (AF>0) --"
  awk -F'\t' 'NR>1 && $9!="." && $9+0>0 {printf "      %-24s n=%-4s %s:%s %s>%s  AC=%s AN=%s AF=%.4f\n",$1,$2,$3,$4,$5,$6,$7,$8,$9}' "$OUTLONG" \
    | sort -k1,1
fi

# ---------- 6. Genotipos individuais dos portadores ----------
echo ">> [6] Genotipos individuais dos portadores"
if [ -s "$R/prlr_fs.vcf.gz" ]; then
  printf "SAMPLE\tCHROM\tPOS\tREF\tALT\tGT\tDP\tGQ\tBREED\n" > "$R/frameshift_carriers.tsv"
  bcftools query -f '[%SAMPLE\t%CHROM\t%POS\t%REF\t%ALT\t%GT\t%DP\t%GQ\n]' "$R/prlr_fs.vcf.gz" 2>/dev/null \
    | awk -F'\t' '$6!="0/0" && $6!="./." && $6!="0|0" && $6!="."' > "$R/_carriers_raw.tsv"
  while IFS=$'\t' read -r s c p rf al gt dp gq; do
    br=$(grep -l -x -- "$s" "$META"/samples_*.txt 2>/dev/null | head -1)
    br=${br##*/samples_}; br=${br%.txt}
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$s" "$c" "$p" "$rf" "$al" "$gt" "$dp" "$gq" "${br:-?}"
  done < "$R/_carriers_raw.tsv" >> "$R/frameshift_carriers.tsv"
  rm -f "$R/_carriers_raw.tsv"
  echo "   portadores: $(( $(wc -l < "$R/frameshift_carriers.tsv") - 1 ))"
  column -t -s$'\t' "$R/frameshift_carriers.tsv" | head -40
fi

echo ""
echo "############ EIXO 1 CONCLUIDO ############"
date
echo "Saidas em $R"
ls -la "$R"
