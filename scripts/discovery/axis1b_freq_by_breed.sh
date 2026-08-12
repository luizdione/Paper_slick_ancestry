#!/bin/bash
# EIXO 1b — frequencia por raca + genotipos dos portadores (corrige bug do 'column')
set -u
META=/data/1000bg_metadata
R=/proj/results/axis1

# Candidatos: frameshift + indels in-frame (todos truncantes/alteradores no CDS)
awk -F'\t' 'NR==1 || $10=="FRAMESHIFT" || $10=="indel_in_frame"' \
  "$R/prlr_cds_variants_classified.tsv" > "$R/prlr_INDEL_candidates.tsv"
echo "=== Candidatos indel no CDS ==="
cat "$R/prlr_INDEL_candidates.tsv"

# regions file para bcftools -R
awk -F'\t' 'NR>1{print $1"\t"$2}' "$R/prlr_INDEL_candidates.tsv" | sort -u -k2,2n > "$R/_cand.rgn"
echo ""
echo "posicoes candidatas: $(wc -l < "$R/_cand.rgn")"

bcftools view -R "$R/_cand.rgn" -Oz -o "$R/prlr_cand.vcf.gz" "$R/prlr_cds.vcf.gz" 2>/dev/null
bcftools index -f "$R/prlr_cand.vcf.gz" 2>/dev/null
echo "variantes extraidas: $(bcftools view -H "$R/prlr_cand.vcf.gz" 2>/dev/null | wc -l)"

# ---------- frequencia por raca ----------
OUT="$R/indel_af_by_breed.tsv"
printf "BREED\tN\tCHROM\tPOS\tREF\tALT\tAC\tAN\tAF\n" > "$OUT"
for f in "$META"/samples_*.txt; do
  b=$(basename "$f"); b=${b#samples_}; b=${b%.txt}
  n=$(grep -c . "$f" 2>/dev/null); [ "${n:-0}" -ge 1 ] || continue
  bcftools view -S "$f" --force-samples "$R/prlr_cand.vcf.gz" 2>/dev/null \
    | bcftools +fill-tags -- -t AF,AC,AN 2>/dev/null \
    | bcftools query -f "${b}\t${n}\t%CHROM\t%POS\t%REF\t%ALT\t%AC\t%AN\t%AF\n" 2>/dev/null >> "$OUT"
done
echo ""
echo "=== RACAS PORTADORAS (AF>0) ==="
awk -F'\t' 'NR>1 && $9!="." && $9+0>0' "$OUT" \
  | sort -k4,4n -k9,9gr \
  | awk -F'\t' '{printf "%-24s n=%-5s %s:%-9s %s>%-14s AC=%-3s AN=%-5s AF=%.4f\n",$1,$2,$3,$4,$5,$6,$7,$8,$9}'

# ---------- genotipos individuais ----------
echo ""
echo "=== PORTADORES INDIVIDUAIS ==="
printf "SAMPLE\tCHROM\tPOS\tREF\tALT\tGT\tDP\tGQ\tBREED\n" > "$R/indel_carriers.tsv"
bcftools query -f '[%SAMPLE\t%CHROM\t%POS\t%REF\t%ALT\t%GT\t%DP\t%GQ\n]' "$R/prlr_cand.vcf.gz" 2>/dev/null \
  | awk -F'\t' '$6!="0/0" && $6!="./." && $6!="0|0" && $6!="."' > "$R/_c.tsv"
while IFS=$'\t' read -r s c p rf al gt dp gq; do
  br=$(grep -l -x -- "$s" "$META"/samples_*.txt 2>/dev/null | head -1)
  br=${br##*/samples_}; br=${br%.txt}
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$s" "$c" "$p" "$rf" "$al" "$gt" "$dp" "$gq" "${br:-NA}"
done < "$R/_c.tsv" >> "$R/indel_carriers.tsv"
rm -f "$R/_c.tsv"
awk -F'\t' 'NR>1{printf "%-16s %s:%-9s %s>%-14s GT=%-4s DP=%-4s GQ=%-4s %s\n",$1,$2,$3,$4,$5,$6,$7,$8,$9}' "$R/indel_carriers.tsv"
echo ""
echo "total portadores: $(( $(wc -l < "$R/indel_carriers.tsv") - 1 ))"
