#!/bin/bash
# 01_sitios.sh -- sitios do Run 8 na janela de 10 Mb, em coordenada GLOBAL (para o merge) e LOCAL
# da isca (para o mpileup dos BAMs). Rodar no container de ferramentas (ver README) com o diretorio de dados em /data.
set -euo pipefail
D=/data/slick_ancestry
OUT=$D/lai10mb
mkdir -p "$OUT"
V=$D/haplo/bta20_win10Mb.vcf.gz
CTG="20:34100000-44100000"; OFF=34099999

bcftools view -m2 -M2 -v snps "$V" -Ou | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' \
  | bgzip -c > "$OUT/sitios_global.tsv.gz"
tabix -f -s1 -b2 -e2 "$OUT/sitios_global.tsv.gz"
zcat "$OUT/sitios_global.tsv.gz" | awk -v c="$CTG" -v o=$OFF 'BEGIN{OFS="\t"}{print c,$2-o,$3,$4}' \
  | bgzip -c > "$OUT/sitios_local.tsv.gz"
tabix -f -s1 -b2 -e2 "$OUT/sitios_local.tsv.gz"
echo "sitios: $(zcat "$OUT/sitios_global.tsv.gz" | wc -l)"
