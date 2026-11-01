#!/bin/bash
# 04d_rfmix_3vias_sensibilidade.sh -- sensibilidade da configuracao de 3 vias ao numero de geracoes
# desde a mistura (-G). Roda depois de 04b_rfmix.sh, cujas entradas reaproveita (consulta principal
# com os controles grandes, referencia e mapa), com a mesma semente e o mesmo numero de threads.
#
# Execucoes:
#   m3_g15, m3_g150 : 3 vias (TAURINO_EUROPEU, TAURINO_AFRICANO, INDICINE), -G 15 e -G 150,
#                     mesma consulta da execucao m3 (Criollo + africanos + controles grandes)
set -euo pipefail
D=/data/slick_ancestry; O=$D/lai10mb; R=$O/rfmix
RFMIX=/t/rfmix2_bin/rfmix
SEED=20261026
for f in $R/q_principal.vcf.gz $R/ref.vcf.gz $R/mapa_rfmix.tsv $O/ref_3vias.tsv; do
  [ -s "$f" ] || { echo "falta $f: rode 04b_rfmix.sh antes"; exit 1; }
done

roda(){  # nome query ref mapa_amostras G
  echo "== $1 (G=$5)"
  $RFMIX -f $2 -r $3 -m $4 -g $R/mapa_rfmix.tsv -o $R/$1 --chromosome=20 -G $5 \
         --n-threads=6 --random-seed=$SEED > $R/$1.log 2>&1
  tail -2 $R/$1.log
}
roda m3_g15  $R/q_principal.vcf.gz $R/ref.vcf.gz $O/ref_3vias.tsv 15
roda m3_g150 $R/q_principal.vcf.gz $R/ref.vcf.gz $O/ref_3vias.tsv 150
echo "feito"
