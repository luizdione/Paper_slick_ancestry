#!/bin/bash
# 04c_rfmix_5vias.sh -- RFMix com um 5o painel, TAURINO ASIATICO (Mongol, Cazaque, Hanwoo, Yanbian,
# Tibetano; fracao indicina <= 0,10). Motivo: o vizinho mais proximo do haplotipo do SLICK2 a
# +-250 kb e Mongol, populacao ausente das 4 referencias; o RFMix de 4 vias o chamou de INDICINO.
# Se a chamada migrar para ASIAN_TAURINE, o "indicino" era painel faltante.
set -euo pipefail
D=/data/slick_ancestry; O=$D/lai10mb; R=$O/rfmix
V=$O/janela10mb_faseado.vcf.gz
cut -f1 $O/ref_5vias.tsv > $R/ref5.ids
bcftools view -S $R/ref5.ids --force-samples $V -Oz -o $R/ref5.vcf.gz && bcftools index -f $R/ref5.vcf.gz
/t/rfmix2_bin/rfmix -f $R/q_principal.vcf.gz -r $R/ref5.vcf.gz -m $O/ref_5vias.tsv -g $R/mapa_rfmix.tsv \
  -o $R/m5 --chromosome=20 -G 50 --n-threads=6 --random-seed=20261026 > $R/m5.log 2>&1
echo "m5 feito"
