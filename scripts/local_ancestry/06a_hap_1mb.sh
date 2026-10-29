#!/bin/bash
# 06a_hap_1mb.sh -- haplotipos faseados em 1 Mb em volta do PRLR, entrada do 06_vizinho.py.
# Rodar no container de ferramentas (ver README) com o diretorio de dados em /data.
#
# Saida: uma linha por sitio, POS seguido do genotipo faseado (0|1) de cada amostra, na ordem do
# VCF; a primeira linha traz os nomes das amostras. Janela 20:38.599.000-39.600.000, que contem
# os +-250 kb que o 06_vizinho.py compara em volta dos sitios focais.
set -euo pipefail
D=/data/slick_ancestry; O=$D/lai10mb; R=$O/rfmix
V=$O/janela10mb_faseado.vcf.gz
{ printf "POS"; bcftools query -l $V | awk '{printf "\t%s", $1}'; echo
  bcftools query -r 20:38599000-39600000 -f '%POS[\t%GT]\n' $V; } | gzip -c > $R/hap_1mb_prlr.tsv.gz
echo "== feito: $(zcat $R/hap_1mb_prlr.tsv.gz | wc -l) linhas, cabecalho incluido"
