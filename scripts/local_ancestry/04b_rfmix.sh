#!/bin/bash
# 04b_rfmix.sh -- ancestralidade local com RFMix v2.03 sobre o faseamento da janela de 10 Mb.
# Rodar no container de ferramentas (ver README) com o diretorio de dados em /data e as
# ferramentas em /t.
#
# Execucoes:
#   m4  : 4 vias (EUR_TAURINE, IBERIAN, NDAMA, INDICINE), consulta = Criollo + africanos + controles
#   m3  : 3 vias (TAURINO_EUROPEU, TAURINO_AFRICANO, INDICINE), mesma consulta
#   c4  : controle dos paineis pequenos: 3 N'Dama e 2 ibericos FORA da referencia, como consulta
#   g15, g150 : sensibilidade ao numero de geracoes desde a mistura (-G), 4 vias, so a consulta
# -G 50 na principal: Criollo formado a partir de ~1500, intervalo de geracao ~5 anos (~100
# geracoes), mas parte da mistura africana e anterior, ja na Iberia; 50 fica no meio e o -G e
# varrido na sensibilidade.
set -euo pipefail
D=/data/slick_ancestry; O=$D/lai10mb; R=$O/rfmix; mkdir -p $R
RFMIX=/t/rfmix2_bin/rfmix
V=$O/janela10mb_faseado.vcf.gz
SEED=20261026

# mapa no formato do RFMix: cromossomo, bp, cM (sem cabecalho)
awk 'BEGIN{OFS="\t"}{print $1,$4,$3}' $O/clarity_bta20.map > $R/mapa_rfmix.tsv

cut -f1 $O/consulta.tsv > $R/consulta.ids
cut -f1 $O/controle_grandes.tsv > $R/ctrl_grandes.ids
cut -f1 $O/controle_pequenos.tsv > $R/ctrl_pequenos.ids
cat $R/consulta.ids $R/ctrl_grandes.ids > $R/q_principal.ids

sub(){ bcftools view -S $1 --force-samples $V -Oz -o $2 && bcftools index -f $2; }
sub $R/q_principal.ids $R/q_principal.vcf.gz
sub $R/consulta.ids    $R/q_consulta.vcf.gz
sub $R/ctrl_pequenos.ids $R/q_ctrl_pequenos.vcf.gz
cut -f1 $O/ref_4vias.tsv > $R/ref.ids
sub $R/ref.ids $R/ref.vcf.gz
grep -vFwf $R/ctrl_pequenos.ids $O/ref_4vias.tsv > $R/ref_4vias_sem_pequenos.tsv
cut -f1 $R/ref_4vias_sem_pequenos.tsv > $R/ref_c4.ids
sub $R/ref_c4.ids $R/ref_c4.vcf.gz

roda(){  # nome query ref mapa_amostras G
  echo "== $1 (G=$5)"
  $RFMIX -f $2 -r $3 -m $4 -g $R/mapa_rfmix.tsv -o $R/$1 --chromosome=20 -G $5 \
         --n-threads=6 --random-seed=$SEED > $R/$1.log 2>&1
  tail -2 $R/$1.log
}
roda m4   $R/q_principal.vcf.gz     $R/ref.vcf.gz    $O/ref_4vias.tsv                50
roda m3   $R/q_principal.vcf.gz     $R/ref.vcf.gz    $O/ref_3vias.tsv                50
roda c4   $R/q_ctrl_pequenos.vcf.gz $R/ref_c4.vcf.gz $R/ref_4vias_sem_pequenos.tsv   50
roda g15  $R/q_consulta.vcf.gz      $R/ref.vcf.gz    $O/ref_4vias.tsv                15
roda g150 $R/q_consulta.vcf.gz      $R/ref.vcf.gz    $O/ref_4vias.tsv                150

# genotipos FASEADOS dos sitios focais e da variante do exon 2, por haplotipo
{ printf "POS"; bcftools query -l $V | awk '{printf "\t%s", $1}'; echo
  bcftools query -r 20:39077999,20:39078000,20:39099213,20:39099226,20:39099228,20:39099267,20:39099321 \
    -f '%POS[\t%GT]\n' $V; } > $R/focais_faseados.tsv
echo "feito"
