#!/bin/bash
# 02_roda.sh -- genotipa os 38 animais novos do painel (31 Criollo + 7 N'Dama) e os WGS ibericos
# nao vazios na janela de 10 Mb. Resumivel: pula saida existente.
set -uo pipefail
D=/data/slick_ancestry; O=$D/lai10mb; mkdir -p $O/gt
FA=$D/criollo/bait/win10Mb.fa
bcftools query -l $D/criollo/incorporado/merged_phased_dedup.vcf.gz | sort > $O/.painel.txt
bcftools query -l $D/haplo/bta20_win10Mb.vcf.gz | sort > $O/.run8.txt
comm -13 $O/.run8.txt $O/.painel.txt > $O/novos.txt
: > $O/.lista
while read s; do
  for dir in bam bam_ndama; do [ -s $D/criollo/$dir/$s.bam ] && echo "$D/criollo/$dir/$s.bam $s" >> $O/.lista && break; done
done < $O/novos.txt
for s in ERR4280259 ERR4280260 ERR4280261 ERR4280262 ERR4280263 ERR4280270 ERR4280271 ERR3305588; do
  echo "$D/criollo/bam_iberian/$s.bam $s" >> $O/.lista
done
echo "BAMs na lista: $(wc -l < $O/.lista) (novos no painel: $(wc -l < $O/novos.txt))"
cat $O/.lista | xargs -P 4 -L 1 bash -c '[ -s '$O'/gt/$1.vcf.gz ] || { [ -s $0.bai ] || samtools index $0; python3 /s/02_genotipa.py '$FA' '$O'/sitios_local.tsv.gz $0 '$O'/gt/$1.vcf.gz $1 >> '$O'/genotipa_stats.tsv; }'
echo "feito: $(ls $O/gt/*.vcf.gz | wc -l) VCFs"
