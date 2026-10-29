#!/bin/bash
# 03_merge_fase.sh -- monta o VCF da janela de 10 Mb (Run 8 + BAMs genotipados + sitios focais) e
# faseia com BEAGLE 5.4 usando o mapa CLARITY (Dm_cM). Rodar no container de ferramentas (ver README) com o diretorio de dados em /data.
#
# Sitios focais vem do VCF da atribuicao por haplotipo (merged_phased_dedup), com as tres
# correcoes de leitura de `tables/derived/suporte_leitura_portadores.tsv`: SLICK3 em ERR12035159 e SLICK6 em ERR12035135
# viram 0/0 (leituras nao sustentam o alternativo); SLICK2 em ERR12035125 vira ausente (4 ref /
# 2 alt, profundidade baixa demais para decidir). ERR12035133 sai inteiro (4% de chamada).
set -euo pipefail
D=/data/slick_ancestry; O=$D/lai10mb; T=$O/tmp; mkdir -p $T
PAINEL=$D/criollo/incorporado/merged_phased_dedup.vcf.gz
RUN8=$D/haplo/bta20_win10Mb.vcf.gz
FOCAIS=20:39099213,20:39099226,20:39099228,20:39099267,20:39099321
MIN_CR=20   # % de chamada minima para entrar no faseamento

# --- amostras genotipadas que entram (ERR12035105 foi o teste, sem linha de estatistica)
awk -F'\t' -v m=$MIN_CR '$4>=m{print $1}' $O/genotipa_stats.tsv > $T/entram.txt
echo ERR12035105 >> $T/entram.txt
sort -u -o $T/entram.txt $T/entram.txt
grep -Fxf $O/novos.txt $T/entram.txt > $T/novos_ok.txt || true
grep -vFxf $O/novos.txt $T/entram.txt > $T/ibericos_ok.txt || true
echo "novos no faseamento: $(wc -l < $T/novos_ok.txt) | ibericos: $(wc -l < $T/ibericos_ok.txt)"

# --- sitios focais, GT nao faseado, com as correcoes de leitura
bcftools view -r $FOCAIS $PAINEL -Ou | bcftools annotate -x INFO,^FORMAT/GT -Ov \
 | awk 'BEGIN{FS=OFS="\t"
          fix["39099226|ERR12035159"]="0/0"; fix["39099267|ERR12035135"]="0/0";
          fix["39099321|ERR12035125"]="./."}
        /^##/{print; next}
        /^#CHROM/{for(i=10;i<=NF;i++) nome[i]=$i; print; next}
        {for(i=10;i<=NF;i++){gsub(/\|/,"/",$i); if($i=="1/0") $i="0/1";
                             k=$2"|"nome[i]; if(k in fix) $i=fix[k]}
         print}' | bgzip -c > $T/focal_todos.vcf.gz
bcftools index -f $T/focal_todos.vcf.gz

# --- bloco A: Run 8 (so GT) + focais nas mesmas amostras e ordem
bcftools annotate -x INFO,^FORMAT/GT $RUN8 -Oz -o $T/A.vcf.gz && bcftools index -f $T/A.vcf.gz
bcftools query -l $T/A.vcf.gz > $T/A.samples
bcftools view -S $T/A.samples $T/focal_todos.vcf.gz -Oz -o $T/A_focal.vcf.gz && bcftools index -f $T/A_focal.vcf.gz
bcftools concat -a $T/A.vcf.gz $T/A_focal.vcf.gz -Ou | bcftools sort -Oz -o $T/A2.vcf.gz && bcftools index -f $T/A2.vcf.gz

# --- bloco B: animais novos do painel (Criollo + N'Dama) + focais
mapfile -t BV < <(sed "s|^|$O/gt/|; s|$|.vcf.gz|" $T/novos_ok.txt)
bcftools merge -m none "${BV[@]}" -Oz -o $T/B.vcf.gz && bcftools index -f $T/B.vcf.gz
bcftools query -l $T/B.vcf.gz > $T/B.samples
bcftools view -S $T/B.samples $T/focal_todos.vcf.gz -Oz -o $T/B_focal.vcf.gz && bcftools index -f $T/B_focal.vcf.gz
bcftools concat -a $T/B.vcf.gz $T/B_focal.vcf.gz -Ou | bcftools sort -Oz -o $T/B2.vcf.gz && bcftools index -f $T/B2.vcf.gz

# --- bloco C: ibericos (sem genotipo focal: ausente, o BEAGLE preenche)
mapfile -t CV < <(sed "s|^|$O/gt/|; s|$|.vcf.gz|" $T/ibericos_ok.txt)
if (( ${#CV[@]} > 1 )); then
  bcftools merge -m none "${CV[@]}" -Oz -o $T/C.vcf.gz
else
  cp "${CV[0]}" $T/C.vcf.gz
fi
bcftools index -f $T/C.vcf.gz

# --- tudo junto
bcftools merge -m none $T/A2.vcf.gz $T/B2.vcf.gz $T/C.vcf.gz -Oz -o $O/janela10mb_nao_faseado.vcf.gz
bcftools index -f $O/janela10mb_nao_faseado.vcf.gz
echo "nao faseado: $(bcftools query -l $O/janela10mb_nao_faseado.vcf.gz | wc -l) amostras x" \
     "$(bcftools view -H $O/janela10mb_nao_faseado.vcf.gz | wc -l) sitios"

# --- mapa genetico CLARITY (Dm_cM), formato PLINK; o BEAGLE interpola entre as ancoras
awk -F'\t' 'NR>1 && $2=="20" && $5!="NA"{print "20\t"$1"\t"$5"\t"$3}' \
  $D/p5_ancestry/clarity/clarity_angus_map.tsv | sort -k4,4n > $O/clarity_bta20.map
echo "ancoras do mapa: $(wc -l < $O/clarity_bta20.map)"

# --- faseamento
java -Xmx10g -jar $D/beagle.jar gt=$O/janela10mb_nao_faseado.vcf.gz map=$O/clarity_bta20.map \
  out=$O/janela10mb_faseado nthreads=6 seed=20261026 > $O/beagle.out 2>&1
bcftools index -f $O/janela10mb_faseado.vcf.gz
echo "faseado: $(bcftools view -H $O/janela10mb_faseado.vcf.gz | wc -l) sitios"
tail -5 $O/janela10mb_faseado.log
