#!/usr/bin/env bash
# incorpora_novos_portadores.sh — junta os criollos e N'Dama baixados do ENA ao
# painel de referencia, refaseia, e refaz a atribuicao de ancestralidade.
#
# Por que refasear e nao concatenar
# ---------------------------------
# Os BAM novos poderiam ser faseados sozinhos e os haplotipos resultantes colados
# ao lado dos do painel. Isso daria um arquivo que abre e uma analise que roda, e
# estaria errado: a fase de um conjunto e resolvida CONTRA os haplotipos daquele
# conjunto. Duas execucoes separadas produzem duas particoes que nao se comparam,
# e a atribuicao passaria a medir a diferenca entre elas. O merge tem de ocorrer
# no nivel de genotipo, com uma unica passada de faseamento sobre tudo.
#
# Coordenadas
# -----------
# Os BAM novos foram alinhados contra a **isca** (`win10Mb.fa`), cujo contig se
# chama `20:34100000-44100000` e cuja posicao 1 corresponde a 34.100.000 em
# ARS-UCD1.2. Antes do merge as posicoes precisam voltar para a referencia; somar
# WIN_START-1 e renomear o contig para `20` faz isso. Sem essa etapa o merge nao
# falha — ele simplesmente nao encontra sitios em comum, e o painel novo entra
# como um bloco de dados ausentes.
set -uo pipefail

D=/data/slick_ancestry/criollo
BAIT=$D/bait/win10Mb.fa
WIN_START=34100000
GENE_INI=39077946
GENE_FIM=39117082
PAINEL=/data/slick_ancestry/haplo/prlr_gene_allmaf.vcf.gz   # genotipos, nao faseado
OUT=$D/incorporado
THREADS=${THREADS:-6}
mkdir -p $OUT

CHR=$(head -1 "$BAIT" | sed 's/^>//' | cut -d' ' -f1)
P_INI=$((GENE_INI - WIN_START + 1))
P_FIM=$((GENE_FIM - WIN_START + 1))

mapfile -t NOVOS < <(ls $D/bam/*.bam $D/bam_ndama/*.bam 2>/dev/null)
(( ${#NOVOS[@]} )) || { echo "ERRO: nenhum BAM novo"; exit 1; }
echo "[1/5] ${#NOVOS[@]} BAM novos | gene na isca: ${P_INI}-${P_FIM}"

# ---------------------------------------------------------------- 1. chamar no gene
printf '%s\n' "${NOVOS[@]}" > $OUT/bamlist.txt
printf '%s\t%d\t%d\n' "$CHR" $((P_INI-1)) "$P_FIM" > $OUT/gene.bed
bcftools mpileup -f "$BAIT" -b $OUT/bamlist.txt -R $OUT/gene.bed \
    -a AD,DP -q 20 -Q 20 --threads $THREADS -Ou \
  | bcftools call -m -Oz -o $OUT/novos_isca.vcf.gz
bcftools index -f $OUT/novos_isca.vcf.gz
N1=$(bcftools view -H $OUT/novos_isca.vcf.gz | wc -l)
echo "[2/5] $N1 sitios chamados nos BAM novos"
(( N1 > 0 )) || { echo "GATE: zero sitios — conferir contig/BED"; exit 3; }

# ---------------------------------------------------------------- 2. isca -> ARS-UCD1.2
# renomeia o contig e desloca a posicao; o nome do BAM vira o nome da amostra
bcftools view -h $OUT/novos_isca.vcf.gz \
  | sed "s|^##contig=<ID=${CHR},|##contig=<ID=20,|" > $OUT/hdr.txt
bcftools view -H $OUT/novos_isca.vcf.gz \
  | awk -v off=$((WIN_START-1)) 'BEGIN{OFS="\t"} {$1="20"; $2=$2+off; print}' \
  > $OUT/body.txt
cat $OUT/hdr.txt $OUT/body.txt | bcftools view -Oz -o $OUT/novos_ars.vcf.gz
bcftools index -f $OUT/novos_ars.vcf.gz
# nomes de amostra legiveis (o VCF guarda o caminho do BAM)
bcftools query -l $OUT/novos_ars.vcf.gz | sed 's|.*/||; s|\.bam$||' > $OUT/novos_nomes.txt
bcftools reheader -s $OUT/novos_nomes.txt -o $OUT/novos.vcf.gz $OUT/novos_ars.vcf.gz
bcftools index -f $OUT/novos.vcf.gz
echo "[3/5] convertido para ARS-UCD1.2 e renomeado"

# ---------------------------------------------------------------- 3. merge com o painel
bcftools merge -m none -Oz -o $OUT/merged.vcf.gz "$PAINEL" $OUT/novos.vcf.gz
bcftools index -f $OUT/merged.vcf.gz
N2=$(bcftools view -H $OUT/merged.vcf.gz | wc -l)
NS=$(bcftools query -l $OUT/merged.vcf.gz | wc -l)
echo "[4/5] merge: $N2 sitios, $NS amostras"
# GATE: o merge tem de somar amostras, nao trocar
NP=$(bcftools query -l "$PAINEL" | wc -l); NN=$(wc -l < $OUT/novos_nomes.txt)
(( NS == NP + NN )) || echo "AVISO: $NS amostras contra $NP+$NN esperadas — conferir duplicatas"

# ---------------------------------------------------------------- 4. refasear TUDO junto
# o jar nao mora em /opt neste container; uma versao anterior do pipeline baixava em
# tempo de execucao e ficou com caminho quebrado — aqui ele e procurado no disco
BEAGLE=$(ls /data/slick_ancestry/beagle.jar /data/tools/beagle.jar 2>/dev/null | head -1)
[[ -s "$BEAGLE" ]] || { echo "ERRO: beagle.jar nao encontrado"; exit 1; }
echo "      beagle: $BEAGLE"
java -Xmx16g -jar "$BEAGLE" gt=$OUT/merged.vcf.gz out=$OUT/merged_phased \
     impute=false nthreads=$THREADS 2>&1 | tail -5
bcftools index -f $OUT/merged_phased.vcf.gz
echo "[5/5] faseado -> $OUT/merged_phased.vcf.gz"

cat <<'FIM'

Proximo passo (fora deste script, porque exige o roster atualizado):
  1. acrescentar os animais novos ao groups_afr_split.tsv (LMN* -> CRIOLLO,
     N'Dama do ENA -> NDAMA)
  2. python gene_ancestry_by_class.py  --vcf .../merged_phased.vcf.gz ...
  3. python bootstrap_ancestry_class.py --vcf .../merged_phased.vcf.gz ...
  4. python haplotype_network.py       --vcf .../merged_phased.vcf.gz ...
  5. atualizar a tabela de atribuicao por classe funcional, as figuras de ranking
     e de controles pareados, e o resumo/discussao dos resultados

O numero a vigiar e a ESTABILIDADE do SLICK2, hoje 0,580/0,643. Se cair com o n
maior, o texto tem de dizer isso — a decisao de incorporar foi tomada com esse
risco declarado.
FIM
