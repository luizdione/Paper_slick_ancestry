#!/usr/bin/env bash
# senepol_no_sitio.sh — o Senepol carrega qual dos dois alelos slick?
#
# Porque existe
# -------------
# Nao ha WGS publico de Senepol, entao a raca nao entra na analise de haplotipo:
# o que se compara ali sao haplotipos criollos. Isso deixava um vao entre a metade
# da ancestralidade, ancorada em Senepol, e a metade do haplotipo, ancorada em
# criollos. O RNA-seq de figado do GSE192471 fecha o vao no nivel do genotipo:
# nao sustenta faseamento, mas responde "a raca carrega o alelo?".
#
# O desenho do GSE192471 e o que torna isto uma verificacao e nao so uma leitura.
# O estudo compara Senepol **heterozigotos e homozigotos para o SLICK1** (coortes S1
# e S2 na tabela de amostras), entao ha genotipo declarado com que confrontar a
# chamada. Cuidado ao ler a coluna `Genotype` do GEO: ela mistura dois vocabularios
# — S1/S2 sao as coortes deste desenho, enquanto SLICK/FREE pertencem a um segundo
# grupo com outro significado. Usar as duas como se fossem a mesma escala produz
# "animais SLICK sem o alelo", que foi o que apareceu na primeira leitura.
#
# Resultado, com o VCF ja produzido no caminho apontado pela variavel $VCF abaixo:
#
#   SLICK1 (20:39,099,213, GCT->GT)
#     SRR17311902  S2 (homozigoto declarado)  1/1   DP 4884   AD 49 ref / 4834 alt
#     SRR17311903  S2 (homozigoto declarado)  1/1   DP   49   AD  0 ref /   49 alt
#   SLICK2 (20:39,099,321)
#     todos os sete animais  0/0, com DP de ate 4271 e zero leituras do alelo derivado
#
# Duas leituras, e a segunda importa tanto quanto a primeira:
#   1. os homozigotos declarados saem 1/1 com a fracao alelica que se espera, o que
#      valida a chamada num indel de 1 bp lido em RNA-seq;
#   2. **nenhum Senepol carrega o SLICK2** — a raca carrega o alelo que a analise de
#      haplotipo atribui ao fundo N'Dama, e nao o de fundo europeu.
#
# Limite honesto: sao dois animais com genotipo declarado, em RNA-seq, sem fase. Isto
# estabelece presenca do alelo na raca; nao estabelece o haplotipo em que ele viaja
# no Senepol, que continua dependendo de WGS proprio.
set -uo pipefail

DATA=${DATA:-/data/slick_ancestry/p4_slick}
VCF=$DATA/results/senepol_validacao/senepol_2sitios.vcf.gz
MAP=$DATA/data/GSE192471/senepol_run_genotype_map.tsv
SLICK1=39099213
SLICK2=39099321

[[ -s "$VCF" ]] || { echo "ERRO: VCF ausente em $VCF"; exit 1; }

echo "=== Senepol nos dois sitios slick (GSE192471, RNA-seq de figado) ==="
printf "%-13s %-10s %-7s %-7s %-16s %s\n" RUN COORTE SLICK1 DP AD_ref/alt SLICK2

# nomes das amostras: o VCF guarda o caminho do BAM, nao o run
mapfile -t RUNS < <(bcftools query -l "$VCF" 2>/dev/null | sed 's|.*/||; s|\.prlr\.bam$||')

# o SLICK1 e um indel: a linha util e a que traz GCT->GT (delecao de 1 bp, c.1382del).
# A posicao tambem aparece como registro so-referencia, que nao carrega genotipo util.
G1=$(bcftools query -r "20:${SLICK1}-${SLICK1}" -f '%REF\t%ALT[\t%GT:%DP:%AD]\n' "$VCF" \
     | awk -F'\t' '$1=="GCT"' | cut -f3-)
G2=$(bcftools query -r "20:${SLICK2}-${SLICK2}" -f '[\t%GT]\n' "$VCF" | head -1)

i=1
for run in "${RUNS[@]}"; do
  coorte=$(tr -d '"' < "$MAP" | awk -F'\t' -v r="$run" '$1==r{print $3}')
  campo=$(printf '%s' "$G1" | cut -f$i)
  gt=${campo%%:*}; resto=${campo#*:}; dp=${resto%%:*}; ad=${resto#*:}
  gt2=$(printf '%s' "$G2" | cut -f$((i+1)))
  printf "%-13s %-10s %-7s %-7s %-16s %s\n" "$run" "${coorte:-?}" "$gt" "$dp" "$ad" "$gt2"
  i=$((i+1))
done

echo
echo "Confronto com o genotipo declarado: as coortes S1/S2 do GSE192471 sao"
echo "heterozigotos e homozigotos para o SLICK1. SLICK/FREE pertencem a outro grupo"
echo "e NAO devem ser lidos na mesma escala."
