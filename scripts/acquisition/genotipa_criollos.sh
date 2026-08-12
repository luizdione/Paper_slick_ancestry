#!/usr/bin/env bash
# genotipa_criollos.sh — genotipa SLICK1/SLICK2 nos BAM produzidos pelo fetch_v2.sh.
#
# O fetch_v2.sh so entrega BAM da janela. A genotipagem vivia colada no fim do v1 e
# nao foi transportada, entao ela e este script — e assim ela pode rodar sobre os BAM
# ja prontos sem esperar as 58 runs terminarem.
#
# O que se procura: um portador do SLICK1 em fundo **menos endogamico**. O unico
# homozigoto conhecido e um Costeno con Cuernos cujo genoma inteiro esta em blocos de
# ~5 Mb; nele o comprimento de trato no locus mede o pedigree, nao a idade do alelo
# (razao ao proprio baseline = 1,0, P = 0,21). Um portador com fundo quebrado
# devolveria o trato como medida util, e so entao a datacao poderia voltar.
#
# ETAPA 0 e um GATE, nao um relatorio. Um filtro por k-mer calibrado errado devolve
# BAM que parecem bons e nao cobrem o sitio; a genotipagem em cima disso produz `0/0`
# em toda amostra, que se le como "nenhum portador novo" quando na verdade e "nao
# olhei". Cobertura no sitio antes de qualquer genotipo.
set -uo pipefail

DATA=/data/slick_ancestry/criollo
WIN_START=34100000
SLICK1=39099213
SLICK2=39099321
MIN_DP=${MIN_DP:-4}          # profundidade minima para o genotipo valer

BAIT="$DATA/bait/win10Mb.fa"
[[ -s "$BAIT" ]] || { echo "ERRO: isca ausente em $BAIT"; exit 1; }

# as coordenadas na isca sao relativas ao inicio da janela
P1=$((SLICK1 - WIN_START + 1))
P2=$((SLICK2 - WIN_START + 1))
CHR=$(head -1 "$BAIT" | sed 's/^>//' | cut -d' ' -f1)

mapfile -t BAMS < <(ls "$DATA"/bam/*.bam 2>/dev/null)
(( ${#BAMS[@]} )) || { echo "ERRO: nenhum BAM em $DATA/bam/"; exit 1; }
echo "[info] ${#BAMS[@]} BAM | isca $CHR | SLICK1 -> $P1 | SLICK2 -> $P2"

# ---------------------------------------------------------------- 0. GATE de cobertura
echo
echo "=== cobertura nos dois sitios (gate) ==="
printf "%-14s %8s %8s %10s %s\n" RUN DP_SLICK1 DP_SLICK2 READS VEREDICTO
OK=0; RUIM=0
for B in "${BAMS[@]}"; do
  RUN=$(basename "$B" .bam)
  # -a força linha mesmo com zero; sem isso um sitio sem cobertura simplesmente
  # some da saida e o campo fica vazio em vez de zero
  D1=$(samtools depth -a -r "${CHR}:${P1}-${P1}" "$B" 2>/dev/null | cut -f3); D1=${D1:-0}
  D2=$(samtools depth -a -r "${CHR}:${P2}-${P2}" "$B" 2>/dev/null | cut -f3); D2=${D2:-0}
  NR=$(samtools view -c "$B" 2>/dev/null)
  if (( D1 >= MIN_DP && D2 >= MIN_DP )); then V="ok"; ((OK++))
  else V="RASO (< $MIN_DP)"; ((RUIM++)); fi
  printf "%-14s %8s %8s %10s %s\n" "$RUN" "$D1" "$D2" "$NR" "$V"
done
echo "-> $OK com cobertura util, $RUIM rasos"
if (( OK == 0 )); then
  echo
  echo "GATE REPROVOU: nenhum BAM cobre os dois sitios."
  echo "Antes de seguir, checar se o filtro por k-mer (mkf=0.5 no fetch_v2.sh) nao"
  echo "esta descartando as reads do locus. Genotipar agora so produziria 0/0."
  exit 2
fi

# ---------------------------------------------------------------- 1. genotipagem
echo
echo "=== genotipando ==="
printf '%s\n' "${BAMS[@]}" > "$DATA/bamlist.txt"

# O contig da isca se chama "20:34100000-44100000" — o nome carrega dois-pontos,
# porque veio do `samtools faidx` sobre uma regiao. Passar isso em `-r` monta
# "20:34100000-44100000:4999014-4999522", que o bcftools nao consegue interpretar:
# ele nao falha, apenas nao seleciona nada, e o VCF sai com **zero linhas**. A
# genotipagem seguinte le um arquivo vazio e reporta "nenhum portador", que se
# confunde com o resultado real. Um BED separa contig e coordenadas por TAB e
# elimina a ambiguidade.
printf '%s\t%d\t%d\n' "$CHR" $((P1-201)) $((P2+200)) > "$DATA/alvo.bed"
bcftools mpileup -f "$BAIT" -b "$DATA/bamlist.txt" \
    -R "$DATA/alvo.bed" -a AD,DP -q 20 -Q 20 -Ou \
  | bcftools call -m -Ov -o "$DATA/slick_sites.vcf"

# GATE: VCF vazio nao e "nenhum portador", e "nao genotipei"
N_VCF=$(grep -vc '^#' "$DATA/slick_sites.vcf" || true)
echo "-> $DATA/slick_sites.vcf ($N_VCF sitios)"
if (( N_VCF == 0 )); then
  echo
  echo "GATE REPROVOU: o VCF saiu vazio. Os BAM tem cobertura (etapa 0 passou),"
  echo "entao o problema esta na selecao de regiao, nao nos dados. Conferir se o"
  echo "nome do contig da isca ($CHR) esta sendo passado por BED e nao por -r."
  exit 3
fi

# ---------------------------------------------------------------- 2. quem carrega o que
echo
echo "=== genotipos nos dois sitios (amostra: GT:DP) ==="
# os nomes de amostra do bcftools sao os caminhos dos BAM; encurtar para o run
bcftools query -l "$DATA/slick_sites.vcf" | sed 's|.*/||; s|\.bam$||' > "$DATA/tmp_samples.txt"
# nomes das amostras, na ordem das colunas do VCF
mapfile -t NOMES < <(bcftools query -l "$DATA/slick_sites.vcf" | sed 's|.*/||; s|\.bam$||')

# Uma linha por (sitio, amostra). O formato antigo empacotava tudo numa linha so e
# era desempacotado por awk com dois separadores; qualquer descasamento de campo
# devolvia "nenhum portador" em silencio, que se le como resultado. Aqui cada
# genotipo sai numa linha propria e o indice da amostra e explicito.
for POS in "$P1" "$P2"; do
  NOME=$([[ "$POS" == "$P1" ]] && echo SLICK1 || echo SLICK2)
  echo "--- $NOME (pos $POS na isca) ---"
  ACHOU=0
  while IFS=$'\t' read -r pp ref alt gts; do
    [[ "$pp" == "$POS" ]] || continue
    echo "    REF=$ref  ALT=$alt"
    i=0
    for campo in $gts; do
      gt=${campo%%:*}; rest=${campo#*:}; dp=${rest%%:*}; ad=${rest#*:}
      if [[ "$gt" != "0/0" && "$gt" != "./." && "$gt" != "." && ${dp:-0} -ge $MIN_DP ]]; then
        printf "    PORTADOR: %-14s %-5s DP=%-4s AD=%s\n" "${NOMES[$i]}" "$gt" "$dp" "$ad"
        ACHOU=$((ACHOU+1))
      fi
      i=$((i+1))
    done
  done < <(bcftools query -f '%POS\t%REF\t%ALT[\t%GT:%DP:%AD]\n' "$DATA/slick_sites.vcf" 2>/dev/null \
           | awk -F'\t' -v p="$POS" '$1==p {printf "%s\t%s\t%s\t", $1,$2,$3; for(i=4;i<=NF;i++) printf "%s ", $i; print ""}')
  (( ACHOU == 0 )) && echo "    nenhum portador com DP >= $MIN_DP"
done

# O SLICK1 e um INDEL (c.1382del), e o `bcftools call` o registra numa linha propria
# cujo POS pode nao coincidir com o do sitio. Procurar so pela posicao devolveria
# "nenhum portador" mesmo havendo deleção chamada dois nucleotideos ao lado.
echo
echo "--- SLICK1: chamadas de indel na vizinhanca (${P1}±25) ---"
bcftools query -f '%POS\t%REF\t%ALT\t%TYPE[\t%GT]\n' "$DATA/slick_sites.vcf" 2>/dev/null \
  | awk -F'\t' -v p="$P1" '$1>=p-25 && $1<=p+25 && $4!="REF"' \
  | head -10 || true
echo "    (vazio = nenhuma variante chamada na vizinhanca do SLICK1)"

echo
echo "Se apareceu portador novo do SLICK1: rodar inbreeding_control.py nesse animal."
echo "So se o trato no locus exceder a linha de base DELE e que a datacao volta."
