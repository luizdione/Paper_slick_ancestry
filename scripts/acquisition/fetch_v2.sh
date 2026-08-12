#!/usr/bin/env bash
# fetch_v2.sh — v2: prioriza Limonero e troca o stream por download paralelo.
#
# A v1 usava `curl | bbduk | bwa` num pipe unico, o que nunca gravava o FASTQ mas
# limitava o download a uma conexao: 2,6 MB/s efetivos, ou ~92 h para os 58 runs.
# O ENA entrega muito mais com conexoes paralelas, e o custo e ter o FASTQ em
# disco por alguns minutos.
#
# Duas mudancas:
#   1. aria2c -x8 -s8 baixa o par, processa, e o FASTQ e APAGADO em seguida.
#      Pico de disco por run: ~17 GB, transitorio.
#   2. Os 30 Limonero (LMN) vem primeiro. E a raca que carrega os dois alelos —
#      os portadores conhecidos do SLICK1 e do SLICK2 sao Limonero e Costeno con
#      Cuernos. Se um portador em fundo menos endogamico existe, esta ali.
#
# CORRECAO 2026-08-04 (a primeira tentativa perdeu dois runs assim):
#   `aria2c URL_1 URL_2` NAO baixa dois arquivos. Multiplas URIs na linha de
#   comando sao tratadas como ESPELHOS DO MESMO arquivo — o aria2 baixa um so,
#   com o nome da primeira, e sai com **codigo 0**. O `_2` nunca chega, o teste
#   de existencia reprova, e o run e pulado depois de 8 GB baixados. Um `-q` no
#   comando deixa o log vazio, entao nao ha nem de onde suspeitar.
#   Cada URI tem de ir numa LINHA do arquivo passado em `-i` (ou usar `-Z`).
#   Somado a isso, o tamanho baixado e conferido contra o `fastq_bytes` que o
#   ENA declara: truncamento tambem nao levanta erro.
set -uo pipefail
D=/data/slick_ancestry/criollo
REF=$D/bait/win10Mb.fa
TSV=/data/slick_ancestry/p4_slick/data/criollo/ena_ERP150979_raw.tsv
TH=${TH:-6}
mkdir -p $D/{bam,tmp,logs}

processa() {
  local RUN=$1 U1=$2 U2=$3 TITLE=$4 BYTES=$5
  [[ -s $D/bam/${RUN}.bam ]] && { echo "[$RUN] ja feito"; return 0; }
  echo "[$RUN] $TITLE — baixando"
  rm -f $D/tmp/${RUN}_* $D/tmp/${RUN}.urls
  # uma URI por linha: e o unico modo em que o aria2 as trata como downloads
  # distintos em vez de espelhos do mesmo arquivo (ver cabecalho)
  printf 'ftp://%s\nftp://%s\n' "$U1" "$U2" > $D/tmp/${RUN}.urls
  aria2c -i $D/tmp/${RUN}.urls -j2 -x8 -s8 --dir=$D/tmp --allow-overwrite=true \
         --summary-interval=0 --console-log-level=warn \
         >> $D/logs/${RUN}.aria.log 2>&1 || { echo "[$RUN] download falhou"; return 1; }
  local F1=$D/tmp/$(basename $U1) F2=$D/tmp/$(basename $U2)
  [[ -s $F1 && -s $F2 ]] || { echo "[$RUN] arquivo ausente"; rm -f $D/tmp/${RUN}_*; return 1; }
  # O aria2 PRE-ALOCA o arquivo no tamanho final: `stat -c%s` devolve o tamanho
  # certo desde o primeiro segundo, entao conferir tamanho nao detecta download
  # incompleto. O marcador honesto e o arquivo de controle `.aria2`, que so some
  # quando o download fecha.
  [[ -e "$F1.aria2" || -e "$F2.aria2" ]] && {
    echo "[$RUN] download incompleto (.aria2 presente)"; rm -f $F1 $F2 $F1.aria2 $F2.aria2; return 1; }
  local E1=${BYTES%%;*} E2=${BYTES##*;}
  local G1=$(stat -c%s "$F1") G2=$(stat -c%s "$F2")
  [[ "$G1" == "$E1" && "$G2" == "$E2" ]] || {
    echo "[$RUN] tamanho divergente do declarado: _1 $G1/$E1  _2 $G2/$E2"; rm -f $F1 $F2; return 1; }
  echo "[$RUN] baixado $(( (G1+G2)/1073741824 )) GB — filtrando + alinhando"
  cat $F1 $F2 \
    | bbduk.sh -Xmx8g in=stdin.fq.gz int=f outm=stdout.fq ref=$REF k=31 hdist=0 mkf=0.5 \
        threads=$TH ow=t 2> $D/logs/${RUN}.bbduk.log \
    | bwa mem -t $TH $REF - 2> $D/logs/${RUN}.bwa.log \
    | samtools view -b -F 4 -q 20 - \
    | samtools sort -@2 -o $D/bam/${RUN}.bam - 2>> $D/logs/${RUN}.bwa.log
  samtools index $D/bam/${RUN}.bam
  rm -f $F1 $F2 $D/tmp/${RUN}.urls    # o FASTQ sai do disco assim que o BAM existe
  echo "[$RUN] OK $(du -h $D/bam/${RUN}.bam | cut -f1) | $(samtools view -c $D/bam/${RUN}.bam) reads"
}

# Limonero primeiro, depois o resto
for FILTRO in '$3 ~ /^LMN/' '$3 !~ /^LMN/'; do
  tail -n +2 "$TSV" | tr -d '\r' | awk -F'\t' "$FILTRO" | while IFS=$'\t' read -r RUN SAMP TITLE ALIAS SCI STRAT BYTES FTP; do
    U1=${FTP%%;*}; U2=${FTP##*;}
    [[ -z "$U1" || -z "$U2" || "$U1" == "$U2" ]] && { echo "[$RUN] sem par de FASTQ — pulado"; continue; }
    processa "$RUN" "$U1" "$U2" "$TITLE" "$BYTES" || echo "[$RUN] pulado"
  done
done
echo "[fim] $(ls $D/bam/*.bam 2>/dev/null | wc -l) BAM prontos"
