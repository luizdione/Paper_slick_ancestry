#!/usr/bin/env bash
# fetch_ndama.sh — amplia o painel N'Dama.
#
# Dez N'Dama excluem so frequencias acima de 13,91%; e nesse limite que hoje mora a
# frase "nao detectado na fonte" no relatorio de resultados. Com 22 individuos o
# limite cai para ~6,5%, e a ausencia passa a dizer mais.
#
# Mesmo pipeline do `fetch_v2.sh`: aria2c baixa o par, bbduk filtra por k-mer contra
# a janela de 10 Mb, bwa alinha, e o FASTQ e apagado. Duas diferencas:
#   - roster proprio (`ndama_ena_roster.tsv`, 12 individuos, 613 GB), montado
#     excluindo os 27 runs que sao **pools de DNA** — um pool nao serve para painel
#     de haplotipo, porque nao ha individuo a que atribuir a copia;
#   - `-x4` em vez de `-x8`, porque roda ao lado do `criollo_fetch` e a banda e a
#     mesma. Dividir as conexoes evita que um sufoque o outro.
set -uo pipefail
D=/data/slick_ancestry/criollo
REF=$D/bait/win10Mb.fa
TSV=$D/ndama_ena_roster.tsv
TH=${TH:-4}
mkdir -p $D/{bam_ndama,tmp_ndama,logs}

processa() {
  local RUN=$1 U1=$2 U2=$3 TITLE=$4 BYTES=$5
  [[ -s $D/bam_ndama/${RUN}.bam ]] && { echo "[$RUN] ja feito"; return 0; }
  echo "[$RUN] $TITLE — baixando"
  rm -f $D/tmp_ndama/${RUN}_* $D/tmp_ndama/${RUN}.urls
  printf 'ftp://%s\nftp://%s\n' "$U1" "$U2" > $D/tmp_ndama/${RUN}.urls
  aria2c -i $D/tmp_ndama/${RUN}.urls -j2 -x4 -s4 --dir=$D/tmp_ndama \
         --allow-overwrite=true --summary-interval=0 --console-log-level=warn \
         >> $D/logs/nd_${RUN}.aria.log 2>&1 || { echo "[$RUN] download falhou"; return 1; }
  local F1=$D/tmp_ndama/$(basename $U1) F2=$D/tmp_ndama/$(basename $U2)
  [[ -s $F1 && -s $F2 ]] || { echo "[$RUN] arquivo ausente"; rm -f $D/tmp_ndama/${RUN}_*; return 1; }
  # o arquivo de controle .aria2 e o unico teste honesto de download completo:
  # o aria2 pre-aloca o arquivo no tamanho final, entao conferir tamanho nao detecta
  # truncamento.
  [[ -e "$F1.aria2" || -e "$F2.aria2" ]] && {
    echo "[$RUN] incompleto"; rm -f $F1 $F2 $F1.aria2 $F2.aria2; return 1; }
  echo "[$RUN] filtrando + alinhando"
  cat $F1 $F2 \
    | bbduk.sh -Xmx8g in=stdin.fq.gz int=f outm=stdout.fq ref=$REF k=31 hdist=0 mkf=0.5 \
        threads=$TH ow=t 2> $D/logs/nd_${RUN}.bbduk.log \
    | bwa mem -t $TH $REF - 2> $D/logs/nd_${RUN}.bwa.log \
    | samtools view -b -F 4 -q 20 - \
    | samtools sort -@2 -o $D/bam_ndama/${RUN}.bam - 2>> $D/logs/nd_${RUN}.bwa.log
  samtools index $D/bam_ndama/${RUN}.bam
  rm -f $F1 $F2 $D/tmp_ndama/${RUN}.urls
  echo "[$RUN] OK $(du -h $D/bam_ndama/${RUN}.bam | cut -f1) | $(samtools view -c $D/bam_ndama/${RUN}.bam) reads"
}

# CRLF: o campo do FTP e o ultimo da linha e carrega o \r; sem o tr a URL nao resolve
# e o run "termina" em segundos com zero leituras.
tail -n +2 "$TSV" | tr -d '\r' | while IFS=$'\t' read -r RUN SAMP TITLE BYTES FTP STUDY; do
  U1=${FTP%%;*}; U2=${FTP##*;}
  [[ -z "$U1" || -z "$U2" || "$U1" == "$U2" ]] && { echo "[$RUN] sem par — pulado"; continue; }
  processa "$RUN" "$U1" "$U2" "$TITLE" "$BYTES" || echo "[$RUN] pulado"
done
echo "[fim] $(ls $D/bam_ndama/*.bam 2>/dev/null | wc -l) BAM de N'Dama prontos"
