#!/usr/bin/env bash
# fetch_criollo_prlr.sh — portadores de SLICK1 em fundos menos consanguineos.
#
# O ENA hospeda 58 runs WGS de criollos (ERP150979, 861 GB de FASTQ). Baixar tudo
# para olhar 10 Mb de BTA20 seria desperdicio de tres ordens de grandeza. A janela
# que interessa e 0,37% do genoma.
#
# O fluxo nunca grava o FASTQ inteiro:
#
#   curl (stream) -> bbduk (k=31 contra a janela) -> so o que casa toca o disco
#                 -> bwa mem -> BAM da janela -> genotipo nos dois sitios slick
#
# Cada run deixa um BAM de dezenas de MB em vez de 15 GB. O gargalo e banda, nao
# disco: o pico de uso transitorio fica em poucos GB.
#
# Por que isso importa cientificamente: o unico homozigoto do SLICK1 e um Costeno
# con Cuernos cujo genoma inteiro esta em blocos de ~5 Mb, e por isso o comprimento
# de trato no locus mede o pedigree dele, nao a idade do alelo. Um portador em
# fundo menos endogamico devolveria o comprimento de trato como medida util.
set -euo pipefail

DATA=/data/slick_ancestry/criollo
REF=/data/slick_ancestry/p4_senepol/idx/bta20.fa
WIN_START=34100000
WIN_END=44100000
SLICK1=39099213
SLICK2=39099321
THREADS=${THREADS:-6}

mkdir -p "$DATA"/{bait,bam,filt,logs}

# ---------------------------------------------------------------- isca de k-mers
# A janela de 10 Mb ao redor do PRLR e o que o bbduk usa para decidir se uma read
# vem daqui. entropy=0.7 descarta k-mers de baixa complexidade, que puxariam reads
# repetitivas do genoma inteiro.
if [[ ! -s "$DATA/bait/win10Mb.fa" ]]; then
  echo "[bait] extraindo BTA20:${WIN_START}-${WIN_END}"
  samtools faidx "$REF" "20:${WIN_START}-${WIN_END}" > "$DATA/bait/win10Mb.fa"
  samtools faidx "$DATA/bait/win10Mb.fa"
  bwa index "$DATA/bait/win10Mb.fa" 2> "$DATA/logs/bwa_index.log"
fi
echo "[bait] $(grep -v '^>' "$DATA/bait/win10Mb.fa" | tr -d '\n' | wc -c) bp de isca"

# ------------------------------------------------------------------ por run
processa() {
  local RUN=$1 URL1=$2 URL2=$3 TITLE=$4
  local BAM="$DATA/bam/${RUN}.bam"
  [[ -s "$BAM" ]] && { echo "[$RUN] ja feito"; return 0; }

  echo "[$RUN] $TITLE — stream + filtro"
  # Tres pontos que custaram uma calibracao:
  #   outm=  guarda o que CASA. O bbduk e um descontaminador: `out=` guarda o
  #          complemento, e usa-lo aqui devolveria o genoma inteiro menos a janela.
  #   mkf=0.5 exige que metade dos k-mers da read casem. Sem isso, um unico k-mer
  #          repetitivo puxa a read, e a janela de 10 Mb captura 24,9% do genoma em
  #          vez dos 0,37% que ela representa. Com mkf=0.5 cai para 6,6%; o
  #          alinhamento descarta o resto. Nao subimos para 0.8 (1,9%) porque uma
  #          read com um SNP ja perde ~26% dos seus k-mers e seria descartada —
  #          justamente as reads informativas.
  #   O FASTQ filtrado nunca e gravado: sai do bbduk direto para o bwa por pipe.
  #          So o BAM da janela toca o disco.
  # Os dois lados do par entram como stream unico e o alinhamento roda em modo
  # single; a informacao de par nao e necessaria para genotipar dois sitios.
  { curl -sL "ftp://$URL1"; curl -sL "ftp://$URL2"; } \
    | bbduk.sh -Xmx8g in=stdin.fq.gz int=f outm=stdout.fq \
        ref="$DATA/bait/win10Mb.fa" k=31 hdist=0 mkf=0.5 \
        threads=$THREADS ow=t 2> "$DATA/logs/${RUN}.bbduk.log" \
    | bwa mem -t $THREADS "$DATA/bait/win10Mb.fa" - 2> "$DATA/logs/${RUN}.bwa.log" \
    | samtools view -b -F 4 -q 20 - \
    | samtools sort -@2 -o "$BAM" - 2>> "$DATA/logs/${RUN}.bwa.log"
  samtools index "$BAM"
  echo "[$RUN] $(du -h "$BAM" | cut -f1) | $(samtools view -c "$BAM") reads"
}

# ------------------------------------------------------------------ fila
TSV=/data/slick_ancestry/p4_slick/data/criollo/ena_ERP150979_raw.tsv
# O TSV do ENA vem com CRLF. Sem o tr, o \r fica colado no fim do campo FTP (que e
# o ultimo da linha), a URL nao resolve, o curl sai calado e o run termina com zero
# reads em poucos segundos. Mesma pegadinha que ja quebrou outro roster do projeto.
tail -n +2 "$TSV" | tr -d '\r' | while IFS=$'\t' read -r RUN SAMP TITLE ALIAS SCI STRAT BYTES FTP; do
  U1=${FTP%%;*}; U2=${FTP##*;}
  [[ -z "$U1" || -z "$U2" ]] && continue
  processa "$RUN" "$U1" "$U2" "$TITLE" || echo "[$RUN] FALHOU, seguindo"
done

# ------------------------------------------------------------- genotipagem
echo "[call] genotipando os dois sitios slick em todos os BAM"
ls "$DATA"/bam/*.bam > "$DATA/bamlist.txt"
# as coordenadas na isca sao relativas ao inicio da janela
P1=$((SLICK1 - WIN_START + 1)); P2=$((SLICK2 - WIN_START + 1))
CHR=$(head -1 "$DATA/bait/win10Mb.fa" | sed 's/^>//' | cut -d' ' -f1)
bcftools mpileup -f "$DATA/bait/win10Mb.fa" -b "$DATA/bamlist.txt" \
    -r "${CHR}:$((P1-200))-$((P2+200))" -a AD,DP -q 20 -Q 20 -Ou \
  | bcftools call -m -Ov -o "$DATA/slick_sites.vcf"
echo "[call] -> $DATA/slick_sites.vcf"
bcftools query -f '%POS\t%REF\t%ALT[\t%GT]\n' "$DATA/slick_sites.vcf" \
  | awk -v a=$P1 -v b=$P2 '$1==a || $1==b'
