#!/bin/bash
# ============================================================================
# SENEPOL — genotipagem do PRLR a partir de RNA-seq de figado (GSE192471)
# Objetivo: localizar o alelo SLICK1 do Senepol nas MESMAS coordenadas
# (ARS-UCD1.3) do frameshift criollo 20:39099213, e testar se sao o MESMO alelo.
# Controle: os genotipos SLICK1 publicados (FREE / SLICK / S1 / S2) sao a verdade.
#
# Estrategia de disco: processa 1 run por vez e APAGA o FASTQ logo apos alinhar.
# Pico de uso ~35 GB mesmo baixando 139 GB no total.
# ============================================================================
set -u
REF=/data/genome_reference/ARS-UCD1.3/Bos_taurus.ARS-UCD1.3.dna.toplevel.fa
WORK=/data/slick_ancestry/p4_senepol   # scratch (grande)
PROJ=/proj
OUTD=$PROJ/results/senepol
RUNS=$PROJ/data/GSE192471/senepol_runs.txt
THREADS=${THREADS:-6}
PRLR_REG="20:39037656-39117082"

mkdir -p "$WORK/fastq" "$WORK/bam" "$WORK/idx" "$OUTD"

echo "########## SENEPOL PRLR GENOTYPING ##########"
date

# ---------- 1. Indice hisat2 do BTA20 ----------
IDX="$WORK/idx/bta20"
if [ ! -f "${IDX}.1.ht2" ]; then
  echo ">> [1] Extraindo BTA20 e construindo indice hisat2"
  [ -f "$WORK/idx/bta20.fa" ] || samtools faidx "$REF" 20 > "$WORK/idx/bta20.fa"
  samtools faidx "$WORK/idx/bta20.fa"
  hisat2-build -p "$THREADS" "$WORK/idx/bta20.fa" "$IDX" > "$WORK/idx/build.log" 2>&1 \
    && echo "   indice OK" || { echo "   ERRO no hisat2-build"; tail -5 "$WORK/idx/build.log"; exit 1; }
else
  echo ">> [1] indice hisat2 ja existe"
fi

# ---------- 2. Loop por run ----------
echo ">> [2] Processando runs"
while IFS= read -r RUN; do
  RUN=$(echo "$RUN" | tr -d '\r' | tr -d '[:space:]')
  [ -n "$RUN" ] || continue
  BAM="$WORK/bam/${RUN}.prlr.bam"
  if [ -s "$BAM" ]; then echo "   [$RUN] ja processado, pulando"; continue; fi

  echo "   ---------------- $RUN ----------------"
  t0=$(date +%s)

  # 2a. download — ENA serve fastq.gz direto: resumivel, paralelo, sem o custo de
  #     conversao do fasterq-dump (que streama do NCBI e e muito mais lento).
  echo "   [$RUN] baixando do ENA..."
  # Verificar tamanho + integridade gzip e OBRIGATORIO: um download truncado
  # faz o hisat2 abortar ("more quality values than read characters") e deixa
  # um BAM parcial que passa por valido.
  META_LINE=$(curl -s --max-time 120 \
    "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${RUN}&result=read_run&fields=fastq_ftp,fastq_bytes&format=tsv" \
    | awk 'NR==2{print}')
  URLS=$(echo "$META_LINE" | cut -f1); SIZES=$(echo "$META_LINE" | cut -f2)
  [ -n "$URLS" ] || { echo "   [$RUN] sem FASTQ no ENA"; continue; }
  IFS=';' read -ra UARR <<< "$URLS"; IFS=';' read -ra SARR <<< "$SIZES"
  DL_OK=1
  for i in "${!UARR[@]}"; do
    u="${UARR[$i]}"; want="${SARR[$i]:-0}"; bn=$(basename "$u"); fn="$WORK/fastq/$bn"
    for attempt in 1 2 3; do
      aria2c -x 8 -s 8 -c -q --max-tries=5 --retry-wait=10 \
             -d "$WORK/fastq" -o "$bn" "https://$u" || true
      got=$(stat -c %s "$fn" 2>/dev/null || echo 0)
      if [ "$want" -gt 0 ] && [ "$got" -ne "$want" ]; then
        echo "   [$RUN] $bn tamanho $got != $want, refazendo"; rm -f "$fn"; continue; fi
      if ! pigz -t "$fn" 2>/dev/null; then
        echo "   [$RUN] $bn gzip corrompido, refazendo"; rm -f "$fn"; continue; fi
      echo "   [$RUN] $bn OK ($got bytes)"; break
    done
    [ -s "$fn" ] || { DL_OK=0; break; }
  done
  [ "$DL_OK" -eq 1 ] || { echo "   [$RUN] FALHA no download - pulando"; rm -f "$WORK/fastq/${RUN}"*; continue; }

  R1="$WORK/fastq/${RUN}_1.fastq.gz"; R2="$WORK/fastq/${RUN}_2.fastq.gz"
  if [ -s "$R1" ] && [ -s "$R2" ]; then PAIRED=1
  elif [ -s "$WORK/fastq/${RUN}.fastq.gz" ]; then PAIRED=0; R1="$WORK/fastq/${RUN}.fastq.gz"
  else echo "   [$RUN] FASTQ ausente apos download"; continue; fi

  # 2b. alinhamento (splice-aware) direto para BAM ordenado, so BTA20
  echo "   [$RUN] hisat2 -> BAM"
  if [ "$PAIRED" -eq 1 ]; then
    hisat2 -p "$THREADS" -x "$IDX" -1 "$R1" -2 "$R2" --no-unal 2> "$WORK/bam/${RUN}.hisat2.log" \
      | samtools sort -@ 2 -o "$WORK/bam/${RUN}.bta20.bam" -
  else
    hisat2 -p "$THREADS" -x "$IDX" -U "$R1" --no-unal 2> "$WORK/bam/${RUN}.hisat2.log" \
      | samtools sort -@ 2 -o "$WORK/bam/${RUN}.bta20.bam" -
  fi
  # o hisat2 pode abortar no meio e deixar BAM parcial: exigir a linha de
  # sumario ("overall alignment rate") antes de aceitar o resultado
  if ! grep -q "overall alignment rate" "$WORK/bam/${RUN}.hisat2.log" 2>/dev/null; then
    echo "   [$RUN] ALINHAMENTO ABORTOU - descartando"
    sed -n '1,3p' "$WORK/bam/${RUN}.hisat2.log" 2>/dev/null | sed 's/^/      /'
    rm -f "$WORK/bam/${RUN}.bta20.bam" "$WORK/fastq/${RUN}"*.fastq.gz
    continue
  fi
  samtools index "$WORK/bam/${RUN}.bta20.bam"

  # 2c. APAGA FASTQ imediatamente (controle de disco)
  rm -f "$WORK/fastq/${RUN}"*.fastq.gz
  echo "   [$RUN] FASTQ removido"

  # 2d. fatia so o PRLR e descarta o BAM do cromossomo inteiro
  samtools view -b -q 30 "$WORK/bam/${RUN}.bta20.bam" "$PRLR_REG" > "$BAM"
  samtools index "$BAM"
  rm -f "$WORK/bam/${RUN}.bta20.bam" "$WORK/bam/${RUN}.bta20.bam.bai"

  d=$(samtools depth -r "$PRLR_REG" "$BAM" 2>/dev/null | awk '{s+=$3;n++} END{if(n)printf "%.1f",s/n; else print 0}')
  t1=$(date +%s)
  echo "   [$RUN] OK — cobertura media PRLR=${d}x  ($(( (t1-t0)/60 )) min)"
done < "$RUNS"

# ---------- 3. Chamada conjunta de variantes no PRLR ----------
echo ""
echo ">> [3] Chamada conjunta de variantes (PRLR)"
ls "$WORK"/bam/*.prlr.bam > "$WORK/bamlist.txt" 2>/dev/null
echo "   BAMs: $(wc -l < "$WORK/bamlist.txt")"
if [ -s "$WORK/bamlist.txt" ]; then
  bcftools mpileup -f "$REF" -r "$PRLR_REG" -b "$WORK/bamlist.txt" \
      -a AD,DP -q 30 -Q 20 --max-depth 5000 -Ou 2>/dev/null \
    | bcftools call -mv -Ov -o "$OUTD/senepol_prlr.vcf" 2>/dev/null
  echo "   variantes chamadas: $(grep -vc '^#' "$OUTD/senepol_prlr.vcf" 2>/dev/null)"
  echo "   -- INDELS (candidatos a SLICK1) --"
  bcftools view -v indels "$OUTD/senepol_prlr.vcf" 2>/dev/null \
    | bcftools query -f '%POS\t%REF\t%ALT\t%QUAL[\t%GT]\n' 2>/dev/null
  echo "   -- cabecalho de amostras --"
  bcftools query -l "$OUTD/senepol_prlr.vcf" 2>/dev/null | tr '\n' ' '; echo
fi

echo ""
echo "########## SENEPOL CONCLUIDO ##########"
date
du -sh "$WORK" 2>/dev/null
