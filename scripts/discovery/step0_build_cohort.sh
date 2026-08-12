#!/bin/bash
# ============================================================================
# STEP 0 — Constroi a matriz genome-wide da coorte de descoberta (pre-requisito eixos 2 e 3)
# Extrai SNPs bialelicos, GT-only, MAF>=0.01 dos 29 autossomos do Run8 TAUIND.
# Saida: /proj/data/cohort/chr{N}.cohort.vcf.gz  + cohort_samples.txt
# ============================================================================
set -u
# pipefail e OBRIGATORIO aqui: sem ele, um "Error: VCF parse error" no meio do
# pipe e mascarado pelo exit 0 do ultimo comando, e o cromossomo fica truncado
# em silencio. Foi o que aconteceu com o chr2 (parou em 134.461.366 de
# 136.231.102 por um registro defeituoso do Run 8 com 343 colunas em vez de 1842).
set -o pipefail
META=/data/1000bg_metadata
VCFDIR=/data/1000bg_vcf
OUT=/proj/data/cohort
LOG=/proj/logs
THREADS=${THREADS:-4}
# Ordem: cromossomos menores primeiro, para dar resultado utilizavel cedo
CHRS=${CHRS:-"25 26 27 28 29 24 23 22 21 20 19 18 17 16 15 14 13 12 11 10 9 8 7 6 5 4 3 2 1"}
mkdir -p "$OUT" "$LOG"

echo "########## STEP 0 — build cohort ##########"
date

# ---------- 1. Monta a coorte ----------
# Raças grandes são subamostradas para equilibrar os grupos e reduzir IO.
build_list() {
  local breed="$1" max="$2"
  local f="$META/samples_${breed}.txt"
  [ -f "$f" ] || { echo "   [!] ausente: $breed" >&2; return; }
  # 10/111 arquivos do roster tem CRLF -> remover \r, senao "X" e "X\r"
  # sobrevivem ao sort -u e o bcftools aborta com "Duplicate sample name"
  tr -d '\r' < "$f" | grep . | head -n "$max"
}

CO="$OUT/cohort_samples.txt"
MAP="$OUT/cohort_groups.tsv"
: > "$CO"; : > "$MAP"

add_group() {
  local grp="$1"; shift
  for spec in "$@"; do
    local b="${spec%%:*}" m="${spec##*:}"
    while IFS= read -r s; do
      [ -n "$s" ] || continue
      echo "$s" >> "$CO"
      printf "%s\t%s\t%s\n" "$s" "$b" "$grp" >> "$MAP"
    done < <(build_list "$b" "$m")
  done
}

add_group TROP_TAURINE   "limonero:99" "costeno_con_cuernos:99" "san_martinero:99" "texas_longhorn:99" "corriente:99"
add_group AFR_TAURINE    "n'dama:99"
add_group AFR_SANGA_ZEBU "ankole:99" "boran:99" "ogaden:99" "nganda:99" "nkiga:99" "nsongora:99" "ntoro:99"
add_group IBERIAN        "pajuna:99" "sayaguesa:99" "maremmana:99" "maronesa:99" "limia:99" "podolica:99" "boskarin:99"
add_group TEMP_TAURINE   "angus:40" "holstein:40" "hereford:40" "jersey:28"
add_group INDICINE       "nelore:40" "gir:36" "brahman:32" "sahiwal:99" "tharparkar:99" "kangayam:99" "vechur:99" "rashoki:99"

# Uma amostra pode constar em >1 arquivo de raca: manter a PRIMEIRA atribuicao
awk -F'\t' '!seen[$1]++' "$MAP" > "$MAP.tmp" && mv "$MAP.tmp" "$MAP"
cut -f1 "$MAP" | sort -u > "$CO"
echo ">> coorte: $(wc -l < "$CO") amostras unicas"
echo "   (amostras em >1 raca, resolvidas p/ 1a: $(cut -f1 "$MAP" | wc -l) linhas no mapa)"
echo "   -- por grupo --"
cut -f3 "$MAP" | sort | uniq -c | sort -rn | sed 's/^/      /'
echo "   -- por raca --"
cut -f2 "$MAP" | sort | uniq -c | sort -rn | sed 's/^/      /'

# ---------- 2. Extrai genome-wide ----------
echo ""
echo ">> Extraindo SNPs bialelicos GT-only (MAF>=0.01, F_MISSING<0.1)"
for c in $CHRS; do
  V="$VCFDIR/Chr${c}-Run8-TAUIND-public.vcf.gz"
  O="$OUT/chr${c}.cohort.vcf.gz"
  if [ ! -f "$V" ]; then echo "   [!] falta $V"; continue; fi
  if [ -s "$O" ] && bcftools index -f "$O" 2>/dev/null; then
    echo "   chr${c}: ja existe, pulando ($(bcftools view -H "$O" 2>/dev/null | wc -l) SNPs)"; continue
  fi
  t0=$(date +%s)
  # --threads acelera a descompressao BGZF, que e o gargalo (1842 genotipos/sitio)
  bcftools view --threads "$THREADS" -S "$CO" --force-samples -m2 -M2 -v snps -Ou "$V" \
    | bcftools +fill-tags -Ou -- -t AF,AC,AN,F_MISSING,MAF \
    | bcftools view -i 'MAF>=0.01 && F_MISSING<0.1' -Ou \
    | bcftools annotate --threads "$THREADS" -x '^FORMAT/GT' -Oz -o "$O" \
    || { echo "   [ERRO] chr${c} falhou"; continue; }
  bcftools index -f --threads "$THREADS" "$O"
  n=$(bcftools view -H "$O" 2>/dev/null | wc -l)
  t1=$(date +%s)
  echo "   chr${c}: ${n} SNPs  ($(( (t1-t0)/60 )) min)  $(du -h "$O" 2>/dev/null | cut -f1)"
done

echo ""
echo ">> TOTAL de SNPs na coorte:"
tot=0
for c in $CHRS; do
  O="$OUT/chr${c}.cohort.vcf.gz"
  [ -s "$O" ] || continue
  n=$(bcftools view -H "$O" 2>/dev/null | wc -l); tot=$((tot+n))
done
echo "   $tot SNPs"
du -sh "$OUT" 2>/dev/null

echo "########## STEP 0 CONCLUIDO ##########"
date
