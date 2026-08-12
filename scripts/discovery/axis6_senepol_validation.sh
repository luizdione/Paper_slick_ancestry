#!/bin/bash
# ============================================================================
# EIXO 6 - VALIDACAO ORTOGONAL dos dois alelos truncantes por RNA-seq de Senepol
#
# Por que isto vale: os BAMs de RNA-seq sao evidencia INDEPENDENTE do call set do
# Run 8. Tecnologia diferente (RNA x DNA), biblioteca diferente, alinhador
# diferente (HISAT2 splice-aware x BWA+GATK), amostras diferentes, e — o que mais
# importa — os animais tem GENOTIPO SLICK1 PUBLICADO (Sosa et al. 2022,
# GSE192471), o que permite medir concordancia contra uma verdade externa.
#
# Reaproveita: os BAMs ja produzidos por senepol_prlr_genotype.sh em
# /data/slick_ancestry/p4_senepol/bam/*.prlr.bam e a imagem p4-tools.
#
# LICOES JA APRENDIDAS, aplicadas aqui:
#  - bcftools mpileup pula chamada de indel acima de DP 250 por amostra
#    (--max-idepth). Sem -L explicito o SLICK1, que e uma delecao, nao e chamado
#    nas amostras de alta cobertura.
#  - A fracao alelica ingenua (ALT/DP) e enganosa em RNA-seq: a coluna DP inclui
#    leituras que atravessam o sitio DENTRO de um operador N do CIGAR (splice),
#    que nao carregam informacao sobre a base. Aqui a fracao e calculada sobre
#    AD, que conta so leituras informativas.
# ============================================================================
set -euo pipefail

REF=/data/genome_reference/ARS-UCD1.3/Bos_taurus.ARS-UCD1.3.dna.toplevel.fa
BAMDIR=/data/slick_ancestry/p4_senepol/bam
MAP=/proj/data/GSE192471/senepol_run_genotype_map.tsv
R=/proj/results/senepol_validacao
THREADS=${THREADS:-3}
mkdir -p "$R"
log() { printf '[%s] %s\n' "$(date '+%H:%M:%S')" "$*"; }

SLICK1_POS=39099213   # GC>G, c.1382del, p.(Ala461ValfsTer2)
SLICK2_POS=39099321   # C>T,  c.1489C>T, p.(Arg497*)
REG="20:${SLICK1_POS}-${SLICK2_POS}"

log "########## EIXO 6 - validacao ortogonal por RNA-seq ##########"

ls "$BAMDIR"/*.prlr.bam > "$R/bams.txt" 2>/dev/null || { echo "ERRO: sem BAMs" >&2; exit 1; }
NB=$(wc -l < "$R/bams.txt")
log "BAMs disponiveis: $NB"
[ "$NB" -ge 1 ] || exit 1

# ---------- 1. chamada conjunta nos dois sitios ----------
if [ ! -s "$R/senepol_2sitios.vcf.gz" ]; then
  log "chamada conjunta (mpileup -L 1000000, sem corte de profundidade)"
  bcftools mpileup -f "$REF" -r "$REG" -b "$R/bams.txt" \
      -a AD,DP -q 20 -Q 20 -L 1000000 -d 100000 --threads "$THREADS" -Ou \
    | bcftools call -m -A -Oz -o "$R/senepol_2sitios.vcf.gz" --threads "$THREADS"
  bcftools index -f "$R/senepol_2sitios.vcf.gz"
fi
log "variantes na janela: $(bcftools view -H "$R/senepol_2sitios.vcf.gz" | wc -l)"

# ---------- 2. evidencia por amostra em cada sitio ----------
log "extraindo evidencia por amostra"
T="$R/evidencia_por_amostra.tsv"
printf 'RUN\tSITIO\tALELO\tGT\tAD_ref\tAD_alt\tDP\tFRACAO_ALT\n' > "$T"
for pos in "$SLICK1_POS" "$SLICK2_POS"; do
  nome=$([ "$pos" = "$SLICK1_POS" ] && echo SLICK1 || echo SLICK2)
  bcftools query -r "20:${pos}-${pos}" \
    -f '[%SAMPLE\t%GT\t%AD\t%DP\n]' "$R/senepol_2sitios.vcf.gz" 2>/dev/null \
  | awk -v OFS='\t' -v p="$pos" -v n="$nome" '{
      # AD vem como "ref,alt[,alt2...]"; a fracao usa SO leituras informativas
      split($3,a,","); r=a[1]+0; alt=0; for(i=2;i<=length(a);i++) alt+=a[i]+0;
      tot=r+alt; f=(tot>0)? alt/tot : 0;
      run=$1; sub(/.*\//,"",run); sub(/\.prlr\.bam$/,"",run);
      printf "%s\t%s\t%s\t%s\t%d\t%d\t%s\t%.4f\n", run, p, n, $2, r, alt, $4, f
    }' >> "$T"
done
log "linhas: $(( $(wc -l < "$T") - 1 ))"

# ---------- 3. concordancia com o genotipo publicado ----------
log "cruzando com o genotipo SLICK1 publicado (Sosa et al. 2022)"
python3 - "$T" "$MAP" "$R" <<'PY'
import csv, sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
ev, mapa, outdir = sys.argv[1], sys.argv[2], sys.argv[3]

# rotulo publicado por run. FREE=0 copias, SLICK=1 (heterozigoto!), S1=1, S2=2.
pub = {}
with open(mapa, encoding='utf-8', errors='replace') as f:
    for r in csv.DictReader(f):
        run = (r.get('Run') or '').strip().strip('"')
        g   = (r.get('Genotype') or '').strip().strip('"')
        if run: pub[run] = g
ESPERADO = {'FREE': 0.00, 'SLICK': 0.50, 'S1': 0.50, 'S2': 1.00}

linhas = list(csv.DictReader(open(ev, encoding='utf-8'), delimiter='\t'))
print(f"\n{'RUN':<13} {'ROTULO':<6} {'ALELO':<7} {'GT':<6} {'ref':>5} {'alt':>5} "
      f"{'frac':>7} {'esperado':>9}  veredito")
print('-'*82)

conc = disc = 0
for l in sorted(linhas, key=lambda x: (x['ALELO'], x['RUN'])):
    run = l['RUN']; rot = pub.get(run, '?')
    f = float(l['FRACAO_ALT']); alt = int(l['AD_alt']); tot = int(l['AD_ref']) + alt
    esp = ESPERADO.get(rot)
    v = ''
    if l['ALELO'] == 'SLICK1' and esp is not None:
        if tot < 5:
            v = 'cobertura insuficiente'
        elif esp == 0.0:
            v = 'CONCORDA' if f < 0.05 else 'DISCORDA (rotulo FREE, alelo presente)'
        elif esp == 1.0:
            v = 'CONCORDA' if f > 0.85 else 'DISCORDA (rotulo S2, nao homozigoto)'
        else:
            v = 'CONCORDA' if 0.15 < f < 0.85 else 'DISCORDA (rotulo het)'
        if v.startswith('CONCORDA'): conc += 1
        elif v.startswith('DISCORDA'): disc += 1
    print(f"{run:<13} {rot:<6} {l['ALELO']:<7} {l['GT']:<6} {l['AD_ref']:>5} "
          f"{l['AD_alt']:>5} {f:>7.3f} {('-' if esp is None else f'{esp:.2f}'):>9}  {v}")

print(f"\nCONCORDANCIA no SLICK1 contra o genotipo publicado: {conc} de {conc+disc}")

# O ponto da validacao ortogonal: o alelo aparece em RNA, tecnologia independente?
s1 = [l for l in linhas if l['ALELO']=='SLICK1' and int(l['AD_alt'])>0]
s2 = [l for l in linhas if l['ALELO']=='SLICK2' and int(l['AD_alt'])>0]
print(f"\nEVIDENCIA ORTOGONAL (leituras de RNA carregando o alelo alternativo):")
print(f"  SLICK1: {len(s1)} amostras com ALT>0; total de leituras ALT = "
      f"{sum(int(l['AD_alt']) for l in s1)}")
print(f"  SLICK2: {len(s2)} amostras com ALT>0; total de leituras ALT = "
      f"{sum(int(l['AD_alt']) for l in s2)}")
print("\n  (o Run 8 sustenta o SLICK1 em 11 leituras ALT somadas nos 4 portadores;")
print("   qualquer numero comparavel aqui vem de tecnologia e alinhador diferentes)")
PY

log "########## EIXO 6 CONCLUIDO ##########"
ls -la "$R"
