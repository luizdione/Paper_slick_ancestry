#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""suporte_leitura.py -- confere cada chamada de portador contra as leituras alinhadas.

O faseamento preenche genotipos ausentes, e um genotipo preenchido num sitio raro pode criar um
portador que nenhuma leitura sustenta. Para cada portador de cada alelo focal no VCF faseado,
conta as leituras de referencia e alternativas no BAM (samtools mpileup, MAPQ e qualidade de
base >= 20) e classifica:

  sem_leitura_imputado        nenhuma leitura no sitio                       -> nao usado
  duvidoso                    < 3 leituras alternativas ou fracao < 0,25     -> nao usado
  heterozigoto                fracao alternativa < 0,85                      -> usado
  homozigoto                  fracao >= 0,85, fase 1|1 e >= 4 leituras       -> usado
  homozigoto_baixa_cobertura  fracao >= 0,85, fase 1|1 e < 4 leituras        -> usado
  homozigoto_por_leitura      fracao >= 0,85 com fase heterozigota           -> usado

Portadores sem BAM (o genotipo publicado e tudo o que existe) nao podem ser conferidos e sao
listados a parte, para que a contagem total fique explicita.

Uso (com samtools e bcftools no PATH):
  suporte_leitura.py <vcf_faseado> <dir_bams> <fasta_isca> <saida.tsv>
       [--contig 20:34100000-44100000] [--offset 34099999]
Os BAMs estao em coordenada LOCAL da isca: global = local + offset.
"""
import argparse
import pathlib
import re
import subprocess
import sys

FOCAIS = {  # alelo: (posicao global ARS-UCD1.2, ref, alt) -- substituicoes de base unica
    "SLICK2": (39099321, "C", "T"),
    "SLICK3": (39099226, "C", "A"),
    "SLICK6": (39099267, "C", "T"),
}


def portadores(vcf, pos):
    out = subprocess.run(["bcftools", "query", "-r", f"20:{pos}", "-f", "[%SAMPLE\t%GT\n]", vcf],
                         capture_output=True, text=True, check=True).stdout
    return [(s, g) for s, g in (l.split("\t") for l in out.splitlines() if l) if "1" in g]


def conta(bam, fasta, contig, local, alt):
    r = subprocess.run(["samtools", "mpileup", "-f", fasta, "-r", f"{contig}:{local}-{local}",
                        "-B", "-q", "20", "-Q", "20", bam], capture_output=True, text=True)
    linhas = [l for l in r.stdout.splitlines() if l.strip()]
    if not linhas:
        return 0, 0
    bases = linhas[0].split("\t")[4]
    bases = re.sub(r"\^.", "", bases).replace("$", "")
    bases = re.sub(r"[+-](\d+)", lambda m: "#" * 0, bases)   # marcador de indel; o texto segue
    limpo, i = [], 0
    while i < len(bases):                                    # remove o texto do indel
        if bases[i] in "+-":
            j = i + 1
            while j < len(bases) and bases[j].isdigit():
                j += 1
            i = j + int(bases[i + 1:j])
            continue
        limpo.append(bases[i])
        i += 1
    s = "".join(limpo)
    return s.count(".") + s.count(","), s.upper().count(alt)


def classifica(ref, alt, gt):
    if ref + alt == 0:
        return "sem_leitura_imputado", "nao", "NA"
    frac = alt / (ref + alt)
    if alt < 3 or frac < 0.25:
        return "duvidoso", "nao", f"{frac:.2f}"
    if frac >= 0.85:
        if gt.replace("/", "|") == "1|1":
            return ("homozigoto" if ref + alt >= 4 else "homozigoto_baixa_cobertura"), "sim", f"{frac:.2f}"
        return "homozigoto_por_leitura", "sim", f"{frac:.2f}"
    return "heterozigoto", "sim", f"{frac:.2f}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("vcf")
    ap.add_argument("bams")
    ap.add_argument("fasta")
    ap.add_argument("saida")
    ap.add_argument("--contig", default="20:34100000-44100000")
    ap.add_argument("--offset", type=int, default=34099999)
    a = ap.parse_args()
    bams = pathlib.Path(a.bams)
    linhas, sem_bam = [], []
    for alelo, (pos, _ref, alt) in FOCAIS.items():
        for amostra, gt in portadores(a.vcf, pos):
            bam = bams / f"{amostra}.bam"
            if not bam.exists():
                sem_bam.append((alelo, amostra))
                continue
            r, x = conta(str(bam), a.fasta, a.contig, pos - a.offset, alt)
            sup, usado, frac = classifica(r, x, gt)
            linhas.append([alelo, amostra, "ERP150979", gt, str(r), str(x), frac, sup, usado])
    cab = ["alelo", "animal", "lote", "genotipo_beagle", "leituras_ref", "leituras_alt",
           "fracao_alt", "suporte", "usado_na_atribuicao"]
    pathlib.Path(a.saida).write_text("\n".join("\t".join(l) for l in [cab] + linhas) + "\n",
                                     encoding="utf-8")
    print(f"{len(linhas)} chamadas conferidas por leitura; {len(sem_bam)} sem BAM (nao conferiveis)")
    print(f"nao usadas: {sum(1 for l in linhas if l[8] == 'nao')}")


if __name__ == "__main__":
    sys.exit(main())
