# -*- coding: utf-8 -*-
"""O sitio codificante informa mais sobre ancestralidade que o intronico? Dois testes.

A atribuicao por classe funcional supoe que o codificante preserva o fundo fundador e que o
intron acumula o fundo local. Este script mede as duas coisas.

TESTE A — informatividade POR SITIO.
  Pergunta: um sitio codificante carrega mais informacao de ancestralidade que um intronico?
  Metrica: informatividade de Rosenberg (In), ja calculada por sitio em aim_*_landscape.tsv.
  Armadilhas controladas:
    - FREQUENCIA. In cresce com a frequencia do alelo menor; as classes tem espectros
      diferentes, entao a media crua compara frequencia, nao classe. Aqui cada sitio
      codificante e pareado a um intronico de MAF proxima (+-`JANELA_MAF`).
    - TAMANHO. 22 codificantes polimorficos contra 625 introns: a media do grupo pequeno
      e muito mais ruidosa. O nulo vem de reamostragem do grupo grande no tamanho do pequeno.
    - CIRCULARIDADE. Os sitios slick focais sao eles proprios nao-sinonimos; entram como
      `nonsyn` e inflariam a classe codificante. Sao removidos.

TESTE B — o intron recupera o fundo LOCAL?
  A hipotese de fundo nao e sobre informatividade por sitio: e sobre
  RECOMBINACAO ter trocado o fundo fundador pelo fundo da raca em que o alelo segrega hoje.
  Isso faz uma previsao falsificavel: a atribuicao por introns deve cair no painel da
  PROPRIA populacao do portador, nao numa fonte ancestral. Aqui essa previsao e conferida
  contra o rotulo de painel de cada animal portador.

Saida: tabela legivel no stdout + `tables/derived/informatividade_classe.json`.
"""
from __future__ import annotations

import pathlib
import io
import json
import os
import subprocess
import sys

import numpy as np

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

RAIZ = pathlib.Path(__file__).resolve().parents[2]   # raiz do repositorio
DATA = pathlib.Path(os.environ.get("SLICK_DATA_DIR", RAIZ / "data"))   # dados brutos (VCF, BAM), fora do repositorio
IMAGEM = os.environ.get("SLICK_TOOLS_IMAGE", "slick-tools:latest")   # container com bcftools
HAPLO = str(DATA / "haplo")
VCF = str(DATA / "criollo/incorporado/merged_phased_dedup.vcf.gz")
GRUPOS = str(RAIZ / "tables" / "groups_afr_split.tsv")
SAIDA = str(RAIZ / "tables" / "derived" / "informatividade_classe.json")

# Posicoes focais slick (OMIA:001372-9913, assembly ARS-UCD1.2/1.3).
FOCAIS = {39099113: "SLICK4", 39099213: "SLICK1", 39099226: "SLICK3",
          39099228: "SLICK5", 39099267: "SLICK6", 39099321: "SLICK2"}
JANELA_MAF = 0.02
N_REP = 10000
SEMENTE = 20261020


def maf_por_posicao() -> dict[int, float]:
    """MAF observada no painel faseado, por posicao.

    (fix 2026-10-20) O cabecalho do VCF DECLARA `AF`, mas os registros gravados pelo BEAGLE
    carregam `AC` e `AN`. Pedir `%AF` devolve "." em toda linha, o float() falha em silencio
    e o dicionario sai VAZIO: o teste pareado some sem erro, dizendo "sitios insuficientes".
    Campo declarado no cabecalho nao e campo presente no registro.
    """
    cmd = ["docker", "run", "--rm", "-v", f"{DATA}:/data",
           IMAGEM, "bash", "-lc",
           "bcftools query -f '%POS\\t%AC\\t%AN\\n' "
           "/data/criollo/incorporado/merged_phased_dedup.vcf.gz 2>/dev/null"]
    r = subprocess.run(cmd, capture_output=True, text=True)
    out: dict[int, float] = {}
    for linha in r.stdout.splitlines():
        c = linha.split("\t")
        if len(c) < 3:
            continue
        try:
            pos, ac, an = int(c[0]), float(c[1].split(",")[0]), float(c[2])
        except ValueError:
            continue
        if an <= 0:
            continue
        f = ac / an
        out[pos] = min(f, 1.0 - f)
    if not out:
        raise SystemExit("MAF vazia: confira os campos do VCF antes de confiar no teste")
    return out


def carrega_landscape(alelo: str) -> list[tuple[int, str, float]]:
    cam = os.path.join(HAPLO, f"aim_{alelo}_coding_landscape.tsv")
    linhas = []
    with io.open(cam, encoding="utf-8") as fh:
        next(fh)
        for l in fh:
            c = l.rstrip("\n").split("\t")
            if len(c) < 4:
                continue
            try:
                linhas.append((int(c[0]), c[1], float(c[3])))
            except ValueError:
                continue
    return linhas


def teste_a(alelo: str, maf: dict[int, float], rng) -> dict:
    dados = carrega_landscape(alelo)
    cod, intr = [], []
    for pos, classe, inn in dados:
        if pos in FOCAIS:                      # circularidade
            continue
        m = maf.get(pos, 0.0)
        if m <= 0.0:                           # monomorfico nao informa nada
            continue
        if classe in ("syn", "nonsyn"):
            cod.append((m, inn))
        elif classe == "intron":
            intr.append((m, inn))
    if len(cod) < 5 or len(intr) < 20:
        return {"alelo": alelo, "erro": "sitios insuficientes"}

    mc = float(np.mean([x[1] for x in cod]))
    mi = float(np.mean([x[1] for x in intr]))

    # (1) nulo por TAMANHO: subamostra introns no n dos codificantes
    iv = np.array([x[1] for x in intr])
    sub = np.array([rng.choice(iv, size=len(cod), replace=False).mean()
                    for _ in range(N_REP)])
    p_tam = float((np.sum(sub >= mc) + 1) / (N_REP + 1))

    # (2) nulo PAREADO POR MAF: para cada codificante, um intron de MAF proxima
    im = np.array([x[0] for x in intr])
    pares = []
    for m, _ in cod:
        cand = np.where(np.abs(im - m) <= JANELA_MAF)[0]
        pares.append(cand if cand.size else np.array([], dtype=int))
    cobertos = sum(1 for c in pares if c.size)
    if cobertos < len(cod) * 0.6:
        pareado = {"aviso": f"so {cobertos}/{len(cod)} codificantes tem par de MAF"}
    else:
        difs = []
        for _ in range(N_REP):
            esc = [iv[rng.choice(c)] for c in pares if c.size]
            difs.append(mc - float(np.mean(esc)))
        difs = np.array(difs)
        pareado = {
            "delta_In_medio": float(difs.mean()),
            "ic95": [float(np.percentile(difs, 2.5)), float(np.percentile(difs, 97.5))],
            "p_bicaudal": float(2 * min((np.sum(difs <= 0) + 1) / (N_REP + 1),
                                        (np.sum(difs >= 0) + 1) / (N_REP + 1))),
            "n_pareados": cobertos,
        }
    return {
        "alelo": alelo,
        "n_codificante_polimorfico": len(cod),
        "n_intron_polimorfico": len(intr),
        "In_medio_codificante": mc,
        "In_medio_intron": mi,
        "p_subamostragem_por_tamanho": p_tam,
        "pareado_por_maf": pareado,
    }


def teste_b() -> dict:
    """A atribuicao por intron cai no painel do PROPRIO portador?"""
    grupo = {}
    with io.open(GRUPOS, encoding="utf-8") as fh:
        for l in fh:
            c = l.rstrip("\n").split("\t")
            if len(c) >= 2:
                grupo[c[0]] = c[1]
    # portadores por alelo, lidos do painel
    port = {}
    for nome, pos in (("SLICK1", 39099213), ("SLICK2", 39099321),
                      ("SLICK3", 39099226), ("SLICK6", 39099267)):
        cmd = ["docker", "run", "--rm", "-v", f"{DATA}:/data",
               IMAGEM, "bash", "-lc",
               f"bcftools query -r 20:{pos}-{pos} -f '[%SAMPLE=%GT\\n]' "
               f"/data/criollo/incorporado/merged_phased_dedup.vcf.gz 2>/dev/null"]
        r = subprocess.run(cmd, capture_output=True, text=True)
        amostras = []
        for l in r.stdout.splitlines():
            s, _, gt = l.partition("=")
            if "1" in gt:
                amostras.append(s)
        port[nome] = amostras
    return {a: {"portadores": v, "painel_do_portador":
                sorted({grupo.get(s, "SEM_ROTULO") for s in v})}
            for a, v in port.items()}


def main() -> None:
    rng = np.random.default_rng(SEMENTE)
    print("=" * 78)
    print("INFORMATIVIDADE DE ANCESTRALIDADE: CODIFICANTE vs INTRON")
    print("=" * 78)
    print("\nlendo MAF do painel faseado...", flush=True)
    maf = maf_por_posicao()
    print(f"  {len(maf)} posicoes com MAF")

    rel = {"janela_maf": JANELA_MAF, "n_rep": N_REP, "semente": SEMENTE, "alelos": []}
    print("\nTESTE A — informatividade de Rosenberg por sitio\n")
    cab = f"  {'alelo':<14}{'n_cod':>7}{'n_int':>7}{'In_cod':>9}{'In_int':>9}{'p_tam':>9}{'delta_pareado':>15}{'p_par':>9}"
    print(cab)
    print("  " + "-" * (len(cab) - 2))
    for alelo in ("SLICK1", "SLICK2", "SLICK3_cand"):
        r = teste_a(alelo, maf, rng)
        rel["alelos"].append(r)
        if "erro" in r:
            print(f"  {alelo:<14}{r['erro']}")
            continue
        pr = r["pareado_por_maf"]
        d = f"{pr['delta_In_medio']:+.5f}" if "delta_In_medio" in pr else "n/d"
        p = f"{pr['p_bicaudal']:.3f}" if "p_bicaudal" in pr else "n/d"
        print(f"  {alelo:<14}{r['n_codificante_polimorfico']:>7}"
              f"{r['n_intron_polimorfico']:>7}{r['In_medio_codificante']:>9.5f}"
              f"{r['In_medio_intron']:>9.5f}{r['p_subamostragem_por_tamanho']:>9.3f}"
              f"{d:>15}{p:>9}")

    print("\nTESTE B — a atribuicao por intron cai no painel do proprio portador?\n")
    b = teste_b()
    rel["teste_b"] = b
    for alelo, v in b.items():
        print(f"  {alelo:<8} portadores em painel: {', '.join(v['painel_do_portador'])}"
              f"   (n={len(v['portadores'])})")
    print("\n  Comparar com a linha 'so introns' da atribuicao por classe: se a atribuicao")
    print("  intronica reproduz o painel do portador, o intron mede fundo LOCAL, nao fonte.")

    io.open(SAIDA, "w", encoding="utf-8").write(
        json.dumps(rel, ensure_ascii=False, indent=1))
    print(f"\n-> {SAIDA}")


if __name__ == "__main__":
    main()
