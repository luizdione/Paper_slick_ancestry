# -*- coding: utf-8 -*-
"""dados_figuras_novas.py -- reune em JSON os dados da Figura 2 (evidencia nas leituras) e das
Figuras S1-S3 do suplementar. Nada e digitado: tudo sai de arquivos de resultado ou do VCF.

  leituras     suporte de leitura de cada chamada conferida (suporte_leitura_portadores.tsv)
  fase         leituras do portador composto ERR12035149 sobre SLICK6 e SLICK2
               (fase_leituras_ERR12035149.tsv, gerado por fase_por_leitura.py --tsv)
  racas        portadores (animais) de cada alelo por raca e animais triados por raca, para as racas
               Criollo; raca do ERP150979 pelo codigo do sample_title no ENA (LMN Limonero, BON Blanco
               Orejinegro, HAR/HAC Harton del Valle -- HAC01 conferido no ENA: "Harton"), do Run 8 pela
               raca declarada no BioSample
  slick5       portadores da troca do SLICK5 (39.099.228) por populacao
  introns      atribuicao intronica por faixa de distancia (resultado_final.json)
  homozigotos  homozigotos do Run 8 por alelo (entrada do controle de endogamia, Figura S3)

Saida: tables/derived/figuras_novas.json e tables/derived/homozigotos_<ALELO>_run8.txt
"""
import collections
import csv
import json
import os
import pathlib
import sys

RAIZ = pathlib.Path(__file__).resolve().parents[2]   # raiz do repositorio
DATA = pathlib.Path(os.environ.get("SLICK_DATA_DIR", RAIZ / "data"))   # dados brutos (VCF, BAM), fora do repositorio
B = RAIZ / "tables" / "derived"   # entradas e saidas desta etapa
sys.path.insert(0, str(RAIZ / "scripts" / "lib"))
from hap_ancestry import load_phased  # noqa: E402

VCF = str(DATA / "criollo/incorporado/merged_phased_dedup.vcf.gz")
ENA = DATA / "p4_slick/data/criollo/ena_ERP150979_raw.tsv"   # filereport do ENA
ROT = RAIZ / "tables" / "rotulos_corrigidos_2026-10-13.tsv"
GRP = RAIZ / "tables" / "groups_iberico_declarado.tsv"
FOCO = {"SLICK1": 39099213, "SLICK3": 39099226, "SLICK6": 39099267, "SLICK2": 39099321}
SLICK5 = 39099228
COD = {"LMN": "Limonero", "BON": "Blanco Orejinegro", "HAR": "Hartón del Valle",
       "HAC": "Hartón del Valle"}
NORMA = {"limonero": "Limonero", "costeno con cuernos": "Costeño con Cuernos",
         "romosinuano": "Romosinuano", "blanco orejinegro": "Blanco Orejinegro",
         "harton del valle": "Hartón del Valle", "san martinero": "San Martinero",
         "sanmartinero": "San Martinero", "caracu": "Caracu", "senepol": "Senepol"}

fin = json.loads((B / "resultado_final.json").read_text(encoding="utf-8"))
sup = list(csv.DictReader(open(B / "suporte_leitura_portadores.tsv", encoding="utf-8"), delimiter="\t"))
fase = list(csv.DictReader(open(B / "fase_leituras_ERR12035149.tsv", encoding="utf-8"), delimiter="\t"))

# rotulos e racas
rot = {r["sample"]: r for r in csv.DictReader(open(ROT, encoding="utf-8"), delimiter="\t")}
ena = {r["run_accession"]: r for r in csv.DictReader(open(ENA, encoding="utf-8"), delimiter="\t")}
painel = {}
for l in open(GRP, encoding="utf-8").read().splitlines():
    c = l.split("\t")
    if len(c) >= 2:
        painel[c[0]] = c[1].strip()


def raca(s):
    if s in ena:
        return COD.get(ena[s]["sample_title"][:3], ena[s]["sample_title"][:3])
    b = (rot.get(s, {}).get("breed_biosample") or "").strip()
    return NORMA.get(b.lower(), b) if b else None


pos, H, samples = load_phased(VCF)
idx = {p: i for i, p in enumerate(pos)}


def gts(p):
    i = idx[p]
    return {s: (int(H[i, 2 * k]), int(H[i, 2 * k + 1])) for k, s in enumerate(samples)}


# portadores retidos (ja sem as chamadas rejeitadas nas leituras)
retidos = {a: sorted({h["amostra"] for h in fin["alelos"][a]["por_haplotipo"]}) for a in FOCO}
criollo_racas = {"Limonero", "Costeño con Cuernos", "Hartón del Valle", "Blanco Orejinegro",
                 "Romosinuano", "San Martinero"}
triados = collections.Counter(raca(s) for s in samples if raca(s) in criollo_racas)
por_raca = {a: dict(collections.Counter(raca(s) for s in retidos[a])) for a in FOCO}

# SLICK5 por populacao
g5 = gts(SLICK5)
car5 = {s: sum(x == 1 for x in g) for s, g in g5.items() if 1 in g}
pop5 = collections.Counter()
hap5 = collections.Counter()
for s, n in car5.items():
    if painel.get(s):
        k = {"EUR_TAURINE": "European taurine", "CRIOLLO": "Criollo", "NDAMA": "N'Dama",
             "OGADEN": "Ogaden", "BORAN": "Boran", "ANKOLE": "Ankole", "IBERIAN": "Iberian"}.get(painel[s], painel[s])
    else:
        r = raca(s)
        k = (r or "no declared breed")
    pop5[k] += 1
    hap5[k] += n
# polo genetico dos sem painel (taurino x indicino), para colorir
polo5 = collections.Counter()
for s in car5:
    if not painel.get(s):
        polo5[(rot.get(s, {}).get("genetic_pole") or "unknown")] += 1

# homozigotos do Run 8 (amostras SAMN/SAMEA/SRS, nao ERR do lote novo)
homo = {}
for a in ("SLICK1", "SLICK2"):
    g = gts(FOCO[a])
    homo[a] = sorted(s for s, x in g.items() if x == (1, 1) and not s.startswith("ERR"))
    (B / f"homozigotos_{a}_run8.txt").write_text("\n".join(homo[a]) + "\n", encoding="utf-8")

introns = {a: fin["alelos"][a]["introns_por_distancia"] for a in FOCO}
out = {"leituras": sup, "fase": fase, "racas": {"triados": dict(triados), "portadores": por_raca},
       "slick5": {"animais": dict(pop5), "haplotipos": dict(hap5), "n_animais": len(car5),
                  "n_haplotipos": int(sum(car5.values())), "polo_sem_painel": dict(polo5)},
       "introns": introns, "homozigotos_run8": homo, "n_paineis": 7}
(B / "figuras_novas.json").write_text(json.dumps(out, ensure_ascii=False, indent=1), encoding="utf-8")
print("triados:", dict(triados))
print("portadores por raca:", por_raca)
print("SLICK5:", len(car5), "animais,", int(sum(car5.values())), "haplotipos;", dict(pop5))
print("homozigotos Run 8:", homo)
