# -*- coding: utf-8 -*-
"""04a_rotulos.py -- paineis de referencia e consultas do RFMix (roda no host, env datascience).

Fontes dos rotulos, nesta ordem:
  - EUR_TAURINE, NDAMA, CRIOLLO, ANKOLE, BORAN, OGADEN e os 4 ibericos do Run 8: exatamente os
    da atribuicao por haplotipo, lidos PELO MOTOR (`load_groups`), com o arquivo como esta
    (ver README: nao converter fim de linha);
  - ibericos novos (Asturiana, `criollo/roster_wgs_iberian.tsv`): so os com >= 40% de chamada;
  - INDICINE: `rotulos_corrigidos_2026-10-13.tsv`, polo indicino, fracao indicina >= 0,95,
    racas asiaticas/americanas (africanas ficam de fora: Boran, Ogaden e Shorthorn Zebu sao
    mistura taurino x indicino).

Paineis grandes sao subamostrados a 40 animais (semente fixa); os que sobram viram CONTROLE:
rodam como consulta e tem de recuperar o proprio painel (regra 1 do projeto). Dos paineis
pequenos (N'Dama, iberico) separa-se um subconjunto fixo para o mesmo controle numa execucao a
parte, sem eles na referencia.
"""
import os
import csv
import pathlib
import random
import sys

RAIZ = RAIZ = pathlib.Path(__file__).resolve().parents[2]   # raiz do repositorio
DATA = pathlib.Path(os.environ.get("SLICK_DATA_DIR", RAIZ / "data"))   # dados brutos (VCF, BAM), fora do repositorio
sys.path.insert(0, str(RAIZ / "scripts" / "lib"))
from hap_ancestry import load_groups  # noqa: E402

O = (DATA / "lai10mb")
GRP = RAIZ / "tables" / "groups_iberico_declarado.tsv"
ROT = (RAIZ / "tables" / "rotulos_corrigidos_2026-10-13.tsv")
MAX_REF, SEMENTE, MIN_CR_REF = 40, 20261026, 40.0
AFRICANAS = {"Boran", "Ogaden", "Shorthorn Zebu"}

no_vcf = set((O / ".run8.txt").read_text().split())
no_vcf |= set((O / "tmp" / "novos_ok.txt").read_text().split())
no_vcf |= set((O / "tmp" / "ibericos_ok.txt").read_text().split())
cr = {l.split("\t")[0]: float(l.split("\t")[3])
      for l in (O / "genotipa_stats.tsv").read_text().splitlines() if l.strip()}
cr.setdefault("ERR12035105", 87.9)

g = {s: p for s, p in load_groups(str(GRP)).items() if s in no_vcf}
for s in ("ERR4280263", "ERR4280261", "ERR4280259", "ERR4280271"):   # Asturiana
    if s in no_vcf and cr.get(s, 0) >= MIN_CR_REF:
        g[s] = "IBERIAN"
with open(ROT, encoding="utf-8") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        try:
            fi = float(r["frac_indicine"])
        except ValueError:
            continue
        if (r["genetic_pole"] == "indicine_pole" and fi >= 0.95 and r["sample"] in no_vcf
                and r["breed_biosample"] not in AFRICANAS and r["sample"] not in g):
            g[r["sample"]] = "INDICINE"

def baixa(s):          # animal novo com chamada baixa nao entra como referencia
    return s in cr and cr[s] < MIN_CR_REF

rng = random.Random(SEMENTE)
ref, ctrl, consulta = {}, {}, {}
for pan in ("EUR_TAURINE", "INDICINE", "NDAMA", "IBERIAN"):
    s = sorted(x for x, p in g.items() if p == pan and not baixa(x))
    rng.shuffle(s)
    if len(s) > MAX_REF:
        ref.update({x: pan for x in s[:MAX_REF]})
        ctrl.update({x: pan for x in s[MAX_REF:MAX_REF + 10]})
    else:
        ref.update({x: pan for x in s})
for x, p in g.items():
    if p in ("CRIOLLO", "ANKOLE", "BORAN", "OGADEN"):
        consulta[x] = p
# controle dos paineis pequenos: 3 N'Dama e 2 ibericos saem da referencia na execucao-controle
peq = {}
for pan, k in (("NDAMA", 3), ("IBERIAN", 2)):
    s = sorted(x for x, p in ref.items() if p == pan)
    peq.update({x: pan for x in rng.sample(s, k)})

def grava(nome, d):
    with open(O / nome, "w", encoding="utf-8", newline="\n") as f:
        for x in sorted(d):
            f.write(f"{x}\t{d[x]}\n")

grava("ref_4vias.tsv", ref)
grava("ref_3vias.tsv", {x: ("TAURINO_EUROPEU" if p in ("EUR_TAURINE", "IBERIAN") else
                            "TAURINO_AFRICANO" if p == "NDAMA" else p) for x, p in ref.items()})
grava("controle_grandes.tsv", ctrl)
grava("controle_pequenos.tsv", peq)
grava("consulta.tsv", consulta)
from collections import Counter
print("referencia:", Counter(ref.values()))
print("controle (sobra dos grandes):", Counter(ctrl.values()))
print("controle (pequenos, saem da referencia):", Counter(peq.values()))
print("consulta:", Counter(consulta.values()))
