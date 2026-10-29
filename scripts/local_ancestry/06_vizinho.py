# -*- coding: utf-8 -*-
"""06_vizinho.py -- vizinho mais proximo de cada haplotipo portador, sem modelo.

O RFMix chama o SLICK2 de INDICINO no PRLR (17 de 22 haplotipos). Uma chamada de ancestralidade
local e sempre "o painel menos distante"; se o haplotipo do SLICK2 nao existir em nenhuma
referencia, o modelo ainda escolhe um. Este script mede a distancia diretamente: para cada
haplotipo portador, o numero de diferencas (Hamming) ate o haplotipo mais parecido de cada
painel e de cada RACA, em janelas centradas no bloco slick, contra TODOS os 1.882 animais
faseados (nao so as referencias do RFMix). Se o vizinho mais proximo for um zebu asiatico, o
fundo e indicino de fato; se for Brahman, pode ser segmento taurino americano no Brahman; se for
outro Criollo nao portador, o haplotipo e local e nao tem par nas referencias.

Saida: `$SLICK_DATA_DIR/lai10mb/rfmix/vizinho.json` + impressao.
"""
import os
import collections
import csv
import gzip
import json
import pathlib

import numpy as np

RAIZ = pathlib.Path(__file__).resolve().parents[2]   # raiz do repositorio
DATA = pathlib.Path(os.environ.get("SLICK_DATA_DIR", RAIZ / "data"))   # dados brutos (VCF, BAM), fora do repositorio

O = (DATA / "lai10mb")
R = O / "rfmix"
ROT = (RAIZ / "tables" / "rotulos_corrigidos_2026-10-13.tsv")
CENTRO = 39099267
JANELAS = (25_000, 50_000, 100_000, 250_000)
FOCO = {"SLICK1": 39099213, "SLICK3": 39099226, "SLICK6": 39099267, "SLICK2": 39099321}
SEM_SUPORTE = {"SLICK2": {"ERR12035125", "ERR12035133"}, "SLICK3": {"ERR12035159"},
               "SLICK6": {"ERR12035135"}}


def le_tsv(p):
    return dict(l.split("\t")[:2] for l in p.read_text().split("\n") if l.strip())


painel = {**le_tsv(O / "consulta.tsv"), **le_tsv(O / "controle_grandes.tsv"),
          **le_tsv(O / "ref_4vias.tsv")}
raca = {}
with open(ROT, encoding="utf-8") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        raca[r["sample"]] = r["breed_biosample"] or r["category_2026_10_13"]
for s, p in painel.items():
    raca.setdefault(s, p)

with gzip.open(R / "hap_1mb_prlr.tsv.gz", "rt") as f:
    amostras = f.readline().rstrip("\n").split("\t")[1:]
    pos, linhas = [], []
    for l in f:
        c = l.rstrip("\n").split("\t")
        pos.append(int(c[0]))
        linhas.append([int(a) for g in c[1:] for a in g.split("|")])
pos = np.array(pos)
H = np.array(linhas, dtype=np.int8)                       # sitios x haplotipos
hap_amostra = [s for s in amostras for _ in (0, 1)]
hap_idx = [i for _ in amostras for i in (0, 1)]
print(f"{H.shape[0]} sitios x {H.shape[1]} haplotipos")

fpos = {p: int(np.where(pos == p)[0][0]) for p in FOCO.values() if p in pos}
portadores = collections.defaultdict(list)
for a, p in FOCO.items():
    for h in np.where(H[fpos[p]] == 1)[0]:
        s = hap_amostra[h]
        if painel.get(s) == "CRIOLLO" and s not in SEM_SUPORTE.get(a, set()):
            portadores[a].append(int(h))
todos_port = set(h for v in portadores.values() for h in v)

rel = {}
for w in JANELAS:
    sel = np.where((pos >= CENTRO - w) & (pos <= CENTRO + w) & ~np.isin(pos, list(FOCO.values())))[0]
    Hs = H[sel]
    rel[w] = {}
    print(f"\n######## janela +-{w // 1000} kb ({len(sel)} SNPs)")
    for a, lst in portadores.items():
        rel[w][a] = []
        for h in lst:
            d = (Hs != Hs[:, [h]]).sum(axis=0)
            d[h] = 10**9
            # o outro haplotipo do mesmo animal e os demais portadores do MESMO alelo saem:
            # queremos o parente mais proximo FORA da linhagem portadora
            for h2 in portadores[a]:
                d[h2] = 10**9
            ordem = np.argsort(d)
            por_painel = {}
            for p in ("EUR_TAURINE", "IBERIAN", "NDAMA", "INDICINE", "CRIOLLO"):
                idx = [k for k in range(len(d)) if painel.get(hap_amostra[k]) == p and d[k] < 10**9]
                if idx:
                    k = min(idx, key=lambda k: d[k])
                    por_painel[p] = {"dist": int(d[k]), "amostra": hap_amostra[k], "raca": raca.get(hap_amostra[k], "?")}
            top = [{"dist": int(d[k]), "amostra": hap_amostra[k], "raca": raca.get(hap_amostra[k], "?"),
                    "painel": painel.get(hap_amostra[k], "-"), "portador_outro_alelo": int(k) in todos_port}
                   for k in ordem[:10]]
            rel[w][a].append({"amostra": hap_amostra[h], "hap": hap_idx[h], "por_painel": por_painel, "top10": top})
            if w in (50_000, 250_000):
                pp = " ".join(f"{p[:4]}:{v['dist']}" for p, v in por_painel.items())
                t3 = ", ".join(f"{x['raca']}({x['dist']})" for x in top[:4])
                print(f"  {a:<6} {hap_amostra[h]:<13} h{hap_idx[h]} | min por painel {pp} | top: {t3}")

(R / "vizinho.json").write_text(json.dumps(rel, ensure_ascii=False, indent=1), encoding="utf-8")
print(f"-> {R / 'vizinho.json'}")
