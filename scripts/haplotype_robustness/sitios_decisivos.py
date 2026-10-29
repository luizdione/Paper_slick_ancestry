# -*- coding: utf-8 -*-
"""Quantos sitios sustentam cada atribuicao? Contribuicao de cada sitio a razao de verossimilhanca.

A reamostragem de sitios (estabilidade_hierarquica.py) mostrou que o "europeu" do SLICK2 vence
em so 59-74% das replicas, e o IC95 da fracao de haplotipos vai de 0 a 1. A explicacao proposta e
que a atribuicao da configuracao codificante dominante (17 dos 22 haplotipos) repousa em poucos
sitios: reamostrar com reposicao deixa cada sitio de fora em ~37% das replicas, e se dois ou
tres sitios decidem tudo, a chamada cai junto com eles.

Para cada configuracao codificante distinta de cada alelo, e para cada classe, lista a
contribuicao de cada sitio ao log da razao de verossimilhanca entre o melhor painel do grupo
vencedor e o melhor painel do outro grupo (africano x europeu), e conta quantos sitios carregam
|log10 RV| > 0,5 (um fator de ~3).
"""
from __future__ import annotations

import os
import collections
import io
import json
import pathlib
import sys

import numpy as np

RAIZ = RAIZ = pathlib.Path(__file__).resolve().parents[2]   # raiz do repositorio
DATA = pathlib.Path(os.environ.get("SLICK_DATA_DIR", RAIZ / "data"))   # dados brutos (VCF, BAM), fora do repositorio
sys.path.insert(0, str(RAIZ / "scripts" / "lib"))
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
from hap_ancestry import (  # noqa: E402
    LocusConfig, load_phased, load_consequences, load_groups, hap_group_labels,
    active_panels, resolve_classes, focal_index, carrier_haplotypes, panel_frequencies,
)
from resultado_final import VCF, CSQ, CFG, GRP, COD, excluidos  # noqa: E402

AFR = ("NDAMA", "ANKOLE", "BORAN", "OGADEN")
EUR = ("EUR_TAURINE", "IBERIAN")
LIMIAR = 0.5
SAIDA = RAIZ / "tables" / "derived" / "sitios_decisivos.json"


def main():
    cfg = LocusConfig.from_yaml(str(CFG))
    pos, H, samples = load_phased(VCF)
    hg = hap_group_labels(samples, load_groups(str(GRP)))
    panels = active_panels(hg, cfg.panels, cfg.min_hap_panel)
    cls = dict(resolve_classes(pos, load_consequences(CSQ, cfg.consequence_classes),
                               cfg.class_report))
    EX = excluidos()
    rel = {"limiar_log10": LIMIAR, "alelos": {}}
    for rotulo, foco in cfg.focal_sites:
        i_foc, exato = focal_index(pos, foco)
        if not exato:
            continue
        todos = np.asarray(carrier_haplotypes(H, i_foc))
        if todos.size == 0:
            continue
        cmask = np.zeros(H.shape[1], bool)
        cmask[todos] = True
        port = np.array([h for h in todos if samples[h // 2] not in EX.get(rotulo, set())])
        F = panel_frequencies(H, hg, panels, cmask, cfg.eps)
        rel["alelos"][rotulo] = {}
        for c in COD:
            sel = np.asarray(cls[c])
            confs = collections.Counter(tuple(H[sel, h].tolist()) for h in port)
            lista = []
            for conf, n in confs.most_common():
                obs = np.array(conf)
                ok = obs >= 0
                # log10 verossimilhanca por sitio, por painel
                L = {p: np.where(obs == 1, np.log10(F[p][sel]), np.log10(1 - F[p][sel]))
                     for p in panels}
                tot = {p: float(L[p][ok].sum()) for p in panels}
                ba = max(AFR, key=lambda p: tot.get(p, -np.inf))
                be = max(EUR, key=lambda p: tot.get(p, -np.inf))
                d = (L[be] - L[ba])
                d[~ok] = 0.0
                ordem = np.argsort(-np.abs(d))
                fortes = [{"pos": int(pos[sel[j]]), "alelo": int(obs[j]),
                           "log10RV_eur_menos_afr": round(float(d[j]), 3),
                           "freq": {p: round(float(F[p][sel[j]]), 3) for p in (be, ba)}}
                          for j in ordem if abs(d[j]) > LIMIAR]
                lista.append({"copias": n, "melhor_europeu": be, "melhor_africano": ba,
                              "log10RV_total": round(float(d.sum()), 3),
                              "n_sitios_chamados": int(ok.sum()),
                              "n_sitios_decisivos": len(fortes), "decisivos": fortes})
                print(f"{rotulo:<7}{c:<18} conf x{n:<3} {be} vs {ba}: log10RV {d.sum():+.2f} "
                      f"| {len(fortes)} sitios |d|>{LIMIAR}: "
                      + ", ".join(f"{f['pos']}({f['log10RV_eur_menos_afr']:+.2f})" for f in fortes))
            rel["alelos"][rotulo][c] = lista
    io.open(SAIDA, "w", encoding="utf-8").write(json.dumps(rel, ensure_ascii=False, indent=1))
    print(f"-> {SAIDA}")


if __name__ == "__main__":
    main()
