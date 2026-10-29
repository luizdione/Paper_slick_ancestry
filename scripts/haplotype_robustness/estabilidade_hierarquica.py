# -*- coding: utf-8 -*-
"""Estabilidade da atribuicao em DOIS niveis: continente e painel.

O `resultado_final.json` mede a estabilidade do painel vencedor exato (2.000 reamostragens de
sitios). Para o SLICK2 ela e baixa (36% e 43%), e o IC95 da fracao vai de 0 a 1: os 17
haplotipos que dividem a mesma configuracao codificante mudam de painel juntos. A pergunta que
isso deixa aberta e ONDE a atribuicao oscila. Se ela troca europeu continental por iberico, a
leitura "fundo europeu" esta firme e so a resolucao dentro da Europa e fraca. Se troca europeu
por N'Dama, a propria afinidade continental cai.

Aqui a mesma reamostragem guarda o vencedor de cada replica (nao so "venceu o observado?") e
agrega em tres grupos: africano (N'Dama, Ankole, Boran, Ogaden), europeu (continental, iberico)
e o fundo local (Criollo nao portador). Mesmos dados, rotulos, portadores e semente do
`resultado_final.py`.
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
    active_panels, resolve_classes, focal_index, carrier_haplotypes,
    panel_frequencies, winners,
)
from hap_ancestry.assign import consensus  # noqa: E402
from resultado_final import VCF, CSQ, CFG, GRP, COD, excluidos  # noqa: E402

GRUPO = {"NDAMA": "africano", "ANKOLE": "africano", "BORAN": "africano", "OGADEN": "africano",
         "EUR_TAURINE": "europeu", "IBERIAN": "europeu", "CRIOLLO": "local"}
REPS, SEMENTE = 2000, 20261027
SAIDA = RAIZ / "tables" / "derived" / "estabilidade_hierarquica.json"


def main():
    rng = np.random.default_rng(SEMENTE)
    cfg = LocusConfig.from_yaml(str(CFG))
    pos, H, samples = load_phased(VCF)
    hg = hap_group_labels(samples, load_groups(str(GRP)))
    panels = active_panels(hg, cfg.panels, cfg.min_hap_panel)
    cls = dict(resolve_classes(pos, load_consequences(CSQ, cfg.consequence_classes),
                               cfg.class_report))
    EX = excluidos()
    rel = {"reps": REPS, "semente": SEMENTE, "grupos": GRUPO, "alelos": {}}
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
            sel = cls[c]
            obs = winners(H, sel, port, F, panels, cfg.min_sites_class)
            v_obs = consensus(obs)
            painel, grupo = collections.Counter(), collections.Counter()
            frac_grupo = []
            for _ in range(REPS):
                s = rng.choice(sel, size=len(sel), replace=True)
                e = winners(H, s, port, F, panels, cfg.min_sites_class)
                if not e:
                    continue
                painel[consensus(e)] += 1
                grupo[GRUPO[consensus(e)]] += 1
                frac_grupo.append(sum(GRUPO[w] == GRUPO[v_obs] for w in e) / len(e))
            n = sum(painel.values())
            lo, hi = np.percentile(frac_grupo, [2.5, 97.5])
            rel["alelos"][rotulo][c] = {
                "vencedor_observado": v_obs, "grupo_observado": GRUPO[v_obs],
                "consenso_por_painel": {k: v / n for k, v in painel.most_common()},
                "consenso_por_grupo": {k: v / n for k, v in grupo.most_common()},
                "fracao_hap_no_grupo_ic95": [float(lo), float(hi)],
                "fracao_hap_no_grupo_mediana": float(np.median(frac_grupo)), "reps": n}
            print(f"{rotulo:<7}{c:<18} obs {v_obs:<12} | painel "
                  + " ".join(f"{k}:{100*v/n:.0f}%" for k, v in painel.most_common())
                  + " | grupo " + " ".join(f"{k}:{100*v/n:.0f}%" for k, v in grupo.most_common())
                  + f" | hap no grupo IC95 {lo:.2f}-{hi:.2f}")
    io.open(SAIDA, "w", encoding="utf-8").write(json.dumps(rel, ensure_ascii=False, indent=1))
    print(f"-> {SAIDA}")


if __name__ == "__main__":
    main()
