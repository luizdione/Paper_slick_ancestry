# -*- coding: utf-8 -*-
"""exon2_dados.py -- grava em JSON tudo o que a figura do exon 2 consome.

Nada aqui e novo: e o `controle_exon2.py` com saida estruturada. Grava:
  - frequencia do alelo alternativo das duas posicoes do exon 2 em cada painel (nao portadores);
  - correlacao entre as duas posicoes;
  - o alelo que cada haplotipo portador carrega nessas posicoes;
  - a atribuicao por classe COM e SEM as duas posicoes (controle por remocao);
  - a contribuicao de CADA sitio codificante ao log10 da razao de verossimilhanca europeu/africano,
    para a configuracao dominante de cada alelo.
Saida: `tables/derived/exon2_dados.json`.
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
    assign_with_posterior,
)
from resultado_final import VCF, CSQ, CFG, GRP, COD, excluidos  # noqa: E402

PAR = (39077999, 39078000)
AFR, EUR = ("NDAMA", "ANKOLE", "BORAN", "OGADEN"), ("EUR_TAURINE", "IBERIAN")
SAIDA = RAIZ / "tables" / "derived" / "exon2_dados.json"

cfg = LocusConfig.from_yaml(str(CFG))
pos, H, samples = load_phased(VCF)
hg = hap_group_labels(samples, load_groups(str(GRP)))
panels = active_panels(hg, cfg.panels, cfg.min_hap_panel)
cls = dict(resolve_classes(pos, load_consequences(CSQ, cfg.consequence_classes), cfg.class_report))
EX = excluidos()
alvo = [i for i, p in enumerate(pos) if int(p) in PAR]

todos_port = np.zeros(H.shape[1], bool)
for _, f in cfg.focal_sites:
    i, ok = focal_index(pos, f)
    if ok:
        todos_port[carrier_haplotypes(H, i)] = True

rel = {"posicoes": list(PAR), "classe": {str(PAR[0]): "nao sinonima", str(PAR[1]): "sinonima"},
       "correlacao": float(np.corrcoef(H[alvo[0]], H[alvo[1]])[0, 1]), "frequencia": {}, "alelos": {}}
for i in alvo:
    f = {}
    for p in panels:
        x = H[i, (hg == p) & ~todos_port]
        x = x[x >= 0]
        f[p] = {"freq_alt": float(x.mean()), "n_hap": int(x.size)}
    rel["frequencia"][str(int(pos[i]))] = f
print(f"correlacao entre as duas posicoes: {rel['correlacao']:.3f}")
for p in panels:
    print(f"  {p:<12} {rel['frequencia'][str(PAR[0])][p]['freq_alt']:.3f} "
          f"(n={rel['frequencia'][str(PAR[0])][p]['n_hap']})")

for rotulo, foco in cfg.focal_sites:
    i_foc, ok = focal_index(pos, foco)
    if not ok:
        continue
    tod = np.asarray(carrier_haplotypes(H, i_foc))
    if tod.size == 0:
        continue
    cm = np.zeros(H.shape[1], bool)
    cm[tod] = True
    port = np.array([h for h in tod if samples[h // 2] not in EX.get(rotulo, set())])
    F = panel_frequencies(H, hg, panels, cm, cfg.eps)
    r = {"n_hap": int(port.size),
         "alelo_exon2": {"ALT": int(sum(H[alvo[0], h] == 1 for h in port)),
                         "REF": int(sum(H[alvo[0], h] == 0 for h in port))},
         "classes": {}}
    for c in COD:
        sel = np.asarray(cls[c])
        sem = np.array([s for s in sel if s not in alvo])
        a, pa = assign_with_posterior(H, sel, port, F, panels, cfg.min_sites_class)
        b, pb = assign_with_posterior(H, sem, port, F, panels, cfg.min_sites_class)
        r["classes"][c] = {"com": dict(a), "sem": dict(b),
                           "posterior_com": float(np.median(pa)), "posterior_sem": float(np.median(pb))}
    # contribuicao por sitio na configuracao dominante
    cod = np.concatenate([np.asarray(cls[c]) for c in COD])
    confs = collections.Counter(tuple(H[cod, h].tolist()) for h in port)
    obs = np.array(confs.most_common(1)[0][0])
    n_dom = confs.most_common(1)[0][1]
    ok_m = obs >= 0
    L = {p: np.where(obs == 1, np.log10(F[p][cod]), np.log10(1 - F[p][cod])) for p in panels}
    tot = {p: float(L[p][ok_m].sum()) for p in panels}
    ba, be = max(AFR, key=lambda p: tot[p]), max(EUR, key=lambda p: tot[p])
    d = L[be] - L[ba]
    d[~ok_m] = 0.0
    r["configuracao_dominante"] = {
        "copias": n_dom, "melhor_europeu": be, "melhor_africano": ba,
        "log10RV_total": float(d.sum()),
        "sitios": [{"pos": int(pos[cod[j]]), "log10RV": float(d[j]),
                    "exon2": int(pos[cod[j]]) in PAR} for j in range(len(cod)) if ok_m[j]]}
    rel["alelos"][rotulo] = r
    print(f"{rotulo}: exon2 ALT {r['alelo_exon2']['ALT']}/{r['n_hap']}; "
          f"log10RV {d.sum():+.2f} ({be} vs {ba})")

io.open(SAIDA, "w", encoding="utf-8").write(json.dumps(rel, ensure_ascii=False, indent=1))
print(f"-> {SAIDA}")
