# -*- coding: utf-8 -*-
"""Controle por remocao: as posicoes 39.077.999 e 39.078.000 do exon 2 decidem a
atribuicao codificante? Imprime a frequencia por painel, a correlacao entre os dois sitios e a
atribuicao de cada alelo com e sem eles."""
import os
import pathlib
import sys, pathlib, collections
import numpy as np
RAIZ = RAIZ = pathlib.Path(__file__).resolve().parents[2]   # raiz do repositorio
DATA = pathlib.Path(os.environ.get("SLICK_DATA_DIR", RAIZ / "data"))   # dados brutos (VCF, BAM), fora do repositorio
sys.path.insert(0, str(RAIZ / "scripts" / "lib")); sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
from hap_ancestry import *
from resultado_final import VCF, CSQ, CFG, GRP, COD, excluidos
cfg = LocusConfig.from_yaml(str(CFG)); pos, H, samples = load_phased(VCF)
hg = hap_group_labels(samples, load_groups(str(GRP))); panels = active_panels(hg, cfg.panels, cfg.min_hap_panel)
cs = load_consequences(CSQ, cfg.consequence_classes)
cls = dict(resolve_classes(pos, cs, cfg.class_report)); EX = excluidos()
alvo = [i for i,p in enumerate(pos) if int(p) in (39077999, 39078000)]
print("indices", alvo, [int(pos[i]) for i in alvo])
for k,v in cs.items():
    hit=[p for p in (39077999,39078000) if p in v]
    if hit: print("classe csq", k, hit)
todos_port = np.zeros(H.shape[1], bool)
for r,f in cfg.focal_sites:
    i,ok = focal_index(pos,f)
    if ok: todos_port[carrier_haplotypes(H,i)] = True
for i in alvo:
    print(f"\n{int(pos[i])}: freq ALT por painel (nao portadores)")
    for p in panels:
        m=(hg==p)&~todos_port; x=H[i,m]; x=x[x>=0]
        print(f"   {p:<12} {x.mean():.3f} (n={x.size})")
    print("   corr com o outro sitio:", np.corrcoef(H[alvo[0]], H[alvo[1]])[0,1].round(4))
for r,f in cfg.focal_sites:
    i,ok = focal_index(pos,f)
    if not ok: continue
    tod = np.asarray(carrier_haplotypes(H,i))
    if tod.size==0: continue
    cm = np.zeros(H.shape[1],bool); cm[tod]=True
    port = np.array([h for h in tod if samples[h//2] not in EX.get(r,set())])
    F = panel_frequencies(H, hg, panels, cm, cfg.eps)
    print(f"\n== {r} ({port.size} hap); alelo nos 2 sitios:", collections.Counter(tuple(H[alvo,h]) for h in port))
    for c in COD:
        sel = np.asarray(cls[c]); sem = np.array([s for s in sel if s not in alvo])
        a,_ = assign_with_posterior(H, sel, port, F, panels, cfg.min_sites_class)
        b,pb = assign_with_posterior(H, sem, port, F, panels, cfg.min_sites_class)
        print(f"   {c:<18} com: {dict(a)} | SEM os 2 sitios: {dict(b)} post.med {np.median(pb):.2f}")
