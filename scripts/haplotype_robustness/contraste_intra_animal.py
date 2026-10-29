# -*- coding: utf-8 -*-
"""Contraste DENTRO do mesmo animal: cada cromossomo de um portador composto em trans.

Se os dois cromossomos de um mesmo animal saem com afinidades diferentes, lote, biblioteca,
profundidade e pipeline ficam constantes; so a ancestralidade do haplotipo muda. E o controle
mais forte que o desenho permite.

Fase: ERR12035149 (SLICK6 / SLICK2) confirmada por LEITURA (`fase_por_leitura.py`: 0 leituras
com os dois alelos em 6). SAMN08166129 (SLICK1 / SLICK2) so por BEAGLE: Run 8, sem BAM local.
"""
import os
import json
import pathlib
import sys

import numpy as np

RAIZ = RAIZ = pathlib.Path(__file__).resolve().parents[2]   # raiz do repositorio
DATA = pathlib.Path(os.environ.get("SLICK_DATA_DIR", RAIZ / "data"))   # dados brutos (VCF, BAM), fora do repositorio
sys.path.insert(0, str(RAIZ / "scripts" / "lib"))
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
from hap_ancestry import (  # noqa: E402
    LocusConfig, load_phased, load_consequences, load_groups, hap_group_labels,
    active_panels, resolve_classes, focal_index, carrier_haplotypes,
    panel_frequencies, winners,
)

CFG = RAIZ / "config" / "slick_iberico_declarado.yaml"
GRP = RAIZ / "tables" / "groups_iberico_declarado.tsv"
VCF = str(DATA / "criollo/incorporado/merged_phased_dedup.vcf.gz")
CSQ = str(DATA / "haplo/prlr_allmaf_csq.tsv")
CASOS = {"ERR12035149": ("SLICK6", "SLICK2"), "SAMN08166129": ("SLICK1", "SLICK2")}

cfg = LocusConfig.from_yaml(str(CFG))
pos, H, samples = load_phased(VCF)
hg = hap_group_labels(samples, load_groups(str(GRP)))
panels = active_panels(hg, cfg.panels, cfg.min_hap_panel)
cls = dict(resolve_classes(pos, load_consequences(CSQ, cfg.consequence_classes), cfg.class_report))
focal = dict(cfg.focal_sites)

# mascara: todos os portadores de qualquer alelo slick fora das frequencias de referencia
todos = np.zeros(H.shape[1], bool)
for _, f in cfg.focal_sites:
    i, ok = focal_index(pos, f)
    if ok:
        todos[carrier_haplotypes(H, i)] = True
F = panel_frequencies(H, hg, panels, todos, cfg.eps)

saida = {}
for amostra, (a1, a2) in CASOS.items():
    k = samples.index(amostra)
    haps = (2 * k, 2 * k + 1)
    print(f"\n=== {amostra} ===")
    saida[amostra] = {}
    for h in haps:
        leva = [a for a in (a1, a2) if H[focal_index(pos, focal[a])[0], h] == 1]
        linha = {}
        for c in ("so sinonimas", "so nao sinonimas"):
            w = winners(H, cls[c], [h], F, panels, cfg.min_sites_class)
            linha[c] = w[0] if w else None
        saida[amostra][f"hap{h % 2}"] = {"carrega": leva, **linha}
        print(f"   cromossomo {h % 2}: carrega {('+'.join(leva) or 'nenhum'):<14} -> "
              f"sin {linha['so sinonimas']:<12} nao-sin {linha['so nao sinonimas']}")

json.dump(saida, open(RAIZ / "tables" / "derived" / "contraste_intra_animal.json", "w", encoding="utf-8"),
          ensure_ascii=False, indent=1)
