#!/usr/bin/env python3
"""
haplotype_network.py — rede dos haplotipos do locus, no espaco dos sitios
codificantes. Wrapper fino sobre o motor `hap_ancestry`.

no = haplotipo distinto; tamanho = numero de copias; fatia = de que painel vem
cada copia; aresta = vizinho mais proximo. Copias identicas colapsam num no — e
esse colapso e informacao: um alelo introduzido uma vez deixa portadores quase
identicos; um recorrente os deixaria dispersos.

So os sitios codificantes (sinonimos + nao sinonimos); os introns ficam de fora de
proposito — carregam o fundo local recombinado, onde todo haplotipo criollo cai
junto de criollo (resultado conhecido, nao a pergunta).
"""
from __future__ import annotations
import argparse, pathlib, sys
import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "lib"))
from hap_ancestry import (                    # noqa: E402
    LocusConfig, load_phased, load_groups, hap_group_labels,
    load_consequences, plot_network,
)

DEFAULT_CONFIG = pathlib.Path(__file__).resolve().parents[2] / "config" / "slick.yaml"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--csq", required=True)
    ap.add_argument("--groups", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--config", default=str(DEFAULT_CONFIG))
    a = ap.parse_args()

    cfg = LocusConfig.from_yaml(a.config)
    pos, H, samples = load_phased(a.vcf)
    hap_grp = hap_group_labels(samples, load_groups(a.groups))

    class_sets = load_consequences(a.csq, cfg.consequence_classes)
    coding: set = set()
    for c in cfg.coding_classes:
        coding |= class_sets.get(c, set())
    coding_idx = np.where(np.isin(pos, list(coding)))[0]

    plot_network(pos, H, hap_grp, coding_idx, cfg, a.out)


if __name__ == "__main__":
    main()
