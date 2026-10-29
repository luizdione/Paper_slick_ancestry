"""hap_ancestry — motor guarda-chuva de ancestralidade de haplotipo por classe
funcional de sitio.

Dado (1) um VCF faseado de um gene/regiao, (2) as consequencias funcionais dos
sitios, (3) grupos/paineis de referencia e (4) uma lista de sitios focais, o motor
diz de qual populacao vem o haplotipo que carrega cada alelo focal, estratificado
por classe funcional, com intervalo por reamostragem de sitios e uma rede de
haplotipos.

O motor NAO conhece nenhum locus: tudo que e especifico (sitios, paineis, cores,
limiares, mapa de consequencias) vem de um `LocusConfig` carregado de YAML. O caso
*slick*/*PRLR* e apenas `configs/slick.yaml`. Para um projeto novo, escreva um YAML
e reuse este mesmo motor.

API publica:
    LocusConfig.from_yaml(path)
    load_phased / load_consequences / load_groups / hap_group_labels
    active_panels / resolve_classes
    focal_index / carrier_haplotypes / panel_frequencies
    winners / assign_with_posterior
    stability (bootstrap)
    plot_network / collapse / mds_layout
"""
from __future__ import annotations

from .config import LocusConfig
from .io import (
    load_phased, load_consequences, load_groups, hap_group_labels,
    active_panels, resolve_classes,
)
from .assign import (
    focal_index, carrier_haplotypes, panel_frequencies,
    winners, assign_with_posterior, consensus,
)
from .bootstrap import stability
from .network import plot_network, collapse, mds_layout, draw_pie

__all__ = [
    "LocusConfig",
    "load_phased", "load_consequences", "load_groups", "hap_group_labels",
    "active_panels", "resolve_classes",
    "focal_index", "carrier_haplotypes", "panel_frequencies",
    "winners", "assign_with_posterior", "consensus",
    "stability",
    "plot_network", "collapse", "mds_layout", "draw_pie",
]
