"""config.py — LocusConfig: descreve UM locus para o motor `hap_ancestry`.

O motor e agnostico ao locus. Tudo que muda entre projetos — sitios focais,
paineis e sua ordem, cores/rotulos, textos das figuras, limiares e o mapa de
consequencias — vem daqui. Um projeto novo escreve um YAML e reusa o mesmo motor,
sem tocar no codigo.

Formato do YAML: ver `configs/slick.yaml` (o caso de referencia) e o README.
"""
from __future__ import annotations
from dataclasses import dataclass
import yaml


@dataclass
class LocusConfig:
    locus_name: str            # ex.: "PRLR" — so para titulo/rotulo
    allele_label: str          # ex.: "slick allele" — texto da legenda da rede
    focal_sites: list          # [(nome, pos), ...]  ORDEM preservada
    panels: list               # paineis candidatos p/ ATRIBUICAO; ORDEM = desempate
    panel_style: dict          # painel -> {"color":..., "label":...} (rede; a ORDEM
                               #   deste dict e a ordem de montagem da rede)
    consequence_classes: dict  # token do csq (minusculo) -> nome de classe
    class_report: list         # [(label, of), ...]; of == "*" (gene inteiro) ou [classes]
    coding_classes: list       # classes que contam como "codificante" (rede)
    eps: float = 1e-4
    min_hap_panel: int = 8
    min_sites_class: int = 3
    title: str = ""            # suptitle da rede; se vazio, montado de locus_name

    @classmethod
    def from_yaml(cls, path):
        d = yaml.safe_load(open(path, encoding="utf-8"))
        th = d.get("thresholds", {})
        focal = [(s["name"], int(s["pos"])) for s in d["focal_sites"]]
        report = [(c["label"], c["of"]) for c in d["class_report"]]
        return cls(
            locus_name=d["locus_name"],
            allele_label=d.get("allele_label", ""),
            focal_sites=focal,
            panels=list(d["panels"]),
            panel_style=dict(d.get("panel_style", {})),
            consequence_classes={str(k).lower(): v
                                 for k, v in d["consequence_classes"].items()},
            class_report=report,
            coding_classes=list(d.get("coding_classes", [])),
            eps=float(th.get("eps", 1e-4)),
            min_hap_panel=int(th.get("min_hap_panel", 8)),
            min_sites_class=int(th.get("min_sites_class", 3)),
            title=d.get("title", ""),
        )
