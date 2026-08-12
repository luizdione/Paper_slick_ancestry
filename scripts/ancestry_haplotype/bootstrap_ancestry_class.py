#!/usr/bin/env python3
"""
bootstrap_ancestry_class.py — intervalo de confianca e estabilidade da atribuicao
do haplotipo portador, por classe funcional de sitio.

Wrapper fino sobre o motor `hap_ancestry` (ver `gene_ancestry_by_class.py`). O
intervalo NAO vem de reamostrar haplotipos, e sim os SITIOS dentro de cada
classe: a pergunta certa e quao estavel e a populacao que cada haplotipo
escolhe, dado este conjunto de sitios informativos.

Duas saidas por classe e por alelo: a fracao de haplotipos no painel vencedor, com
IC95 percentil; e a fracao de replicas em que aquele painel venceu o consenso, que
e uma medida direta de estabilidade da resposta.
"""
from __future__ import annotations
import argparse, pathlib, sys
import numpy as np
import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "lib"))
from hap_ancestry import (                    # noqa: E402
    LocusConfig, load_phased, load_groups, hap_group_labels, active_panels,
    load_consequences, resolve_classes, focal_index, carrier_haplotypes,
    panel_frequencies, stability,
)

DEFAULT_CONFIG = pathlib.Path(__file__).resolve().parents[2] / "config" / "slick.yaml"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--csq", required=True)
    ap.add_argument("--groups", required=True)
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--config", default=str(DEFAULT_CONFIG))
    ap.add_argument("--reps", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=20260804)
    a = ap.parse_args()

    cfg = LocusConfig.from_yaml(a.config)
    rng = np.random.default_rng(a.seed)
    pos, H, samples = load_phased(a.vcf)
    hap_grp = hap_group_labels(samples, load_groups(a.groups))
    panels = active_panels(hap_grp, cfg.panels, cfg.min_hap_panel)
    class_sets = load_consequences(a.csq, cfg.consequence_classes)
    classes = resolve_classes(pos, class_sets, cfg.class_report)

    linhas = []
    for rotulo, foco in cfg.focal_sites:
        i_foc, exato = focal_index(pos, foco)
        if not exato:
            print(f"[{rotulo}] sitio focal ausente — pulando")
            continue
        port = carrier_haplotypes(H, i_foc)

        carrier_mask = np.zeros(H.shape[1], bool)
        carrier_mask[port] = True
        F = panel_frequencies(H, hap_grp, panels, carrier_mask, cfg.eps)

        print(f"\n=== {rotulo} — {len(port)} haplotipos portadores ===")
        for nome, sel in classes:
            if len(sel) < cfg.min_sites_class:
                continue
            r = stability(H, sel, port, F, panels, cfg.min_sites_class, a.reps, rng)
            if r is None:
                continue
            print(f"   {nome:18s} {len(sel):5d} sitios | {r['vencedor']:15s} "
                  f"{r['frac_obs']:.2f} [IC95 {r['ic_lo']:.2f}–{r['ic_hi']:.2f}] | "
                  f"estabilidade {r['estabilidade']:.3f}")
            linhas.append({"alelo": rotulo, "classe": nome, "n_sitios": len(sel),
                           "n_hap": r["n_hap"], "painel": r["vencedor"],
                           "fracao": r["frac_obs"], "ic_lo": r["ic_lo"],
                           "ic_hi": r["ic_hi"], "estabilidade": r["estabilidade"],
                           "reps": r["reps"]})

    out = pathlib.Path(f"{a.out_prefix}_bootstrap_classe.tsv")
    pd.DataFrame(linhas).to_csv(out, sep="\t", index=False, float_format="%.4f")
    print(f"\n  -> {out}")


if __name__ == "__main__":
    main()
