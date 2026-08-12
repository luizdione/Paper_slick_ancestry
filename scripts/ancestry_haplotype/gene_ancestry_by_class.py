#!/usr/bin/env python3
"""
gene_ancestry_by_class.py — de qual populacao vem o haplotipo que carrega cada
alelo focal, medido por classe funcional de sitio. Produz a tabela de atribuicao
por classe funcional e os dados de entrada para a figura de ancestralidade por
classe.

Este arquivo e um WRAPPER FINO sobre o motor guarda-chuva `hap_ancestry`. Toda a
logica (que antes vivia aqui, cravando SLICK1/SLICK2) foi extraida para o motor; o
que era especifico do *slick* virou `config/slick.yaml`. Rodar sem `--config`
reproduz exatamente o comportamento historico (validado contra baseline). Para
outro gene/projeto, escreva um YAML e aponte `--config` para ele — o motor e o
mesmo.

A estratificacao por classe e o teste. Sitios sinonimos sao neutros por
construcao: se a atribuicao por sinonimas concorda com a por nao sinonimas, ela
nao pode ser explicada por selecao atuando sobre o haplotipo. Os introns entram
como contraste — sao os flancos, onde a recombinacao ja teve um seculo para trocar
o fundo fundador pelo local.
"""
from __future__ import annotations
import argparse, pathlib, sys
import numpy as np
import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "lib"))
from hap_ancestry import (                    # noqa: E402
    LocusConfig, load_phased, load_consequences, load_groups,
    hap_group_labels, active_panels, resolve_classes,
    focal_index, carrier_haplotypes, panel_frequencies, assign_with_posterior,
)
# re-exporta o carregamento p/ quem importava daqui (compat), agora vindo do motor
from hap_ancestry.io import load_phased, load_consequences  # noqa: E402,F811

DEFAULT_CONFIG = pathlib.Path(__file__).resolve().parents[2] / "config" / "slick.yaml"

# constantes historicas, mantidas para retrocompatibilidade de imports externos
EPS = 1e-4
SLICK1_POS = 39_099_213
SLICK2_POS = 39_099_321
MIN_HAP_PANEL = 8
MIN_SITES_CLASS = 3


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--vcf", required=True, help="gene faseado, sem filtro de MAF")
    ap.add_argument("--csq", required=True, help="TSV do bcftools csq")
    ap.add_argument("--groups", required=True, help="amostra<TAB>painel")
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--config", default=str(DEFAULT_CONFIG),
                    help="YAML do locus (default: slick.yaml)")
    a = ap.parse_args()

    cfg = LocusConfig.from_yaml(a.config)
    pos, H, samples = load_phased(a.vcf)
    hap_grp = hap_group_labels(samples, load_groups(a.groups))

    panels = active_panels(hap_grp, cfg.panels, cfg.min_hap_panel)
    print(f"[gene] {len(pos):,} sitios | paineis: " +
          " ".join(f"{p}({int((hap_grp == p).sum())})" for p in panels))

    class_sets = load_consequences(a.csq, cfg.consequence_classes)
    classes = resolve_classes(pos, class_sets, cfg.class_report)
    counts = {nome: len(s) for nome, s in class_sets.items()}
    fora = len(pos) - sum(counts.values())
    print("[classes] " + " | ".join(f"{k} {v}" for k, v in counts.items()) +
          f" | fora (UTR/splice/nao anotados) {fora}")

    linhas = []
    for rotulo, foco in cfg.focal_sites:
        i_foc, exato = focal_index(pos, foco)
        if not exato:
            # controle positivo: o sitio focal TEM de estar no VCF faseado; se nao,
            # a fase nao esta resolvida ali e o resultado seria inferido.
            print(f"\n[{rotulo}] SITIO FOCAL AUSENTE (mais proximo: {pos[i_foc]:,}) — pulando")
            continue

        port = carrier_haplotypes(H, i_foc)
        animais = sorted({samples[h // 2] for h in port})
        print(f"\n=== {rotulo} — {len(port)} haplotipos portadores em {len(animais)} animais ===")
        print(f"    {', '.join(animais)}")

        carrier_mask = np.zeros(H.shape[1], bool)
        carrier_mask[port] = True
        F = panel_frequencies(H, hap_grp, panels, carrier_mask, cfg.eps)

        for nome, sel in classes:
            if len(sel) < cfg.min_sites_class:
                continue
            cont, post_v = assign_with_posterior(H, sel, port, F, panels, cfg.min_sites_class)
            ordem = sorted(cont.items(), key=lambda kv: -kv[1])
            med = float(np.median(post_v)) if post_v else float("nan")
            print(f"   {nome:18s} {len(sel):5d} sitios | " +
                  " ".join(f"{k}:{v}" for k, v in ordem) + f"   (post. mediana {med:.2f})")
            linhas.append({"alelo": rotulo, "classe": nome, "n_sitios": len(sel),
                           "n_hap": sum(cont.values()),
                           "painel": ordem[0][0] if ordem else "NA",
                           "n_painel": ordem[0][1] if ordem else 0,
                           "post_mediana": med})

    out = pathlib.Path(f"{a.out_prefix}_por_classe.tsv")
    pd.DataFrame(linhas).to_csv(out, sep="\t", index=False, float_format="%.4f")
    print(f"\n  -> {out}")


if __name__ == "__main__":
    main()
