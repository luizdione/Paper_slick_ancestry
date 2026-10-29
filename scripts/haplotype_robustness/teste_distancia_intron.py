# -*- coding: utf-8 -*-
"""O fundo fundador se apaga com a DISTANCIA ao sitio causal? Teste direto.

Contexto. A atribuicao ancora a ancestralidade no nucleo codificante e trata os introns como
ruido local. O teste de informatividade (`teste_informatividade_classe.py`) mostrou que,
pareados por frequencia alelica, sitios codificantes NAO carregam mais informacao de
ancestralidade que introns: a vantagem bruta de In era artefato do espectro de frequencias.

Resta uma hipotese que nao depende de informatividade: PROXIMIDADE. Os tres sitios slick ficam no exon terminal, dentro de ~110 pb
um do outro. Se um seculo de recombinacao trocou o fundo fundador pelo fundo local, a troca
tem de ser GRADUAL NA DISTANCIA: introns vizinhos ao sitio causal ainda devem acompanhar a
atribuicao codificante, e introns distantes devem convergir para o painel da raca em que o
alelo segrega hoje.

Essa previsao e falsificavel, e e o que este script mede: a mesma atribuicao do motor, rodada
sobre introns binados por |pos - sitio focal|.

Leitura do resultado:
  - atribuicao que MUDA com a distancia, do painel codificante para o painel do portador,
    sustenta o argumento de recombinacao;
  - atribuicao CONSTANTE em toda a distancia nao sustenta: aponta para diferenca de poder
    entre classes, nao para erosao do fundo.
"""
from __future__ import annotations

import os
import io
import json
import pathlib
import sys

import numpy as np

RAIZ = RAIZ = pathlib.Path(__file__).resolve().parents[2]   # raiz do repositorio
DATA = pathlib.Path(os.environ.get("SLICK_DATA_DIR", RAIZ / "data"))   # dados brutos (VCF, BAM), fora do repositorio
sys.path.insert(0, str(RAIZ / "scripts" / "lib"))
sys.stdout.reconfigure(encoding="utf-8", errors="replace")

from hap_ancestry import (                                        # noqa: E402
    LocusConfig, load_phased, load_consequences, load_groups,
    hap_group_labels, active_panels,
    focal_index, carrier_haplotypes, panel_frequencies, assign_with_posterior,
)

VCF = str(DATA / "criollo/incorporado/merged_phased_dedup.vcf.gz")
CSQ = str(DATA / "haplo/prlr_allmaf_csq.tsv")
GRUPOS = RAIZ / "tables" / "groups_afr_split.tsv"
CONFIG = RAIZ / "config" / "slick_com_slick6.yaml"
SAIDA = RAIZ / "tables" / "derived" / "distancia_intron.json"

# Bordas em pares de base. A ultima e aberta.
BORDAS = [0, 1_000, 5_000, 10_000, 20_000, 10**9]
ROTULOS = ["<1 kb", "1-5 kb", "5-10 kb", "10-20 kb", ">20 kb"]


def main() -> None:
    cfg = LocusConfig.from_yaml(str(CONFIG))
    pos, H, samples = load_phased(VCF)
    hap_grp = hap_group_labels(samples, load_groups(str(GRUPOS)))
    panels = active_panels(hap_grp, cfg.panels, cfg.min_hap_panel)
    class_sets = load_consequences(CSQ, cfg.consequence_classes)

    intron_pos = class_sets.get("intron", set())
    cod_pos = class_sets.get("syn", set()) | class_sets.get("nonsyn", set())
    idx_intron = np.array([i for i, p in enumerate(pos) if int(p) in intron_pos])
    idx_cod = np.array([i for i, p in enumerate(pos) if int(p) in cod_pos])
    print(f"[painel] {len(pos):,} sitios | {len(idx_intron)} introns | "
          f"{len(idx_cod)} codificantes | paineis: {', '.join(panels)}")

    rel = {"bordas_pb": BORDAS, "alelos": []}
    for rotulo, foco in cfg.focal_sites:
        i_foc, exato = focal_index(pos, foco)
        if not exato:
            print(f"\n[{rotulo}] sitio focal ausente, pulando")
            continue
        port = carrier_haplotypes(H, i_foc)
        if len(port) == 0:
            print(f"\n[{rotulo}] sem portadores, pulando")
            continue
        animais = sorted({samples[h // 2] for h in port})
        cmask = np.zeros(H.shape[1], bool)
        cmask[port] = True
        F = panel_frequencies(H, hap_grp, panels, cmask, cfg.eps)

        print(f"\n=== {rotulo} ({foco:,}) — {len(port)} haplotipos em {len(animais)} animais ===")

        # referencia: o nucleo codificante, ancora da atribuicao
        cont, postv = assign_with_posterior(H, idx_cod, port, F, panels, cfg.min_sites_class)
        ref = " ".join(f"{k}:{v}" for k, v in sorted(cont.items(), key=lambda kv: -kv[1]))
        print(f"   {'codificante (ancora)':<22}{len(idx_cod):>6} sitios | {ref}"
              f"   (post. mediana {np.median(postv):.2f})")

        dist = np.abs(pos[idx_intron].astype(np.int64) - int(foco))
        linhas = []
        for k, rot in enumerate(ROTULOS):
            sel = idx_intron[(dist >= BORDAS[k]) & (dist < BORDAS[k + 1])]
            if len(sel) < cfg.min_sites_class:
                print(f"   {('intron ' + rot):<22}{len(sel):>6} sitios | (poucos sitios)")
                linhas.append({"faixa": rot, "n_sitios": int(len(sel)), "atribuicao": None})
                continue
            cont, postv = assign_with_posterior(H, sel, port, F, panels, cfg.min_sites_class)
            ordem = sorted(cont.items(), key=lambda kv: -kv[1])
            txt = " ".join(f"{k2}:{v}" for k2, v in ordem)
            med = float(np.median(postv)) if postv else float("nan")
            print(f"   {('intron ' + rot):<22}{len(sel):>6} sitios | {txt}"
                  f"   (post. mediana {med:.2f})")
            linhas.append({"faixa": rot, "n_sitios": int(len(sel)),
                           "atribuicao": {k2: int(v) for k2, v in ordem},
                           "posterior_mediana": med})
        rel["alelos"].append({"alelo": rotulo, "pos": int(foco),
                              "n_haplotipos": int(len(port)),
                              "ancora_codificante": ref, "introns_por_distancia": linhas})

    io.open(SAIDA, "w", encoding="utf-8").write(json.dumps(rel, ensure_ascii=False, indent=1))
    print(f"\n-> {SAIDA}")


if __name__ == "__main__":
    main()
