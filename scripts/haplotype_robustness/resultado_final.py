# -*- coding: utf-8 -*-
"""Fonte UNICA dos numeros de haplotipo: tabelas e figuras leem o JSON que este script grava.

Configuracao final, decidida em 2026-10-26:
  - paineis: os seis originais + IBERIAN (4 animais com raca DECLARADA no BioSample);
  - rotulos: `groups_iberico_declarado.tsv`, base lida PELO MOTOR (ver README);
  - portadores: so as chamadas com suporte de leitura (`suporte_leitura_portadores.tsv`);
  - alelos: SLICK1, SLICK2, SLICK3, SLICK6 (SLICK4 monomorfico; SLICK5 candidato nao confirmado).

Produz, num JSON so, tudo o que a Tabela 2 e as figuras 2-4 consomem:
  por alelo x classe: contagem, posterior mediana, estabilidade de bootstrap por sitio (2.000),
  painel modal sob rarefacao de todos os paineis ao menor (500), separacao por lote, n de
  configuracoes codificantes distintas, atribuicao por faixa de distancia ao sitio causal (introns),
  e a atribuicao de CADA haplotipo portador (para a rede e o contraste intra-animal).
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
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
from hap_ancestry import (  # noqa: E402
    LocusConfig, load_phased, load_consequences, load_groups, hap_group_labels,
    active_panels, resolve_classes, focal_index, carrier_haplotypes,
    panel_frequencies, assign_with_posterior, winners,
)
from hap_ancestry.bootstrap import stability  # noqa: E402

VCF = str(DATA / "criollo/incorporado/merged_phased_dedup.vcf.gz")
CSQ = str(DATA / "haplo/prlr_allmaf_csq.tsv")
CFG = RAIZ / "config" / "slick_iberico_declarado.yaml"
GRP = RAIZ / "tables" / "groups_iberico_declarado.tsv"
SUP = RAIZ / "tables" / "derived" / "suporte_leitura_portadores.tsv"
SAIDA = RAIZ / "tables" / "derived" / "resultado_final.json"
CLASSES = ("so sinonimas", "so nao sinonimas", "so introns", "gene inteiro")
COD = ("so sinonimas", "so nao sinonimas")
BORDAS = [0, 1_000, 5_000, 10_000, 20_000, 10**9]
ROT_FAIXA = ["<1 kb", "1-5 kb", "5-10 kb", "10-20 kb", ">20 kb"]
REPS_BOOT, REPS_RAR, SEMENTE = 2000, 500, 20261026


def excluidos():
    """Chamadas de portador SEM suporte de leitura, por alelo."""
    ex = collections.defaultdict(set)
    for l in io.open(SUP, encoding="utf-8").read().splitlines()[1:]:
        c = l.split("\t")
        if len(c) > 8 and c[8].strip() == "nao":
            ex[c[0]].add(c[1])
    return ex


def main():
    rng = np.random.default_rng(SEMENTE)
    cfg = LocusConfig.from_yaml(str(CFG))
    pos, H, samples = load_phased(VCF)
    hg = hap_group_labels(samples, load_groups(str(GRP)))
    panels = active_panels(hg, cfg.panels, cfg.min_hap_panel)
    cls = dict(resolve_classes(pos, load_consequences(CSQ, cfg.consequence_classes),
                               cfg.class_report))
    intron_pos = set(load_consequences(CSQ, cfg.consequence_classes).get("intron", set()))
    idx_intron = np.array([i for i, p in enumerate(pos) if int(p) in intron_pos])
    EX = excluidos()
    tam = {p: int((hg == p).sum()) for p in panels}
    print(f"[paineis] {tam}\n[excluidos sem leitura] { {k: sorted(v) for k, v in EX.items()} }")

    rel = {"config": str(CFG.name), "rotulos": str(GRP.name), "paineis_haplotipos": tam,
           "excluidos_sem_suporte": {k: sorted(v) for k, v in EX.items()},
           "reps_bootstrap": REPS_BOOT, "reps_rarefacao": REPS_RAR, "semente": SEMENTE,
           "alelos": {}}

    for rotulo, foco in cfg.focal_sites:
        i_foc, exato = focal_index(pos, foco)
        if not exato:
            continue
        todos = np.asarray(carrier_haplotypes(H, i_foc))
        if todos.size == 0:
            continue
        cmask = np.zeros(H.shape[1], bool)
        cmask[todos] = True                       # conservador: chamados todos fora da referencia
        port = np.array([h for h in todos if samples[h // 2] not in EX.get(rotulo, set())])
        F = panel_frequencies(H, hg, panels, cmask, cfg.eps)
        animais = sorted({samples[h // 2] for h in port})
        r = {"pos": int(foco), "n_haplotipos": int(port.size), "n_animais": len(animais),
             "animais": animais, "classes": {}}
        print(f"\n=== {rotulo}: {port.size} haplotipos em {len(animais)} animais ===")

        for c in CLASSES:
            sel = cls[c]
            cont, post = assign_with_posterior(H, sel, port, F, panels, cfg.min_sites_class)
            st = stability(H, sel, port, F, panels, cfg.min_sites_class, REPS_BOOT, rng) \
                if c in COD else None
            r["classes"][c] = {
                "n_sitios": int(len(sel)), "contagem": dict(cont),
                "posterior_mediana": float(np.median(post)) if post else None,
                "estabilidade_bootstrap": float(st["estabilidade"]) if st else None,
                # IC95 percentil da FRACAO de haplotipos no painel vencedor, sobre as replicas:
                # a incerteza da atribuicao, e nao so o painel vencedor.
                "vencedor_bootstrap": st["vencedor"] if st else None,
                "fracao_ic95": [float(st["ic_lo"]), float(st["ic_hi"])] if st else None,
            }
            stab = r["classes"][c]["estabilidade_bootstrap"]
            print(f"   {c:<18} {len(sel):>6} sitios | "
                  + " ".join(f"{k}:{v}" for k, v in sorted(cont.items(), key=lambda kv: -kv[1]))
                  + f"  (post. {np.median(post) if post else float('nan'):.2f}"
                  + (f"; estab. {100*stab:.1f}%" if stab is not None else "") + ")")

        # rarefacao: todos os paineis ao tamanho efetivo do menor
        efetivo = {p: int(((hg == p) & ~cmask).sum()) for p in panels}
        k = min(efetivo.values())
        venc = {c: collections.Counter() for c in COD}
        for _ in range(REPS_RAR):
            g2 = np.array(hg, dtype=object)
            for p in panels:
                idx = np.where((hg == p) & ~cmask)[0]
                g2[rng.choice(idx, size=idx.size - k, replace=False)] = "RAREFEITO"
            F2 = panel_frequencies(H, g2, panels, cmask, cfg.eps)
            for c in COD:
                cont, _ = assign_with_posterior(H, cls[c], port, F2, panels, cfg.min_sites_class)
                if cont:
                    venc[c][max(cont.items(), key=lambda kv: kv[1])[0]] += 1
        r["rarefacao"] = {"k": k, "modal_por_replica": {c: dict(venc[c]) for c in COD}}

        # lote
        lote = np.array(["RUN8" if samples[h // 2].startswith("SAMN") else "ERP150979" for h in port])
        r["lote"] = {}
        for L in sorted(set(lote)):
            r["lote"][L] = {c: dict(assign_with_posterior(H, cls[c], port[lote == L], F, panels,
                                                          cfg.min_sites_class)[0]) for c in COD}

        # configuracoes codificantes
        cod_idx = np.concatenate([cls[c] for c in COD])
        confs = collections.Counter(tuple(H[cod_idx, h].tolist()) for h in port)
        r["configuracoes"] = {"n": len(confs), "copias": sorted(confs.values(), reverse=True)}

        # introns por distancia ao sitio causal
        dist = np.abs(pos[idx_intron].astype(np.int64) - int(foco))
        r["introns_por_distancia"] = []
        for j, rot in enumerate(ROT_FAIXA):
            sel = idx_intron[(dist >= BORDAS[j]) & (dist < BORDAS[j + 1])]
            if len(sel) < cfg.min_sites_class:
                r["introns_por_distancia"].append({"faixa": rot, "n_sitios": int(len(sel))})
                continue
            cont, post = assign_with_posterior(H, sel, port, F, panels, cfg.min_sites_class)
            r["introns_por_distancia"].append({"faixa": rot, "n_sitios": int(len(sel)),
                                               "contagem": dict(cont),
                                               "posterior_mediana": float(np.median(post))})

        # cada haplotipo portador, nas classes codificantes
        r["por_haplotipo"] = []
        for h in port:
            d = {"amostra": samples[h // 2], "cromossomo": int(h % 2)}
            for c in COD:
                w = winners(H, cls[c], [h], F, panels, cfg.min_sites_class)
                d[c] = w[0] if w else None
            r["por_haplotipo"].append(d)
        rel["alelos"][rotulo] = r

    io.open(SAIDA, "w", encoding="utf-8").write(json.dumps(rel, ensure_ascii=False, indent=1))
    print(f"\n-> {SAIDA}")


if __name__ == "__main__":
    main()
