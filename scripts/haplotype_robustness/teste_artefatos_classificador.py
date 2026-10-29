# -*- coding: utf-8 -*-
"""A afinidade europeia de SLICK2, SLICK3 e SLICK6 e ancestralidade ou artefato?

O contraste africano x europeu existe SO no nucleo codificante, e nenhuma das duas
justificativas para privilegiar essa classe se sustenta: a informatividade por sitio foi
refutada (teste pareado por MAF, `teste_informatividade_classe.py`) e a proximidade fisica
tambem (os codificantes estao, em mediana, a ~1,1 kb do sitio causal,
`teste_distancia_intron.py`).

Restam tres explicacoes de ARTEFATO para o "europeu", cada uma testavel com o dado em disco:

TESTE 1 - TAMANHO DE PAINEL.
  O painel europeu e o maior (162 haplotipos; N'Dama 34; Ankole 16). Com piso de frequencia
  `eps`, um painel pequeno e penalizado sempre que o haplotipo traz alelo ausente nele, entao o
  maior tende a vencer por construcao. Teste: rarefazer TODOS os paineis ao mesmo tamanho (o do
  menor) N vezes e ver se o vencedor se mantem.

TESTE 2 - LOTE / PIPELINE.
  Os portadores de SLICK1 sao todos do Run 8 (SAMN); os de SLICK3 e SLICK6 sao todos do lote novo
  ERP150979 (ERR), com captura de janela, bwa e 4,25% de genotipos ausentes imputados contra um
  painel dominado por europeus do Run 8. O SLICK2 mistura os dois lotes, e e o unico teste
  interno possivel: atribuir separadamente os portadores SAMN e os ERR.

TESTE 3 - INDEPENDENCIA.
  Quatro animais carregam dois alelos. Se os tres "europeus" vivem nos mesmos animais, nao sao
  tres observacoes. O teste lista as sobreposicoes e as configuracoes codificantes distintas.

IMPORTANTE: os rotulos vem de `groups_afr_split.tsv` lido pelo motor como sempre. NAO
"corrigir" o carriage return: o conteudo derrubado e contaminado por Brahman e Hanwoo (ver README).
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

from hap_ancestry import (                                         # noqa: E402
    LocusConfig, load_phased, load_consequences, load_groups,
    hap_group_labels, active_panels, resolve_classes,
    focal_index, carrier_haplotypes, panel_frequencies, assign_with_posterior,
)

VCF = str(DATA / "criollo/incorporado/merged_phased_dedup.vcf.gz")
CSQ = str(DATA / "haplo/prlr_allmaf_csq.tsv")
# (2026-10-26) Rotulos, config e saida por argumento, para rodar a mesma bateria com o painel
# IBERIAN. Sem argumentos, reproduz exatamente a execucao original.
#   uso: teste_artefatos_classificador.py [rotulos.tsv] [config.yaml] [saida.json]
_A = sys.argv[1:]
GRUPOS = pathlib.Path(_A[0]) if len(_A) > 0 else RAIZ / "tables" / "groups_afr_split.tsv"
CONFIG = pathlib.Path(_A[1]) if len(_A) > 1 else RAIZ / "config" / "slick_com_slick6.yaml"
SAIDA = pathlib.Path(_A[2]) if len(_A) > 2 else RAIZ / "tables" / "derived" / "teste_artefatos_classificador.json"
# (2026-10-26) 4o argumento opcional: JSON {alelo: [amostras]} com chamadas de portador SEM
# suporte de leitura (ver tables/derived/suporte_leitura_portadores.tsv). Essas amostras saem da
# lista de haplotipos ATRIBUIDOS, mas continuam mascaradas nas frequencias dos paineis, o que e
# o lado conservador: um portador possivelmente real nao contamina a referencia.
EXCLUIR = json.loads(_A[3]) if len(_A) > 3 else {}
N_REP = 500
SEMENTE = 20261026
CLASSES_COD = ("so sinonimas", "so nao sinonimas")


def resumo(cont):
    return " ".join(f"{k}:{v}" for k, v in sorted(cont.items(), key=lambda kv: -kv[1]))


def main() -> None:
    rng = np.random.default_rng(SEMENTE)
    cfg = LocusConfig.from_yaml(str(CONFIG))
    pos, H, samples = load_phased(VCF)
    hap_grp = hap_group_labels(samples, load_groups(str(GRUPOS)))
    panels = active_panels(hap_grp, cfg.panels, cfg.min_hap_panel)
    class_sets = load_consequences(CSQ, cfg.consequence_classes)
    classes = dict(resolve_classes(pos, class_sets, cfg.class_report))
    tam = {p: int((hap_grp == p).sum()) for p in panels}
    print(f"[paineis] {tam}")

    rel = {"semente": SEMENTE, "n_rep": N_REP, "paineis": tam, "alelos": {}}
    portadores = {}

    for rotulo, foco in cfg.focal_sites:
        i_foc, exato = focal_index(pos, foco)
        if not exato:
            continue
        port = np.asarray(carrier_haplotypes(H, i_foc))
        if port.size == 0:
            continue
        cmask = np.zeros(H.shape[1], bool)
        cmask[port] = True            # mascara com TODOS os chamados (lado conservador)
        fora = set(EXCLUIR.get(rotulo, []))
        if fora:
            antes = port.size
            port = np.array([h for h in port if samples[h // 2] not in fora])
            print(f"\n[{rotulo}] {antes - port.size} haplotipo(s) sem suporte de leitura removidos "
                  f"da atribuicao: {sorted(fora)}")
            if port.size == 0:
                continue
        portadores[rotulo] = port
        print(f"\n=== {rotulo}: {port.size} haplotipos ===")
        r = {"n_hap": int(port.size)}

        # ---------- linha de base (paineis completos, sem rarefacao)
        F = panel_frequencies(H, hap_grp, panels, cmask, cfg.eps)
        base = {}
        for c in CLASSES_COD:
            cont, _ = assign_with_posterior(H, classes[c], port, F, panels, cfg.min_sites_class)
            base[c] = dict(cont)
            print(f"   base        {c:<18} {resumo(cont)}")
        r["base"] = base

        # ---------- TESTE 1: rarefacao de todos os paineis ao tamanho do menor
        # tamanho efetivo por painel, ja excluidos os portadores deste alelo
        efetivo = {p: int(((hap_grp == p) & ~cmask).sum()) for p in panels}
        k = min(efetivo.values())
        venc = {c: collections.Counter() for c in CLASSES_COD}
        for _ in range(N_REP):
            g2 = np.array(hap_grp, dtype=object)
            for p in panels:
                idx = np.where((hap_grp == p) & ~cmask)[0]
                fora = rng.choice(idx, size=idx.size - k, replace=False)
                g2[fora] = "RAREFEITO"
            F2 = panel_frequencies(H, g2, panels, cmask, cfg.eps)
            for c in CLASSES_COD:
                cont, _ = assign_with_posterior(H, classes[c], port, F2, panels,
                                                cfg.min_sites_class)
                if cont:
                    venc[c][max(cont.items(), key=lambda kv: kv[1])[0]] += 1
        r["rarefacao"] = {"k_por_painel": k, "efetivo": efetivo,
                          "painel_modal_por_replica": {c: dict(venc[c]) for c in CLASSES_COD}}
        for c in CLASSES_COD:
            tot = sum(venc[c].values()) or 1
            txt = " ".join(f"{p}:{100*n/tot:.0f}%" for p, n in venc[c].most_common())
            print(f"   rarefeito k={k:<3} {c:<18} {txt}")

        # ---------- TESTE 2: separar por lote (so faz sentido onde ha os dois)
        lote = np.array(["RUN8" if samples[h // 2].startswith("SAMN") else "NOVO" for h in port])
        r["lote"] = dict(collections.Counter(lote.tolist()))
        if len(set(lote)) > 1:
            r["por_lote"] = {}
            for L in ("RUN8", "NOVO"):
                sub = port[lote == L]
                r["por_lote"][L] = {}
                for c in CLASSES_COD:
                    cont, _ = assign_with_posterior(H, classes[c], sub, F, panels,
                                                    cfg.min_sites_class)
                    r["por_lote"][L][c] = dict(cont)
                    print(f"   lote {L:<6} n={sub.size:<3} {c:<18} {resumo(cont)}")
        else:
            print(f"   lote: todos {lote[0]} ({lote.size} haplotipos) - sem contraste interno")

        # ---------- configuracoes codificantes distintas (o n efetivo)
        cod = np.concatenate([classes[c] for c in CLASSES_COD])
        confs = collections.Counter(tuple(H[cod, h].tolist()) for h in port)
        r["configuracoes_codificantes"] = len(confs)
        r["copias_por_configuracao"] = sorted(confs.values(), reverse=True)
        print(f"   configuracoes codificantes distintas: {len(confs)} "
              f"(copias {sorted(confs.values(), reverse=True)})")
        rel["alelos"][rotulo] = r

    # ---------- TESTE 3: sobreposicao de portadores entre alelos
    animais = {a: {samples[h // 2] for h in p} for a, p in portadores.items()}
    sobre = []
    nomes = list(animais)
    for i, a in enumerate(nomes):
        for b in nomes[i + 1:]:
            comum = sorted(animais[a] & animais[b])
            if comum:
                sobre.append({"par": f"{a}+{b}", "animais": comum})
                print(f"\n[sobreposicao] {a}+{b}: {', '.join(comum)}")
    rel["sobreposicao"] = sobre

    io.open(SAIDA, "w", encoding="utf-8").write(json.dumps(rel, ensure_ascii=False, indent=1))
    print(f"\n-> {SAIDA}")


if __name__ == "__main__":
    main()
