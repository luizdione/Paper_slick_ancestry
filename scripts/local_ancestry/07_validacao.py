# -*- coding: utf-8 -*-
"""07_validacao.py -- consolida a VALIDACAO do RFMix num JSON so, fonte unica das figuras.

Tres camadas de validacao, nesta ordem (regra 1 do projeto: controle antes de conclusao):

  1. SIMULACAO INTERNA do proprio RFMix: ele cria haplotipos misturados a partir das referencias
     e mede quanto acerta. E a matriz impressa em `*.log` sob "Simulation results": COLUNA =
     ancestralidade verdadeira, LINHA = chamada. Le-se pela diagonal.
  2. CONTROLE POR OMISSAO: animais de referencia deixados de fora e submetidos como consulta.
  3. CONTROLE POSITIVO DE MISTURA: Ankole, Boran e Ogaden sao taurino africano x zebu conhecidos;
     tem de sair misturados, nao puros.

E a camada que decide a leitura: so a configuracao de 3 vias passa na camada 1 para TODOS os
paineis. Saida: `$SLICK_DATA_DIR/lai10mb/rfmix/validacao.json`.
"""
import os
import json
import pathlib
import re

RAIZ = pathlib.Path(__file__).resolve().parents[2]   # raiz do repositorio
DATA = pathlib.Path(os.environ.get("SLICK_DATA_DIR", RAIZ / "data"))   # dados brutos (VCF, BAM), fora do repositorio

R = (DATA / "lai10mb/rfmix")
RUNS = {"m4": "4 vias", "m3": "3 vias", "m5": "5 vias", "c4": "4 vias, paineis pequenos reduzidos"}


def matriz(run):
    txt = (R / f"{run}.log").read_text()
    bloco = txt.split("Simulation results")[1].splitlines()
    pops = [x for x in bloco[1].split("\t") if x.strip()]
    linhas = []
    for l in bloco[3:]:
        c = l.split("\t")
        if not c[0].strip().isdigit():
            break
        linhas.append([float(x) for x in c[1:] if x.strip()])
    # m[i][j]: fracao chamada como pops[i] quando a verdade e pops[j]
    sim = {"populacoes": pops,
           "acerto": {p: linhas[i][i] for i, p in enumerate(pops)},
           "matriz": linhas}
    for i, p in enumerate(pops):          # para onde vai o erro de cada populacao
        col = [(pops[k], linhas[k][i]) for k in range(len(pops)) if k != i]
        sim.setdefault("principal_confusao", {})[p] = max(col, key=lambda x: x[1])
    m = re.search(r"Maximum scoring weight is (\d+)", txt)
    sim["peso_crf"] = int(m.group(1)) if m else None
    m = re.search(r"Internally simulated (\d+) samples from (\d+)", txt)
    if m:
        sim["n_simulados"], sim["n_parentes"] = int(m.group(1)), int(m.group(2))
    n = re.search(r"setting up CRF points", txt)
    return sim


rel = {"runs": {}}
for run, desc in RUNS.items():
    sim = matriz(run)
    rel["runs"][run] = {"descricao": desc, "simulacao_interna": sim}
    piso = min(sim["acerto"].values())
    pior = min(sim["acerto"], key=sim["acerto"].get)
    rel["runs"][run]["passa_todos_os_paineis"] = piso >= 70.0
    rel["runs"][run]["pior_painel"] = {"painel": pior, "acerto": piso,
                                       "confundido_com": sim["principal_confusao"][pior][0]}
    print(f"{run} ({desc}): " + " ".join(f"{p} {a:.0f}%" for p, a in sim["acerto"].items())
          + f" | peso CRF {sim['peso_crf']} | pior {pior} {piso:.0f}% "
          f"-> {sim['principal_confusao'][pior][0]} | "
          + ("PASSA" if piso >= 70 else "REPROVA"))

# camadas 2 e 3 vem do resumo ja calculado
res = json.loads((R / "resumo_rfmix.json").read_text(encoding="utf-8"))
for run in ("m4", "m3", "c4"):
    r = res.get(run, {})
    rel["runs"][run]["controle_por_omissao"] = {
        k: {"n_hap": v["n_hap"], "certo_no_prlr": v["certo_no_prlr"],
            "fracao_media_do_proprio_painel": v["fracao_media_do_proprio_painel"]}
        for k, v in r.get("controles", {}).items()}
    if "africanos_misturados" in r:
        rel["runs"][run]["controle_positivo_mistura"] = r["africanos_misturados"]
    rel["runs"][run]["janela_prlr"] = r.get("janela_prlr")
    rel["runs"][run]["n_janelas"] = r.get("n_janelas")

(R / "validacao.json").write_text(json.dumps(rel, ensure_ascii=False, indent=1), encoding="utf-8")
print(f"-> {R / 'validacao.json'}")
