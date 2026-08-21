#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""apply_run9_codebook.py -- resolve o grupo de ancestralidade de cada Run9 a partir de
TODAS as evidencias disponiveis (NCBI breed + ENA sample_title/study_title/alias + EBI breed),
aplicando o codebook curado com casamento por FRONTEIRA DE PALAVRA (regra #3: nada de substring
ingenuo tipo 'red' casando 'Red Fulani').

Prioridade de evidencia por amostra (a 1a que resolver para grupo != unresolved/ambiguous vence):
  1. NCBI breed (run9_breed_map.tsv col2), se != MISSING/empty
  2. ENA sample_title (ena_run9_evidence col sample_title), se nao-generico
  3. EBI breed characteristic (ena_run9_evidence col ebi_breed)
  4. ENA study_title (procura nomes de raca do codebook no titulo do estudo)
  5. alias prefix (KHO->Kholmogory, YAK->Yakutian) quando o estudo cita 2 racas
  6. organism (Bos indicus -> indicine coarse; Bos taurus -> unresolved taurino)

Saidas:
  run9_labels_v2.tsv       sample_id, resolved_breed, ancestry_group, ancestry_coarse,
                           confidence, evidence_source, evidence_string
  panel_labels_v2.tsv      painel completo (run9 v2 + linhas nao-run9 ja confiaveis, remapeadas)
  apply_run9_codebook.log  resumo (quantos resolvidos, por fonte, residuo)
REGRA #1: rotulo so vira REFERENCIA depois da sanidade de painel. Aqui so resolvemos o texto.
"""
import re, csv, sys, time
from pathlib import Path

BASE = Path(__file__).resolve().parent
LAB  = BASE / "sample_breed_labels.tsv"
EV   = BASE / "ena_run9_evidence.tsv"
MAP  = BASE / "run9_breed_map.tsv"
CB   = BASE / "run9_breed_codebook.tsv"
OUT  = BASE / "run9_labels_v2.tsv"
PANEL= BASE / "panel_labels_v2.tsv"
LOG  = BASE / "apply_run9_codebook.log"

COARSE = {
    "african_taurine": "african", "african_zebu": "african", "african_sanga": "african",
    "european_taurine": "european", "east_asian_taurine": "east_asian",
    "asian_indicine": "indicine", "iberian_criollo": "iberian_criollo",
    "iberian_taurine": "iberian",
    "ambiguous_geographic": "unresolved", "unresolved": "unresolved",
    "exclude": "exclude",
}
ALIAS_PREFIX = {"KHO": "Kholmogory", "YAK": "Yakutian"}

# mapeia os rotulos do painel PRIOR (nao-run9, ja confiaveis) para o esquema coarse novo
PRIOR_COARSE = {
    "EUR_TAURINE": "european", "HOLSTEIN": "european", "HOLSTEIN_introgr": "european",
    "IBERIAN": "iberian",
    "AFR_unresolved": "african", "AFR_SANGA": "african", "AFR_COMPOSITE": "african",
    "AFR_north": "african", "AFR_TAURINE_cand": "african",
    "SENEPOL": "senepol",
}
# amostras query (portadoras/alvo), NAO referencias, mesmo com rotulo coarse
PRIOR_QUERY_ONLY = {"SENEPOL", "HOLSTEIN_introgr", "AFR_TAURINE_cand"}

def log(m):
    with open(LOG, "a", encoding="utf-8") as f:
        f.write(f"[{time.strftime('%F %T')}] {m}\n")
    print(m, flush=True)

def load_codebook():
    rows = []
    for ln in open(CB, encoding="utf-8"):
        if ln.startswith("#") or ln.startswith("raw_pattern") or not ln.strip():
            continue
        p = ln.rstrip("\n").split("\t")
        if len(p) < 7:
            continue
        rows.append(dict(raw=p[0], mtype=p[1], breed=p[2], group=p[3],
                         subtype=p[4], conf=p[5], note=p[6]))
    return rows

def norm(s):
    """Normaliza para casamento: separa CamelCase (EthiopianBoran->Ethiopian Boran) e
    fronteiras letra<->digito, colapsa espaco, minusculas. Evita token colado escapar do \\b."""
    if not s:
        return ""
    s = re.sub(r"(?<=[a-z])(?=[A-Z])", " ", s)      # camelCase
    s = re.sub(r"(?<=[A-Za-z])(?=[0-9])|(?<=[0-9])(?=[A-Za-z])", " ", s)  # letra<->digito
    s = re.sub(r"[_\-/]+", " ", s)
    return re.sub(r"\s+", " ", s).strip().lower()

def match_string(s, cb):
    """Casa uma string livre contra o codebook. Retorna a MELHOR regra (grupo resolvido tem
    prioridade sobre ambiguous/unresolved), ou None."""
    if not s:
        return None
    low = norm(s)
    best = None
    for r in cb:
        pat = r["raw"].lower()
        hit = False
        if r["mtype"] == "exact":
            hit = (low.strip() == pat)
        elif r["mtype"] == "word":
            hit = re.search(r"(?<![a-z])" + re.escape(pat) + r"(?![a-z])", low) is not None
        elif r["mtype"] == "prefix":
            hit = low.startswith(pat)
        if not hit:
            continue
        resolved = r["group"] not in ("unresolved", "ambiguous_geographic")
        rank = (1 if resolved else 0, {"high": 3, "medium": 2, "low": 1}.get(r["conf"], 0))
        if best is None or rank > best[0]:
            best = (rank, r)
    return best[1] if best else None

def study_title_breed(title, cb):
    """Procura nomes de raca do codebook dentro do study_title (word match); ignora ambiguous."""
    if not title:
        return None
    low = norm(title)
    hits = []
    for r in cb:
        if r["group"] in ("unresolved", "ambiguous_geographic"):
            continue
        pat = r["raw"].lower()
        if r["mtype"] in ("word", "exact") and re.search(r"(?<![a-z])" + re.escape(pat) + r"(?![a-z])", low):
            hits.append(r)
        # tambem tenta o nome canonico
        bn = r["breed"].lower()
        if bn and re.search(r"(?<![a-z])" + re.escape(bn) + r"(?![a-z])", low) and r not in hits:
            hits.append(r)
    return hits or None

def main():
    for f in (OUT, PANEL, LOG):
        try: f.unlink()
        except FileNotFoundError: pass
    cb = load_codebook()

    # NCBI raw breed per input_acc
    ncbi = {}
    if MAP.exists():
        for ln in open(MAP, encoding="utf-8"):
            p = ln.rstrip("\n").split("\t")
            if len(p) >= 2:
                ncbi[p[0]] = p[1]

    # ENA evidence per input_acc
    ev = {}
    if EV.exists():
        rdr = csv.DictReader(open(EV, encoding="utf-8"), delimiter="\t")
        for row in rdr:
            ev[row["input_acc"]] = row

    # panel rows
    rows = [ln.rstrip("\n").split("\t") for ln in open(LAB, encoding="utf-8")]
    header = rows[0]
    run9 = [r for r in rows[1:] if len(r) >= 2 and r[1] == "run9_chr20"]
    other = [r for r in rows[1:] if len(r) >= 2 and r[1] != "run9_chr20"]

    def resolve(sample_id):
        e = ev.get(sample_id, {})
        title = (e.get("sample_title") or "").strip()
        study = (e.get("study_title") or "").strip()
        alias = (e.get("sample_alias") or "").strip()
        ebi   = (e.get("ebi_breed") or "").strip()
        organism = (e.get("scientific_name") or "").strip()
        nb = (ncbi.get(sample_id) or "").strip()
        gen_title = (not title) or title.lower().startswith("model organism") or title.lower() == "bos taurus" or title.lower() == "bos indicus"

        # 0. EXCLUSAO: amostra nao-individual (mistura de leite de tanque) ou especie fora do alvo
        blob = f"{title} {study}".lower()
        if "tank milk" in blob or re.match(r"^milksample", title.strip(), re.I) or "biological mixtures" in blob:
            return dict(breed="(tank-milk mixture)", group="exclude", subtype="pooled/non-individual", conf="high", note="drop"), "exclude_mixture", (title or study)
        if organism and organism.lower() not in ("bos taurus", "bos indicus"):
            return dict(breed=f"({organism})", group="exclude", subtype="non-target species", conf="high", note="drop"), "exclude_species", organism

        # 1. NCBI breed
        if nb and nb.upper() not in ("MISSING",):
            r = match_string(nb, cb)
            if r and r["group"] not in ("unresolved", "ambiguous_geographic"):
                return r, "ncbi_breed", nb
        # 2. ENA sample_title (nao-generico)
        if not gen_title:
            r = match_string(title, cb)
            if r and r["group"] not in ("unresolved", "ambiguous_geographic"):
                return r, "ena_title", title
        # 3. EBI breed
        if ebi:
            r = match_string(ebi, cb)
            if r and r["group"] not in ("unresolved", "ambiguous_geographic"):
                return r, "ebi_breed", ebi
        # 4. study_title
        hits = study_title_breed(study, cb)
        if hits:
            if len(hits) == 1:
                return hits[0], "study_title", study
            # varias racas no estudo -> tenta desambiguar por alias prefix
            pref = alias[:3].upper()
            want = ALIAS_PREFIX.get(pref)
            if want:
                for r in hits:
                    if r["breed"].lower() == want.lower():
                        return r, "study+alias", f"{study} | alias={alias}"
            # ambiguo entre varias: devolve o 1o como low
            r = dict(hits[0]); r["conf"] = "low"
            return r, "study_title_multi", study
        # 5. alias prefix isolado
        pref = alias[:3].upper()
        if pref in ALIAS_PREFIX:
            r = match_string(ALIAS_PREFIX[pref], cb)
            if r:
                rr = dict(r); rr["conf"] = "low"
                return rr, "alias_prefix", alias
        # 6. organism coarse
        if organism.lower() == "bos indicus" or (nb == "" and title.lower() == "model organism or animal sample from bos indicus"):
            return dict(breed="(Bos indicus)", group="asian_indicine", subtype="organism-level", conf="low", note="coarse"), "organism_indicus", title or organism
        # unresolved
        base = nb if nb and nb.upper() != "MISSING" else (title or "(none)")
        return dict(breed="(unresolved)", group="unresolved", subtype="", conf="low", note=""), "none", base

    src_count, grp_count, coarse_count, conf_count = {}, {}, {}, {}
    with open(OUT, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["sample_id", "resolved_breed", "ancestry_group", "ancestry_coarse",
                    "confidence", "evidence_source", "evidence_string"])
        for r in run9:
            sid = r[0]
            rule, source, estr = resolve(sid)
            group = rule["group"]
            coarse = COARSE.get(group, "unresolved")
            conf = rule.get("conf", "low")
            w.writerow([sid, rule["breed"], group, coarse, conf, source, estr[:120]])
            src_count[source] = src_count.get(source, 0) + 1
            grp_count[group] = grp_count.get(group, 0) + 1
            coarse_count[coarse] = coarse_count.get(coarse, 0) + 1
            conf_count[conf] = conf_count.get(conf, 0) + 1

    # painel completo v2: run9 v2 (coarse->grupo) + linhas nao-run9 remapeadas para o schema novo
    with open(OUT, encoding="utf-8") as f:
        run9v2 = {row["sample_id"]: row for row in csv.DictReader(f, delimiter="\t")}
    with open(PANEL, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["sample_id", "panel", "resolved_breed", "ancestry_group",
                    "ancestry_coarse", "confidence", "evidence_source"])
        for r in run9:
            v = run9v2[r[0]]
            w.writerow([r[0], "run9_chr20", v["resolved_breed"], v["ancestry_group"],
                        v["ancestry_coarse"], v["confidence"], v["evidence_source"]])
        for r in other:
            # r = sample, panel, breed, ancestry_group(orig), confidence
            breed = r[2] if len(r) > 2 else ""
            ag = r[3] if len(r) > 3 else ""
            conf = r[4] if len(r) > 4 else ""
            coarse = PRIOR_COARSE.get(ag, "unresolved")
            src = "prior_panel_query" if ag in PRIOR_QUERY_ONLY else "prior_panel"
            w.writerow([r[0], r[1], breed, ag, coarse, conf, src])

    log("== RESOLUCAO RUN9 v2 ==")
    log(f"total run9: {len(run9)}")
    log(f"por FONTE: {dict(sorted(src_count.items(), key=lambda x:-x[1]))}")
    log(f"por GRUPO fino: {dict(sorted(grp_count.items(), key=lambda x:-x[1]))}")
    log(f"por GRUPO coarse: {dict(sorted(coarse_count.items(), key=lambda x:-x[1]))}")
    log(f"por CONFIANCA: {dict(sorted(conf_count.items(), key=lambda x:-x[1]))}")
    resolved = len(run9) - grp_count.get("unresolved", 0)
    log(f"RESOLVIDOS (grupo != unresolved): {resolved}/{len(run9)}  ({100*resolved/len(run9):.1f}%)")
    log(f"RESIDUO unresolved: {grp_count.get('unresolved',0)}")
    log(f"-> {OUT.name}, {PANEL.name}")

if __name__ == "__main__":
    main()
