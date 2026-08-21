#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""resolver_run9_ena.py -- resolucao PROFUNDA da raca dos 1038 Run9 via ENA/EBI (nao so NCBI).

Motivo: o resolver_run9_biosample.py consultou apenas NCBI db=biosample lendo <Attribute>.
Isso PERDE a raca de:
  * SAMEA* (BioSamples do EBI): a raca vive no `sample_title` da ENA (ex. "Brown Swiss"),
    nao num Attribute breed harmonizado.
  * SRX*   (experimentos SRA): nao sao BioSample -> efetch biosample devolve nada. A raca
    esta no `study_title` (ex. "...Yakut and Kholmogory cattle breeds...") + prefixo do alias.
  * SAMN*  (NCBI reais) que vieram 'MISSING': muitos genericos (U. Missouri) -> so restam
    para inferencia genetica / BioProject.

Este script consulta a ENA portal API (filereport, result=read_run) para CADA id da coluna 1
(a chave de amostra do VCF, que pode ser SAMEA/SAMN/SRX), capturando o registro RICO:
  input_acc, run_accession, sample_accession, sample_title, sample_alias,
  study_accession, study_title, scientific_name, ebi_breed(from BioSamples json)

NAO decide o grupo aqui (isso e o codebook, etapa seguinte). So COLETA a evidencia bruta,
de forma resumivel (cache), com retry e rate-limit gentil. REGRA #1: rotulo so vira
referencia depois da sanidade de painel.

Saidas em /data/slick_ancestry/haplo/expansao/:
  ena_run9_evidence.tsv   (uma linha por input_acc; append/resume)
  resolver_run9_ena.log
"""
import sys, time, io, csv, json, urllib.request, urllib.parse
from pathlib import Path

import os
BASE = Path(os.environ.get("RUN9_BASE", Path(__file__).resolve().parent))
LAB  = BASE / "sample_breed_labels.tsv"
EV   = BASE / "ena_run9_evidence.tsv"
LOG  = BASE / "resolver_run9_ena.log"

ENA_FILEREPORT = "https://www.ebi.ac.uk/ena/portal/api/filereport"
EBI_BIOSAMPLE  = "https://www.ebi.ac.uk/biosamples/samples/{}.json"
FIELDS = ["run_accession", "sample_accession", "sample_title", "sample_alias",
          "study_accession", "study_title", "scientific_name", "broker_name", "center_name"]
SLEEP = 0.12
UA = {"User-Agent": "slick_ancestry/2.0 (luiz.melo@ifrj.edu.br)"}

def log(m):
    with open(LOG, "a", encoding="utf-8") as f:
        f.write(f"[{time.strftime('%F %T')}] {m}\n")
    print(m, flush=True)

def http_get(url, tries=4, timeout=60):
    for t in range(tries):
        try:
            req = urllib.request.Request(url, headers=UA)
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return r.read().decode("utf-8", "replace")
        except Exception as e:
            if t == tries - 1:
                return None
            time.sleep(1.5 + 2.0 * t)
    return None

def ena_filereport(acc):
    """Retorna dict do 1o run row, ou {} se nada."""
    q = urllib.parse.urlencode({"accession": acc, "result": "read_run",
                                "fields": ",".join(FIELDS), "format": "tsv"})
    txt = http_get(f"{ENA_FILEREPORT}?{q}")
    if not txt:
        return {}
    rows = list(csv.DictReader(io.StringIO(txt), delimiter="\t"))
    if not rows:
        return {}
    # prefere um run com sample_title nao-generico; senao o primeiro
    def score(r):
        t = (r.get("sample_title") or "").strip().lower()
        generic = (not t) or t.startswith("model organism") or t.startswith("bos ")
        return 0 if generic else 1
    rows.sort(key=score, reverse=True)
    return rows[0]

def ebi_breed(sample_acc):
    """Puxa characteristic 'breed'/'organism part' do BioSamples do EBI, se houver."""
    if not sample_acc:
        return ""
    txt = http_get(EBI_BIOSAMPLE.format(sample_acc), tries=2, timeout=40)
    if not txt:
        return ""
    try:
        ch = json.loads(txt).get("characteristics", {})
    except Exception:
        return ""
    for key in ("breed", "Breed", "organism part", "ecotype", "strain", "population"):
        if key in ch:
            vals = [x.get("text", "") for x in ch[key] if x.get("text")]
            if vals:
                return vals[0]
    return ""

def load_done():
    done = set()
    if EV.exists():
        for ln in open(EV, encoding="utf-8"):
            p = ln.split("\t", 1)
            if p and p[0] and p[0] != "input_acc":
                done.add(p[0])
    return done

def main():
    rows = [ln.rstrip("\n").split("\t") for ln in open(LAB, encoding="utf-8")]
    targets = [r[0] for r in rows if len(r) >= 2 and r[1] == "run9_chr20"]
    done = load_done()
    todo = [a for a in targets if a not in done]
    log(f"Run9 alvo: {len(targets)}; ja no cache: {len(done)}; a resolver agora: {len(todo)}")

    header = ["input_acc", "run_accession", "sample_accession", "sample_title",
              "sample_alias", "study_accession", "study_title", "scientific_name",
              "broker_name", "center_name", "ebi_breed"]
    new = not EV.exists()
    with open(EV, "a", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        if new:
            w.writerow(header)
        for i, acc in enumerate(todo, 1):
            rec = ena_filereport(acc)
            sample_acc = rec.get("sample_accession", "") or ""
            st = (rec.get("sample_title") or "").strip()
            eb = ""
            # so gasta a 2a chamada (EBI) quando o titulo esta vazio/generico
            gen = (not st) or st.lower().startswith("model organism") or st.lower().startswith("bos ")
            if gen:
                eb = ebi_breed(sample_acc or acc)
            w.writerow([acc, rec.get("run_accession", ""), sample_acc, st,
                        rec.get("sample_alias", ""), rec.get("study_accession", ""),
                        (rec.get("study_title") or "").strip(),
                        rec.get("scientific_name", ""), rec.get("broker_name", ""),
                        rec.get("center_name", ""), eb])
            if i % 50 == 0:
                f.flush()
                log(f"  {i}/{len(todo)} resolvidos (ultimo: {acc} -> title='{st[:40]}' ebi_breed='{eb[:30]}')")
            time.sleep(SLEEP)
    log(f"CONCLUIDO. Evidencia em {EV.name} ({len(targets)} alvos, {len(todo)} novos nesta corrida).")

if __name__ == "__main__":
    main()
