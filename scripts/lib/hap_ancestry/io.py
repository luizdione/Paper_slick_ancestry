"""io.py — leitura de VCF faseado, consequencias e grupos, e resolucao de classes.

Generico: nada aqui conhece *slick* ou *PRLR*. O que muda entre projetos (mapa de
consequencias, classes reportadas) entra por parametro, vindo do `LocusConfig`.
"""
from __future__ import annotations
import gzip
import numpy as np
import pandas as pd


def load_phased(path):
    """VCF faseado -> (posicoes, matriz haplotipo x sitio, amostras).

    Uma coluna por HAPLOTIPO (duas por amostra, na ordem das colunas do VCF).
    Alelo ausente vira -1. Le tanto `|` quanto `/` como separador de fase.
    """
    pos, rows, samples = [], [], []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            f = line.rstrip("\n").split("\t")
            if line.startswith("#CHROM"):
                samples = f[9:]
                continue
            r = np.empty(2 * len(samples), dtype=np.int8)
            for i, g in enumerate(f[9:]):
                a = g.split(":")[0]
                sep = "|" if "|" in a else "/"
                x = a.split(sep)
                r[2 * i] = int(x[0]) if x[0].isdigit() else -1
                r[2 * i + 1] = int(x[1]) if len(x) > 1 and x[1].isdigit() else -1
            pos.append(int(f[1]))
            rows.append(r)
    return np.array(pos), np.vstack(rows), samples


def load_consequences(path, consequence_classes):
    """POS -> classe funcional, a partir do TSV do bcftools csq.

    `consequence_classes` mapeia o token de consequencia (em minusculas) para o
    nome da classe: p.ex. {'synonymous':'syn', 'missense':'nonsyn', 'intron':
    'intron'}. Tokens fora do mapa sao IGNORADOS de proposito — UTR, splice e nao
    anotados nao sao codificantes nem intronicos, e agrupa-los com uma classe de
    contraste a inflaria com sitios de outra natureza.

    O arquivo tem QUATRO colunas (POS, REF, ALT, consequencia) e a consequencia e
    a ultima, nao a oitava — um parser que procure em `f[7]` nao acha nada, nao
    levanta erro, e devolve tudo na classe default. Cada linha traz a consequencia
    repetida por transcrito, separada por virgula; basta o primeiro campo (os
    transcritos concordam na classe, divergindo so na numeracao do residuo).

    Retorna dict {nome_da_classe: set(POS)}.
    """
    classes: dict[str, set] = {}
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if len(f) < 4:
            continue
        try:
            p = int(f[0])
        except ValueError:
            continue
        c = f[3].split(",")[0].split("|")[0].lower()
        if c.startswith("@") or c == ".":   # continuacao do csq / sem anotacao
            continue
        nome = consequence_classes.get(c)
        if nome is None:
            continue
        classes.setdefault(nome, set()).add(p)
    return classes


def load_groups(path):
    """TSV `amostra<TAB>grupo` (sem cabecalho) -> dict {amostra: grupo}."""
    grp = pd.read_csv(path, sep="\t", header=None, names=["sample", "grupo"])
    return dict(zip(grp["sample"], grp["grupo"]))


def hap_group_labels(samples, group_of, missing="NA"):
    """Rotulo de grupo por HAPLOTIPO (dois por amostra, na ordem das colunas)."""
    return np.array([group_of.get(s, missing) for s in samples for _ in (0, 1)])


def active_panels(hap_grp, panels, min_hap_panel):
    """Filtra os paineis candidatos aos que tem >= min_hap_panel haplotipos,
    PRESERVANDO a ordem passada — o desempate do argmax na atribuicao depende dela.
    """
    return [p for p in panels if (hap_grp == p).sum() >= min_hap_panel]


def resolve_classes(pos, class_sets, class_report):
    """Traduz o `class_report` da config em selecoes de indices no array `pos`.

    Cada item e (label, of): `of == "*"` significa o gene inteiro (todos os
    sitios, inclusive UTR/splice); `of == [nomes]` significa a UNIAO dos POS das
    classes nomeadas. Retorna lista de (label, idx_array), na ordem dada.
    """
    out = []
    for label, of in class_report:
        if of == "*" or of == ["*"]:
            idx = np.arange(len(pos))
        else:
            wanted: set = set()
            for nome in of:
                wanted |= class_sets.get(nome, set())
            idx = np.where(np.isin(pos, list(wanted)))[0]
        out.append((label, idx))
    return out
