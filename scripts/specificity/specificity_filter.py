"""specificity_filter.py — Filtro Estatistico de Restricao de Grupo Fenotipico.

Generaliza, como funcao reutilizavel, um teste de especificidade populacional:
uma mutacao LoF (loss-of-function) so e candidata a alelo causal de um
fenotipo restrito (aqui, o pelame *slick* do gado criollo tropical) se ela
propria estiver restrita — ou ao menos concentrada — nas populacoes onde o
fenotipo existe. Uma variante presente em multiplos paineis TAXONOMICAMENTE
DISTANTES do grupo fenotipico (taurino temperado, sanga, zebu) nao pode ser a
variante responsavel, qualquer que seja a frequencia publicada numa unica
amostra isolada.

Duas leituras da mesma pergunta, aplicadas em conjunto:

1. TESTE BINOMIAL POR PAINEL DE CONTROLE. Sob H0, a variante e rara/ruido
   no painel de controle k (frequencia <= f0). Observar x_k copias em n_k
   alelos chamados so e compativel com H0 se P(X >= x_k | n_k, f0) nao for
   pequena. Quando essa probabilidade cai abaixo de alpha, H0 e rejeitada:
   a variante e "comum demais" fora do grupo fenotipico para ser rara ali,
   e o painel de controle "vota" pela exclusao do candidato.
   E o mesmo calculo binomial usado para estimar o limite de deteccao de um
   painel de referencia ((1-p)^n >= 0,05), so que perguntado ao contrario —
   aqui ha portadores observados, e a pergunta e se a frequencia deles e alta
   demais para ser ruido, nao se a ausencia de portadores exclui uma
   frequencia.

2. REGRA DE DECISAO AGREGADA. Conta em quantos paineis de controle DISTINTOS
   a variante aparece com pelo menos um portador. Um candidato exploratorio
   real testado com esta rotina (coordenada c.1397A>T) apareceu em
   EUR_TAURINE, BORAN, NDAMA e OGADEN — 4 paineis taxonomicamente distantes —
   e foi descartado sem meio-termo. O limiar `min_control_panels_hit`
   codifica esse criterio como numero, nao como julgamento caso a caso.

Uso tipico (linha de comando):

    python specificity_filter.py \
        --vcf candidatos_lof.vcf.gz \
        --groups groups.tsv \
        --phenotype-groups CRIOLLO \
        --control-groups EUR_TAURINE NDAMA ANKOLE BORAN OGADEN ZEBU \
        --f0 0.01 --alpha 0.05 --min-control-panels-hit 2 \
        --out specificity_report.tsv

Uso como biblioteca (dentro de um notebook ou de outro script do projeto):

    from specificity_filter import run_specificity_filter
    report = run_specificity_filter(vcf_path, groups_path,
                                     phenotype_groups=["CRIOLLO"],
                                     control_groups=["EUR_TAURINE","NDAMA",
                                                      "ANKOLE","BORAN","OGADEN"])

Depende so de numpy/pandas/scipy e do modulo `io` do motor generico de
ancestralidade de haplotipo (reaproveitado, nao duplicado) para ler o VCF
faseado e o TSV de grupos no mesmo formato que o resto do projeto ja usa.
"""
from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import binomtest

# Reusa o parser de VCF faseado e de grupos do motor existente em vez de
# duplicar. O motor fica em scripts/lib/hap_ancestry; adiciona ao path se
# este script for chamado de fora da arvore do projeto.
_LIB = Path(__file__).resolve().parents[1] / "lib"
if str(_LIB) not in sys.path:
    sys.path.insert(0, str(_LIB))
from hap_ancestry import io as hio  # noqa: E402


@dataclass
class SpecificityVerdict:
    pos: int
    ref: str
    alt: str
    n_phenotype_carriers: int
    n_phenotype_haplotypes: int
    freq_phenotype: float
    control_hits: dict = field(default_factory=dict)   # painel -> (x, n, p_value)
    n_control_panels_hit: int = 0
    decision: str = ""     # "mantido" | "descartado"
    reason: str = ""


def binomial_specificity_pvalue(x: int, n: int, f0: float) -> float | None:
    """P(X >= x | n, f0) sob H0: frequencia no painel de controle <= f0.

    p pequeno => rejeita H0 => a frequencia observada no controle e alta
    demais para ser explicada por ruido/erro de genotipagem numa variante
    rara — o painel de controle "vota" pela exclusao do candidato.
    Retorna None se n == 0 (painel sem chamadas, sem evidencia nenhuma).
    """
    if n == 0:
        return None
    return binomtest(x, n, f0, alternative="greater").pvalue


def load_lof_vcf_sites(vcf_path: str) -> tuple[np.ndarray, np.ndarray, np.ndarray, list]:
    """Le um VCF (faseado, bgzip) de variantes LoF candidatas.

    Reaproveita hap_ancestry.io.load_phased, que ja devolve posicoes,
    matriz haplotipo x sitio (uma coluna por haplotipo, -1 = ausente) e a
    lista de amostras na ordem das colunas do VCF. Alelo REF/ALT por sitio
    e lido a parte porque load_phased descarta essas colunas (o motor
    original nao precisava delas).
    """
    pos, hap_matrix, samples = hio.load_phased(vcf_path)
    ref, alt = [], []
    import gzip
    with gzip.open(vcf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            ref.append(f[3])
            alt.append(f[4])
    return pos, hap_matrix, samples, list(zip(ref, alt))


def run_specificity_filter(
    vcf_path: str,
    groups_path: str,
    phenotype_groups: list[str],
    control_groups: list[str],
    f0: float = 0.01,
    alpha: float = 0.05,
    min_control_panels_hit: int = 2,
    min_control_pvalue_panels: int = 1,
) -> pd.DataFrame:
    """Aplica o filtro de especificidade a toda variante LoF do VCF.

    Parametros
    ----------
    vcf_path : VCF faseado (bgzip) com as variantes LoF candidatas —
        tipicamente a saida de `bcftools csq` filtrada para stop_gained/
        stop_lost/start_lost/frameshift, na mesma janela do gene, faseada
        sem filtro de frequencia.
    groups_path : TSV `amostra<TAB>grupo` sem cabecalho (mesmo formato de
        hap_ancestry.io.load_groups) — o rotulo populacional de cada amostra
        (CRIOLLO, EUR_TAURINE, NDAMA, ANKOLE, BORAN, OGADEN, ZEBU, ...).
    phenotype_groups : grupos onde o fenotipo *slick* existe (tipicamente
        so CRIOLLO — os demais nao tem o pelame liso documentado).
    control_groups : grupos "nao-tropicais/nao-criollos" contra os quais a
        especificidade e testada (taurino temperado, sanga, zebu — NUNCA
        o proprio grupo fenotipico).
    f0 : frequencia nula por painel de controle — o teto de frequencia
        compativel com "rara/ruido" nesse painel. Default 1%, o mesmo
        patamar de MAF minimo tipicamente usado no controle de qualidade
        de paineis de ancestralidade.
    alpha : nivel de significancia do teste binomial por painel.
    min_control_panels_hit : minimo de paineis de controle DISTINTOS com
        pelo menos 1 portador para descartar por contagem bruta (o criterio
        que fechou o candidato c.1397A>T: 4 paineis, muito acima deste
        limiar).
    min_control_pvalue_panels : minimo de paineis cujo teste binomial
        rejeita H0 (p < alpha) para descartar por evidencia estatistica,
        mesmo que so 1 painel de controle tenha sido amostrado.

    Retorna
    -------
    DataFrame, uma linha por variante, com a decisao e a razao — pronto
    para virar tabela suplementar ou para alimentar o motor hap_ancestry
    (so os candidatos "mantido" devem entrar em focal_sites no YAML).
    """
    pos, hap_matrix, samples, ref_alt = load_lof_vcf_sites(vcf_path)
    group_of = hio.load_groups(groups_path)
    hap_grp = hio.hap_group_labels(samples, group_of)

    phen_mask = np.isin(hap_grp, phenotype_groups)
    rows = []
    for i, p in enumerate(pos):
        site = hap_matrix[i]
        called = site >= 0
        alt_carrier = site == 1   # ALT indexado 1; sitios multialelicos nao
                                    # esperados aqui (LoF ja filtrado biallelico)

        n_phen = int((called & phen_mask).sum())
        x_phen = int((alt_carrier & phen_mask).sum())
        freq_phen = x_phen / n_phen if n_phen else float("nan")

        control_hits = {}
        n_panels_hit = 0
        n_panels_reject_h0 = 0
        for panel in control_groups:
            panel_mask = hap_grp == panel
            n_k = int((called & panel_mask).sum())
            x_k = int((alt_carrier & panel_mask).sum())
            pval = binomial_specificity_pvalue(x_k, n_k, f0)
            control_hits[panel] = (x_k, n_k, pval)
            if x_k > 0:
                n_panels_hit += 1
            if pval is not None and pval < alpha:
                n_panels_reject_h0 += 1

        descartado = (
            n_panels_hit >= min_control_panels_hit
            or n_panels_reject_h0 >= min_control_pvalue_panels
        )
        if descartado:
            razao_partes = []
            if n_panels_hit >= min_control_panels_hit:
                razao_partes.append(
                    f"presente em {n_panels_hit} paineis de controle "
                    f"(limiar {min_control_panels_hit})"
                )
            if n_panels_reject_h0 >= min_control_pvalue_panels:
                razao_partes.append(
                    f"{n_panels_reject_h0} painel(is) rejeitam H0 f<={f0:.1%} "
                    f"a alpha={alpha}"
                )
            razao = "; ".join(razao_partes)
            decisao = "descartado"
        else:
            razao = "restrito (ou concentrado) ao grupo fenotipico sob os limiares dados"
            decisao = "mantido"

        r, a = ref_alt[i]
        rows.append({
            "pos": int(p), "ref": r, "alt": a,
            "n_haplotipos_fenotipo": n_phen,
            "portadores_fenotipo": x_phen,
            "freq_fenotipo": freq_phen,
            **{f"{panel}_x": control_hits[panel][0] for panel in control_groups},
            **{f"{panel}_n": control_hits[panel][1] for panel in control_groups},
            **{f"{panel}_pval": control_hits[panel][2] for panel in control_groups},
            "n_paineis_controle_com_portador": n_panels_hit,
            "n_paineis_rejeitam_H0": n_panels_reject_h0,
            "decisao": decisao,
            "razao": razao,
        })

    df = pd.DataFrame(rows).sort_values("pos").reset_index(drop=True)
    return df


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--vcf", required=True, help="VCF LoF faseado, bgzip")
    ap.add_argument("--groups", required=True, help="TSV amostra<TAB>grupo, sem cabecalho")
    ap.add_argument("--phenotype-groups", nargs="+", required=True)
    ap.add_argument("--control-groups", nargs="+", required=True)
    ap.add_argument("--f0", type=float, default=0.01)
    ap.add_argument("--alpha", type=float, default=0.05)
    ap.add_argument("--min-control-panels-hit", type=int, default=2)
    ap.add_argument("--min-control-pvalue-panels", type=int, default=1)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    df = run_specificity_filter(
        args.vcf, args.groups, args.phenotype_groups, args.control_groups,
        f0=args.f0, alpha=args.alpha,
        min_control_panels_hit=args.min_control_panels_hit,
        min_control_pvalue_panels=args.min_control_pvalue_panels,
    )
    df.to_csv(args.out, sep="\t", index=False)
    n_mantidos = (df["decisao"] == "mantido").sum()
    n_descartados = (df["decisao"] == "descartado").sum()
    print(f"[specificity_filter] {len(df)} variante(s) avaliada(s): "
          f"{n_mantidos} mantida(s), {n_descartados} descartada(s). "
          f"Relatorio: {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
