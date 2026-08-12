#!/usr/bin/env python3
"""
fig_class_ancestry.py — qual painel ancestral cada haplotipo portador recebe,
por classe funcional de sitio.

A figura existe para mostrar duas coisas ao mesmo tempo:

  1. As classes codificantes (sinonimas e nao sinonimas) tendem a concordar
     dentro de cada alelo. As sinonimas sao neutras por construcao, entao a
     concordancia com as nao sinonimas remove a objecao de que a atribuicao
     refletiria selecao.
  2. Os introns, muito mais numerosos, ficam nos flancos, onde a recombinacao
     teve um seculo para trocar o fundo fundador pelo local. O nucleo
     funcional guarda a procedencia melhor que o flanco.

O tamanho do marcador carrega o numero de sitios, que e a informacao que
separa uma atribuicao de 24 sitios de uma de 1.360.

Le o TSV real gravado pelo script de calculo (nao numero copiado a mao para
dentro deste arquivo).
"""
from __future__ import annotations
import argparse
import pathlib
import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

OUT = pathlib.Path(__file__).parent / "figuras"
OUT.mkdir(exist_ok=True)

# Ordem de exibicao das classes (rotulo do TSV -> rotulo em ingles no eixo).
# "gene inteiro" fica fora do grafico por padrao (ver --com-gene-inteiro);
# o objetivo central da figura e o contraste codificante x intron.
CLASSES_EIXO = [
    ("so sinonimas", "Synonymous"),
    ("so nao sinonimas", "Non-synonymous"),
    ("so introns", "Introns only"),
]
CLASSE_GENE = "gene inteiro"

# Cores fixas por painel (mesma paleta da config declarativa do locus).
# Paineis que aparecerem no TSV mas nao estiverem aqui recebem uma cor cinza
# de reserva, para a figura nunca quebrar por um painel novo.
COR = {
    "NDAMA": "#1b7837", "N'Dama": "#1b7837",
    "EUR_TAURINE": "#2166ac", "European taurine": "#2166ac",
    "CRIOLLO": "#b2abd2", "Criollo": "#b2abd2",
    "ANKOLE": "#8c6d31", "Ankole": "#8c6d31",
    "BORAN": "#d6604d", "Boran": "#d6604d",
    "OGADEN": "#4d4d4d", "Ogaden": "#4d4d4d",
}
COR_RESERVA = "#999999"
# Rotulo de exibicao para paineis que vem do TSV em maiusculo/codigo interno.
ROTULO_PAINEL = {
    "NDAMA": "N'Dama", "EUR_TAURINE": "European taurine", "CRIOLLO": "Criollo",
    "ANKOLE": "Ankole", "BORAN": "Boran", "OGADEN": "Ogaden",
}


def rotulo(painel: str) -> str:
    return ROTULO_PAINEL.get(painel, painel)


def cor(painel: str) -> str:
    return COR.get(painel, COR_RESERVA)


def ler_tsv(caminho: pathlib.Path, incluir_gene_inteiro: bool = False):
    """Le gene_ancestry_*_por_classe.tsv e devolve {alelo: [linhas na ordem do eixo]}.

    Cada linha: (rotulo_classe, n_sitios, painel, n_painel, n_hap, post_mediana).
    """
    linhas_por_classe = {}
    with open(caminho, encoding="utf-8") as f:
        header = f.readline().rstrip("\n").split("\t")
        idx = {nome: i for i, nome in enumerate(header)}
        for linha in f:
            campos = linha.rstrip("\n").split("\t")
            alelo = campos[idx["alelo"]]
            classe = campos[idx["classe"]]
            n_sit = int(campos[idx["n_sitios"]])
            n_hap = int(campos[idx["n_hap"]])
            painel = campos[idx["painel"]]
            n_painel = int(campos[idx["n_painel"]])
            post = float(campos[idx["post_mediana"]])
            linhas_por_classe.setdefault(alelo, {})[classe] = (
                n_sit, painel, n_painel, n_hap, post
            )

    dados = {}
    classes_disp = list(CLASSES_EIXO)
    if incluir_gene_inteiro:
        classes_disp = classes_disp + [(CLASSE_GENE, "Whole gene")]
    for alelo, por_classe in linhas_por_classe.items():
        linhas_alelo = []
        for chave_tsv, rotulo_eixo in classes_disp:
            if chave_tsv not in por_classe:
                continue
            n_sit, painel, n_painel, n_hap, post = por_classe[chave_tsv]
            linhas_alelo.append((rotulo_eixo, n_sit, painel, n_painel, n_hap, post))
        # ordem de exibicao no eixo Y: flancos embaixo, codificantes em cima
        # (mesma convencao usada nas versoes anteriores — inverte a ordem de leitura)
        dados[alelo] = list(reversed(linhas_alelo))
    return dados


def montar_figura(dados: dict, out_pdf: pathlib.Path):
    n_alelos = len(dados)
    fig, axes = plt.subplots(1, n_alelos, figsize=(5.6 * n_alelos, 4.3), sharey=False)
    if n_alelos == 1:
        axes = [axes]

    paineis_usados = set()

    for ax, (alelo, linhas) in zip(axes, dados.items()):
        n_classes_codificantes = sum(
            1 for rotulo_eixo, *_ in linhas if rotulo_eixo in ("Synonymous", "Non-synonymous")
        )
        for y, (classe, n_sit, painel, n_painel, n_hap, post) in enumerate(linhas):
            paineis_usados.add(painel)
            frac = n_painel / n_hap if n_hap else 0.0
            ax.barh(y, frac, height=0.52, color=cor(painel), alpha=0.85,
                    edgecolor="black", linewidth=0.6, zorder=2)
            ax.text(frac + 0.02, y, f"{rotulo(painel)}  {n_painel}/{n_hap}",
                    va="center", fontsize=8.5, zorder=3)
            ax.text(-0.03, y, f"{n_sit:,}", va="center", ha="right",
                    fontsize=7.6, color="#555", zorder=3)
        if n_classes_codificantes and n_classes_codificantes < len(linhas):
            y_div = len(linhas) - n_classes_codificantes - 0.5
            ax.axhline(y_div, color="#999", lw=0.9, ls=(0, (4, 3)), zorder=1)
            ax.text(0.985, y_div + 0.12, "flanks", fontsize=7.2, color="#777",
                    ha="right", va="bottom", style="italic")
            ax.text(0.985, y_div - 0.12, "coding core", fontsize=7.2, color="#777",
                    ha="right", va="top", style="italic")
        ax.set_yticks(range(len(linhas)))
        ax.set_yticklabels([c for c, *_ in linhas], fontsize=9)
        ax.set_xlim(0, 1.52)
        ax.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
        ax.set_xticklabels(["0", "", "50%", "", "100%"], fontsize=8)
        ax.set_xlabel("carrier haplotypes assigned", fontsize=9)
        ax.set_title(alelo, fontsize=11, fontweight="bold", pad=8)
        ax.spines[["top", "right"]].set_visible(False)
        ax.grid(axis="x", alpha=0.25, zorder=0)

    axes[0].text(-0.145, 1.045, "sites", transform=axes[0].transAxes,
                 fontsize=7.6, color="#555", ha="center")

    ordem_legenda = [p for p in ("NDAMA", "EUR_TAURINE", "CRIOLLO", "ANKOLE", "BORAN", "OGADEN")
                      if p in paineis_usados]
    ordem_legenda += sorted(p for p in paineis_usados if p not in ordem_legenda)
    leg = [Line2D([], [], marker="s", ls="", ms=8.5, mfc=cor(p), mec="black",
                  mew=0.5, label=rotulo(p)) for p in ordem_legenda]
    fig.legend(handles=leg, loc="lower center", ncol=min(len(leg), 6), frameon=False,
               fontsize=8.8, bbox_to_anchor=(0.5, -0.035))
    fig.tight_layout(rect=(0, 0.045, 1, 1))
    for ext in ("pdf", "png"):
        fig.savefig(out_pdf.with_suffix(f".{ext}"), dpi=300, bbox_inches="tight")
    print(f"-> {out_pdf.with_suffix('.pdf')}")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--tsv", type=pathlib.Path,
                     default=pathlib.Path(__file__).parents[1] /
                     "analises/03_slick/results/gene_ancestry_dedup_por_classe.tsv",
                     help="saida de gene_ancestry_by_class.py")
    ap.add_argument("--out", type=pathlib.Path, default=OUT / "Fig6_class_ancestry.pdf")
    ap.add_argument("--com-gene-inteiro", action="store_true",
                     help="inclui a linha 'gene inteiro' alem das 3 classes padrao")
    args = ap.parse_args()

    if not args.tsv.exists():
        sys.exit(f"ERRO: TSV nao encontrado: {args.tsv}\n"
                  f"Rode antes: python analises/03_slick/scripts/gene_ancestry_by_class.py ...")

    dados = ler_tsv(args.tsv, incluir_gene_inteiro=args.com_gene_inteiro)
    if not dados:
        sys.exit(f"ERRO: nenhum alelo lido de {args.tsv}")
    montar_figura(dados, args.out)


if __name__ == "__main__":
    main()
