"""network.py — rede dos haplotipos no espaco dos sitios codificantes.

no = haplotipo distinto; tamanho = numero de copias; fatia = painel de origem de
cada copia; aresta = vizinho mais proximo, rotulavel pelo numero de diferencas.
Copias identicas colapsam num unico no, e esse colapso e informacao: um alelo
introduzido uma vez deixa portadores quase identicos; um alelo recorrente os
deixaria dispersos. Por isso uma rede, e nao um MDS de copias (que empilha copias
identicas no mesmo pixel).

Generico: um subplot por sitio focal do `LocusConfig`; cores, rotulos e textos
saem do `panel_style`/`allele_label`.
"""
from __future__ import annotations
import pathlib
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.patheffects as pe
from .assign import focal_index, carrier_haplotypes


def collapse(Hs, labels):
    """Agrupa haplotipos identicos (colunas de `Hs`). Devolve (tipos, composicao),
    onde composicao[k] e {painel: n_copias} do tipo k."""
    tipos, comp, chave = [], [], {}
    for k in range(Hs.shape[1]):
        t = tuple(Hs[:, k])
        if t not in chave:
            chave[t] = len(tipos)
            tipos.append(np.array(t))
            comp.append({})
        c = comp[chave[t]]
        c[labels[k]] = c.get(labels[k], 0) + 1
    return tipos, comp


def mds_layout(tipos, seed=20260804):
    """Posiciona os nos por MDS classico sobre as diferencas entre TIPOS (nao
    copias). Devolve (XY, D) — sempre os dois, mesmo com um unico tipo."""
    n = len(tipos)
    if n == 1:
        return np.zeros((1, 2)), np.zeros((1, 1))
    D = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            ok = (tipos[i] >= 0) & (tipos[j] >= 0)
            D[i, j] = D[j, i] = int((tipos[i][ok] != tipos[j][ok]).sum())
    J = np.eye(n) - np.ones((n, n)) / n
    B = -0.5 * J @ (D ** 2) @ J
    w, V = np.linalg.eigh(B)
    idx = np.argsort(w)[::-1][:2]
    XY = V[:, idx] * np.sqrt(np.clip(w[idx], 0, None))
    if XY.shape[1] == 1:
        XY = np.column_stack([XY, np.zeros(n)])
    return XY, D


def draw_pie(ax, x, y, raio, comp, colors):
    """Desenha o no como pizza da composicao por painel."""
    tot = sum(comp.values())
    ini = 90.0
    for p, k in sorted(comp.items(), key=lambda kv: -kv[1]):
        ext = 360.0 * k / tot
        ax.add_patch(plt.matplotlib.patches.Wedge(
            (x, y), raio, ini - ext, ini, facecolor=colors.get(p, "#cccccc"),
            edgecolor="black", linewidth=0.7, zorder=4))
        ini -= ext


def plot_network(pos, H, hap_grp, coding_idx, cfg, out):
    """Desenha a rede (um subplot por sitio focal) e grava PDF+PNG e o TSV de stats.

    `coding_idx` sao os indices dos sitios codificantes (uniao das classes de
    `cfg.coding_classes`), computados pelo chamador a partir das consequencias.
    A ordem de montagem dos paineis segue a ordem do `cfg.panel_style`.
    """
    colors = {p: cfg.panel_style[p]["color"] for p in cfg.panel_style}
    rotulos = {p: cfg.panel_style[p]["label"] for p in cfg.panel_style}
    panels = [p for p in cfg.panel_style if (hap_grp == p).sum() >= cfg.min_hap_panel]

    nfoc = len(cfg.focal_sites)
    # Layout VERTICAL, um alelo por linha: cabe numa coluna de pagina; na horizontal so
    # cabia em pagina inteira e os rotulos encolhiam ate ficarem ilegiveis.
    fig, axes = plt.subplots(nfoc, 1, figsize=(3.4, 3.1 * nfoc))
    if nfoc == 1:
        axes = [axes]
    axes = np.atleast_1d(axes).ravel()

    linhas = []
    for ax, (rotulo, foco) in zip(axes, cfg.focal_sites):
        i_foc, _ = focal_index(pos, foco)
        port = carrier_haplotypes(H, i_foc)

        # haplotipos dos paineis (excluindo os portadores) + os portadores
        idx, rot = [], []
        for p in panels:
            cand = np.setdiff1d(np.where(hap_grp == p)[0], port)
            idx += list(cand); rot += [p] * len(cand)
        idx += list(port); rot += ["CARRIER"] * len(port)

        Hs = H[np.ix_(coding_idx, np.array(idx))]
        tipos, comp = collapse(Hs, rot)
        XY, D = mds_layout(tipos)

        n_port_tipos = sum(1 for c in comp if "CARRIER" in c)
        print(f"[{rotulo}] {len(port)} haplotipos portadores -> {n_port_tipos} tipos "
              f"distintos | {len(tipos)} tipos no total")
        linhas.append({"alelo": rotulo, "n_portadores": len(port),
                       "tipos_portadores": n_port_tipos, "tipos_total": len(tipos)})

        # arestas: cada no ao vizinho mais proximo
        for i in range(len(tipos)):
            j = int(np.argmin(np.where(np.arange(len(tipos)) == i, np.inf, D[i])))
            ax.plot([XY[i, 0], XY[j, 0]], [XY[i, 1], XY[j, 1]],
                    color="#999999", lw=0.6, zorder=1)

        span = max(np.ptp(XY[:, 0]), np.ptp(XY[:, 1]), 1.0)
        for i, c in enumerate(comp):
            tot = sum(c.values())
            r = span * 0.018 * (1 + 0.55 * np.log1p(tot))
            eh_port = "CARRIER" in c
            comp_pizza = {k: v for k, v in c.items() if k != "CARRIER"}
            if eh_port and not comp_pizza:
                ax.add_patch(plt.Circle((XY[i, 0], XY[i, 1]), r, facecolor="white",
                                        edgecolor="black", lw=1.6, zorder=4))
            else:
                draw_pie(ax, XY[i, 0], XY[i, 1], r, comp_pizza or c, colors)
            if eh_port:
                ax.add_patch(plt.Circle((XY[i, 0], XY[i, 1]), r * 1.55, fill=False,
                                        edgecolor="#d73027", lw=2.0, zorder=6))
                # o rotulo caia sobre o anel do no vizinho, na mesma cor, e sumia.
                # Contorno branco e afastamento maior separam os dois sem mover o no.
                ax.annotate(f"{c['CARRIER']}", (XY[i, 0], XY[i, 1] + r * 2.25),
                            ha="center", fontsize=7.5, color="#d73027",
                            fontweight="bold", zorder=8,
                            path_effects=[pe.withStroke(linewidth=2.4, foreground="white")])

        ax.set_title(f"{rotulo}: {len(port)} carrier haplotypes, "
                     f"{n_port_tipos} distinct", fontsize=9, fontweight="bold", pad=4)
        ax.set_aspect("equal"); ax.axis("off")

    leg = [Line2D([], [], marker="o", ls="", ms=8, mfc=colors[p], mec="black",
                  mew=0.5, label=rotulos[p]) for p in panels]
    leg.append(Line2D([], [], marker="o", ls="", ms=11, mfc="none", mec="#d73027",
                      mew=2, label=f"carries the {cfg.allele_label}"))
    fig.legend(handles=leg, loc="lower center", ncol=2, frameon=False, fontsize=7.5,
               bbox_to_anchor=(0.5, -0.012))
    fig.tight_layout(rect=(0, 0.055, 1, 0.995))

    out = pathlib.Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(out.with_suffix("." + ext), dpi=300, bbox_inches="tight")
    pd.DataFrame(linhas).to_csv(out.with_name("haplotype_network_stats.tsv"),
                                sep="\t", index=False)
    print(f"  -> {out}")
    return linhas
