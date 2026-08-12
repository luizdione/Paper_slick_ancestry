"""aim_informativeness.py — mapa de informatividade de ancestralidade ao longo
do locus SLICK, e reatribuicao concentrada nos SNPs mais discriminantes.

Motivacao: as mutacoes slick sao recentes demais para relogio molecular; a
ancestralidade ja e lida dos SNPs co-herdados vizinhos (motor hap_ancestry,
classes funcionais). ESTE script vai um passo alem: em vez de pontuar todos
os sitios de uma classe por igual, mede a INFORMATIVIDADE de cada sitio para
distinguir as populacoes-fonte candidatas (delta de frequencia e In de
Rosenberg 2003), acha as posicoes de alta variabilidade/alto poder
discriminante (AIMs locais), e recalcula a atribuicao (a) so nos top-K AIMs e
(b) ponderada por informatividade — para ver se a confianca (posterior/
estabilidade) sobe frente a versao "todos os sitios iguais". Nao data nada;
apenas concentra o sinal de PROVENIENCIA nos marcadores que carregam sinal.

Anti-circularidade: a informatividade e as frequencias de painel SEMPRE excluem
os haplotipos portadores (mesma regra do assign.py). Selecionar AIMs a partir
das REFERENCIAS (nao dos portadores) evita escolher sitios que so discriminam
por acaso nesta amostra pequena.

Uso:
    python aim_informativeness.py \
        --vcf phased_prlr_10mb.vcf.gz --groups groups.tsv \
        --config ../../config/slick.yaml \
        --consequences prlr_csq.tsv \
        --allele SLICK1 --contrast NDAMA EUR_TAURINE \
        --top-k 50 --out-prefix aim_slick1

Saidas:
  <prefix>_landscape.tsv  — por sitio: pos, classe, delta, In, rank
  <prefix>_reassign.tsv   — atribuicao all-sites vs top-K AIM vs In-ponderada
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

_LIB = Path(__file__).resolve().parents[1] / "lib"
if str(_LIB) not in sys.path:
    sys.path.insert(0, str(_LIB))
from hap_ancestry import io as hio            # noqa: E402
from hap_ancestry.config import LocusConfig    # noqa: E402
from hap_ancestry.assign import focal_index, carrier_haplotypes, panel_frequencies  # noqa: E402


def rosenberg_in(freqs_by_panel: np.ndarray) -> np.ndarray:
    """Informativeness for assignment I_n (Rosenberg et al. 2003), por sitio,
    para K populacoes. freqs_by_panel: matriz (n_sitios x K) com a freq do alelo
    ALT em cada painel. Alelo bialelico: p e 1-p.

    I_n = -sum_a pbar_a ln pbar_a + (1/K) sum_k sum_a p_ka ln p_ka
    onde a percorre os 2 alelos {ALT, REF}, pbar = media entre paineis.
    Retorna vetor (n_sitios,). Maior = mais informativo.
    """
    K = freqs_by_panel.shape[1]
    p = np.clip(freqs_by_panel, 1e-9, 1 - 1e-9)     # ALT
    q = 1.0 - p                                       # REF
    pbar = p.mean(axis=1)
    qbar = q.mean(axis=1)
    term_bar = -(pbar * np.log(pbar) + qbar * np.log(qbar))
    term_k = (p * np.log(p) + q * np.log(q)).mean(axis=1)
    return term_bar + term_k


def delta_contrast(F: dict, a: str, b: str) -> np.ndarray:
    """|f_a - f_b| por sitio entre dois paineis (o discriminante classico de AIM)."""
    return np.abs(F[a] - F[b])


def weighted_loglik(obs, ok, F, panels, sel, w):
    """Log-verossimilhanca PONDERADA por informatividade w (mesmos sitios sel/ok)."""
    ll = {}
    for p in panels:
        fr = np.clip(F[p][sel][ok], 1e-9, 1 - 1e-9)
        wj = w[ok]
        ll[p] = float(np.sum(wj * np.where(obs[ok] == 1, np.log(fr), np.log(1 - fr))))
    return ll


def assign_one(obs, ok, F, panels, sel, w=None, min_sites=3):
    """Atribuicao de UM haplotipo em sel; w=None -> nao ponderada. Retorna
    (painel_vencedor, posterior_do_vencedor, n_sitios_usados) ou (None,...)."""
    if ok.sum() < min_sites:
        return None, np.nan, int(ok.sum())
    if w is None:
        w = np.ones(len(sel))
    ll = weighted_loglik(obs, ok, F, panels, sel, w)
    mx = max(ll.values())
    post = {p: np.exp(v - mx) for p, v in ll.items()}
    tot = sum(post.values())
    post = {p: v / tot for p, v in post.items()}
    best = max(post, key=post.get)
    return best, post[best], int(ok.sum())


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--groups", required=True)
    ap.add_argument("--config", required=True)
    ap.add_argument("--consequences", required=True)
    ap.add_argument("--allele", required=True, help="nome do sitio focal no YAML")
    ap.add_argument("--contrast", nargs=2, default=["NDAMA", "EUR_TAURINE"],
                     help="par de paineis para o delta de AIM (o contraste de interesse)")
    ap.add_argument("--top-k", type=int, default=50)
    ap.add_argument("--out-prefix", required=True)
    args = ap.parse_args()

    cfg = LocusConfig.from_yaml(args.config)
    pos, H, samples = hio.load_phased(args.vcf)
    class_sets = hio.load_consequences(args.consequences, cfg.consequence_classes)
    group_of = hio.load_groups(args.groups)
    hap_grp = hio.hap_group_labels(samples, group_of)
    panels = hio.active_panels(hap_grp, cfg.panels, cfg.min_hap_panel)

    focal_by_name = dict(cfg.focal_sites)
    if args.allele not in focal_by_name:
        sys.exit(f"'{args.allele}' nao esta em focal_sites do YAML")
    i_foc, exato = focal_index(pos, focal_by_name[args.allele])
    if not exato:
        sys.exit("sitio focal ausente do VCF faseado (fase nao resolvida)")
    port = carrier_haplotypes(H, i_foc)
    carrier_mask = np.zeros(H.shape[1], dtype=bool)
    carrier_mask[port] = True

    # Frequencias de painel EXCLUINDO portadores (anti-circularidade)
    F = panel_frequencies(H, hap_grp, panels, carrier_mask, cfg.eps)
    Fmat = np.column_stack([F[p] for p in panels])

    # --- 1) Mapa de informatividade por sitio ---
    In = rosenberg_in(Fmat)
    a, b = args.contrast
    if a not in F or b not in F:
        sys.exit(f"contraste {a}/{b} nao esta entre os paineis ativos: {panels}")
    delta = delta_contrast(F, a, b)

    # classe funcional por sitio (para o landscape)
    pos_class = np.array(["outro"] * len(pos), dtype=object)
    for classe, s in class_sets.items():
        pos_class[np.isin(pos, list(s))] = classe

    land = pd.DataFrame({
        "pos": pos, "classe": pos_class,
        f"delta_{a}_vs_{b}": delta, "In_rosenberg": In,
    }).sort_values("In_rosenberg", ascending=False).reset_index(drop=True)
    land["rank_In"] = np.arange(1, len(land) + 1)
    land.to_csv(f"{args.out_prefix}_landscape.tsv", sep="\t", index=False,
                float_format="%.5f")

    # --- 2) Reatribuicao: all-sites vs top-K AIM vs In-ponderada ---
    # Restringe aos sitios CODIFICANTES (as classes que sustentam o claim) OU
    # gene inteiro? Usamos gene inteiro aqui: a ideia e concentrar por
    # informatividade, nao por classe. (Rode com --consequences do gene.)
    order = np.argsort(-In)                    # sitios por In decrescente
    topk_idx = np.sort(order[:args.top_k])     # indices dos top-K AIMs
    sel_all = np.arange(len(pos))
    w_in = In.copy()

    rows = []
    for h in port:
        obs_all = H[:, h]
        for tag, sel, w in (
            ("all_sites", sel_all, None),
            (f"top{args.top_k}_AIM", topk_idx, None),
            ("In_weighted", sel_all, w_in),
        ):
            obs = H[sel, h]
            ok = obs >= 0
            wsel = None if w is None else w[sel]
            best, post, nused = assign_one(obs, ok, F, panels, sel, wsel)
            rows.append({"haplotipo_col": int(h), "metodo": tag,
                         "atribuicao": best, "posterior": post,
                         "n_sitios": nused})
    reassign = pd.DataFrame(rows)
    reassign.to_csv(f"{args.out_prefix}_reassign.tsv", sep="\t", index=False,
                    float_format="%.4f")

    # Resumo: posterior mediano por metodo (concentrar nos AIMs sobe a confianca?)
    print("=== posterior mediano por metodo (maior = atribuicao mais confiante) ===")
    print(reassign.groupby("metodo")["posterior"].median().round(4).to_string())
    print(f"\n[ok] landscape: {args.out_prefix}_landscape.tsv "
          f"(top sitio In={land['In_rosenberg'].iloc[0]:.3f} @ pos "
          f"{int(land['pos'].iloc[0])})")
    print(f"[ok] reassign: {args.out_prefix}_reassign.tsv")


if __name__ == "__main__":
    main()
