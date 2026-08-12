"""bootstrap.py — estabilidade da atribuicao por REAMOSTRAGEM DE SITIOS.

Nao se reamostra haplotipos: reamostrar cinco portadores com reposicao mede o
tamanho da amostra de portadores, que e fixo e conhecido, e devolveria um
intervalo largo por construcao sem informar nada sobre a atribuicao. O que se
reamostra sao os SITIOS dentro de cada classe. A pergunta que isso responde e a
certa: dado este conjunto de sitios informativos, quao estavel e a populacao que
cada haplotipo escolhe? Se a atribuicao depende de um punhado de sitios sortudos,
a reamostragem a desmancha; se esta distribuida, ela se mantem.

O consumo do RNG e sensivel a ordem: o chamador deve iterar (alelo, classe) na
mesma ordem em toda execucao para os intervalos serem reproduziveis com uma dada
semente.
"""
from __future__ import annotations
import numpy as np
from .assign import winners, consensus


def stability(H, sel, port, F, panels, min_sites_class, reps, rng):
    """Estabilidade da atribuicao numa classe de sitios.

    Devolve dict com: o painel vencedor observado; a fracao de haplotipos nele,
    com IC95 percentil sobre as reamostragens; e a fracao de replicas em que esse
    mesmo painel venceu o consenso (a `estabilidade`, medida direta de robustez).
    Retorna None se nenhum haplotipo tem sitios chamados suficientes.
    """
    obs = winners(H, sel, port, F, panels, min_sites_class)
    if not obs:
        return None
    vencedor = consensus(obs)
    frac_obs = obs.count(vencedor) / len(obs)

    fracs, venceu = [], 0
    for _ in range(reps):
        s = rng.choice(sel, size=len(sel), replace=True)   # reamostra SITIOS
        e = winners(H, s, port, F, panels, min_sites_class)
        if not e:
            continue
        fracs.append(e.count(vencedor) / len(e))
        if consensus(e) == vencedor:
            venceu += 1
    lo, hi = np.percentile(fracs, [2.5, 97.5]) if fracs else (np.nan, np.nan)
    estab = venceu / len(fracs) if fracs else float("nan")
    return {"vencedor": vencedor, "frac_obs": frac_obs, "ic_lo": lo, "ic_hi": hi,
            "estabilidade": estab, "n_hap": len(obs), "reps": len(fracs)}
