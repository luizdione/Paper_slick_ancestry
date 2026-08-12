"""assign.py — de qual populacao vem o haplotipo que carrega um alelo focal.

A atribuicao e por maxima verossimilhanca binomial: em cada sitio informativo, a
probabilidade de observar o alelo do haplotipo portador sob a frequencia de cada
painel de referencia. O painel que maximiza a soma das log-verossimilhancas e a
origem mais provavel daquele haplotipo, naquela classe de sitios.

Generico: o motor nao conhece *slick*. Sitios focais, paineis e limiares vem do
`LocusConfig`.
"""
from __future__ import annotations
import numpy as np


def focal_index(pos, focal_pos):
    """(indice do sitio focal no array `pos`, exato?).

    `exato=False` quando o sitio focal NAO esta no VCF faseado — ali a fase nao
    esta resolvida e a leitura do haplotipo portador seria inferida, nao lida.
    O chamador deve tratar isso como controle positivo (pular, nao adivinhar).
    """
    i = int(np.argmin(np.abs(pos - focal_pos)))
    return i, bool(pos[i] == focal_pos)


def carrier_haplotypes(H, i_foc):
    """Colunas (haplotipos) que carregam o alelo derivado no sitio focal."""
    return np.where(H[i_foc] == 1)[0]


def panel_frequencies(H, hap_grp, panels, carrier_mask, eps):
    """Frequencia do alelo derivado por painel, EXCLUINDO os haplotipos portadores
    (para nao se auto-avaliarem). Clipada em [eps, 1-eps] para o log nao estourar.
    """
    F = {}
    for p in panels:
        m = (hap_grp == p) & ~carrier_mask
        sub = H[:, m].astype(float)
        sub[sub < 0] = np.nan
        F[p] = np.clip(np.nanmean(sub, axis=1), eps, 1 - eps)
    return F


def _loglik(obs, ok, F, panels, sel):
    """Log-verossimilhanca de um haplotipo sob cada painel, nos sitios `sel`
    efetivamente chamados (`ok`)."""
    ll = {}
    for p in panels:
        fr = F[p][sel][ok]
        ll[p] = float(np.sum(np.where(obs[ok] == 1, np.log(fr), np.log(1 - fr))))
    return ll


def consensus(choices):
    """Painel mais escolhido entre `choices`; empate resolvido pela ordem de
    PRIMEIRA APARICAO em `choices` (deterministica).

    Substitui o antigo `max(set(choices), key=choices.count)`, cujo desempate
    seguia a ordem de iteracao de um `set` — randomizada por processo pelo
    PYTHONHASHSEED. Aquele desempate tornava a `estabilidade` do bootstrap NAO
    reproduzivel (oscilava ~0,03 entre execucoes nas classes de baixo consenso).

    Usa-se `dict.fromkeys` (ordem de insercao, estavel) e nao a ordem dos paineis:
    desempatar por prioridade de painel enviesaria o consenso para o primeiro
    painel da lista. A ordem de aparicao vem dos dados (a ordem dos haplotipos
    portadores), que e neutra e fixa.
    """
    return max(dict.fromkeys(choices), key=choices.count)


def winners(H, sel, port, F, panels, min_sites_class):
    """Painel de maxima verossimilhanca de CADA haplotipo portador, no conjunto de
    sitios `sel`. Haplotipos com menos de `min_sites_class` sitios chamados sao
    omitidos. Devolve a lista de paineis vencedores (um por haplotipo elegivel).

    E a base tanto do bootstrap (reamostra `sel`) quanto do argmax do relatorio.
    """
    out = []
    for h in port:
        obs = H[sel, h]
        ok = obs >= 0
        if ok.sum() < min_sites_class:
            continue
        ll = _loglik(obs, ok, F, panels, sel)
        out.append(max(ll, key=ll.get))
    return out


def assign_with_posterior(H, sel, port, F, panels, min_sites_class):
    """Como `winners`, mas devolve tambem o posterior do vencedor de cada
    haplotipo (softmax das log-verossimilhancas, prior uniforme). Para a tabela e
    a figura do relatorio.

    Retorna (contagem_por_painel, lista_do_posterior_do_vencedor).
    """
    cont: dict = {}
    post_vencedor = []
    for h in port:
        obs = H[sel, h]
        ok = obs >= 0
        if ok.sum() < min_sites_class:
            continue
        ll = _loglik(obs, ok, F, panels, sel)
        mx = max(ll.values())
        post = {p: np.exp(v - mx) for p, v in ll.items()}
        tot = sum(post.values())
        post = {p: v / tot for p, v in post.items()}
        melhor = max(post, key=post.get)
        cont[melhor] = cont.get(melhor, 0) + 1
        post_vencedor.append(post[melhor])
    return cont, post_vencedor
