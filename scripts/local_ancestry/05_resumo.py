# -*- coding: utf-8 -*-
"""05_resumo.py -- le as saidas do RFMix e resume, NESTA ORDEM (regra 1 do projeto):

  1. CONTROLES: animais de referencia deixados de fora recuperam o proprio painel?
     (sobra dos paineis grandes na execucao m4; 3 N'Dama + 2 ibericos na execucao c4)
     e os africanos misturados (Ankole, Boran, Ogaden) saem como mistura indicino x N'Dama?
  2. FUNDO: ancestralidade dos haplotipos Criollo NAO portadores na janela do PRLR.
  3. PORTADORES: chamada e posterior na janela do PRLR, comprimento do trato de mesma
     ancestralidade em torno do locus, e o alelo da variante do exon 2 (39.077.999) ao lado.
  4. SENSIBILIDADE: a chamada dos portadores com -G 15 e -G 150.

Saida: `$SLICK_DATA_DIR/lai10mb/rfmix/resumo_rfmix.json` + impressao.
"""
import os
import collections
import json
import pathlib

RAIZ = pathlib.Path(__file__).resolve().parents[2]   # raiz do repositorio
DATA = pathlib.Path(os.environ.get("SLICK_DATA_DIR", RAIZ / "data"))   # dados brutos (VCF, BAM), fora do repositorio

R = (DATA / "lai10mb/rfmix")
O = R.parent
FOCO = {"SLICK1": 39099213, "SLICK3": 39099226, "SLICK6": 39099267, "SLICK2": 39099321}
EXON2 = 39077999
PRLR = 39099267                     # meio do bloco dos sitios slick (39.099.213-39.099.321)
SEM_SUPORTE = {"SLICK2": {"ERR12035125", "ERR12035133"}, "SLICK3": {"ERR12035159"},
               "SLICK6": {"ERR12035135"}}


def le_tsv(p):
    return dict(l.split("\t")[:2] for l in p.read_text().split("\n") if l.strip())


def le_msp(nome):
    linhas = (R / f"{nome}.msp.tsv").read_text().splitlines()
    codigos = dict(x.split("=") for x in linhas[0].split(":", 1)[1].split())
    cod = {int(v): k for k, v in codigos.items()}
    cab = linhas[1].lstrip("#").split("\t")
    haps = cab[6:]
    jan = []
    for l in linhas[2:]:
        c = l.split("\t")
        jan.append({"spos": int(c[1]), "epos": int(c[2]), "sg": float(c[3]), "eg": float(c[4]),
                    "call": [cod[int(x)] for x in c[6:]]})
    return haps, jan


def le_fb_no_prlr(nome):
    with open(R / f"{nome}.fb.tsv") as f:
        f.readline()
        cab = f.readline().rstrip("\n").split("\t")
        melhor, dist = None, None
        for l in f:
            c = l.rstrip("\n").split("\t")
            d = abs(int(c[1]) - PRLR)
            if dist is None or d < dist:
                melhor, dist = c, d
            elif d > dist:
                break
    post = collections.defaultdict(dict)
    for k, v in zip(cab[4:], melhor[4:]):
        s, h, pop = k.split(":::")
        post[(s, int(h.replace("hap", "")) - 1)][pop] = float(v)
    return post, int(melhor[1])


def janela_prlr(jan):
    return next(i for i, j in enumerate(jan) if j["spos"] <= PRLR <= j["epos"])


def trato(jan, k, i0):
    alvo = jan[i0]["call"][k]
    a = i0
    while a > 0 and jan[a - 1]["call"][k] == alvo:
        a -= 1
    b = i0
    while b < len(jan) - 1 and jan[b + 1]["call"][k] == alvo:
        b += 1
    return {"ancestralidade": alvo, "mb": round((jan[b]["epos"] - jan[a]["spos"]) / 1e6, 2),
            "cm": round(jan[b]["eg"] - jan[a]["sg"], 2), "inicio": jan[a]["spos"],
            "fim": jan[b]["epos"], "na_borda": a == 0 or b == len(jan) - 1}


def frac_genoma(jan, k):
    tot = collections.Counter()
    for j in jan:
        tot[j["call"][k]] += j["eg"] - j["sg"]
    s = sum(tot.values())
    return {p: round(v / s, 3) for p, v in tot.most_common()}


def hapid(h):
    s, i = h.rsplit(".", 1)
    return s, int(i)


def main():
    rot = {**le_tsv(O / "consulta.tsv"), **le_tsv(O / "controle_grandes.tsv"),
           **le_tsv(O / "controle_pequenos.tsv")}
    # genotipos faseados dos sitios focais
    lin = (R / "focais_faseados.tsv").read_text().splitlines()
    amostras = lin[0].split("\t")[1:]
    gt = {}
    for l in lin[1:]:
        c = l.split("\t")
        gt[int(c[0])] = {s: g.split("|") for s, g in zip(amostras, c[1:])}
    rel = {}

    for nome in ("m4", "m3", "c4"):
        haps, jan = le_msp(nome)
        post, pos_fb = le_fb_no_prlr(nome)
        i0 = janela_prlr(jan)
        J = jan[i0]
        r = {"janela_prlr": [J["spos"], J["epos"], J["sg"], J["eg"]], "n_janelas": len(jan),
             "fb_posicao": pos_fb}
        print(f"\n######## {nome}: {len(jan)} janelas; janela do PRLR {J['spos']:,}-{J['epos']:,} "
              f"({J['eg'] - J['sg']:.2f} cM)")

        # 1. controles
        ctrl = {}
        for k, h in enumerate(haps):
            s, i = hapid(h)
            if nome != "c4" and rot.get(s) not in ("EUR_TAURINE", "INDICINE"):
                continue
            if nome == "c4" or s in le_tsv(O / "controle_grandes.tsv"):
                ctrl.setdefault(rot[s], []).append(
                    {"hap": h, "prlr": J["call"][k], "genoma": frac_genoma(jan, k)})
        if ctrl:
            r["controles"] = {}
            for pan, lst in ctrl.items():
                certo_prlr = sum(x["prlr"] == pan or (nome == "m3" and pan == "EUR_TAURINE"
                                 and x["prlr"] == "TAURINO_EUROPEU") for x in lst)
                med = sum(x["genoma"].get(pan if nome != "m3" else
                          {"EUR_TAURINE": "TAURINO_EUROPEU", "IBERIAN": "TAURINO_EUROPEU",
                           "NDAMA": "TAURINO_AFRICANO"}.get(pan, pan), 0) for x in lst) / len(lst)
                r["controles"][pan] = {"n_hap": len(lst), "certo_no_prlr": certo_prlr,
                                       "fracao_media_do_proprio_painel": round(med, 3),
                                       "detalhe": lst}
                print(f"  CONTROLE {pan:<12} {len(lst)} hap: proprio painel no PRLR "
                      f"{certo_prlr}/{len(lst)}; fracao media do proprio painel na janela {med:.2f}")
        if nome == "c4":
            rel[nome] = r
            continue

        # africanos misturados: controle positivo de mistura
        mist = collections.defaultdict(collections.Counter)
        for k, h in enumerate(haps):
            s, _ = hapid(h)
            if rot.get(s) in ("ANKOLE", "BORAN", "OGADEN"):
                for p, v in frac_genoma(jan, k).items():
                    mist[rot[s]][p] += v
        r["africanos_misturados"] = {}
        for pop, c in mist.items():
            n = sum(c.values())
            r["africanos_misturados"][pop] = {p: round(v / n, 3) for p, v in c.most_common()}
            print(f"  MISTURA {pop:<7} " + " ".join(f"{p}:{v / n:.2f}" for p, v in c.most_common()))

        # 2 e 3. Criollo: fundo e portadores
        portador = {}
        for a, p in FOCO.items():
            for s in amostras:
                if rot.get(s) != "CRIOLLO" or s in SEM_SUPORTE.get(a, set()):
                    continue
                for i in (0, 1):
                    if gt[p][s][i] == "1":
                        portador.setdefault((s, i), []).append(a)
        fundo = collections.Counter()
        port = []
        for k, h in enumerate(haps):
            s, i = hapid(h)
            if rot.get(s) != "CRIOLLO":
                continue
            if (s, i) in portador:
                port.append({"amostra": s, "hap": i, "alelos": portador[(s, i)],
                             "prlr": J["call"][k], "posterior": post[(s, i)],
                             "exon2_alt": gt[EXON2][s][i] == "1", "trato": trato(jan, k, i0)})
            else:
                fundo[J["call"][k]] += 1
        n = sum(fundo.values())
        r["fundo_criollo_nao_portador"] = {"n_hap": n, **{p: v for p, v in fundo.most_common()}}
        print(f"  FUNDO Criollo nao portador ({n} hap) no PRLR: "
              + " ".join(f"{p}:{v}" for p, v in fundo.most_common()))
        r["portadores"] = port
        for a in FOCO:
            lst = [x for x in port if a in x["alelos"]]
            c = collections.Counter(x["prlr"] for x in lst)
            print(f"  {a}: {len(lst)} hap -> " + " ".join(f"{p}:{v}" for p, v in c.most_common()))
            for x in lst:
                pp = max(x["posterior"].values())
                t = x["trato"]
                print(f"     {x['amostra']:<13} h{x['hap']} {'+'.join(x['alelos']):<14} "
                      f"{x['prlr']:<16} post {pp:.2f} | trato {t['mb']:>5} Mb {t['cm']:>5} cM"
                      f"{' (borda)' if t['na_borda'] else ''} | exon2 {'ALT' if x['exon2_alt'] else 'REF'}")
        rel[nome] = r

    # 4. sensibilidade ao -G
    rel["sensibilidade_G"] = {}
    for nome in ("g15", "g150"):
        haps, jan = le_msp(nome)
        i0 = janela_prlr(jan)
        sens = {}
        for k, h in enumerate(haps):
            s, i = hapid(h)
            for a, p in FOCO.items():
                if rot.get(s) == "CRIOLLO" and s not in SEM_SUPORTE.get(a, set()) and gt[p][s][i] == "1":
                    sens.setdefault(a, collections.Counter())[jan[i0]["call"][k]] += 1
        rel["sensibilidade_G"][nome] = {a: dict(c) for a, c in sens.items()}
        print(f"\n  SENSIBILIDADE {nome}: " + " | ".join(
            f"{a} " + " ".join(f"{p}:{v}" for p, v in c.most_common()) for a, c in sens.items()))

    (R / "resumo_rfmix.json").write_text(json.dumps(rel, ensure_ascii=False, indent=1),
                                         encoding="utf-8")
    print(f"\n-> {R / 'resumo_rfmix.json'}")


if __name__ == "__main__":
    main()
