"""titration_admixture.py — titulacao/limite de deteccao do ensaio de
ancestralidade usando ADMIXTURE supervisionado (substitui a titration.py do SQP).

Desenho: hospedeiro Jersey (JER), doador N'Dama (NDAM), introgressao em doses
conhecidas, 40 individuos sinteticos por dose, ambos os alelos de cada
individuo sorteados preservando o espectro de frequencia empirico dos pools
reais. A diferenca e o ESTIMADOR: aqui os sinteticos entram num
`.bed` junto das referencias EUT/AFT/ZEB pinadas e o ADMIXTURE supervisionado K=3
le a fracao africana (AFT) de cada sintetico. Curva dose-resposta + limite de
deteccao (menor dose cujo IC bootstrap separa do controle dose-0) + checagem de
monotonicidade (dose conhecida deve produzir resposta monotonica crescente).

Circularidade (declarada): o doador NDAM tambem compoe o painel AFT — mesmo
desenho da titulacao SQP que este script substitui; a titulacao mede
RESOLUCAO do ensaio, nao predicao com doador retido. Rode `--holdout-donor` para
a variante mais conservadora (NDAM fora do AFT).

Fluxo (2 passos, com Docker no meio):
  1) python titration_admixture.py gen  ... -> escreve titr.vcf.gz, titr.pop, dose_index.tsv
  2) [Docker] plink --vcf titr.vcf.gz --make-bed --out titr ; admixture --supervised titr.bed 3
  3) python titration_admixture.py analyze --q titr.3.Q --fam titr.fam --dose-index dose_index.tsv
O passo 2 e impresso pronto por `gen` (mesma imagem rocker/rstudio:4.2 usada no fluxo).
"""
from __future__ import annotations

import argparse
import gzip
import os
import pathlib
import sys

import numpy as np
import pandas as pd

DATA_DIR = pathlib.Path(os.environ.get("SLICK_DATA_DIR", "data"))
ZDIR = str(DATA_DIR / "senepol_widde" / "decker2014" / "plink" / "zebu")
DOSES = [0.0, 0.005, 0.01, 0.02, 0.03, 0.05, 0.10, 0.25, 0.50, 1.0]
N_PER_DOSE = 40
SEED = 20260811


def load_raw(zdir):
    """Le all.raw (dosagens 0/1/2, NA) + panel.pop alinhado. Retorna
    (dos int8 [N x M] com -1=missing, fam DataFrame[FID,IID,pop], marker_names)."""
    raw = pd.read_csv(zdir + r"\all.raw", sep=r"\s+")
    meta = raw.iloc[:, :6]
    markers = raw.columns[6:]
    dos = raw.iloc[:, 6:].to_numpy(dtype=float)
    dos[np.isnan(dos)] = -1
    dos = dos.astype(np.int8)
    pop = pd.read_csv(zdir + r"\panel.pop", header=None, names=["pop"])
    fam = pd.DataFrame({"FID": meta.iloc[:, 0].values,
                        "IID": meta.iloc[:, 1].values,
                        "pop": pop["pop"].values})
    return dos, fam, list(markers)


def pool_freq(dos_rows):
    """Frequencia do alelo contado por marcador (media/2, nan-aware sobre -1)."""
    x = dos_rows.astype(float)
    x[x < 0] = np.nan
    return np.clip(np.nanmean(x, axis=0) / 2.0, 0.0, 1.0)


def make_synthetics(host_freq, donor_freq, dose, n, rng):
    """MODO INDEPENDENTE (obsoleto p/ ADMIXTURE): cada alelo ~ Bernoulli(freq do
    pool). NAO preserva LD -> ADMIXTURE supervisionado nao ve a dose (o sinal
    africano espalhado e absorvido como ruido europeu). Mantido so como controle
    negativo documentado; use make_synthetics_blocks para ADMIXTURE."""
    M = len(host_freq)
    out = np.empty((n, M), dtype=np.int8)
    for i in range(n):
        a_src = rng.random(M) < dose
        a_fr = np.where(a_src, donor_freq, host_freq)
        allele_a = (rng.random(M) < a_fr).astype(np.int8)
        b_src = rng.random(M) < dose
        b_fr = np.where(b_src, donor_freq, host_freq)
        allele_b = (rng.random(M) < b_fr).astype(np.int8)
        out[i] = allele_a + allele_b
    return out


def make_synthetics_blocks(host_dos, donor_dos, dose, n, rng, block):
    """MODO BLOCO (preserva LD): o genoma (colunas em ordem de cromossomo) e
    tiled em blocos de `block` marcadores; cada bloco e COPIADO inteiro do
    genotipo de um individuo real — do pool doador (N'Dama) com prob `dose`,
    senao do hospedeiro (Jersey). Injeta blocos de haplotipo africano REAL, com
    o LD real, que o ADMIXTURE detecta (ao contrario de alelos independentes).
    Copia genotipo (dosagem 0/1/2, -1=missing); e uma aproximacao do mosaico de
    ancestralidade sem exigir fase — suficiente para calibrar a resolucao."""
    M = host_dos.shape[1]
    nblocks = int(np.ceil(M / block))
    out = np.empty((n, M), dtype=np.int8)
    for i in range(n):
        row = np.empty(M, dtype=np.int8)
        for b in range(nblocks):
            s, e = b * block, min((b + 1) * block, M)
            pool = donor_dos if rng.random() < dose else host_dos
            a = int(rng.integers(len(pool)))
            row[s:e] = pool[a, s:e]
        out[i] = row
    return out


def cmd_gen(args):
    rng = np.random.default_rng(SEED)
    dos, fam, markers = load_raw(args.zdir)
    M = len(markers)

    host = dos[fam["FID"].values == args.host]
    donor = dos[fam["FID"].values == args.donor]
    if len(host) == 0 or len(donor) == 0:
        sys.exit(f"host {args.host} (n={len(host)}) ou donor {args.donor} "
                 f"(n={len(donor)}) ausente")
    host_freq = pool_freq(host)
    donor_freq = pool_freq(donor)

    def gen_dose(d):
        if args.ld_block > 0:
            return make_synthetics_blocks(host, donor, d, args.n_per_dose, rng,
                                          args.ld_block)
        return make_synthetics(host_freq, donor_freq, d, args.n_per_dose, rng)

    # Referencias a pinar: EUT/AFT/ZEB. Opcao holdout: tira o doador do AFT.
    is_ref = fam["pop"].isin(["EUT", "AFT", "ZEB"]).values
    if args.holdout_donor:
        is_ref = is_ref & ~(fam["FID"].values == args.donor)
    ref_idx = np.where(is_ref)[0]
    ref_dos = dos[ref_idx]
    ref_pop = fam["pop"].values[ref_idx]
    ref_iid = fam["IID"].values[ref_idx]

    # Sinteticos por dose
    synth_blocks, dose_rows = [], []
    for d in DOSES:
        s = gen_dose(d)
        synth_blocks.append(s)
        for k in range(args.n_per_dose):
            dose_rows.append((f"SYN_d{d:.4f}_{k:03d}", d))
    synth = np.vstack(synth_blocks)
    synth_iid = [r[0] for r in dose_rows]

    # Matriz combinada [refs; sinteticos] e rotulos .pop
    combined = np.vstack([ref_dos, synth])           # (N x M)
    sample_iids = [str(x) for x in ref_iid] + synth_iid
    pop_labels = list(ref_pop) + ["-"] * len(synth_iid)
    N = combined.shape[0]

    # Escreve VCF (posicao/aleles placeholder — ADMIXTURE ignora mapa)
    gt_map = {0: "0/0", 1: "0/1", 2: "1/1", -1: "./."}
    with gzip.open(args.out_vcf, "wt") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
                 + "\t".join(sample_iids) + "\n")
        col = combined.T  # (M x N)
        for m in range(M):
            gts = [gt_map[int(v)] for v in col[m]]
            fh.write(f"20\t{m+1}\t{markers[m].split('_')[0] or ('s'+str(m))}"
                     f"\tA\tC\t99\tPASS\t.\tGT\t" + "\t".join(gts) + "\n")

    with open(args.out_pop, "w") as fh:
        fh.write("\n".join(pop_labels) + "\n")
    pd.DataFrame(dose_rows, columns=["IID", "dose"]).to_csv(
        args.dose_index, sep="\t", index=False)

    print(f"[gen] {N} amostras ({len(ref_idx)} refs + {len(synth_iid)} sinteticos), "
          f"{M} marcadores. holdout_donor={args.holdout_donor}", file=sys.stderr)
    print(f"[gen] VCF={args.out_vcf}  POP={args.out_pop}  DOSE_IDX={args.dose_index}",
          file=sys.stderr)
    base = args.out_vcf.replace(".vcf.gz", "")
    host_mount = os.environ.get("SLICK_DATA_DIR", "$PWD/data")
    print("\n# ---- rode no Docker (mesma imagem do fluxo) ----")
    print(f"# docker run --rm -v {host_mount}:/data rocker/rstudio:4.2 bash -c '"
          f"cd /tmp && tar xzf /data/tools/admixture_linux-1.3.0.tar.gz; "
          f"ADM=$(find /tmp -name admixture -type f|head -1); "
          f"OUT=<DIR_DO_VCF_EM_/data>; cd $OUT; "
          f"plink --vcf {args.out_vcf.split(chr(92))[-1]} --make-bed --out titr --allow-extra-chr --double-id; "
          f"cp {args.out_pop.split(chr(92))[-1]} titr.pop; "
          f'"$ADM" --supervised titr.bed 3 -j4 -s 42 > titr.log 2>&1; echo EXIT $?\'')


def cmd_analyze(args):
    Q = pd.read_csv(args.q, sep=r"\s+", header=None)
    fam = pd.read_csv(args.fam, sep=r"\s+", header=None, usecols=[1], names=["IID"])
    di = pd.read_csv(args.dose_index, sep="\t")
    df = pd.concat([fam.reset_index(drop=True), Q], axis=1)

    # Identifica a coluna AFT pelas referencias AFT (as que estao no .pop original).
    # Aqui detectamos pela coluna que e ~1 nas linhas cujo IID NAO e sintetico e
    # que o .pop marcou AFT — mas o .fam nao tem o rotulo; usamos o dose_index:
    # sinteticos tem prefixo SYN_. Para achar a coluna AFT, usamos o fato de que
    # sinteticos de dose alta MAXIMIZAM a coluna africana.
    syn = df[df["IID"].astype(str).str.startswith("SYN_")].merge(di, on="IID")
    qcols = list(range(Q.shape[1]))
    # coluna africana = a que mais cresce com a dose (correlacao dose vs coluna)
    corrs = {c: np.corrcoef(syn["dose"], syn[c])[0, 1] for c in qcols}
    AFT = max(corrs, key=corrs.get)
    print(f"[analyze] coluna AFT = {AFT} (corr dose~col: "
          f"{ {c: round(v,3) for c,v in corrs.items()} })")

    syn["afr"] = syn[AFT] * 100.0
    rng = np.random.default_rng(SEED)
    rows = []
    for d, g in syn.groupby("dose"):
        vals = g["afr"].to_numpy()
        boots = np.array([rng.choice(vals, len(vals), replace=True).mean()
                          for _ in range(2000)])
        lo, hi = np.percentile(boots, [2.5, 97.5])
        rows.append({"dose_pct": d * 100, "afr_mean": vals.mean(),
                     "ic95_lo": lo, "ic95_hi": hi, "n": len(vals)})
    res = pd.DataFrame(rows).sort_values("dose_pct")

    base_hi = res.iloc[0]["ic95_hi"]  # teto do IC do controle dose-0
    res["separa_do_controle"] = res["ic95_lo"] > base_hi
    det = res[res["separa_do_controle"] & (res["dose_pct"] > 0)]
    limite = det["dose_pct"].min() if len(det) else np.nan
    mono = bool(np.all(np.diff(res["afr_mean"].to_numpy()) >= -1e-9))

    print("\n=== TITULACAO ADMIXTURE supervisionado (curva dose-resposta) ===")
    print(res.to_string(index=False, float_format=lambda x: f"{x:.3f}"))
    print(f"\n  monotonica crescente: {mono}  (dose conhecida deve produzir resposta monotonica)")
    print(f"  limite de deteccao (menor dose que separa do controle): {limite}%")
    res.to_csv(args.out, sep="\t", index=False, float_format="%.4f")
    print(f"[ok] {args.out}")


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    g = sub.add_parser("gen")
    g.add_argument("--zdir", default=ZDIR)
    g.add_argument("--host", default="JER")
    g.add_argument("--donor", default="NDAM")
    g.add_argument("--n-per-dose", type=int, default=N_PER_DOSE)
    g.add_argument("--ld-block", type=int, default=300,
                   help="tamanho do bloco (marcadores) p/ preservar LD; 0=alelos independentes")
    g.add_argument("--holdout-donor", action="store_true")
    g.add_argument("--out-vcf", required=True)
    g.add_argument("--out-pop", required=True)
    g.add_argument("--dose-index", required=True)
    g.set_defaults(func=cmd_gen)
    a = sub.add_parser("analyze")
    a.add_argument("--q", required=True)
    a.add_argument("--fam", required=True)
    a.add_argument("--dose-index", required=True)
    a.add_argument("--out", required=True)
    a.set_defaults(func=cmd_analyze)
    args = ap.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
