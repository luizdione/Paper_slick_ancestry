# Geographic origin of *PRLR* slick alleles in tropical cattle — calculation code

This repository holds the calculation code, curated tables and figures for a study
of the biogeographic origin of truncating (*slick*) alleles of the bovine prolactin
receptor gene (*PRLR*, BTA20), combining two complementary analyses:

1. **Haplotype ancestry.** For each *slick* allele (SLICK1 = c.1382del, SLICK2 =
   c.1489C>T, and the candidate c.1393C>A), the coding haplotype carrying it is
   assigned to a reference population by the log-likelihood of its alleles under
   each panel's frequencies, stratified by functional class of site (synonymous,
   non-synonymous, intronic). A haplotype network and resampling-based stability
   accompany the assignment.
2. **Global genomic ancestry.** Supervised ADMIXTURE (K=3: European taurine,
   African taurine, zebu) with an unsupervised sNMF replication, using a
   lineage-paired control design, quantifies the excess of West African ancestry in
   a tropically adapted composite breed relative to its founding British lineage.

The two analyses use different animals and share no method; that they point to the
same source population is the study's central result.

## Repository layout

```
scripts/
  lib/hap_ancestry/     generic haplotype-ancestry engine (config-driven; io,
                        assign, bootstrap, network). Reusable for any locus.
  acquisition/          streaming download + on-the-fly filtering/alignment of
                        public reads to per-window BAMs; genotype merge & re-phasing.
  discovery/            unbiased discovery of PRLR truncating variants, per-breed
                        frequency, cohort matrix, RNA-seq site validation, HGVS
                        c. position, shared-haplotype length.
  ancestry_genome/      supervised maximum-likelihood ancestry, ADMIXTURE titration,
                        lineage-paired and leave-one-out controls.
  ancestry_haplotype/   thin wrappers over the engine: assignment by functional
                        class, bootstrap stability, haplotype network, inbreeding
                        control, tract sharing, phase-by-sharing.
  specificity/          panel-specificity test that excludes loss-of-function
                        candidates not restricted to the phenotype group.
  aim/                  ancestry-informative-marker informativeness map.
  figures/              figure and table generators.
config/slick.yaml       the only locus-specific file: focal sites, panels, colours,
                        consequence map, thresholds.
tables/                 curated result tables (+ LEGENDS.md).
figures/                final figures (+ LEGENDS.md).
data/                   placeholder; heavy inputs are reproduced, not versioned.
```

## Requirements

- Python 3.10+ with the packages in `requirements.txt` (numpy, scipy, pandas,
  matplotlib, pyyaml).
- External command-line tools for the shell pipelines (install via conda-forge /
  bioconda): `bcftools`, `samtools`, `bwa`, `bbmap` (bbduk), `aria2c`, `plink`/
  `plink2`, `beagle` (Java), `admixture`.

The Python scripts read their large inputs from a data root given by the
`SLICK_DATA_DIR` environment variable (default `./data`). Point it at your local
data directory; the shell pipelines take the same variable.

## Data provenance (all public)

- **1000 Bull Genomes Project, Run 8** — WGS reference panels.
- **Decker et al. 2014** (Dryad `10.5061/dryad.th092`) — BovineSNP50 array, 134 breeds.
- **ENA `ERP150979`** — Criollo (Limonero, Hartón del Valle, Blanco Orejinegro,
  Costeño con Cuernos) sequencing reads.
- **GEO `GSE192471`** — Senepol liver RNA-seq (site-level genotype validation only).
- Reference genome **ARS-UCD1.2** (chromosome BTA20 identical to ARS-UCD1.3 in the
  *PRLR* window); genetic map **CLARITY** (already in ARS-UCD1.2).

Access to controlled datasets follows each source's own terms.

## How to run (pipeline order)

```bash
export SLICK_DATA_DIR=/path/to/your/data

# 1) Acquisition — stream public reads to per-window BAMs and genotype the sites
bash scripts/acquisition/fetch_v2.sh
bash scripts/acquisition/genotipa_criollos.sh
bash scripts/acquisition/incorpora_novos_portadores.sh   # merge + joint re-phasing

# 2) Discovery — find PRLR truncating variants and their per-breed frequency
bash scripts/discovery/axis1_prlr_discovery.sh
bash scripts/discovery/axis1b_freq_by_breed.sh
python scripts/discovery/calc_c_pos.py

# 3) Global genomic ancestry (supervised ADMIXTURE + controls)
python scripts/ancestry_genome/paired_control.py
python scripts/ancestry_genome/titration_admixture.py

# 4) Haplotype ancestry by functional class (uses config/slick.yaml by default)
python scripts/ancestry_haplotype/gene_ancestry_by_class.py  --vcf <phased.vcf.gz> --csq <csq.tsv> --groups <groups.tsv> --out-prefix out/slick
python scripts/ancestry_haplotype/bootstrap_ancestry_class.py --vcf <phased.vcf.gz> --csq <csq.tsv> --groups <groups.tsv> --out-prefix out/slick
python scripts/ancestry_haplotype/haplotype_network.py        --vcf <phased.vcf.gz> --csq <csq.tsv> --groups <groups.tsv> --out out/slick_net.pdf
python scripts/ancestry_haplotype/inbreeding_control.py

# 5) Specificity test and figures
python scripts/specificity/specificity_filter.py
python scripts/figures/fig_class_ancestry.py
python scripts/figures/fig_genomic_ancestry.py
```

The haplotype-ancestry wrappers accept `--config another_locus.yaml` to run the
engine on a different gene/project without changing code (see
`scripts/lib/hap_ancestry/README.md`).

## License

MIT — see `LICENSE`.
