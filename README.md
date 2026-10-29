# Slick *PRLR* alleles in Criollo cattle — calculation code

Calculation code, input tables, result tables and figures for a study of the truncating
(*slick*) alleles of the bovine prolactin receptor gene (*PRLR*, BTA20) in Criollo cattle.
The six published slick positions (SLICK1–SLICK6) are screened in phased whole genomes,
carrier calls are checked against the aligned reads, and the haplotype background of each
segregating allele is estimated by three independent routes:

1. **Coding-site classifier, 39 kb.** Each carrier haplotype is assigned to a reference panel
   by the log-likelihood of its alleles under each panel's frequencies, stratified by
   functional class of site, with site-resampling stability, panel rarefaction, batch
   separation and removal of the most informative sites.
2. **Local ancestry, 10 Mb.** RFMix v2.03 over a 10 Mb window of BTA20 phased with BEAGLE 5.4
   on the CLARITY genetic map, in three reference configurations, each validated by internal
   simulation, held-out reference animals and an admixed positive control before any carrier
   is read.
3. **Nearest neighbour, model-free.** The number of differences between each carrier haplotype
   and the closest haplotype of every panel around the locus.

A genome-wide ancestry estimate for the Senepol breed (supervised maximum likelihood with a
lineage-paired control and sNMF replication) is kept in `scripts/ancestry_genome/`.

## Repository layout

```
scripts/
  lib/hap_ancestry/       haplotype-ancestry engine (config-driven: io, assign, bootstrap,
                          network); reusable for any locus.
  acquisition/            streaming download and per-window alignment of public reads;
                          genotype merge and joint re-phasing.
  discovery/              PRLR variant calls, per-breed frequency, Senepol site genotypes.
  haplotype_robustness/   read support of carrier calls, read-based phase, the 39 kb
                          classifier result and its robustness tests.
  local_ancestry/         10 Mb local ancestry: forced genotyping of BAMs, phasing, reference
                          panels, RFMix runs, validation and nearest-neighbour distances.
  ancestry_haplotype/     thin wrappers over the engine (assignment by class, bootstrap,
                          network, inbreeding control for tract length).
  ancestry_genome/        supervised genome-wide ancestry, titration and paired controls.
  aim/                    per-site ancestry informativeness (Rosenberg's In).
  specificity/            panel-specificity test for candidate variants.
  panel_labeling/         breed labels for the Run 8 panel.
  figures/                make_figures.py, Figures 1-4.
config/                   locus configurations (focal sites, panels, colours, consequences).
tables/                   input label tables, Tables 1-3 and S1-S6 as TSV (+ LEGENDS.md).
tables/derived/           every intermediate result the tables and figures are built from.
figures/                  Figures 1-4, PDF and PNG (+ LEGENDS.md).
data/                     placeholder; raw data are public and not versioned.
```

## Requirements

- Python 3.11 with the packages in `requirements.txt`.
- A container image with `bcftools` >= 1.17, `samtools`, `bwa`, Java with BEAGLE 5.4
  (`beagle.jar`), RFMix v2.03 (built from source, https://github.com/slowkoni/rfmix) and
  `python3` with numpy. `scripts/haplotype_robustness/teste_informatividade_classe.py` calls it
  by the name in `SLICK_TOOLS_IMAGE` (default `slick-tools:latest`).

Python scripts read raw data from `SLICK_DATA_DIR` (default `./data`). The shell steps run
inside the container and expect three mounts:

```bash
docker run --rm \
  -v "$SLICK_DATA_DIR":/data/slick_ancestry \
  -v /path/to/tools:/t \
  -v "$PWD/scripts/local_ancestry":/s \
  slick-tools:latest bash /s/01_sitios.sh
```

`/t` holds `rfmix2_bin/rfmix`; `beagle.jar` and the CLARITY map are read from the data
directory (`beagle.jar`, `p5_ancestry/clarity/clarity_angus_map.tsv`).

## Data (all public)

| Source | Accession | Used for |
|---|---|---|
| 1000 Bull Genomes, Run 8 | ENA PRJEB42783 | phased panel at *PRLR*; 10 Mb window |
| Criollo whole genomes | ENA ERP150979 | carriers; read-level validation |
| N'Dama whole genomes | ENA PRJEB36894, PRJEB39353, PRJEB39924; SRA PRJNA604048 | N'Dama panel |
| Asturiana whole genomes | ENA PRJEB38981 | Iberian reference (local ancestry) |
| Senepol liver RNA-seq | GEO GSE192471 | Senepol genotypes at the focal sites |
| BovineSNP50 genotypes, 134 breeds | Dryad doi:10.5061/dryad.th092 | Senepol genome-wide ancestry |
| CLARITY genetic map | doi:10.5281/zenodo.17909700 | phasing and RFMix |

Reference genome ARS-UCD1.2. The capture BAMs of the 10 Mb window use a local contig
(`20:34100000-44100000`): global position = local position + 34,099,999.

Files expected under `SLICK_DATA_DIR`:

```
criollo/incorporado/merged_phased_dedup.vcf.gz   phased 39 kb panel with the focal sites
criollo/bam/, criollo/bam_ndama/, criollo/bam_iberian/   capture BAMs of the 10 Mb window
haplo/prlr_allmaf_csq.tsv                        bcftools csq consequences, all MAF
haplo/aim_*_landscape.tsv                        per-site informativeness (scripts/aim)
haplo/bta20_win10Mb.vcf.gz                       Run 8, 10 Mb window, unphased
lai10mb/                                         written by scripts/local_ancestry
```

## How to run

```bash
export SLICK_DATA_DIR=/path/to/data

# 1) Read support of every carrier call, then the 39 kb classifier and its robustness tests
python scripts/haplotype_robustness/suporte_leitura.py <phased.vcf.gz> <bam_dir> <bait.fa> \
       tables/derived/suporte_leitura_portadores.tsv
python scripts/haplotype_robustness/resultado_final.py           # -> resultado_final.json
python scripts/haplotype_robustness/estabilidade_hierarquica.py
python scripts/haplotype_robustness/sitios_decisivos.py
python scripts/haplotype_robustness/exon2_dados.py
python scripts/haplotype_robustness/controle_exon2.py            # printed only
python scripts/haplotype_robustness/teste_informatividade_classe.py
python scripts/haplotype_robustness/teste_distancia_intron.py
python scripts/haplotype_robustness/contraste_intra_animal.py
python scripts/haplotype_robustness/teste_artefatos_classificador.py   # default panels
python scripts/haplotype_robustness/teste_artefatos_classificador.py \
       tables/groups_iberico_declarado.tsv config/slick_iberico_declarado.yaml \
       tables/derived/teste_artefatos_iberico_declarado.json
python scripts/haplotype_robustness/teste_artefatos_classificador.py \
       tables/groups_iberico_declarado.tsv config/slick_iberico_declarado.yaml \
       tables/derived/teste_artefatos_iberico_suportado.json '<{allele: [samples]}>'
python scripts/haplotype_robustness/fase_por_leitura.py <bam> <nameA:pos:REF:ALT> <nameB:pos:REF:ALT>

# 2) Local ancestry in the 10 Mb window (container steps marked [c])
bash   /s/01_sitios.sh            # [c] site list, global and local coordinates
bash   /s/02_roda.sh              # [c] forced genotyping of every BAM (02_genotipa.py)
bash   /s/03_merge_fase.sh        # [c] merge with Run 8 and the focal sites; BEAGLE phasing
python scripts/local_ancestry/04a_rotulos.py   # reference panels, queries, held-out controls
bash   /s/04b_rfmix.sh            # [c] 3-way and 4-way runs, controls, -G sensitivity
bash   /s/04c_rfmix_5vias.sh      # [c] 5-way run with an Asian taurine panel
python scripts/local_ancestry/05_resumo.py     # controls first, then background and carriers
bash   /s/06a_hap_1mb.sh          # [c] phased haplotypes, 1 Mb around PRLR
python scripts/local_ancestry/06_vizinho.py    # nearest-neighbour differences
python scripts/local_ancestry/07_validacao.py  # internal simulation accuracy per panel

# 3) Figures (read tables/derived/ only)
python scripts/figures/make_figures.py
```

The fourth argument of the last `teste_artefatos_classificador.py` run lists the calls marked
`usado_na_atribuicao = nao` in `tables/derived/suporte_leitura_portadores.tsv`, as a JSON object
`{"SLICK2": [...], "SLICK3": [...], "SLICK6": [...]}`. The local ancestry steps write to
`$SLICK_DATA_DIR/lai10mb/`; the copies the figures read are in `tables/derived/local_ancestry/`.
Read counts in `suporte_leitura.py` are taken with base alignment quality recalibration
disabled (`samtools mpileup -B`). All resampling uses fixed seeds.

## Label tables

`tables/groups_afr_split.tsv` is read exactly as stored. 807 of its lines contain a stray
carriage return before the tab; pandas splits those lines and the affected samples are left
without a panel label. The labels that survive this reading agree with the provenance-checked
breed records in `tables/rotulos_corrigidos_2026-10-13.tsv`, whereas normalising the line
endings places mislabelled animals, Brahman and Hanwoo among them, in the European taurine
panel. The two `groups_iberico_*.tsv` files used by the main analysis were written from the
surviving labels and contain no stray carriage returns. `.gitattributes` keeps every label file
byte-exact; do not convert their line endings.

## Archive and licence

Archived at Zenodo, doi:10.5281/zenodo.22184899. MIT licence, see `LICENSE`.
