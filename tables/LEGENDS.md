# Table legends

Tables 1–3 and S1–S6 are exported as TSV, one row per table row, with markdown formatting
removed. Each is built from the files in `derived/`, listed under every legend.

## Result tables

### `Table1_datasets.tsv`
Source of each dataset, data type, what it contributes and the analysis it supports.

### `Table2_alleles.tsv`
The four segregating slick alleles: position, HGVS, carrier haplotypes and animals, carrier
calls checked against the reads out of all calls and how many of those were rejected (Run 8
calls have no alignments to check), distinct coding configurations, and the background assigned
by the 39 kb coding classifier, by local ancestry over 10 Mb (three-way configuration) and by the
closest panel within ±50 kb, decided haplotype by haplotype.
From `derived/resultado_final.json`, `derived/suporte_leitura_portadores.tsv`,
`derived/local_ancestry/resumo_rfmix.json`, `validacao.json` and `vizinho.json`.

### `Table3_local_ancestry_validation.tsv`
Validation of each RFMix reference configuration: internal simulation accuracy per panel, the
weakest panel and what it is mistaken for, held-out reference haplotypes recovered, and mean
window size.
From `derived/local_ancestry/validacao.json`.

### `TableS1_read_support.tsv`
Read support for each carrier call with an alignment available: genotype after phasing,
reference and alternative read counts (`samtools mpileup -B`, mapping and base quality ≥ 20),
alternative fraction, and whether the call is kept.
From `derived/suporte_leitura_portadores.tsv` (`scripts/haplotype_robustness/suporte_leitura.py`).

### `TableS2_carrier_haplotypes.tsv`
Every carrier haplotype: its allele at the exon 2 site chr20:39,077,999, the coding-classifier
assignment at synonymous and non-synonymous sites, and the local ancestry call at the locus in
the three-way, four-way and five-way configurations.
From `derived/resultado_final.json`, `derived/exon2_dados.json`,
`derived/local_ancestry/m3.msp.tsv`, `m4.msp.tsv`, `m5.msp.tsv` and `focais_faseados.tsv`.

### `TableS3_rarefaction_batch.tsv`
By allele and class of site: the modal panel when every reference panel is drawn down to the
size of the smallest (500 replicates), and the assignment of carriers split by sequencing batch
(Run 8 and ERP150979).
From `derived/resultado_final.json`.

### `TableS4_intronic_by_distance.tsv`
Intronic assignment by distance band to the causal site: number of sites, modal panel, median
posterior and whether the band discriminates between panels.
From `derived/resultado_final.json`.

### `TableS5_ndama_detection_limit.tsv`
Highest SLICK1 frequency in N'Dama compatible with observing no copy, p = 1 − 0.05^(1/n) for n
sampled alleles.

### `TableS6_inbreeding_control.tsv`
Length of the homozygous tract at the locus against the baseline of 150 random positions in the
same 10 Mb window, each at least 1 Mb from the receptor, with the empirical *P*.
From `scripts/ancestry_haplotype/inbreeding_control.py`.

## Input tables

### `groups_iberico_declarado.tsv`, `groups_iberico_mais_alentejana.tsv`
Sample-to-panel labels used by the classifier and by the local ancestry panels: the six original
panels plus IBERIAN, the four Run 8 animals whose BioSample declares an Iberian breed (Sayaguesa,
Maronesa, Limia, Pajuna). The second file adds one Alentejana animal identified by its sample
alias only.

### `groups_afr_split.tsv`
Base label file the two above were written from; read as stored (see README, "Label tables").

### `rotulos_corrigidos_2026-10-13.tsv`
Provenance-checked breed record of every Run 8 sample: declared BioSample breed, study,
principal component, indicine fraction, category and genetic pole.

### `run9_breed_codebook.tsv`
Codebook mapping BioSample breed strings to ancestry groups, used by
`scripts/panel_labeling/apply_run9_codebook.py`.

## `derived/`

| File | Written by |
|---|---|
| `suporte_leitura_portadores.tsv` | `haplotype_robustness/suporte_leitura.py` |
| `resultado_final.json` | `haplotype_robustness/resultado_final.py` |
| `estabilidade_hierarquica.json` | `haplotype_robustness/estabilidade_hierarquica.py` |
| `sitios_decisivos.json` | `haplotype_robustness/sitios_decisivos.py` |
| `exon2_dados.json` | `haplotype_robustness/exon2_dados.py` |
| `informatividade_classe.json` | `haplotype_robustness/teste_informatividade_classe.py` |
| `distancia_intron.json` | `haplotype_robustness/teste_distancia_intron.py` |
| `contraste_intra_animal.json` | `haplotype_robustness/contraste_intra_animal.py` |
| `teste_artefatos_*.json` | `haplotype_robustness/teste_artefatos_classificador.py` (three runs, see README) |
| `fase_leituras_ERR12035149.tsv` | `haplotype_robustness/fase_por_leitura.py --tsv` |
| `figuras_novas.json`, `homozigotos_SLICK{1,2}_run8.txt` | `haplotype_robustness/dados_figuras_novas.py` |
| `endogamia_SLICK{1,2}.json` | `ancestry_haplotype/inbreeding_control.py --json` |
| `local_ancestry/consulta.tsv`, `ref_*vias.tsv`, `controle_*.tsv` | `local_ancestry/04a_rotulos.py` |
| `local_ancestry/genotipa_stats.tsv` | `local_ancestry/02_roda.sh` |
| `local_ancestry/m3.msp.tsv`, `m4.msp.tsv`, `focais_faseados.tsv` | `local_ancestry/04b_rfmix.sh` |
| `local_ancestry/m5.msp.tsv` | `local_ancestry/04c_rfmix_5vias.sh` |
| `local_ancestry/resumo_rfmix.json` | `local_ancestry/05_resumo.py` |
| `local_ancestry/vizinho.json` | `local_ancestry/06_vizinho.py` |
| `local_ancestry/validacao.json` | `local_ancestry/07_validacao.py` |

Column names inside `derived/` are as the scripts write them (Portuguese: `alelo` allele,
`classe` class, `painel` panel, `sitios` sites, `suporte` support, `usado_na_atribuicao` used in
the assignment).
