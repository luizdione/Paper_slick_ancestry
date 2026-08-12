# Table legends

Curated result tables. Column headers are as emitted by the scripts (Portuguese
tokens: `alelo`=allele, `classe`=class, `n_sitios`=n sites, `n_hap`=n haplotypes,
`painel`=panel, `fracao`=fraction, `estabilidade`=stability, `so sinonimas`=
synonymous only, `so nao sinonimas`=non-synonymous only, `so introns`=introns only,
`gene inteiro`=whole gene). Other tables in the article that are computed inline
(inbreeding control; N'Dama non-detection limit; the Senepol-vs-paired-control
composition) are reproduced by running the corresponding scripts.

### `gene_ancestry_dedup_por_classe.tsv`
Ancestral panel assigned to each *slick* allele's carrier haplotypes by functional
class of site, with the number of carrier haplotypes matching that panel and the
median posterior. Produced by `scripts/ancestry_haplotype/gene_ancestry_by_class.py`
on the deduplicated panel. Note the coding-site counts: 24 synonymous + 32
non-synonymous = 56 coding sites.

### `bootstrap_dedup_bootstrap_classe.tsv`
Resampling stability of the assignment above: fraction assigned to the winning
panel, 95% interval, and stability across 2000 site-resamples, per allele and class.
Produced by `scripts/ancestry_haplotype/bootstrap_ancestry_class.py`.

### `haplotype_network_stats.tsv`
Haplotype-network summary per allele: number of carriers, number of distinct carrier
haplotype types, and total distinct types. Produced by
`scripts/ancestry_haplotype/haplotype_network.py`.

### `Table1_PRLR_indels.tsv`
Truncating/indel variants discovered in the *PRLR* CDS by the unbiased scan
(position, alleles, quality, allele count, frame effect, class). The SLICK1
frameshift is at 20:39,099,213. Produced by `scripts/discovery/axis1_prlr_discovery.sh`.

### `Table2_carriers.tsv`
Genotypes of the carriers of each discovered variant (sample, site, genotype, depth,
genotype quality, breed). Produced by `scripts/discovery/axis1b_freq_by_breed.sh`.

### `TableS1_cohort.tsv`
Composition of the reference cohort used throughout: ancestral group, breed, and
number of animals. Produced by `scripts/discovery/step0_build_cohort.sh`.
