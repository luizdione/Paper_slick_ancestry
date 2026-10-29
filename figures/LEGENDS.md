# Figure legends

All four figures are drawn by `scripts/figures/make_figures.py` (`python
scripts/figures/make_figures.py [1 2 3 4]`) from the files in `tables/derived/`; no number is
typed into the script except the coordinates and protein domains of Figure 1, which are
annotation. Each figure is saved as vector PDF and as PNG.

### Figure 1 — The *PRLR* locus and the six published slick alleles
(a) Gene model of ENSBTAT00000014437 on BTA20 (ARS-UCD1.2), with the final coding exon
highlighted. (b) The six published truncating alleles (SLICK1–SLICK6) in the final exon, with
HGVS nomenclature and the number of carrier haplotypes found in the phased panel. (c) The
581-residue receptor with the transmembrane segment, the Box1 motif and the distal STAT5
docking tyrosine; bars give the cytoplasmic tail each allele retains. Box1 and the STAT5
tyrosine are placed by homology with the human receptor, not taken from UniProt.
Input: `derived/resultado_final.json`.

### Figure 2 — Validation of the local ancestry inference
(a, b) RFMix internal simulation in the three-way and four-way reference configurations:
haplotypes simulated from each panel (columns) and the panel they are returned to (rows).
(c) Reference animals held out of the panels and run as queries: share of their haplotypes
returned to their own panel at the locus. (d) Positive control: three African breeds of known
taurine × indicine admixture run as queries. (e) Mean window size of each configuration over the
10 Mb interval.
Input: `derived/local_ancestry/validacao.json`, `m3.msp.tsv`, `m4.msp.tsv`, `m5.msp.tsv`,
`consulta.tsv`.

### Figure 3 — Dependence of the 39 kb coding classifier on the exon 2 site pair
(a) Frequency of the alternative allele at chr20:39,077,999 in each reference panel, carriers
excluded, with the number of carrier haplotypes of each slick allele that carry it. (b) Coding
sites ranked by their contribution to the European-over-African log-likelihood ratio, as a
cumulative share of the total. (c) Panel assigned to each allele at non-synonymous sites with
all coding sites and with the exon 2 pair (39,077,999/39,078,000) removed. (d) Share of 2,000
site-resampling replicates that keep the observed continental group.
Input: `derived/exon2_dados.json`, `derived/estabilidade_hierarquica.json`,
`derived/resultado_final.json`.

### Figure 4 — Haplotype background over the 10 Mb window
(a) Local ancestry along the 10 Mb interval, three-way configuration, for every carrier
haplotype and for twelve Criollo haplotypes carrying no slick allele. (b) Differences between
each carrier haplotype and the closest haplotype of each panel within 50 kb of the locus; bars
give the median. (c) Ancestry called at the locus for each carrier under the three reference
configurations.
Input: `derived/local_ancestry/m3.msp.tsv`, `focais_faseados.tsv`, `vizinho.json`,
`consulta.tsv`, `m4.msp.tsv`, `m5.msp.tsv`.
