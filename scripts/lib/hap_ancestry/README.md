# hap_ancestry — motor de ancestralidade de haplótipo por classe funcional

Motor **guarda-chuva**, reutilizável por qualquer projeto de ancestralidade de
haplótipo. Dado um gene/região faseado, ele responde: **de qual população vem o
haplótipo que carrega um alelo focal**, estratificado por classe funcional de
sítio, com intervalo por reamostragem de sítios e uma rede de haplótipos.

O motor **não conhece nenhum locus**. Tudo que é específico — sítios focais,
painéis de referência e sua ordem, cores/rótulos, limiares, mapa de consequências
— vem de um `LocusConfig` carregado de um YAML. O caso *slick*/*PRLR* é apenas um
arquivo `slick.yaml` de exemplo. **Para um projeto novo, escreva um YAML e reuse
este mesmo motor** — sem tocar no código.

## Estrutura

```
scripts/lib/hap_ancestry/
    io.py         load_phased, load_consequences, load_groups, hap_group_labels,
                  active_panels, resolve_classes
    config.py     LocusConfig (dataclass) + from_yaml
    assign.py     atribuição por verossimilhança: panel_frequencies, winners,
                  assign_with_posterior, consensus
    bootstrap.py  stability (reamostragem de SÍTIOS, não de haplótipos)
    network.py    plot_network, collapse, mds_layout (rede de haplótipos)
    __init__.py   API pública
```

Os pontos de entrada são scripts wrappers finos sobre o motor (retrocompatíveis:
rodar sem `--config` usa o caso slick):

| Script | Produz |
|---|---|
| `gene_ancestry_by_class.py` | tabela de atribuição por classe funcional (a saída central do cálculo) |
| `bootstrap_ancestry_class.py` | IC e estabilidade da atribuição |
| `haplotype_network.py` | figura de rede de haplótipos |

## Como rodar (caso slick, default)

```bash
export SLICK_DATA_DIR=/path/to/data   # default: ./data
PY=python
SC=scripts/ancestry_haplotype
VCF=$SLICK_DATA_DIR/haplo/prlr_gene_phased.vcf.gz    # gene faseado, sem filtro de MAF
CSQ=$SLICK_DATA_DIR/haplo/prlr_allmaf_csq.tsv        # bcftools csq (4 colunas)
GRP=$SLICK_DATA_DIR/groups_afr_split.tsv             # amostra<TAB>painel (fornecido)

$PY $SC/gene_ancestry_by_class.py  --vcf $VCF --csq $CSQ --groups $GRP --out-prefix out/slick
$PY $SC/bootstrap_ancestry_class.py --vcf $VCF --csq $CSQ --groups $GRP --out-prefix out/slick
$PY $SC/haplotype_network.py        --vcf $VCF --csq $CSQ --groups $GRP --out out/slick_net.pdf
```

Todos aceitam `--config outro_locus.yaml` para rodar sobre outro gene/projeto.

## Como usar em OUTRO projeto/locus

1. Copie `config/slick.yaml` e ajuste: `focal_sites` (nome+posição na assembly do
   VCF), `panels` (ordem = desempate do argmax), `panel_style` (cores/rótulos da
   rede), `consequence_classes` (mapa do seu anotador), `class_report`, limiares.
2. Rode qualquer wrapper com `--config seu_locus.yaml`, ou importe o motor direto:

```python
import sys; sys.path.insert(0, "scripts/lib")
from hap_ancestry import (LocusConfig, load_phased, load_consequences, load_groups,
                          hap_group_labels, active_panels, resolve_classes,
                          focal_index, carrier_haplotypes, panel_frequencies,
                          assign_with_posterior, stability, plot_network)
cfg = LocusConfig.from_yaml("seu_locus.yaml")
```

## Formato do YAML

Ver `config/slick.yaml` — está inteiramente comentado. Campos: `locus_name`,
`allele_label`, `focal_sites[]`, `panels[]`, `panel_style{}`, `thresholds{}`,
`consequence_classes{}`, `class_report[]`, `coding_classes[]`.

## Notas de reprodutibilidade

- **Bootstrap determinístico.** O consenso desempata pela ordem de *primeira
  aparição* nos dados (`dict.fromkeys`), não pela ordem de um `set` (que o Python
  randomiza por processo via `PYTHONHASHSEED`). Um desempate por ordem de `set`
  faz a `estabilidade` do bootstrap oscilar ~0,03 entre execuções nas classes de
  baixo consenso; este método não. O desempate **não** usa a ordem dos painéis,
  para não enviesar o consenso ao primeiro painel da lista.
- **Fase lida no próprio sítio.** O sítio focal tem de estar no VCF faseado; se não
  estiver, `focal_index` devolve `exato=False` e o wrapper pula (não infere).
- **Painéis excluem os portadores** ao estimar frequência (não se auto-avaliam).

## Proveniência

Motor extraído de scripts de cálculo específicos de um locus (que hoje são
wrappers finos sobre ele). Equivalência validada contra a versão anterior:
atribuição idêntica (diff 0); estabilidades codificantes idênticas a 2 casas
decimais; as de contraste (introns/gene inteiro) passaram de não-reprodutíveis a
fixas — dois valores de intron mudaram de leitura (SLICK1 0.843→0.792, SLICK2
0.585→0.649), sinalizados automaticamente por verificação de integridade sempre
que a análise é reexecutada.
