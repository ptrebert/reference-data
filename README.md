# refdata
Collection of scripts and converters assembled as pipeline to process annotation data

# Portal links
- ENCODE annotations:   https://www.encodeproject.org/data/annotations

## Glossary

### ENCODE promoters (v3)

- based on H3K4me3 / DNaseI
- human: 107 cell types
- mouse: 14 cell types
- proximal: (DNase peak less than 2kb away from annotated TSS) *AND* (peak in top 10000 peaks,
 ranking based on all nearby expressed transcripts)
- distal: higher ranked according to the above schema, but not TSS-proximal; could be unannotated TSS or transcribed enhancers

### ENCODE enhancers (v3)

- based on H3K27ac / DNaseI
- human: 47 cell types
- mouse: 14 cell types
- distal: (region more than 2kb away from any TSS) *AND* (region in top 20000 based on (i) anchoring with DNase peak and (ii) ranking
 by H3K27ac signal); notably, a region in this set can contain several DNase peaks since the boundaries are defined based on H3K27ac peaks
- proximal: higher ranked according to the above schema, but too close to a TSS; could be promoter or enhancer with promoter-like activty


# Download log
- 2016-02-10	http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit
- 2016-02-10	http://hgdownload.soe.ucsc.edu/goldenPath/mm9/bigZips/mm9.2bit
- 2016-02-10	http://hgdownload.soe.ucsc.edu/goldenPath/bosTau7/bigZips/bosTau7.2bit
- 2016-02-10	http://hgdownload.soe.ucsc.edu/goldenPath/mm9/liftOver/mm9ToBosTau7.over.chain.gz
- 2016-02-10	http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToBosTau7.over.chain.gz
- 2016-02-10	http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToPanTro3.over.chain.gz
- 2016-02-10	http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToMm9.over.chain.gz
- 2016-02-10	http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToDanRer7.over.chain.gz
- 2016-02-10	http://hgdownload.soe.ucsc.edu/goldenPath/mm9/liftOver/mm9ToDanRer7.over.chain.gz
- 2016-02-10	http://hgdownload.soe.ucsc.edu/goldenPath/mm9/liftOver/mm9ToPanTro3.over.chain.gz
- 2016-02-10	http://hgdownload.soe.ucsc.edu/goldenPath/mm9/liftOver/mm9ToHg19.over.chain.gz
- 2016-02-10	http://hgdownload.soe.ucsc.edu/goldenPath/panTro3/bigZips/panTro3.2bit
- 2016-02-15	http://hgdownload.soe.ucsc.edu/goldenPath/susScr2/bigZips/susScr2.2bit
- 2016-04-04	http://hgdownload.soe.ucsc.edu/goldenPath/canFam3/bigZips/canFam3.2bit
- 2016-05-10	http://fantom.gsc.riken.jp/5/datafiles/latest/extra/Enhancers/human_permissive_enhancers_phase_1_and_2.bed.gz
- 2016-05-10	http://fantom.gsc.riken.jp/5/datafiles/latest/extra/Enhancers/mouse_permissive_enhancers_phase_1_and_2.bed.gz
- 2016-05-24	ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M1/gencode.vM1.annotation.gtf.gz
- 2016-05-24	ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz
- 2016-06-06	ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M1/gencode.vM1.pc_transcripts.fa.gz
- 2016-06-06	ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.pc_transcripts.fa.gz
- 2016-07-07	http://128.206.12.216/drupal/sites/bovinegenome.org/files/data/btau_4.6.1/Ensembl75_liftOver_Btau_4.6.1_genes.gff3.gz
- 2016-07-11	http://hgdownload.soe.ucsc.edu/goldenPath/canFam3/database/ensGene.txt.gz
- 2016-07-11	http://hgdownload.soe.ucsc.edu/goldenPath/canFam3/database/ensemblSource.txt.gz
- 2016-07-11	http://hgdownload.soe.ucsc.edu/goldenPath/canFam3/database/ensemblToGeneName.txt.gz
- 2016-07-11	http://hgdownload.soe.ucsc.edu/goldenPath/susScr2/database/ensemblSource.txt.gz
- 2016-07-11	http://hgdownload.soe.ucsc.edu/goldenPath/susScr2/database/ensGene.txt.gz
- 2016-07-11	http://hgdownload.soe.ucsc.edu/goldenPath/susScr2/database/ensemblToGeneName.txt.gz
- 2016-08-23	https://groups.google.com/a/soe.ucsc.edu/d/msg/genome/diLb_BKSr9Q/SchnwWBa7loJ [track version - info about Ensembl version]
- 2016-09-15	http://enhancer.lbl.gov/frnt_page_n.shtml [no stable URLs for VISTA dataset available]
- 2016-09-15	http://bioinfo.au.tsinghua.edu.cn/dbsuper/data/bed/mm9/all_mm9_bed.bed [super enhancer mm9 - doi: 10.1093/nar/gkv1002 ]
- 2016-09-15	http://bioinfo.au.tsinghua.edu.cn/dbsuper/data/bed/hg19/all_hg19_bed.bed [super enhancer hg19 - doi: 10.1093/nar/gkv1002 ]
- 2016-09-15	http://genome.ucsc.edu/cgi-bin/hgTables [UCSC table browser; all CpG islands]
- 2016-10-03    http://genome.ucsc.edu/cgi-bin/hgTables [UCSC table browser; all assembly gaps]
- 2016-10-14    ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
- 2016-11-30    http://zlab-annotations.umassmed.edu/promoters/20161130-062049-promoter-like-hg19-BothDNaseAndH3K4me3.v3.zip
- 2016-11-30    http://zlab-annotations.umassmed.edu/promoters/20161130-062104-promoter-like-mm10-BothDNaseAndH3K4me3.v3.zip
- 2016-11-30    http://zlab-annotations.umassmed.edu/enhancers/20161130-062139-enhancer-like-hg19-BothDNaseAndH3K27ac.v3.zip
- 2016-11-30    http://zlab-annotations.umassmed.edu/enhancers/20161130-062201-enhancer-like-mm10-BothDNaseAndH3K27ac.v3.zip
- 2016-12-02    http://genome.ucsc.edu/cgi-bin/hgTables [UCSC table browser; phastCons conserved elements] hg19_phastConsElem_vert_46way.tsv.gz
- 2016-12-02    http://genome.ucsc.edu/cgi-bin/hgTables [UCSC table browser; phastCons conserved elements] mm9_phastConsElem_vert_30way.tsv.gz
- 2016-12-13    http://www.bio-bigdata.com/SEA/download/SEA00201.bed [super enhancer mm9 - doi: 10.1093/nar/gkv1243 ]
- 2017-01-04    ftp://ftp.ensembl.org/pub/grch37/release-87/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna_sm.primary_assembly.fa.gz
- 2017-01-17    ftp://ftp.ebi.ac.uk/pub/databases/genenames/hcop/human_* [human orthologs for mouse, pig, cow etc]
- 2017-01-27    http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/hg19-human/wgEncodeHg19ConsensusSignalArtifactRegions.bed.gz
- 2017-01-27    http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/mm10-mouse/mm10.blacklist.bed.gz
- 2017-01-27    http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/mm9-mouse/mm9-blacklist.bed.gz
- 2017-01-27    http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/hg38-human/hg38.blacklist.bed.gz
- 2017-02-28    http://www.orthomam.univ-montp2.fr/orthomam/html/ [selection: OrthoMaM v9: CDS; homo, mus, sus, canis, bos: 10010 results]
- 2017-03-03    http://www.orthodb.org/v9/download/odb9_OGs.tab.gz [http://www.orthodb.org/v9/ - PMID:25428351]
- 2017-03-03    http://www.orthodb.org/v9/download/odb9_OG2genes.tab.gz [http://www.orthodb.org/v9/ - PMID:25428351]
- 2017-03-03    http://www.orthodb.org/v9/download/odb9_OGs.tab.gz [http://www.orthodb.org/v9/ - PMID:25428351]


## OrthoDB flat files
odb9_OGs.tab.gz - Ortho DB orthologous groups

1.	OG unique id (not stable between releases)
2.	level tax_id on which the cluster was built
3.	OG name (the group's most common gene name)

odb9_OG2genes.tab.gz - OGs to genes correspondence

1.	OG unique id
2.	Ortho DB unique gene id

odb9_genes.tab.gz - genes with some info

odb9_genes.tab
1.	Ortho DB unique gene id (not stable between releases)
2.	organism tax id
3.	protein sequence id, as downloaded together with the sequence
4.	Uniprot id, evaluated by mapping
5.	ENSEMBL gene name, evaluated by mapping
6.	NCBI gid, evaluated by mapping
7.	description, evaluated by mapping


## UCSC Chainfiles

### Chains and Nets
Following the descriptions in the [GenomeWiki](http://genomewiki.ucsc.edu/index.php/Chains_Nets) *(last modified status: 16 April 2015, 19:10)*,
the process to generate reciprocal best chains/nets generates the following output:

> The reciprocal best process uses that output:
> the query-referenced (but target-centric / target single-cov) net is turned back into component chains,
> and then those are netted to get single coverage in the query too;
> the two outputs of that netting are reciprocal-best in query and target coords.
> Reciprocal-best nets are symmetrical again.

Additionally, the following naming convention is used by UCSC for pre-computed chain/net files:

> In chain and net lingo, the target is the reference genome sequence and the query is some other genome sequence.
> For example, if you are viewing Human-Mouse alignments in the Human genome browser, human is the target and mouse is the query.

> The file names reflect the assembly conversion data contained within
> in the format <db1>To<Db2>.over.chain.gz. For example, a file named
> hg15ToHg16.over.chain.gz file contains the liftOver data needed to
> convert hg15 (Human Build 33) coordinates to hg16 (Human Build 34).

> a net is single-coverage for target but not for query,
> unless it has been filtered to be single-coverage on both target and query.
> By convention we add "rbest" to the net filename in that case.

The generic process to generate reciprocal best/symmetric single coverage files is also documented in the
[Genome Wiki](http://genomewiki.ucsc.edu/index.php/HowTo:_Syntenic_Net_or_Reciprocal_Best) *(last modified status: 12 January 2016, 19:00)*
The documentation is incomplete or imprecise or just not error-free - see following paragraph.

The *error* report on the non-symmetric output can be found in this [Google groups thread](https://groups.google.com/a/soe.ucsc.edu/forum/#!topic/genome/vwgVtaXvyug)
and includes the following statement by the UCSC support:

> Chaining and netting are not simple operations. They may not be symmetrical operations,
> there may be some slight difference in each direction.
> I tried taking these results and running them around in another cycle to get comparable bed files
> in the same coordinate system to see what might be missing, but this led to even more missing bases.
> Evidently the cycle itself does something to cause bases to go missing.

Notably, the tool *chainNet* has to be executed with non-default parameters. The defaults are:

1) *-minSpace=N - minimum gap size to fill, default 25*

2) *-minScore=N - minimum chain score to consider, default 2000.0*

In the script provided by the UCSC support, these values are set to:

1) *-minSpace=1*

2) *-minScore=0*

Setting these parameters to non-default values dramatically increases the number of missing bases,
for unclear reasons.


### Chain scoring

Apparently, there is no precise/formal description of how chain scores are derived (truth may be in the code).
In this [UCSC google groups thread](https://groups.google.com/a/soe.ucsc.edu/d/msg/genome/C7Wxn01IzQk/wHpOkr9mWAoJ), Angie (UCSC) states that

> In a nutshell, the chain scoring scheme is somewhat complicated,
> but ballpark estimates of scores can be made from expected size of aligned blocks, 
> percent identity, gap size and gap frequency. [...]
> Gaps in chains are penalized with a piecewise linear function that penalizes gap openings the most, 
> with less harsh penalties as gaps are larger 
> (we expect large gaps in cross-species alignments due to insertions, rearrangements etc.) [...]
> Another approach to making sense of chain scores is to look at chain score histograms, 
> e.g. for all chains or for chains of a particular length and gap size extracted using chainFilter.

Besides (manual) score filtering using the above mentioned histogram method, in this
[UCSC google groups thread](https://groups.google.com/a/soe.ucsc.edu/d/msg/genome/zcHuWtmw-LE/5hClwQ4furIJ), Rachel (UCSC) mentions that

> Regarding filtering, if there are a lot of low scoring chains e.g. those
> with score < 5000 then often we filter these out using the minScore option
> for axtChain. Chains can also be filtered after they are made using the
> chainFilter program.  chainPreNet does remove chains that do not have a
> chance of being netted. Then chainNet makes the alignments nets from
> chains using the highest scoring chains in the top level. Gaps are
> filled in with other chains at level 2 and then gaps in the
> level 2 chains can be filled in with chains in level 3 etc. In the net,
> chains are trimmed to fit into these sections that are not covered by a
> higher-scoring chain. We also use netFilter with the minGap option set to
> 12 before loading the net into the database. This restricts the nets to
> those with a gap size >= 12 bp.

However, there seems not to be an overall good strategy as stated by Angie (UCSC) in this
[UCSC google groups thread](https://groups.google.com/a/soe.ucsc.edu/d/msg/genome/aoOw3z7u0Qk/cMzhkAdjeuwJ)
 
> Picking a score threshold for chains is a tricky business... scores
> vary hugely with length as well as conservation.  This scoring scheme
> allows us to recognize long chains in syntenic regions, but it also
> retains almost anything from blastz.  That's why we also have the
> "net" tracks -- to keep the best chains and ignore most of the
> "fluff".