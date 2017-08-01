# mutationProfiles

## Installation

Install from github:

```
git clone https://github.com/nriddiford/mutationProfiles.git

cd mutationProfiles
```
## Dependencies

`trinucs.pl` requires [BioPerl](http://bioperl.org/INSTALL.html), which can be installed using cpanm:

```
brew install cpanm
cpanm Bio::Perl
```

## Extracting SNV calls from Mutect2 or Freebayes vcf files or Varscan2 native format

Move all `.vcf` files into `data/` and run `bash run_trinucs.sh -g <path to genome.fasta>`  
For Varscan native data run: `bash run_trinucs.sh -v -g <path to genome.fasta>`

This will run `script/trinucs.pl` on each `.vcf` file in `data/`, and write data from all samples to `data/combined_snvs.txt` in the following format:

```
[sample] [chrom] [pos] [ref] [alt] [tri context] [ref>alt] [decomposed trinuc context] [decomposed ref>alt] [type]
```
```
A512R17	2L	229832	A	C	CAG	A>C	CTG	T>G	somatic
A512R17	2L	1819239	T	A	TTC	T>A	TTC	T>A	somatic
A512R17	2L	2439881	C	T	GCC	C>T	GCC	C>T	somatic
A512R17	2L	3154318	C	G	GCC	C>G	GCC	C>G	somatic
A512R17	2L	3511198	G	A	CGA	G>A	TCG	C>T	somatic
A512R17	2L	4565784	C	G	CCT	C>G	CCT	C>G	somatic
A512R17	2L	5233457	T	G	TTA	T>G	TTA	T>G	somatic
A512R17	2L	6478473	G	C	GGT	G>C	ACC	C>G	somatic
A512R17	2L	9792284	C	T	GCC	C>T	GCC	C>T	somatic
```

## Annotate SNVs with gene and feature it's contained within

Run `perl script/snv2gene.pl -i data/combined_snvs.txt` to annotate the gene and feature hit by each SNV

e.g.:
```
[sample] [chrom] [pos] [ref] [alt] [tri context] [ref>alt] [decomposed trinuc context] [decomposed ref>alt] [type] [feature] [gene]
```

```
A512R17 2L      229832  A       C       CAG     A>C     CTG     T>G     somatic intron  kis
A512R17 2L      1819239 T       A       TTC     T>A     TTC     T>A     somatic intergenic      intergenic
A512R17 2L      2439881 C       T       GCC     C>T     GCC     C>T     somatic intron  dpp
A512R17 2L      3154318 C       G       GCC     C>G     GCC     C>G     somatic intron  Mad
A512R17 2L      3511198 G       A       CGA     G>A     TCG     C>T     somatic exon_5  LeuRS
A512R17 2L      4565784 C       G       CCT     C>G     CCT     C>G     somatic intron  dpy
A512R17 2L      5233457 T       G       TTA     T>G     TTA     T>G     somatic intron  tkv
```


## Explore snv data

Start an R session, and install package:

```{R}
library(devtools)
install_github("nriddiford/mutationProfiles")
library(mutationProfiles)
setwd('mutationProfiles')
```

## mutationProfiles

The following functions are included:

```{R}
chromDist : function (object = NA, notch = 0)  
cleanTheme : function (base_size = 12)  
featuresHit : function ()  
geneHit : function (n = 10)  
genomeSnvs : function ()  
genTris : function ()  
getData : function (infile = "data/annotated_snvs.txt")  
mutSigs : function (samples = NA, pie = NA)  
mutSpectrum : function ()  
notchSnvs : function ()  
samplesPlot : function (count = NA)  
setCols : function (df, col)  
snvStats : function ()  
triFreq : function (genome = NA, count = NA)  
```
### See some stats

```{R}
snvStats()
```

### Plot mutations per sample

```{R}
samplesPlot()
```

### Plot mutation spectrum for all samples combined

```{R}
mutSpectrum()
```

### Plot mutational signatures in data

This plots the output of the package [deconstructSigs](https://github.com/raerose01/deconstructSigs/tree/master/R)

```{R}
mutSigs()
```

### Plot distribution of snvs across chromosomes

```{R}
chromDist()
```

### Plot number of times a feature type has been hit_ref

```{R}
featuresHit()
```

### Show the 20 most hit genes

```{R}
geneHit(n=20)
```
