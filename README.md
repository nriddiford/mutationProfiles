# mutationProfiles

## Installation

Install from github:

```
git clone https://github.com/nriddiford/mutationProfiles.git

cd mutationProfiles
```
## Dependencies

`trinucs.pl` requires [BioPerl](http://bioperl.org/INSTALL.html), which can be installed from cpanm:

```
brew install cpanm
cpanm Bio::Perl
```

## Extracting SNV calls from Mutect2 or Freebayes vcf files or Varscan2 native format:

Move all `.vcf` files into `data/` and run `bash run_trinucs.sh -g <path to genome.fasta>`  
For Varscan native data run: `bash run_trinucs.sh -v -g <path to genome.fasta>`

This will run `script/trinucs.pl` for each `.vcf` file in `data/`, and write data from all samples to `data/GW.snv.dist.txt`

e.g.:
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
cleanTheme : function (base_size = 12)  
genomeSnvs : function ()  
getData : function (infile = "data/GW.snv.dist.txt")  
mutSigs : function (samples = NA, pie = NA)  
mutSpectrum : function ()  
notchSnvs : function ()  
genTris : function ()  
samplesPlot : function (count = NA)  
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
