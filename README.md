# sigsim

*sigsim* (aka *ssssim*) is a Stupidly Simple Signal Simulator for simulating nanopore raw signal data. It is currently in early development stage and therefore anticipate bugs and interface changes. Do not hesitate to open an [issue](https://github.com/hasindu2008/sigsim) if you found a bug, something is not clear or for any feature requests.

*sigsim* uses traditional pore models and gaussian noise for simulation. Due to simplicity, simulation would not be perfect, but takes miniscule effort to setup and run. You may checkout sophisticated deep-learning-based simulators such as [Nanosim](https://github.com/bcgsc/NanoSim) and [DeepSimulator](https://github.com/liyu95/DeepSimulator), if you wish.

Currently, reads directly extracted from the reference genome are simulated without any mutations/variants. If you want to have variants in your simulated data, you can first apply a set of variants to the reference using [bcftools](http://www.htslib.org/download/) and use that as the input to *sigsim*.

*sigsim* at the moment is single threaded and not optimised for efficient runtime performance or RAM usage. Generating 4000 reads from human genome takes like 2-3 minutes with around 3 GB of RAM. Depending on the demand, I will optimise this.

## Building

```
sudo apt-get install zlib1g-dev   #install zlib development libraries
git clone https://github.com/hasindu2008/sigsim
cd sigsim
make
```

The commands to install zlib __development libraries__ on some popular distributions :
```sh
On Debian/Ubuntu : sudo apt-get install zlib1g-dev
On Fedora/CentOS : sudo dnf/yum install zlib-devel
On OS X : brew install zlib
```

## Usage

The simplest command to generate reads:
```
sigsim [OPTIONS] ref_genome.fa -o out_signal.blow5 -n NUM_READS
```

Bu default, DNA PromethION reads will be simulated. Specify the `-x STR` option to set a different profile from the following available pre-sets.
- `dna-r9-min`: genomic DNA on MinION R9.4.1 flowcells
- `dna-r9-prom`: genomic DNA on PromethION R9.4.1 flowcells
- `rna-r9-min`: direct RNA on MinION R9.4.1 flowcells
- `rna-r9-prom`: direct RNA on PromethION R9.4.1 flowcells

If a genomic DNA profile is selected, the input reference must be the reference genome in fasta format. *sigsim* will randomly sample the genome under a uniform distribution and generate reads whose lengths are from a gamma distribution (based on `-r`). If a direct RNA profile is selected, the input reference must be the transcriptome is fasta format. For RNA, *sigsim* will randomly pick transcripts under a uniform distribution and the whole transcript length is simulated.

You can basecall the generated raw signal directly from the [BLOW format](https://www.nature.com/articles/s41587-021-01147-4) using our fork of [dorado basecaller](https://github.com/hiruna72/dorado/releases/tag/v0.0.1) or using the SLOW5 Guppy wrapper called [buttery-eel](https://github.com/Psy-Fer/buttery-eel).  Alternatively, if you love FAST5 that much, use [slow5tools](https://github.com/hasindu2008/slow5tools) to convert the BLOW5 to FAST5 and then use original Guppy basecaller.

Generated read IDs encodes the true mapping positions in a format like `S1_33!chr1!225258409!225267761!-` which is the format compatible with [*mapeval* command in *paftools.js* under Minimap2 repository](https://github.com/lh3/minimap2/blob/master/misc/README.md#evaluation).

Basic options are as below:
- `-o FILE`: SLOW5/BLOW5 file to write.
- `-x STR`: Parameter profile (always applied before other options). Available profiles are: dna-r9-min, dna-r9-prom, rna-r9-min, rna-r9-prom.
- `-n INT`: Number of reads to simulate.
- `-q FILE`: Save the original reads directly taken from the reference genome (without any basecalling errors) in FASTA format. Note that these are perfect reads from the reference and for representative nanopore reads you must basecall the SLOW5/BLOW5 file.
- `-r` median read length (estimated median only, unused for RNA).
- `--ideal`: To generate perfect signals with no noise.

Advanced options are as below:
- `--full-contigs`: generate a complete raw signal per each contig in the input reference genome (incompatible with `-n`, only effective for DNA).
-  `--prefix`: generate prefixes such as adaptor (and polya for RNA).
-  `--kmer-model FILE`: custom nucleotide k-mer model file (format similar to [here](https://github.com/hasindu2008/f5c/blob/master/test/r9-models/r9.4_450bps.nucleotide.6mer.template.model))
-  `--seed INT`: seed for random generators (if 0, will be autogenerated). Giving the same seed will produce same results.
-  `--ideal-time`: Generate signals with no noise in the time domain. Each k-mer will have the same number of signal samples equal to the mean dwell.
-  `--ideal-amp `: Generate signals with no noise in the amplitude domain. All samples for a given k-mer/base will have same signal values.
-  `--dwell-mean FLOAT`: Mean of number of signal samples per k-mer/base. This is usually the sampling rate (4000Hz for DNA and 3000Hz for RNA) divided by translocation speed in bases per second (450 for R9.4.1 pore for DNA and 70 for RNA).
-  `--dwell-std FLOAT`: Standard deviation of number of signal samples per k-mer/base. Increasing this will increase time-domain noise. Setting this to 0 is same as `--ideal-time`.

## Examples

DNA examples:

```
# generate 150,000 PromethION DNA reads from a reference genome
sigsim hg38noAlt.fa -x dna-r9-prom -o reads.blow5 -n 150000

# generate 30,000 MinION ultra-long DNA reads with mean readlength of around 50,000 bases
sigsim hg38noAlt.fa -x dna-r9-min -o reads.blow5 -n 30000 -r 50000

# generate 1000 PromethION DNA reads with perfect signals with no noise
sigsim hg38noAlt.fa -x dna-r9-prom -o reads.blow5 -n 1000 --ideal

```

RNA examples:
```
# generate 4000 PromethION direct RNA reads from a transcriptome while including the adaptor and polyA tail
sigsim gencode.v40.transcripts.fa -x rna-r9-prom -o reads.blow5 -n 4000 --prefix
```

DNA example with variants that requires [bcftools](http://www.htslib.org/download/):

```
# ploidy 1; coronavirus (reference ~30,000 bases) at ~500X depth with mean readlength of around 300 bases (approximately 30,000*500/300=50,000 reads); apply some variants
bcftools consensus -f nCoV-2019.reference.fasta alpha.vcf -o alpha.fa
sigsim alpha.fa -x dna-r9-prom -o reads.blow5 -n 50000 -r 300

# ploidy 2; chr22 (reference ~50,000,000 bases) at ~30X coverage with mean readlength of around 10,000 bases (approximately 50,000,000*30/10,000=150,000 reads); apply na12878 truthset from genome in a bottle consortium

bcftools consensus -H 1 -f hg38noAlt_chr22.fa na12878_chr22.vcf.gz -o na12878_chr22_1.fa
bcftools consensus -H 2 -f hg38noAlt_chr22.fa na12878_chr22.vcf.gz -o na12878_chr22_2.fa
cat na12878_chr22_1.fa na12878_chr22_2.fa > na12878_chr22.fa
sigsim na12878_chr22.fa -x dna-r9-prom -o reads.blow5 -n 150000 -r 10000
```

## Acknowledgement

The pore-models are from [Nanopolish](https://github.com/jts/nanopolish).
Some code snippets have been taken from [Minimap2](https://github.com/lh3/minimap2), [Samtools](http://samtools.sourceforge.net/).
Kseq from [klib](https://github.com/attractivechaos/klib) is used.

