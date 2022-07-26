# sigsim

*sigsim* (aka *ssssim*) is a Stupidly Simple Signal Simulator for nanopore data. Uses traditional pore models and gaussian noise for simulation. Due to simplicity, simulation would not be perfect, but takes miniscule effort to setup and run. You may checkout sophisticated deep-learning-based simulators such as [Nanosim](https://github.com/bcgsc/NanoSim) and [DeepSimulator](https://github.com/liyu95/DeepSimulator) if you wish. Currently, reads directly extracted from the reference genome are simulated without any mutations/variants.

## Building

```
sudo apt-get install zlib1g-dev   #install zlib development libraries
git clone https://github.com/hasindu2008/sigsim
cd sigsim
make
```

The commands to zlib __development libraries__ on some popular distributions :
```sh
On Debian/Ubuntu : sudo apt-get install zlib1g-dev
On Fedora/CentOS : sudo dnf/yum install zlib-devel
On OS X : brew install zlib
```

## Usage

To generate reads:
```
sigsim [OPTIONS] ref_genome.fa -o out_signal.blow5 -n NUM_READS
```

You can basecall the generated raw signal directly from the BLOW format using our fork of [dorado basecaller](https://github.com/hiruna72/dorado/releases/tag/v0.0.1). Alternatively, use [slow5tools](https://github.com/hasindu2008/slow5tools) to convert the BLOW5 to FAST5 and then use Guppy basecaller.

Generated read IDs encodes the true mapping positions in a format like `S1_33!chr1!225258409!225267761!-` which is the format compatible with [*mapeval* command in *paftools.js* under Minimap2 repository](https://github.com/lh3/minimap2/blob/master/misc/README.md#evaluation).


You may use `-q` option in *sigsim* if you want to save the original reads directly taken from the reference genome (without basecalling errors) in FASTA format. The `--ideal` flag can be used to generate perfect signals with no noise.  `--full-contigs` option can be used for generating a complete raw signal per each contig in the input reference genome (incompatible with `-n`). `-r` option can be used for setting the median read length (estimated median only).

## Acknowledgement

The pore-models are from [Nanopolish](https://github.com/jts/nanopolish).
Some code snippets have been taken from [Minimap2](https://github.com/lh3/minimap2), [Samtools](http://samtools.sourceforge.net/).
Kseq from [klib](https://github.com/attractivechaos/klib) is used.