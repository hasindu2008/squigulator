# slow5lib

*slow5lib* is a software library for reading & writing SLOW5 files. *slow5lib* is designed to facilitate use of data in SLOW5 format by third-party software packages. Existing packages that read/write data in FAST5 format can be easily modified to support SLOW5.

**About SLOW5 format:**<br/>
SLOW5 is a new file format for storing signal data from Oxford Nanopore Technologies (ONT) devices. SLOW5 was developed to overcome inherent limitations in the standard FAST5 signal data format that prevent efficient, scalable analysis and cause many headaches for developers. SLOW5 can be encoded in human-readable ASCII format, or a more compact and efficient binary format (BLOW5) - this is analogous to the seminal SAM/BAM format for storing DNA sequence alignments. The BLOW5 binary format supports  *zlib* (DEFLATE) compression, or other compression methods (see [notes](https://github.com/hasindu2008/slow5lib#notes)), thereby minimising the data storage footprint while still permitting efficient parallel access. Detailed benchmarking experiments have shown that SLOW5 format is an order of magnitude faster and significantly smaller than FAST5.


Full documentation: https://hasindu2008.github.io/slow5lib<br/>
Pre-print: https://www.biorxiv.org/content/10.1101/2021.06.29.450255v1<br/>
Publication: https://www.nature.com/articles/s41587-021-01147-4<br/>
SLOW5 specification: https://hasindu2008.github.io/slow5specs<br/>

[![BioConda Install](https://img.shields.io/conda/dn/bioconda/pyslow5.svg?style=flag&label=BioConda%20install)](https://anaconda.org/bioconda/pyslow5)
[![PyPI](https://img.shields.io/pypi/v/pyslow5.svg?style=flat)](https://pypi.python.org/pypi/pyslow5)
[![C/C++ CI](https://github.com/hasindu2008/slow5lib/actions/workflows/c-cpp.yml/badge.svg)](https://github.com/hasindu2008/slow5lib/actions/workflows/c-cpp.yml)
[![Python CI](https://github.com/hasindu2008/slow5lib/actions/workflows/python.yml/badge.svg)](https://github.com/hasindu2008/slow5lib/actions/workflows/python.yml)


## Building

To build the C/C++ library :

```sh
sudo apt-get install zlib1g-dev   #install zlib development libraries
git clone https://github.com/hasindu2008/slow5lib
cd slow5lib
make
```

This will generate *lib/libslow5.a* for static linking and *libslow5.so* for dynamic linking.

The commands to zlib __development libraries__ on some popular distributions :
```sh
On Debian/Ubuntu : sudo apt-get install zlib1g-dev
On Fedora/CentOS : sudo dnf/yum install zlib-devel
On OS X : brew install zlib
```


#### Optional zstd compression

You can optionally enable [*zstd* compression](https://facebook.github.io/zstd) support when building *slow5lib* by invoking `make zstd=1`. This requires __zstd 1.3 or higher development libraries__ installed on your system:

```sh
On Debian/Ubuntu : sudo apt-get libzstd1-dev # libzstd-dev on newer distributions if libzstd1-dev is unavailable
On Fedora/CentOS : sudo yum libzstd-devel
On OS X : brew install zstd
```

SLOW5 files compressed with *zstd* offer smaller file size and better performance compared to the default *zlib*. However, *zlib* runtime library is available by default on almost all distributions unlike *zstd* and thus files compressed with *zlib* will be more 'portable'.

#### Without SIMD

*slow5lib* from version 0.3.0 onwards uses code from [StreamVByte](https://github.com/lemire/streamvbyte) and by default requires vector instructions (SSSE3 or higher for Intel/AMD and neon for ARM). If your processor is an ancient processor with no such vector instructions, invoke make as `make no_simd=1`.


## Usage

Simply include `<slow5/slow5.h>` in your C program and call the API functions. To compile your program and statically link against slow5lib:

```
gcc [OPTIONS] -I path/to/slow5lib/include your_program.c path/to/slow5lib/lib/libslow5.a -lm -lz
```
*path/to/slow5lib/* is the absolute or relative path to the *slow5lib* repository cloned above. To dynamically link:
```
gcc [OPTIONS] -I path/to/slow5lib/include your_program.c -L path/to/slow5lib/lib/ -lslow5 -lm -lz
```

If you compiled *slow5lib* with *zstd* support enabled, make sure you append `-lzstd` to the above two commands.


For the documentation of the C API visit [here](https://hasindu2008.github.io/slow5lib/slow5_api/slow5.html) and for the Python API visit [here](https://hasindu2008.github.io/slow5lib/pyslow5_api/pyslow5.html).


### Examples

Examples are provided under [examples](https://github.com/hasindu2008/slow5lib/tree/master/examples).
- *sequential_read.c* demonstrates how to read a slow5/blow5 file, sequentially from start to end.
- *random_read.c* demonstrates how to fetch a given read ID from a slow5/blow5 file.
- *header_attribute.c* demonstrates how to fetch a header data attribute from a slow5/blow5 file.
- *auxiliary_field.c* demonstrates how to fetch a auxiliary field from a slow5/blow5 file.
- *random_read_pthreads.c* demonstrates how to fetch given read IDs in parallel from a slow5/blow5 file using *pthreads*.
- *random_read_openmp.c* demonstrates how to fetch given read IDs in parallel from a slow5/blow5 file using openMP.

Following examples will be added upon request. If you are interested, open a GitHub issue rather than spending your time trying to figure out - will try to provide within a day or two.
- using multiple threads for parsing and decompression when performing sequential reading.
- accessing header attributes when multiple read groups are present.
- an efficient input, processing, output pipeline for processing raw signals using slow5lib
- any other use case.

You can invoke `examples/build.sh` to compile the example programmes. Have a look at the script to see the commands used for compiling and linking. If you compiled *slow5lib* with *zstd* support enabled, make sure you append `-lzstd` to the compilation commands.


### pyslow5

Python wrapper for slow5lib or *pyslow5* can be installed using conda as `conda install pyslow5 -c bioconda -c conda-forge` or pypi as `pip install pyslow5`.
The instructions to build *pyslow5* and the usage instructions are [here](https://hasindu2008.github.io/slow5lib/pyslow5_api/pyslow5.html).


### Current limitations & future work

slow5lib is a reference implementation for SLOW5 format. Depending on the interest from the community, the following limitations could be overcome and more performance optimisations can be performed. Open a GitHub issue if you are interested. Contributions are welcome.

- No native windows support: slow5lib works well on Windows through WSL, in fact, this is my primary development environment. I am not aware of anyone using native Windows for nanopore bioinformatics. However, if needed, [methods used for minimap2](https://github.com/lh3/minimap2/issues/19) can be adopted.
- Does not support big-endian systems: Big-endian systems are rare nowadays and I do not have access to one to test. If necessary, it is a matter of writing a layer that swaps the bytes before/after writing to disk. Note: Not to be confused with big.LITTLE architecture which is something else on which slow5lib already works.
- When running with >64 threads, malloc() calls could reduce the thread efficiency. If that is the case, frequent mallocs could be replaced with kalloc in [klib](https://github.com/attractivechaos/klib).
- Aggressive compiler optimisations (e.g.,  -O3) and architecture-specific compiler optimisations (e.g., -march=native) are not used in the makefile. These flags will improve performance at the cost of limited portability. These could be provided in a separate make target.
- multi-threaded decompression and parsing could be implemented inside the C library, but for efficiency and portability it is preferred to be done at user-level. For python bindings in-built multi-threading is already implemented.


### Notes

*slow5lib* from version 0.3.0 onwards has built in [StreamVByte](https://github.com/lemire/streamvbyte) compression support to enable even smaller file sizes, which is applied to the raw signal by default when producing BLOW5 files.  *zlib* compression is then applied by default to each SLOW5 record. If *zstd* is used instead of *zlib* on top of *StreamVByte*, it is similar to ONT's latest [vbz](https://github.com/nanoporetech/vbz_compression) compression. BLOW5 files with *zstd+StreamVByte* are still about 25% smaller than vbz compressed FAST5 files.


## Acknowledgement
slow5lib uses [klib](https://github.com/attractivechaos/klib) and [StreamVByte](https://github.com/lemire/streamvbyte). Some code snippets have been taken from [Minimap2](https://github.com/lh3/minimap2) and [Samtools](http://samtools.sourceforge.net/).
