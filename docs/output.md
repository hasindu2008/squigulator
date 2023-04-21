# Outputs

1. SLOW5/BLOW5 file (-o) containing the simulated raw signal data
2. FASTA file (-q) containing the perfect simulated reads directly from the reference. Note that these reads do not contain sequencing errors. You must basecall the S/BLOW5 file to produce a FASTQ file with sequencing errors.
3. PAF FILE (-c) containing signal to read/reference alignment information
3. SAM FILE (-a) containing signal to reference alignment information

## PAF output

When -c is specified, a PAF output that contains signal to read alignment is generated. Please refer [here](https://hasindu2008.github.io/f5c/docs/output#resquiggle-paf-output-format) for a detailed explanation of the specification, along with examples. Note that in Squigulator's context the target sequence (column 6 to 9 in PAF) refers to the perfect sequences output through -q option (rather than basecalled reads in [here](https://hasindu2008.github.io/f5c/docs/output#resquiggle-paf-output-format)). So when reading that specification, assume that the basecalled sequences refers to output generated through -q.

When --paf-ref is specified along with -c, the ouput PAF file will contain signal to reference alignment information (rather than signal to read). The target reference (column 6 to 9 in PAF) will contain the reference sequence information, rather than the read in this case.

## SAM output

When -a is specified, a SAM file that contains signal to reference alignment information is created.

|Col|Type  |Name |Description                               |
|--:|:----:|:----|:-----------------------------------------|
|1  |string|QNAME|Read identifier name                       |
|2  |int   |FLAG|Bitwise flag   (0 if "+" strand and 16 if '-')                 |
|3  |string   |RNAME|Reference sequence name |
|4  |int   |POS|Reference sequence start index for the mapping (1-based; open)       |
|5 |int   |mapq|Mapping quality (always 255)  |
|6 |string|CIGAR| CIGAR string (read to reference)|
|7  |string|RNEXT| always "*" |
|8  |int   |PNEXT| always 0 |
|9  |int   |TLEN| always 0  |
|10|string|SEQ| the simulated perfect read sequence directly extracted from the reference |
|11|string|QUAL| always "*"|

Following optional tags are present:

|Tag|Type  |Description                               |
|--:|:----:|:-----------------------------------------|
|si  |Z   |coordinates associated with the ss tag below (explained below)                     |
|ss  |Z   |signal alignment string in format described [here](https://hasindu2008.github.io/f5c/docs/output#resquiggle-paf-output-format)   |

*si* tag contains four comma separated values *start_raw*, *end_raw*, *start_kmer* and *end_kmer*, respectively. Those values have the same  as the columns 3,4,8 and 9 in the PAF format explained above when --paf-ref is specified along with -c
