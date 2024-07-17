# Parameter profiles

Inspired by parameter presets in Minimap2, Squigulator has parameter profiles/presets that can be used via the `-x option`. This provides the user with a convenient way to simulate different chemistries and flowcells, without having to use advanced parameters.

Parameter profiles are currently available for R9.4.1, R10.4.1 amd RNA004. More profiles will be added in future when new nanopore chemistries are introduced.

## R9 profiles

Following profiles are available for R9.4.1.

- `dna-r9-min`: DNA R9.4.1 MinION (or GridION)
- `dna-r9-prom`: DNA R9.4.1 promethION (or P2solo)
- `rna-r9-min`: RNA R9.4.1 MinION (or GridION)
- `rna-r9-prom`:  RNA R9.4.1 promethION (or P2solo)

These profile sets the following advanced parameters in squigulator.

| Parameter            | dna-r9-min   | dna-r9-prom   | rna-r9-min   | rna-r9-prom   |
|----------------------|--------------|---------------|--------------|---------------|
| digitisation         | 8192         | 2048          | 8192         | 2048          |
| sample-rate          | 4000         | 4000          | 3012         | 3000          |
| bps                  | 450          | 450           | 70           | 70            |
| range                | 1443.030273  | 748.5801      | 1126.47      | 548.788269    |
| offset-mean          | 13.7222605   | -237.4102     | 4.65491888   | -231.9440589  |
| Offset-std           | 10.25279688  | 14.1575       | 4.115262472  | 12.87185278   |
| median-before-mean   | 200.815801   | 214.2890337   | 242.6584118  | 238.5286796   |
| median-before-std    | 20.48933762  | 18.0127916    | 10.60230888  | 21.1871794    |
| dwell-mean           | 9.0          | 9.0           | 43.0         | 43.0          |
| dwell-std            | 4.0          | 4.0           | 35.0         | 35.0          |


## R10 and RNA004 profiles

Following profiles are available for R10.4.1 and RNA004.

- `dna-r10-min`:  DNA R10.4.1 MinION (or GridION)
- `dna-r10-prom`: DNA R10.4.1 promethION (or P2solo)
- `rna004-min`: RNA004 MinION (or GridION)
- `rna004-prom`: RNA004 promethION (or P2solo)

These profile sets the following advanced parameters in squigulator.


| Parameter           | dna-r10-min   | dna-r10-prom   | rna004-min   | rna004-prom   |
|---------------------|---------------|-----------------|--------------|---------------|
| digitisation        | 8192          | 2048            | 8192         | 2048          |
| sample-rate         | 4000          | 4000            | 4000         | 4000          |
| bps                 | 400           | 400             | 130          | 130           |
| range               | 1536.598389   | 281.345551      | 1437.976685  | 299.432068    |
| offset-mean         | 13.380569389  | -127.5655735    | 12.47686423863 | -259.421128  |
| offset-std         | 16.311471649  | 19.377283387665 | 10.442126577137 | 16.010841823643 |
| median-before-mean  | 202.154074388 | 189.87607393756 | 205.08496731088 | 189.87607393756 |
| median-before-std   | 13.406139242  | 15.788097978713 | 8.6671292866233 | 15.788097978713 |
| dwell-mean          | 10.0          | 10.0            | 31.0         | 31.0          |
| dwell-std           | 4.0           | 4.0             | 0.0          | 0.0           |

## Determining parameters for a profile

This section briefly describes how these parameter profiles are created for a given flow-cell chemistry. Following requirements and assumptions must be met.

- A high quality pore model must be available. ONT usually provides k-mer models for their chemistries [here](https://github.com/nanoporetech/kmer_models). For R9.4.1 both the level means and the standard deviation values were provided, however, for R10.4.1 and RNA004 only level means have been provided. The method described [here](https://hasindu2008.github.io/f5c/docs/r10train) can be used for deducing the standard deviation values.

- You need some sample signal data in BLOW5 format and [slow5tools](https://github.com/hasindu2008/slow5tools) set up. If the data is in FAST5, use [slow5tools f2s](https://github.com/hasindu2008/slow5tools) to convert. If data is in POD5, use [blue-crab](https://github.com/Psy-Fer/blue-crab).

- Install `datamash` (e.g., apt-get install datamash)


Now let us see how the the parameters are deduced.

- `digitisation`. This is the digitisation field in the BLOW5 file as explained [here](https://hasindu2008.github.io/slow5specs/summary). The digitisation value must be the same across the whole dataset. Infact, so far as far as we know, MinIONs/GridIONs always has a digitisation of 8192 and promethION/P2 has 2048. You can deduce the digitisation using the dataset by using the example command below (you should only see one value):

  ```
  slow5tools skim -t40 reads_500k.blow5 | cut -f 3 | tail -n+2 | sort -u
  2048
  ```

- `sample-rate`. This is the sample-rate field in the BLOW5 file as explained [here](https://hasindu2008.github.io/slow5specs/summary). This value must be the same across the whole dataset. To deduce this value, you can use the following example command (you should only see one value):

  ```
  slow5tools skim -t40 reads_500k.blow5 | cut -f 6 | tail -n+2 | sort -u
  4000
  ```

- `bps`. This is the translocation speed which can be found on the relevant Guppy/Dorado basecalling model. For example, the Dorado model for rna004 is `rna004_130bps_sup@v3.0.1`. The bps is 130.

- `range`. This is the range field in the BLOW5 file as explained [here](https://hasindu2008.github.io/slow5specs/summary). This value also must be the same across the whole dataset. Just as before,   example command (you should only see one value):

  ```
  slow5tools skim -t40 reads_500k.blow5 | cut -f 5 | tail -n+2 |  sort -u
  299.432068
  ```

- `offset-mean` and `offset-std`. This is the mean and standard deviation of the [offset field in the BLOW5 file](https://hasindu2008.github.io/slow5specs/summary).  Example command to get the two parameters using datamash:

  ```
   slow5tools skim -t40 500k.blow5 | cut -f 4 | tail -n+2 |  datamash mean 1 sstdev 1
   -259.421128     16.010841823643
  ```

- `median-before-mean` and `median-before-std`. This is the mean and standard deviation of the [median_before field in the BLOW5 file](https://hasindu2008.github.io/slow5specs/summary). Example command to get the two parameters:

 ```
  slow5tools skim -t40 reads_500k.blow5 | awk -v c="median_before" 'NR==1{for (i=1; i<=NF; i++) if ($i==c){p=i; break};} {if ($p!=".") print $p}' | tail -n+2 | datamash mean 1 sstdev 1
  205.63935594369 8.3994882799157
 ```

- `dwell-mean`. This is the mean of the number of signal samples per base. This must be equal to the `sample_rate`/`bps`, and is just used as a sanity check.

- `dwell-std`. This is the standard deviation of number of signal samples per base. Setting this to 0 gives the highest basecalling accuracy. Increasing the value will reduce the accuracy.
The best value for this parameter can be empirically determined by basecalling the data simulated with different values for this parameter and choosing the value that gives the closest accuracy to basecalls of the real data.
