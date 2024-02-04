# Parameter profiles

-x option. Sets the parameters below.

## R9 profiles

- dna-r9-min: DNA R9.4.1 MinION (or GridION)
- dna-r9-prom: DNA R9.4.1 promethION (or P2solo)
- rna-r9-min: RNA R9.4.1 MinION (or GridION)
- rna-r9-prom:  RNA R9.4.1 promethION (or P2solo)

| Profile              | dna-r9-min   | dna-r9-prom   | rna-r9-min   | rna-r9-prom   |
|----------------------|--------------|---------------|--------------|---------------|
| digitisation         | 8192         | 2048          | 8192         | 2048          |
| sample-rate          | 4000         | 4000          | 3012         | 3000          |
| bps                  | 450          | 450           | 70           | 70            |
| range                | 1443.030273  | 748.5801      | 1126.47      | 548.788269    |
| offset-mean          | 13.7222605   | -237.4102     | 4.65491888   | -231.9440589  |
| Offset STD           | 10.25279688  | 14.1575       | 4.115262472  | 12.87185278   |
| Median Before Mean   | 200.815801   | 214.2890337   | 242.6584118  | 238.5286796   |
| Median Before STD    | 20.48933762  | 18.0127916    | 10.60230888  | 21.1871794    |
| Dwell Mean           | 9.0          | 9.0           | 43.0         | 43.0          |
| Dwell STD            | 4.0          | 4.0           | 35.0         | 35.0          |


## R10 and RNA004 profiles

- dna-r10-min:  DNA R10.4.1 MinION (or GridION)
- dna-r10-prom: DNA R10.4.1 promethION (or P2solo)
- rna004-min: RNA004 MinION (or GridION)
- rna004-prom: RNA004 promethION (or P2solo)

| Profile             | dna-r10-min   | dna-r10-prom   | rna004-min   | rna004-prom   |
|---------------------|---------------|-----------------|--------------|---------------|
| digitisation        | 8192          | 2048            | 2048         | 2048          |
| sample-rate         | 4000          | 4000            | 4000         | 4000          |
| bps                 | 400           | 400             | 130          | 130           |
| range               | 1536.598389   | 281.345551      | TBD  | 299.432068    |
| offset-mean         | 13.380569389  | -127.5655735    | TBD5 | -259.421128  |
| offset-std         | 16.311471649  | 19.377283387665 | TBD | 16.010841823643 |
| median-before-mean  | 202.154074388 | 189.87607393756 | TBD | 189.87607393756 |
| median-before-std   | 13.406139242  | 15.788097978713 | TBD | 15.788097978713 |
| dwell-mean          | 10.0          | 10.0            | TBD         | 31.0          |
| dwell-std           | 4.0           | 4.0             | TBD          | 4.0           |

## Determining parameters for a profile

Assume S/BLOW5. Need slow5tools and datamash.
Convert using X and Y. Following methods:
Assume you have a pore model.


- digitisation. This is the [digitisation field in the BLOW5 file](https://hasindu2008.github.io/slow5specs/summary). Observe that this is the same across the whole dataset (Infact, MinIONs/GridIONs so far has 8192 and promethION/P2 has 2048).

  Example command (you should only see one value):
  ```
  slow5tools skim -t40 PNXRXX240011_reads_500k.blow5 | cut -f 3 | tail -n+2 | sort -u
  2048
  ```

- sample-rate. This is the [sample-rate field in the BLOW5 file](https://hasindu2008.github.io/slow5specs/summary). Observe that this is the same across the whole dataset.

  Example command (you should only see one value):
  ```
  slow5tools skim -t40 PNXRXX240011_reads_500k.blow5 | cut -f 6 | tail -n+2 | sort -u
  4000
  ```

- bps. This is the translocation speed which can be found on the relevant Guppy/Dorado model. For example, the Dorado model for rna004 is `rna004_130bps_sup@v3.0.1`. The bps is 130.

- range. This is the [digitisation field in the BLOW5 file](https://hasindu2008.github.io/slow5specs/summary).

  Example command (you should only see one value):
  ```
  slow5tools skim -t40 PNXRXX240011_reads_500k.blow5 | cut -f 5 | tail -n+2 |  sort -u
  299.432068
  ```

- offset-mean and offset-std. This is the mean and standard deviation of the [offset field in the BLOW5 file](https://hasindu2008.github.io/slow5specs/summary).

  Example command to get the two parameters:
  ```
   slow5tools skim -t40 PNXRXX240011_reads_500k.blow5 | cut -f 4 | tail -n+2 |  datamash mean 1 sstdev 1
   -259.421128     16.010841823643
  ```

-
- median-before-mean and median-before-std. This is the mean and standard deviation of the [median_before in the BLOW5 file](https://hasindu2008.github.io/slow5specs/summary).

Example command to get the two parameters:
 ```
  slow5tools skim -t40 PNXRXX240011_reads_500k.blow5 | awk -v c="median_before" 'NR==1{for (i=1; i<=NF; i++) if ($i==c){p=i; break};} {if ($p!=".") print $p}' | tail -n+2 | datamash mean 1 sstdev 1
  205.63935594369 8.3994882799157
 ```

- dwell-mean. This must be equal to the sample_rate/bps, and acts as a sanity check currently.

- dwell-std.
