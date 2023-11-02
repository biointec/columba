# Result Reproduction
This README provides instructions on how to reproduce the results described in our publications on Columba. [Section 1](#reproducing-results-1) covers our paper titled "Dynamic partitioning of search patterns for approximate pattern matching using search schemes" published in iScience. 
Next, [section 2](#reproducing-results-2) details the experiments laid out in our second paper "Approximate pattern matching using search schemes and in-text verification". 
Finally, [section 3](reproducing-results-3) gives instructions how to reproduce the results as outlined in "Automated design of efficient
search schemes for lossless approximate pattern matching".

# Reproducing results 1
The results form our paper: [Dynamic partitioning of search patterns for approximate pattern matching using search schemes](https://doi.org/10.1016/j.isci.2021.102687) can be reproduced by using the following instructions.

First download and compile the source code of [Columba 1.0](https://github.com/biointec/columba/releases/tag/v1.0).

## Dataset

### Reference
We used the [human genome](https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.26/) as our reference genome. From this text we removed all N's and substituted them with a random nucleotide.
Then we removed the first line, concatenated the chromosomes, removed newline characters and added a dollar to the end of file by executing the following command
```bash
tail -n +2 hs.grch38.fasta.non | sed "s/>.*//g" | tr -d '\n' | cat - dollar > hs.grch38.txt
```
Where `hs.grch38.fasta.non` is the result of the substitution of the N's and `hs.grch38.txt` is the reference text used to build the index. `dollar` is a text file containing a single character `$`. This file should be present in your directory before executing the above command. 

For this reference text the index needs to be build according to the build instructions described above.

### Reads
We sampled a 100 000 reads from [an Illumina experiment dataset](ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR194/ERR194147/ERR194147_1.fastq.gz), which only contain `A`, `C`, `G` and `T` characters. One way to do this is to select the first 100 000 such reads. For example using the following python3 script

```python
import sys

read_file = sys.argv[1]
write_file = sys.argv[2]
f = open(read_file, "r")
w = open(write_file, "w")

reads_written = 0

while(reads_written < 100000):
    fastq_record = [next(f) for x in range(4)]
    read = fastq_record[1]
    if(read.count('A') + read.count('C') + read.count('G') + read.count('T') != len(read.strip())):
        # contains a non-ACGT character
        continue
    for line in fastq_record:
        w.write(line)
    reads_written += 1

f.close()
w.close()


```

This python script takes an input name and requested output name. It samples 100 000 reads and writes them to a fastq file [output].fastq, it is assumed that the input file is a fastq file. Of course sampling of the reads can be achieved in a myriad of ways.

The sampled dataset we used for the results in the paper is available [here](https://github.com/biointec/columba/releases/download/v1.0/sampled_illumina_reads.fastq).

## Partitioning strategies
To reproduce the results from table 1, run Columba with `kuch1` as option for search scheme and the different partitioning strategies `[part]` and maximal edit distances `[k]`, from the `build` folder.
```bash
./columba -ss kuch1 -e [k] -partitioning [part] basefile readsfile
```
Another option is to use the custom search schemes:
```bash
./columba -ss custom ../search_schemes/kuch_k+1 -e [k] -partitioning [part] basefile readsfile
```


## Non-interleaved bitvectors
Table 2 contains the results on interleaved and non-interleaved bitvectors. The results for interleaved bitvectors are the same as those aquired from the previous section for k = 4. For non-interleaved bitvectors some steps need to be taken.
Navigate to file `src/bwtrepr.h`, on line 40 replace:
```c++
BitvecIntl<S - 1> bv; // bitvector representation of the BWT
```
by 
```c++
BitvecNonIntl<S - 1> bv; // bitvector representation of the BWT
```
Save this file and navigate to the build folder.
Here run the following command
```bash
make
```
You will need to rebuild the index (see [Building the index](#Building-the-index)).
After rebuilding the index you can get the results with non-interleaved bitvectors for k = 4 with:
```bash
./columba -ss kuch1 -e 4 -partitioning [part] basefile readsfile
```

---
**NOTE**

Make sure to undo these changes before using the tool, as interleaved bitvectors are superior to non-interleaved bitvectors. To undo the changes, simply revert `src/bwtrepr.h` back to the original, compile again with the `make` command and rebuild the index. 

---

## Effect of reducing the redundancy 
Table 3 compares the effect of a naive and optimized implementation of the edit distance metric. The results for the optimized implementation are the same as achieved [in the section partitioning strategies](#Partitioning-strategies) with dynamic partitioning. The results for a naive implementation can be reproduced for different values of k by running:

```bash
./columba -ss kuch1 -e [k] -p dynamic -m editnaive basefile readsfile
```

## Other search schemes
Tables 4, 6, 7, 8, 9 and 10 contains results for other search schemes. To reproduce these results run:
```bash
./columba -ss [search scheme] -e [k] -partitioning [p] basefile readsfile
```

where `[search scheme]` is one of: `kuch2`, `kianfar`, `manbest`, `pigeon`, `01*0`, and `[p]` is one of: `uniform`, `static`, `dynamic`

---
**NOTE**

`manbest` only accepts k = 4

---

## Results on PacBio Data
Table 5 contains results on 100 000 PacBio seeds. These seeds are sampled from  [this PacBio experiment](https://www.ebi.ac.uk/ena/browser/view/SRR1304331).

The sampled dataset we used for the results in the paper is available [here](https://github.com/biointec/columba/releases/download/v1.0/sampled_pacbio_seeds.fastq)
To reproduce the results run 
```
./columba -ss kuch1 -e [k] -partitioning [part] basefile readsfile
```
from the `build` folder. Where readsfile is the file containing the PacBio seeds.

 

# Reproducing results 2

To reproduce the results from our paper: "Approximate pattern matching using search schemes and in-text verification" follow these instructions.
Download and compile the source code of [Columba 1.1](https://github.com/biointec/columba/releases/tag/v1.1).

## Dataset

### Reference
We reuse the reference as in the experiments of Columba 1.0 ([see above](#reference)).  Note that we need to rebuild the index for this newer version of Columba.

### Reads
We sampled a 100 000 reads from [an Illumina experiment dataset](ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR909/009/SRR9091899/SRR9091899_1.fastq.gz), which only contain `A`, `C`, `G` and `T` characters. The average length of the reads is 150.
Our sampled dataset is available [here](https://github.com/biointec/columba/releases/download/v1.1/sampled_illumina_reads_150.fastq)

## Adapted Search schemes
To reproduce the comparison of the adapted and original search schemes run:

```bash
./columba -ss custom [path] -e [k] basefile readsfile
```

Where `[path]` is the path to the folder containing the search scheme and `k` is the number of allowed errors.

## In-text verification
To reproduce the results for different tipping points and a maximal edit distance of 4, run:


```bash
./columba -ss custom ../search_schemes/kuch_k+1_adapted/ -e 4 -i [t] basefile readsfile
```
with `[t]` the tipping point.

## SA space-time trade-off
To reproduce the results for different tipping points and different sparseness factors, for `k` = 4 run:

```bash
./columba -ss custom ../search_schemes/kuch_k+1_adapted/ -e 4 -i [t] -s [s] basefile readsfile
```
with `[t]` the tipping point and `[s]` the sparseness factor.

# Reproducing results 3

To reproduce the results from our paper: "Automated design of efficient search schemes for lossless approximate pattern matching" follow these instructions.
Download and compile the source code of [Columba 1.2](https://github.com/biointec/columba/releases/tag/v1.2).

## Dataset

### Reference
We reuse the reference as in the experiments of Columba 1.0 ([see above](#reference)). Note that we need to rebuild the index for this newer version of Columba.

### Reads
We sampled a 100 000 reads from [an Illumina experiment dataset](ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR909/009/SRR9091899/SRR9091899_1.fastq.gz), which only contain `A`, `C`, `G` and `T` characters. The average length of the reads is 150.
Our sampled dataset is available [here](https://github.com/biointec/columba/releases/download/v1.1/sampled_illumina_reads_150.fastq)

## Comparing Search Schemes Across Search Space Size
To reproduce the comparison of the different search schemes run:

```bash
./columba -ss custom [path] -e [k] -p uniform [base] [readsfile]
```

Where `[path]` is the path to the folder containing the search scheme and `k` is the number of allowed errors. 

Below, a table is provided that links the search scheme name, as used in the paper, to the name of the directory in [the search schemes folder](../search_schemes).

| Name in paper             | Name of directory                 | 
| --------                  | --------                          | 
| Kucherov k + 1            | kuch_k+1/                         | 
| Kucherov k + 2            | kuch_k+2/                         | 
| Pigeonhole principle      | pigeon/                           | 
| 01*0                      | 01star0/                          | 
| Man<sub>best</sub>        | manbest/                          | 
| Kianfar                   | kianfar/                          | 
| Greedy Kucherov k + 1     | kuch_k+1_adapted/                 | 
| Greedy Pigeonhole princ.  | pigeon_adapted/                   | 
| Greedy 01*0               | 01star0_adapted/                  | 
| minU (ILP)                | multiple_opt/individual_schemes/  | 

Note that for the minU schemes, the average value was reported over the co-optimal schemes in the respective directory.

## Dynamic Selection

To get the results with dynamic selection run:

```bash
./columba -ss mutliple ../search_schemes/multiple_opt -e [k] -p uniform [base] [readsfile]
```

## Comparison with the State-of-the-art

To get the timings for Columba 1.2 used in the comparison with the state of the art run:


```bash
./columba -ss multiple ../search_schemes/multiple_opt -e [k] -i 5 [base] [readsfile]
```

for even `k` and:

```bash
./columba -ss multiple ../search_schemes/multiple_opt/individual_schemes/scheme1 -e [k] -i 5 [base] [readsfile]
```
for odd `k`.

Note that to get the timings for Columba 1.1 for `k >= 5` the Columba 1.1's source code needs to be slightly adapted.

