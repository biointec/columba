# Columba
Fast Approximate Pattern Matching using Search Schemes

Columba was introduced in our [paper](https://doi.org/10.1016/j.isci.2021.102687). If you find this code useful in your research, please cite: 
```
@article{RENDERS2021102687,
title = {Dynamic partitioning of search patterns for approximate pattern matching using search schemes},
journal = {iScience},
volume = {24},
number = {7},
pages = {102687},
year = {2021},
issn = {2589-0042},
doi = {https://doi.org/10.1016/j.isci.2021.102687},
url = {https://www.sciencedirect.com/science/article/pii/S2589004221006556},
author = {Luca Renders and Kathleen Marchal and Jan Fostier},
keywords = {Algorithms, Bioinformatics, Computer science, High-performance computing in bioinformatics},
```

These instructions will get you a copy of the project up and running on your local machine.

#  Prerequisites

This package requires a number of packages to be install on your system. Required: CMake (3.0 or higher); Google's Sparsehash; gcc (GCC 6.4 or a more recent version) 

How to install these packages:

As a root, execute the following commands:

on Redhat / Fedora distributions
```bash
yum install cmake
yum install sparsehash-devel
``` 

on Ubuntu / Debian distributions
```bash
apt-get install cmake
apt-get install libsparsehash-dev
```  
# Installing Columba:

The installation is now simple. First, clone columba from the Github address

    git clone "https://github.com/biointec/columba.git"

From this directory, run the following commands:
```bash
mkdir build
cd build
cmake ..
make 
```
# Usage:
Columba aligns reads to a bidirectional FM-index. To do this you need to build the FM-index based on the input data. Currenly we only support input data with an alphabet of length 5 (e.g. for DNA A, C, G, T + $).

## Building the index
To build the bidirectional FM-index the text and the suffix arrays of the text and reverse text are required. 
To reverse a text you can use the commando `rev`:
```bash
rev [basefile].txt > [basefile].rev.txt
```
This commando will create the reverse text and store it in `[basefile].rev.txt`
We recommend the use of [radixSA64](https://github.com/mariusmni/radixSA64) for building the suffix arrays.
```bash
# make SA for original text
[pathToRadixSA]/radixSA [basefile].txt [basefile].sa
# make SA for reversed text
[pathToRadixSA]/radixSA [basefile].rev.txt [basefile].rev.sa
```
To build the FM-index run the following command in the `build` folder. 
```bash
./columba-build [basefile]
```
### Example 1
After installing columba, the columba directory should look like this:

    .
    ├── cmake
    ├── build
    ├── search_schemes
    └── src
    
In this example we will build the FM index for the 21st chromosome of the human genome. In order to do this we will create an `example` folder.
    To create this folder navigate to the columba folder. Here enter the following command
```bash
mkdir example
```
To this new directoy, copy the example file found [here](https://github.com/biointec/columba/releases/download/v1.0/genome.hs.chr_21.txt). This is the 21st chromosome of the human genome where all non-ACGT character were removed. A sentinel character was also appended to this file.

Reverse the text using the `rev` command:
```bash
rev genome.hs.chr_21.txt > genome.hs.chr_21.rev.txt
```
After reversing the text your directory structure should look like:

    .
    ├── cmake
    ├── build
    ├── example
    |   ├── genome.hs.chr_21.rev.txt
    |   └── genome.hs.chr_21.txt
    ├── search_schemes
    └── src

Now we need to create the suffix arrays. To do this enter the following commands:
```bash
# make SA for original text
[pathToRadixSA]/radixSA genome.hs.chr_21.txt`genome.hs.chr_21.sa
# make SA for reversed text
[pathToRadixSA]/radixSA genome.hs.chr_21.rev.txt genome.hs.chr_21.rev.sa
```

Where `[pathToRadixSA]` is the path to where you installed radixSA.

After this operation your directory structure will look like:

    .
    ├── cmake
    ├── build
    ├── example
    |   ├── genome.hs.chr_21.rev.sa
    |   ├── genome.hs.chr_21.rev.txt
    |   ├── genome.hs.chr_21.sa    
    |   └── genome.hs.chr_21.txt
    ├── search_schemes
    └── src

Finally, everything is in place to build the FM-index!
To build the FM-index, navigate to the `build` folder and run the following command:
```bash
./columba-build ../example/genome.hs.chr_21
```
The index files are then written to the same folder. Your directory structure will now look like:
 ```
    .
    ├── cmake
    ├── build
    ├── example 
    |   ├── genome.hs.chr_21.brt
    |   ├── genome.hs.chr_21.bwt
    |   ├── genome.hs.chr_21.cct
    |   ├── genome.hs.chr_21.rev.brt   
    |   ├── genome.hs.chr_21.rev.sa
    |   ├── genome.hs.chr_21.rev.txt
    |   ├── genome.hs.chr_21.sa
    |   ├── genome.hs.chr_21.sa.1 
    |   ├── genome.hs.chr_21.sa.128 
    |   ├── genome.hs.chr_21.sa.16      
    |   ├── genome.hs.chr_21.sa.2
    |   ├── genome.hs.chr_21.sa.32
    |   ├── genome.hs.chr_21.sa.4   
    |   ├── genome.hs.chr_21.sa.64
    |   ├── genome.hs.chr_21.sa.8
    |   ├── genome.hs.chr_21.sa.bv.1  
    |   ├── genome.hs.chr_21.sa.bv.128
    |   ├── genome.hs.chr_21.sa.bv.16
    |   ├── genome.hs.chr_21.sa.bv.2
    |   ├── genome.hs.chr_21.sa.bv.32
    |   ├── genome.hs.chr_21.sa.bv.4
    |   ├── genome.hs.chr_21.sa.bv.64
    |   ├── genome.hs.chr_21.sa.bv.8
    |   └── genome.hs.chr_21.txt
    ├── search_schemes
    └── src
```

Congratulations! You have used columba to build the FM-index of the 21st chormosome of the human genome!

---
**CAUTION!**

If you have used Columba 1.0 to build your index files you will have to rerun the build process to work with Columba 1.1 (and vice versa)!
---


## Using the index
Columba can align reads in a fasta (`.FASTA`, `.fasta`, `.fa`) or fastq (`.fq`, `.fastq`) format. 
To align your reads use the following format 

```bash
./columba [options] basefile readfile.[ext]
```

options:

```
[options]
  -e  --max-ed          maximum edit distance [default = 0]
  -s  --sa-sparseness   suffix array sparseness factor [default = 1]
  -p  --partitioning    Add flag to do uniform/static/dynamic partitioning [default = dynamic]
  -m  --metric          Add flag to set distance metric (editnaive/editopt/hamming) [default = editopt];
  -i  --in-text	The tipping point for in-text verification [default = 5]
  -ss --search-scheme   Choose the search scheme
  options:
        kuch1   Kucherov k + 1
        kuch2   Kucherov k + 2
        kianfar  Optimal Kianfar scheme
        manbest  Manual best improvement for kianfar scheme (only for ed = 4)
        pigeon   Pigeonhole scheme
        01*0    01*0 search scheme
        custom  custom search scheme, takes one parameter which is the path to the folder containing this search scheme

[ext]
        one of the following: fq, fastq, FASTA, fasta, fa
Following input files are required:
        <base filename>.txt: input text T
        <base filename>.cct: character counts table
        <base filename>.sa.[saSF]: suffix array sample every [saSF] elements
        <base filename>.bwt: BWT of T
        <base filename>.brt: Prefix occurrence table of T
        <base filename>.rev.brt: Prefix occurrence table of the reverse of T
    
```

The number of nodes, duration, and number of reported/unique matches will be printed to stdout, as well as the number of matches found entirely in the index, the number of unique matches found via in-text verification, the number of started and failed in-text verification procedures and the number of searches that started in the index.
The matches will be written to a custom output file in the folder where your readfile was. This output file will be a tab- separated file with the fields: `identifier`, `position`, `length`, `ED`, `CIGAR` and `reverse strand`. For each optimal alignment under the maximal given edit distance a line will be present. This output file will be called `readfile_output.txt`.






### Example 2
Consider the final directory structure from [example 1](##Example-1). 
Copy this [file](https://github.com/biointec/columba/releases/download/v1.0/genome.hs.chr_21.reads.fasta) to this directory. 
This file contains 100 000 reads of length 100 all sampled from the reference text. Thus, each read will have at least one exact occurrence.
If you want to align these reads using the Pigeonhole scheme with k = 3 and using the edit distance and optimal static partitioning to our reference text, run the following command in the `build` folder:
```bash
./columba -e 3 -ss pigeon -m editopt -p static ../example/genome.hs.chr_21 ../example/genome.hs.chr_21.reads.fasta
```

After this operation your directory structure will look like:
``` .
    ├── cmake
    ├── build
    ├── example 
    |   ├── genome.hs.chr_21.brt
    |   ├── genome.hs.chr_21.bwt
    |   ├── genome.hs.chr_21.cct    
    |   ├── genome.hs.chr_21.reads.fasta
    |   ├── genome.hs.chr_21.reads.fasta_output.txt
    |   ├── genome.hs.chr_21.rev.brt   
    |   ├── genome.hs.chr_21.rev.sa
    |   ├── genome.hs.chr_21.rev.txt
    |   ├── genome.hs.chr_21.sa
    |   ├── genome.hs.chr_21.sa.1    
    |   ├── genome.hs.chr_21.sa.2
    |   ├── genome.hs.chr_21.sa.4   
    |   ├── genome.hs.chr_21.sa.8
    |   ├── genome.hs.chr_21.sa.16    
    |   ├── genome.hs.chr_21.sa.32
    |   ├── genome.hs.chr_21.sa.64
    |   ├── genome.hs.chr_21.sa.128
    |   └── genome.hs.chr_21.txt
    └── src
```

The results can be found in `genome.hs.chr_21.reads.fasta_output.txt`.

To align the reads with maximal hamming distance 2, dynamic partitioning and the Kucherov K + 1 search schem run
```bash
./columba -m hamming -e 2 -ss kuch1 -p dynamic ../example/genome.hs.chr_21 ../example/genome.hs.chr_21.reads.fasta
```

---
**NOTE**

This second alignment of reads will overwrite the `genome.hs.chr_21.reads.fasta_output.txt`. Before running a second time, make sure to back-up the original file if you would like to keep it stored.

---

Congratulations! You are now able to use Columba to align reads to the 21st chromosome of the human genome!

## Custom Search Schemes

The search scheme can either be one of the hardcoded search schemes present in Columba or you can provide a custom search scheme. In the `search_schemes` folder a number of search schemes is already present. 

To make your own search scheme you need to create a folder containing at least a file called `name.txt`, which contains the name of your scheme on the first line. 
For every maximum edit/hamming distance a subfolder should be present, which contains at least the file `searches.txt`. In this file the searches of your scheme are written line per line. Each line contains of three space-separated arrays: pi, L and U. Each array is written between curly braces {} and the values are comma-separated.

---
**NOTE**

The pi array should be zero-based! The connectivity property should always be satifsfied. The L and U array cannot decrease.

---
### Static Partitioning
If you want to provide optimal static partitioning you can create a file named `static_partitioning.txt` in the folder of the maximum edit/hamming distance this partitioning is for. This file should contain one line with percentages (values between 0 and 1) separated by spaces. The ith percentage corresponds to the starting position (relative to the size of the pattern) of the (i + 1)th part (again this is zero based). The starting position of the first part is always zero and should **not** be provided.

### Dynamic Partitionining
Similarly, to provide values for dynamic partitioning you can create a file called `dynamic_partitioning.txt`. This file should contain two lines. The first line are percentages (again between 0 and 1) that correspond to the seeding positions, relative to the size of the pattern, of all parts, except the first and last part. 
The second line should contain space-separated integers corresponding to the weights of each part.

### Folder Structure Example
Consider a search scheme which supports maximal edit/hamming distances 1, 2 and 4. For distance 1 no static or dynamic partitioning values are known. For distance 2 only static partitioning values are known and for distance 4 both static and dynamic partitioning values are known. The folder structure of this search scheme should look like this:

``` .
    ├── 1
    |   ├── searches.txt
    ├── 2 
    |   ├── searches.txt
    |   ├── static_partitioning.txt
    ├── 4
    |   ├── dynamic_partitioning.txt
    |   ├── searches.txt
    |   ├── static_partitioning.txt
    └── name.txt
```

### Example `searches.txt`
Consider the pigeon hole search scheme for maximum edit distance 4. The `searches.txt` file should look like:

```
{0,1,2,3,4} {0,0,0,0,0} {0,4,4,4,4}
{1,2,3,4,0} {0,0,0,0,0} {0,4,4,4,4}
{2,3,4,1,0} {0,0,0,0,0} {0,4,4,4,4}
{3,4,2,1,0} {0,0,0,0,0} {0,4,4,4,4}
{4,3,2,1,0} {0,0,0,0,0} {0,4,4,4,4}
```

### Example 2
The adapted search schemes based on those by Kucherov can be found in the directories:
`search_schemes/kuch_k+1_adapted` and `search_schemes/kuch_k+2_adapted`.

### Other examples
In the `search_schemes` folder the hardcoded search schemes of Columba are available as custom search schemes. 

## In-text verification
Columba 1.1 introduces the ability to switch to in-text verification if the number of occurrences in the reference text is lower then some tipping point t. This tipping point can be set via the parameter `-i` or `--in-text`.


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
We reuse the reference as in the experiments of Columba 1.0

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

## SA space-time tradeoff
To reproduce the results for different tipping points and different sparseness factors, for `k` = 4 run:

```bash
./columba -ss custom ../search_schemes/kuch_k+1_adapted/ -e 4 -i [t] -s [s] basefile readsfile
```
with `[t]` the tipping point and `[s]` the sparseness factor.