# Columba
Fast Approximate Pattern Matching using Search Schemes

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
    ├── build
    └── src
    
In this example we will build the FM index for the 21st chromosome of the human genome. In order to do this we will create an `example` folder.
    To create this folder navigate to the columba folder. Here enter the following command
```bash
mkdir example
```
To this new directoy, copy the example file found [here](https://github.com/biointec/columba/releases/download/example/genome.hs.chr_21.txt). This is the 21st chromosome of the human genome where all non-ACGT character were removed. A sentinel character was also appended to this file.

Reverse the text using the `rev` command:
```bash
rev genome.hs.chr_21.txt > genome.hs.chr_21.rev.txt
```
After reversing the text your directory structure should look like:

    .
    ├── build
    ├── example
    |   ├── genome.hs.chr_21.rev.txt
    |   └── genome.hs.chr_21.txt
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
    ├── build
    ├── example
    |   ├── genome.hs.chr_21.rev.sa
    |   ├── genome.hs.chr_21.rev.txt
    |   ├── genome.hs.chr_21.sa    
    |   └── genome.hs.chr_21.txt
    └── src

Finally, everything is in place to build the FM-index!
To build the FM-index, navigate to the `build` folder and run the following command:
```bash
./columba-build ../example/genome.hs.chr_21
```
The index files are then written to the same folder. Your directory structure will look like:
 ```
    .
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

Congratulations! You have used columba to build the FM-index of the 21st chormosome of the human genome!
## Using the index
Columba can align reads in a fasta (`.FASTA`, `.fasta`, `.fa`) or fastq (`.fq`, `.fastq`) format. 
To align your reads use the following format 

```bash
./columba [options] basefile readfile
```

options:

```
[options]
  -e  --max-ed          maximum edit distance [default = 0]
  -s  --sa-sparseness   suffix array sparseness factor [default = 1]
  -p  --partitioning    Add flag to do uniform/static/dynamic partitioning [default = dynamic]
  -h   --hamming        Add flag to use hamming distance [default = false]
  -ss --search-scheme   Choose the search scheme
  options:
        kuch1   Kucherov k + 1
        kuch2   Kucherov k + 2
        kianfar  Optimal Kianfar scheme
        manbest  Manual best improvement for kianfar scheme (only for ed = 4)
        pigeon   Pigeonhole scheme
        01*0    01*0 search scheme

[ext]
        one of the following: fq, fastq, FASTA, fasta, faFollowing input files are required:
        <base filename>.txt: input text T
        <base filename>.cct: charachter counts table
        <base filename>.sa.[saSF]: suffix array sample every [saSF] elements
        <base filename>.bwt: BWT of T
        <base filename>.brt: Prefix occurrence table of T
        <base filename>.rev.brt: Prefix occurrence table of the reverse of T
    
```

The number of nodes, number of matrix elements, duration, and number of matches will be printed to stdout. 
The matches will be written to a custom output file in the folder where your readfile was. This output file will be a tab-seperated file with the fields: `identifier`, `position`, `length` and `ED`. For each optimal alignment under the maximal given edit distance a line will be present. This output file will be called `readfile_output.txt`.

### Example 2
Consider the final directory structure from [example 1](##Example-1). 
Copy this [file](https://github.com/biointec/columba/releases/download/example/genome.hs.chr_21.reads.fasta) to this directory. 
This file contains 100 000 reads of length 100 all sampled from the reference text. Thus, each read will have at least one exact occurrence.
If you want to align these reads using the Pigeonhole scheme with k = 3 and using the edit distance and optimal static partitioning to our refercen text, run the following command in the `build` folder:
```bash
./columba -e 3 -ss pigeon -p static ../example/genome.hs.chr_21 ../example/genome.hs.chr_21.reads.fasta
```

After this operation your directory structure will look like:
```    .
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
./columba -h -e 2 -ss kuch1 -p dynamic ../example/genome.hs.chr_21 ../example/genome.hs.chr_21.reads.fasta
```

---
**NOTE**

This second alignment of reads will overwrite the `genome.hs.chr_21.reads.fasta_output.txt`. Before running a second time, make sure to back-up the original file if you would like to keep it stored.

---

Congratulations! You are now able to use Columba to align reads to the 21st chromosome of the human genome!

# Reproducing results

## Dataset

### Reference
We used the [human genome](https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.26/) as our reference genome. From this text we removed all N's and substituted them with a random nucleotide.
Then we removed the first line, concatenated the chromosomes, removed newline characters and added a dollar to the end of file by executing the following command
```bash
tail -n +2 hs.grch38.fasta.non | sed "s/>.*//g" | tr -d '\n' | cat - dollar > hs.grch38.txt
```
Where `hs.grch38.fasta.non` is the result of the substitution of the N's and `hs.grch38.txt` is the reference text used to build the index. `dollar` is a text file containing a single character `$`. This file should be present in your directory before executing the above command. 

### Reads
Sample 100000 reads of length 100 from the acquired human genome with only {A, C, T, G} characters. 
The following python3 script is one way to  achieve this. 
```python
import sys
import numpy as np


filename = sys.argv[1]
output = sys.argv[2]
f = open(filename , "r")
o = open(output + ".fasta", "w")
content = f.read()
f.close()

for i in range(100000):
    # get postion and string
    p = np.random.randint(len(content) - 100)
    s = content[p: p + 100]
    # write identifier line (where id = position)
    o.write("> " + str(p) + "\n")
    # write string lines of max length 80
    o.write(s[0:80] + "\n")
    o.write(s[81:] + "\n")
   
o.close()

```

This script takes a filename and requested output name. It samples 100 000 reads of length 100 and writes them to a fasta file `[output].fasta`. Of course sampling of the reads can be achieved in a myriad of ways.

## Partitioning strategies
To reproduce the results from table 1, run Columba with `kuch1` as option for search scheme and the different partitioning strategies `[part]` and maximal edit distances `[k]`.
```bash
./columba -ss kuch1 -e [k] -partitioning [part] basefile readsfile
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

Make sure to undo these changes before using the tool, as interleaved bitvectors are superior to non-interleaved bitvectors. To undo the changes, simply revert `src/bwtrepr.h` back to the original and rebuild the index. 

---

## Other search schemes
Table 3 contains results for other search schemes for dynamic partitioning. To reproduce these results run:
```bash
./columba -ss [search scheme] -e [k] -partitioning dynamic basefile readsfile
```

where `[search scheme]` is one of: `kuch2`, `kianfar`, `manbest`, `pigeon`, `01*0`

---
**NOTE**

`manbest` only accepts k = 4

---
