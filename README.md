# Columba

Fast Approximate Pattern Matching using Search Schemes

Columba is a powerful open-source read-mapper developed to significantly enhance the performance of lossless approximate pattern matching. This README provides an overview of the [key features and benefits](#key-features-and-benefits) of Columba, along with instructions for [installation](#installation), [usage](#usage), [result reproduction](#result-reproduction) and [citation](#citation).
Finally, [contact information](#contact) is provided.

## Key features and benefits

- Columba is a **lossless** read-mapper, meaning that **all** occurrences up to the provided distance are found. Both the edit distance and the hamming distance are supported.
- Columba can handle any valid search scheme, so if new search schemes are discovered no update is needed. See [custom search schemes](#custom-search-schemes) below.
- Columba supports a dynamic selection of search schemes, where the search scheme best suited for the current read and reference text is chosen, resulting in faster runtimes! More info [here](#dynamic-selection---multiple-search-schemes).
- Columba supports dynamic partitioning of the read, boosting the performance!
- Columba comes with tunable in-text verification, which improves runtimes significantly. More info [here](#in-text-verification).
- Columba is **fast**, outperforming other lossless aligners, thanks to all kinds of algorithmic tricks like an interleaved bit-vector representation and redundancy avoidance for the edit distance metric.

## Installation

The following instructions will get you a copy of the project up and running on your local machine.

### Prerequisites

This package requires a number of packages to be install on your system. Required: CMake (3.0 or higher); Google's Sparsehash; gcc (tested on 8.3.0 and more recent compilers).

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

### Installing Columba

The installation is now simple. First, clone columba from the Github address

```terminal
    git clone "https://github.com/biointec/columba.git"
```

From this directory, run the following commands to install columba:

```bash
mkdir build
cd build
cmake ..
make 
```

---
**NOTE!**

If your reference genome is longer than 4.29M characters, you need to compile Columba in 64-bit mode. To do so use these commands instead:

```bash
mkdir build
cd build
cmake -DTHIRTY_TWO=OFF ..
make 
```

---

## Usage

Columba aligns reads to a bidirectional FM-index. To do this you need to build the FM-index based on the input data. Currently we only support input data with an alphabet of length 5 (e.g. for DNA A, C, G, T + $).

### Building the index

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

#### Example 1

After installing columba, the columba directory should look like this:

```plaintext
    .
    ├── cmake
    ├── build
    ├── search_schemes
    └── src
```

In this example we will build the FM index for the 21st chromosome of the human genome. In order to do this we will create an `example` folder.
    To create this folder navigate to the columba folder. Here enter the following command

```bash
mkdir example
```

To this new directory, copy the example file found [here](https://github.com/biointec/columba/releases/download/v1.0/genome.hs.chr_21.txt). This is the 21st chromosome of the human genome where all non-ACGT character were removed. A sentinel character was also appended to this file.

Reverse the text using the `rev` command:

```bash
rev genome.hs.chr_21.txt > genome.hs.chr_21.rev.txt
```

After reversing the text your directory structure should look like:

```plaintext
    .
    ├── cmake
    ├── build
    ├── example
    |   ├── genome.hs.chr_21.rev.txt
    |   └── genome.hs.chr_21.txt
    ├── search_schemes
    └── src
```

Now we need to create the suffix arrays. To do this enter the following commands:

```bash
# make SA for original text
[pathToRadixSA64]/radixSA64 genome.hs.chr_21.txt genome.hs.chr_21.sa
# make SA for reversed text
[pathToRadixSA64]/radixSA64 genome.hs.chr_21.rev.txt genome.hs.chr_21.rev.sa
```

Where `[pathToRadixSA64]` is the path to where you installed radixSA64.

After this operation your directory structure will look like:

```plaintext
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
```

Finally, everything is in place to build the FM-index!
To build the FM-index, navigate to the `build` folder and run the following command:

```bash
./columba-build ../example/genome.hs.chr_21
```

The index files are then written to the same folder. Your directory structure will now look like:

```plaintext
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

Congratulations! You have used columba to build the FM-index of the 21st chromosome of the human genome!

Note that the suffix array is created with different sparseness factors. This gives the user the option to choose any sparseness factor, depending on the amount of available RAM.

---
**CAUTION!**

If you have used an older version of Columba to build your index files you will have to rerun the build process to work with the newest version (and vice versa)!

---

### Using the index

This section explains how to use a built index for a reference genome. First a general overview is given.
Finally, the details on how to use [custom search schemes](#custom-search-schemes) and [dynamic selection of search schemes](#dynamic-selection---multiple-search-schemes) are provided.

Columba can align reads in a fasta (`.FASTA`, `.fasta`, `.fa`) or fastq (`.fq`, `.fastq`) format.
To align your reads use the following format:

```bash
./columba [options] basefile readfile.[ext]
```

options:

```plaintext
 [options]
  -e  --max-ed          maximum edit distance [default = 0]
  -s  --sa-sparseness   suffix array sparseness factor [default = 1]
  -p  --partitioning    Add flag to do uniform/static/dynamic partitioning [default = dynamic]
  -m   --metric         Add flag to set distance metric (editnaive/editopt/hamming) [default = editopt]
  -i  --in-text         The tipping point for in-text verification [default = 5]
  -ks --kmer-size       The size of the seeds for dynamic partitioning [default = 10]
  -o  --output          The name of the outputfile. This file will be in .sam format. [default = ColumbaOutput.sam]
  -ss --search-scheme   Choose the search scheme
  options:
        kuch1           Kucherov k + 1
        kuch2           Kucherov k + 2
        kianfar         Optimal Kianfar scheme
        manbest         Manual best improvement for kianfar scheme (only for ed = 4)
        pigeon          Pigeon hole scheme
        01*0            01*0 search scheme
        custom          custom search scheme, the next parameter should be a path to the folder containing this search scheme
        multiple        multiple search scheme, the next parameter should be a path to the folder containing the different search schemes to choose from with dynamic selection.

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
The matches will be written in SAM format to the given output file (default = ColumbaOutput.sam).

#### Example 2

Consider the final directory structure from [example 1](#example-1).
Copy this [file](https://github.com/biointec/columba/releases/download/v1.0/genome.hs.chr_21.reads.fasta) to this directory.
This file contains 100 000 reads of length 100 all sampled from the reference text. Thus, each read will have at least one exact occurrence.
If you want to align these reads using the Pigeonhole scheme with k = 3 and using the edit distance and optimal static partitioning to our reference text, run the following command in the `build` folder:

```bash
./columba -e 3 -ss pigeon -m editopt -p static -o ../example/output.sam ../example/genome.hs.chr_21 ../example/genome.hs.chr_21.reads.fasta
```

After this operation your directory structure will look like:

```plaintext
    .
    ├── cmake
    ├── build
    ├── example 
    |   ├── genome.hs.chr_21.brt
    |   ├── genome.hs.chr_21.bwt
    |   ├── genome.hs.chr_21.cct    
    |   ├── genome.hs.chr_21.reads.fasta
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
    |   ├── genome.hs.chr_21.txt
    |   └── output.sam
    └── src
```

The results can be found in `output.sam`.

To align the reads with maximal hamming distance 2, dynamic partitioning, in-text verification at 5, and the Kucherov K + 1 search scheme run

```bash
./columba -m hamming -e 2 -ss kuch1 -p dynamic -i 5 -o output2.sam ../example/genome.hs.chr_21 ../example/genome.hs.chr_21.reads.fasta
```

The results are in file `output2.Sam`.

Congratulations! You are now able to use Columba to align reads to the 21st chromosome of the human genome!

### In-text verification

To speed up the read-mapper, we introduced in-text verification in Columba 1.1. Using this technique the potential in-index matches are looked up using the suffix array and verified directly in the text if the number of potential in-index matches falls below or equals some threshold `t`. This can reduce runtimes by  up to 50%!

We advise to use a threshold of 5.
However, if you need to use a very sparse suffix array due to memory constraints, a lower value could perform better.

### Custom Search Schemes

The search scheme can either be one of the hardcoded search schemes present in Columba or you can provide a custom search scheme. In the [search_schemes](./search_schemes/) folder a number of search schemes is already present. Including the search schemes created by Hato.

To make your own search scheme you need to create a folder containing at least a file called `name.txt`, which contains the name of your search scheme on the first line.
For every maximum edit/hamming distance a subfolder should be present, which contains at least the file `searches.txt`. In this file the searches of your scheme are written line per line. Each line contains of three space-separated arrays: pi, L and U. Each array is written between curly braces {} and the values are comma-separated.

---
**NOTE!**

The pi array should be zero-based! The connectivity property should always be satisfied. The L and U array cannot decrease. All error distributions should be covered! To check if your search scheme is valid you can use a python script provided [in this directory](./validitychecker/).

---

#### Static Partitioning

If you want to provide optimal static partitioning you can create a file named `static_partitioning.txt` in the folder of the maximum edit/hamming distance this partitioning is for. This file should contain one line with percentages (values between 0 and 1) separated by spaces. The ith percentage corresponds to the starting position (relative to the size of the pattern) of the (i + 1)th part (again this is zero based). The starting position of the first part is always zero and should **not** be provided.

#### Dynamic Partitioning

Similarly, to provide values for dynamic partitioning you can create a file called `dynamic_partitioning.txt`. This file should contain two lines. The first line are percentages (again between 0 and 1) that correspond to the seeding positions, relative to the size of the pattern, of all parts, except the first and last part.
The second line should contain space-separated integers corresponding to the weights of each part.

#### Folder Structure Example

Consider a search scheme which supports maximal edit/hamming distances 1, 2 and 4. For distance 1 no static or dynamic partitioning values are known. For distance 2 only static partitioning values are known and for distance 4 both static and dynamic partitioning values are known. The folder structure of this search scheme should look like this:

```plaintext
    .
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

#### Example `searches.txt`

Consider the pigeon hole search scheme for maximum edit distance 4. The `searches.txt` file should look like:

```plaintext
{0,1,2,3,4} {0,0,0,0,0} {0,4,4,4,4}
{1,2,3,4,0} {0,0,0,0,0} {0,4,4,4,4}
{2,3,4,1,0} {0,0,0,0,0} {0,4,4,4,4}
{3,4,2,1,0} {0,0,0,0,0} {0,4,4,4,4}
{4,3,2,1,0} {0,0,0,0,0} {0,4,4,4,4}
```

#### Other examples

In the [`search_schemes` directory](./search_schemes), a number of search schemes are available as custom search schemes.

### Dynamic Selection - Multiple Search Schemes

To use dynamic selection you can provide a set of search schemes to be used.
For this you need to create a folder with a file called `name.txt` in it and for all values of `k`, where you want to provide multiple search schemes you need to create a directory.
The files in these subdirectory are named `scheme<x>.txt`, where `<x>` is the number of the search scheme (starting from 1).
Each file contains the searches from one of the search schemes and the file needs to be structured like the `searches.txt` file of a custom search strategy (see [above](#custom-search-schemes)).

#### Example Folder Structure

Consider a strategy where multiple search schemes are used which supports maximal edit/hamming distances  2, 4 and 6. For these values of `k` there exist respectively 2, 3 and 4 search schemes. The folder structure then should look like this:

```plaintext
    .
    ├── 2
    |   ├── scheme1.txt
    |   ├── scheme2.txt
    ├── 4 
    |   ├── scheme1.txt
    |   ├── scheme2.txt
    |   ├── scheme3.txt
    ├── 6
    |   ├── scheme1.txt
    |   ├── scheme2.txt
    |   ├── scheme3.txt
    |   ├── scheme4.txt
    └── name.txt
```

#### MinU search schemes

The minU search schemes, introduced in our paper "Automated design of efficient search schemes for lossless approximate pattern matching" are available [here](./search_schemes/multiple_opt).
In this directory there are folders made for dynamic selection for even values of `k`. It also contains a subdirectory `individual_schemes`, where each (co-)optimal scheme is presented.

## Result reproduction

To reproduce the results presented in our papers please refer to [these instructions](./result_reproduction/README.md).

## Citation

Columba was first introduced in our [paper](https://doi.org/10.1016/j.isci.2021.102687). If you find this code useful in your research, please cite:

```bibtex
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

Columba 1.1 was introduced in a [conference paper](https://doi.org/10.1007/978-3-031-07802-6_36). If you use version 1.1 in your research please cite both versions:

```bibtex
@InProceedings{10.1007/978-3-031-07802-6_36,
author="Renders, Luca
and Depuydt, Lore
and Fostier, Jan",
editor="Rojas, Ignacio
and Valenzuela, Olga
and Rojas, Fernando
and Herrera, Luis Javier
and Ortu{\~{n}}o, Francisco",
title="Approximate Pattern Matching Using Search Schemes and In-Text Verification",
booktitle="Bioinformatics and Biomedical Engineering",
year="2022",
publisher="Springer International Publishing",
address="Cham",
pages="419--435",
isbn="978-3-031-07802-6"
}
```

Columba 1.2 was introduced in our newest paper "Automated design of efficient search schemes for lossless approximate pattern matching", accepted at RECOMB 2024.

```bibtex
@unpublished{renders2023,
  author = {Renders, L. and Depuydt, L. and Rahmann, S. and Fostier, J.},
  title = {Automated design of efficient search schemes for lossless approximate pattern matching},
  year = {2023},
  note = {Submitted to RECOMB 2024}
}
```

## Contact

Questions and suggestions can be directed to:

- [Luca.Renders@UGent.be](mailto:Luca.Renders@UGent.be)
- [Jan.Fostier@UGent.be](mailto:Jan.Fostier@UGent.be)
