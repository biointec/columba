# Columba
Fast Approximate Pattern Matching using Search Schemes

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. 

#  Prerequisites

This package requires a number of packages to be install on your system. Required: CMake; Google's Sparsehash; gcc (GCC **TODO check gcc version** or a more recent version) 

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

The installation is now simple. First, clone browniealigner from the Github address

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
We recommend the use of [radixSA64](https://github.com/mariusmni/radixSA64) for building these arrays.
Name these arrays respectively [basefile].sa and [basefile].rev.sa, where your text is named [basefile].txt
To build the FM-index run the following command in the `build` folder. 
```bash
./columba-build [basefile]
```
### Example
Consider the following directory structure. Here, `example.txt` is the file you want to build the index for. You already created `example.sa` and `example.rev.sa`.

    .
    ├── build
    ├── example
    |   ├── example.sa
    |   ├── example.rev.sa
    |   └── example.txt
    └── src

To build the index, navigate to the `build` folder and run the following command:
```bash
./columba-build ../example/example
```
The index files are then written to the same folder.

## Using the index
Columba can align reads in a fasta (`.FASTA`, `.fasta`, `.fa`) or fastq (`.fq`, `.fastq`) format. 
To align your reads use the following format 

```bash
./columba [options] basefile readfile
```

options:

```
-e  --max-ed		maximum edit distance [default = 0]
-s  --sa-sparseness	suffix array sparseness factor [default = 1]
-st  --static	Add flag to do static partitioning [default = false]
-ss --search-scheme   Choose the search scheme [default = kuch1]
  options:
        kuch1   Kucherov k + 1
        kuch2   Kucherov k + 2
        kianfar  Optimal Kianfar scheme
        manbest  Manual best improvement for kianfar scheme (only for ed = 4)
        pigeon   Pigeon hole scheme
        01*0    01*0 search scheme
    
```

The number of nodes, number of matrix elements, duration, and number of matches will be printed to stdout. 
The matches will be written to a custom output file in the folder where your readfile was. This output file will be a tab-seperated file with the fields: `identifier`, `position`, `length` and `ED`. For each optimal alignment under the maximal given edit distance a line will be present. This output file will be called `readfile_output.txt`.







