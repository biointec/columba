# Example instructions

This document quickly explains how to run Columba 2.0 on a test dataset. The test dataset includes the following files:

- **Reference Files:**  
  `HG00096.1.19.fasta`
  `HG00097.1.19.fasta`

- **Paired-End Reads Files:**  
  `SRR17981962_1_sampled.fastq`  
  `SRR17981962_2_sampled.fastq`

You can download a zip archive containing all four files from the following link:  
[Download Test Dataset (zip)](https://zenodo.org/records/15090246/files/test_data.zip)

The reference files are sequences that are derived from the **1000 Genomes Project** and represent two haplotypes of chromosome 19 from human individuals HG00096 and HG00097.

The reads files contain 100k paired-end reads sampled from a larger run (accession no. SRR17981962).

---

## Setup Instructions

### 1. Clone the Repository

```bash
git clone "https://github.com/biointec/columba.git"
cd columba
```

### 2. Build Columba

Columba has two build options: **Vanilla** and **RLC**. The interface is identical for both, but the executables are stored in different directories.
For this example you can try both options or [pick the one you want to use](../README.md/#choosing-the-right-columba-flavor).  
You can follow their install instructions to install it on your system.

#### **Vanilla Build**

Columba Vanilla uses 32-bit numbers by default.

```bash
bash build_script.sh Vanilla
```

Build directory: `build_Vanilla`

#### **RLC Build** (requires SDSL-lite, not supported on Windows)

As our reference file is less than 4.29 billion characters long, we will let Columba RLC know to use 32-bit numbers.
Note that Columba RLC needs the [SDSL-lite library](https://github.com/simongog/sdsl-lite) installed on your system.
You can install this library using [their install instructions](https://github.com/simongog/sdsl-lite//blob/master/README.md#installation).

```bash
bash build_script.sh RLC 32
```

Build directory: `build_RLC_32`

---

## Indexing the Reference Genome

First go the respective build directory (`build_Vanilla` or `build_RLC_32`):

Download and unzip the test data in this build directory.

```bash
wget https://zenodo.org/records/15090246/files/test_data.zip
unzip test_data.zip
```

Run the following command from the build directory.

```bash
./columba_build -r chrom19_index -f reference_files/HG00096.1.19.fasta reference_files/HG00097.1.19.fasta
```

This creates index files for the two reference files with the base name `chrom19_index`in our current directory.
You can opt to create the index files in a different existing directory by providing /path/to/directory/name. This will create index files with base name `name`in `/path/to/directory`.

---

## Running Columba on the Test Dataset

The following commands are all executed from the respective build directory (`build_Vanilla` or `build_RLC_32`):
By default the alignments are written to `ColumbaOutput.sam` in this directory. You can redirect the output using the `-o` flag.

### **Single-End All-Best Matches (Default, Identity ≥ 95%)**

```bash
./columba -r chrom19_index -f read_files/SRR17981962_1_sampled.fastq
```

Reports all best matches with at least 95% identity for each read in the reads file.

#### **Single-End All Matches (-a all flag, Allowing up to 2 Edit Distance Errors with -e flag)**

```bash
./columba -r chrom19_index -f read_files/SRR17981962_1_sampled.fastq -a all -e 2
```

Reports all occurrences of each read with up to 2 errors

### **Paired-End All-Best Matches (Default, Identity ≥ 95%)**

```bash
./columba -r chrom19_index -f read_files/SRR17981962_1_sampled.fastq -F read_files/SRR17981962_2_sampled.fastq
```

Reports all best paired matches with at least 95% identity for paired-end reads. If no paired match could be found the reads are treated as single end reads. )

#### **Paired-End All Matches (-a all flag, Allowing up to 2 Edit Distance Errors with -e flag)**

```bash
./columba -r chrom19_index -f SRR17981962_1_sampled.fastq -F SRR17981962_2_sampled.fastq -a all -e 2
```

Reports all paired matches with at most two errors. If no valid paired match could be found the reads are treated as single-end reads.

## Other settings

Feel free to adjust any parameters (e.g., thread count using the `-t`flag) according to your specific needs.
For more information on the parameters use:

```bash
./columba_build --help
./columba --help
```
