# Columba 2.0

<!-- TOC depthFrom:4 -->

Fast Approximate Pattern Matching using Search Schemes

Columba is a powerful open-source read-mapper developed to significantly enhance the performance of lossless approximate pattern matching. This README provides an overview of the [key features and benefits](#key-features-and-benefits) of Columba, along with instructions for [installation](#installation), [usage](#usage), [result reproduction](#result-reproduction) and [citation](#citation).
Finally, [contact information](#contact) and [license details](#license-and-dependencies) are provided.

## Key features and benefits

- Columba is a **lossless** read-mapper, meaning that **all** occurrences up to the provided distance or **all-best** occurrences with the minimum identity are found. Both the edit distance and the hamming distance are supported.
- Columba can handle any valid search scheme, so if new search schemes are discovered no update is needed. See [custom search schemes](./further_info/advanced_options/README.md#custom-search-schemes) below.
- Columba supports a dynamic selection of search schemes, where the search scheme best suited for the current read and reference text is chosen, resulting in faster runtimes!
- Columba supports dynamic partitioning of the read, boosting the performance!
- Columba comes in two flavors: Vanilla and RLC. Vanilla is based on the bidirectional FM Index and is faster. RLC is based on the b-move structure and is optimized for memory usage by run-length-compressing the index.
- Columba Vanilla comes with tunable in-text verification, which improves runtimes significantly. More info [here](#advanced-options).
- Columba is **fast**, outperforming other lossless aligners, thanks to all kinds of algorithmic tricks like an interleaved bit-vector representation and redundancy avoidance for the edit distance metric.
- Columba outputs to either SAM format or a custom read hit summary (or a gzipped variation) if the position in the genome is not relevant.
- Columba Vanilla is supported on both Unix and Windows systems!

## Choosing the Right Columba Flavor

When selecting the appropriate Columba flavor for your needs, consider the following:

- **Operating System Compatibility**: If you are using a Windows system, Columba Vanilla is your only option.
- **Memory and Reference Text Considerations**:
  - *Columba Vanilla* loads the entire reference text (or combined reference texts) and the full FM Index into memory. It offers superior speed, making it ideal if memory usage is not a concern or if your reference text is not too large. Loading the reference text allows for the calculation of CIGAR strings.
  - *Columba RLC*, however, employs a run-length compressed index and does not require loading the reference text into memory. This makes it a better choice for handling large pan-genomes or when memory is limited but it comes with slower execution time and no calculation of CIGAR strings.
- **Accuracy**: Columba RLC does not perform in-text verification, which may lead to an overestimation of the edit distance at the edge of certain reference sequences in rare cases.

## Installation

The following instructions will get you a copy of the project up and running on your local machine.

### Prerequisites

This package requires CMake (3.14 or higher) and a compiler.
For Unix systems, the recommended compiler is GCC (tested on 8.3.0 and more recent versions).
For Windows systems, we recommend the use of [Clang installed via MSYS2](https://packages.msys2.org/package/mingw-w64-x86_64-clang).
Click [here](https://www.msys2.org/) for more info on MSYS2.
It is recommended that you also install CMake and Ninja via MSYS2.

**Note for Windows users:** Do not use the mingw64 compiler as long as [this bug](https://github.com/msys2/MINGW-packages/issues/2519) related to thread_local variables is not fixed. Using this compiler will result in a segmentation fault when the processing threads are destroyed.

#### Extra Prerequisite RLC

Columba RLC needs the [SDSL library](https://github.com/simongog/sdsl-lite).
You can follow their install instructions to install it on your system.
If you have installed SDSL to a non-standard location, you can point CMake to the installation location by adding `-DSDSL_INCLUDE_DIR=<path-to-sdsl>/` and `-DSDSL_LIBRARY=<path-to-sdsl-lib>` to the CMake command.

### Installing Columba Vanilla

The installation is now simple. First, clone Columba from the GitHub address and checkout this beta-branch.

```bash
    git clone "https://github.com/biointec/columba.git"
    git checkout columba-2.0-beta
```

From this directory, run the following commands to install Columba:

```bash
mkdir build
cd build
cmake ..
```

For Unix-like systems (Linux, macOS):

```bash
make
```

For Windows:

```bash
ninja
```

### Installing Columba RLC

The installation for Columba RLC is completely analogous to Columba Vanilla ([provided SDSL is already installed](#extra-prerequisite-rlc)). The only difference is that you must add `-DRUN_LENGTH_COMPRESSION=ON`to the CMake command.

### Installing Both Flavors

We provide a script that installs both flavors to two separate directories.
To do this, you must first clone the repository, as detailed above, and then, from the top directory in the repository run `bash buildbothversions.sh`.
This will install columba and columba_build in directories build_vanilla and build_rlc.

### Compiler options

#### Word-size

By default, Columba Vanilla uses 32-bit integers for the suffix array and positions in the text.
However, if your reference genome (or concatenation of genomes) is longer than 4.29 billion characters this will no longer fit.
In that case, you need to let the compiler know that you want 64-bit integers.
You can do this by changing your CMake command.

```bash
mkdir build
cd build
cmake -DTHIRTY_TWO=OFF ..
```

Followed by `make` or `ninja`.

Columba RLC uses 64-bit integers by default.
You can turn this off by adding `-DTHIRTY_TWO=ON`to the CMake command.

---

## Usage

Columba aligns reads to a bidirectional FM Index or the bidirectional move structure (b-move).
To do this you need to build the index based on the input data.
Currently, we only support input data with an alphabet of length 5 (e.g. for DNA A, C, G, T + $).
If your input data has other characters they will be replaced by a random character chosen from A, C, G and T.

### Building the index

To build the bidirectional FM-index from your reference genome (or collection of genomes), you need FASTA files, where each sequence has a unique identifier.
You can use more than one FASTA file, but the identifiers must be unique over all files.
To build the index use this command from the directory where you built Columba (normally this is the build/ directory):

```bash
./columba_build -r <reference_base_name> -f <list_of_FASTA_files>
```

or 

```bash
./columba_build -r <reference_base_name> -F <file_with_paths>
```

where `<list_of_FASTA_files>` is a space separated list of FASTA input files, `<reference_base_name>` is the basename of the index files to write and `<file_with_paths>` is a text file where each line is a path to a FASTA input file.
Note that you can also combine the `-f` and `-F` flags.

**Warning** All paths cannot contain any spaces!

This will create all necessary files in the same directory.
Note that the build process may take some time.

To find info on all options use:

```console
./columba_build --help
```

#### Backward compatibility

Columba's alignment mode is not guaranteed to be compatible with indexes that were built using an older version.
Consider rebuilding your index.

#### Vanilla option: suffix array sparseness

By default in Columba Vanilla, the suffix array is sampled with a sparseness factor of 4.
This has proven to be a good mix between memory and runtime requirements.
However, in Columba Vanilla you can choose which sparseness factor you would like to use (as long as it is a power of two), by adding `-s <factor>`to the build command.

You can also use `-a` instead of `-s` to create all sampled versions with sparseness factors 1 to 128.
If you choose to use another sparseness factor, you should pass it to the aligner with the `-s` flag.
For more information see the info on [Advanced options](./further_info/advanced_options/README.md).

#### RLC options

Columba RLC has two build modes.
The first mode is analogous to Vanilla mode and uses `n*8` more memory than the second mode, with `n` the length of the (concatenation of) reference text(s).
The second mode uses prefix-free parsing and makes use of [Big-BWT](https://gitlab.com/manzai/Big-BWT).
Note that this mode requires python 3.8 (or greater) and package psutil to be installed on your system.
To use the second mode, run

```bash
bash columba_build_pfp.sh -r <reference_base_name> -f <list_of_FASTA_files> [-w <ws>] [-p <mod>]
```

or

```bash
bash columba_build_pfp.sh -r <reference_base_name> -F <file_with_paths> [-w <ws>] [-p <mod>]
```


from your build folder.

- `-w <ws>`: (Optional) Window size for Big-BWT. If this option is **unset**, Big-BWT will use its **default window size**.  
- `-p <mod>`: (Optional) Mod value for Big-BWT. If this option is **unset**, Big-BWT will use its **default mod value**.

**Note**:  
If you encounter hash collisions during Big-BWT processing, try increasing the window size (`-w`) and/or the mod value (`-p`) to resolve the issue.

**Warning**: The PLCP phase of the build algorithm can take a long time without updates on your screen.
This does not mean the program has halted.

##### Seeded replacement

Columba RLC achieves better memory usage by run-length compressing the index.
However, if your text contains blocks of many non-ACGT characters, they will be randomly replaced.
This introduces a lot of randomness which significantly hurts the efficiency of run-length compressing.
To combat this, seeded replacement can be used.
This ensures that a replacement is always taken from a seed.
You can set the seed length by adding `-l <seed_length>` to the build command.

Columba Vanilla does not use seeded replacement by default (`-l 0`).
Columba RLC does use seeded replacement by default (`-l 100`).

### Using the index

This section explains how to use a built index for a reference genome. First, a general overview is given with a minimal example.
Later the options for the aligner are explored.

Columba can align reads in a FASTA (`.FASTA`, `.fasta`, `.fa`, `.fna`) or FASTQ (`.fq`, `.fastq`) format (or .gz variations thereof).
A minimal working example to align your reads is shown below, where `reference_base_name` is the reference base name used while constructing the index, and `read_file`is the path to a file with the reads you want to align.
 
```bash
./columba <options> -r <reference_base_name> -f <read_file>
```

This will align the reads to the reference genome(s) and output the alignments in SAM format to ColumbaOutput.sam.

In the following sections, the optional arguments for Columba are discussed.
After each section, a summary of the options is listed.

To find info on all options use:

```console
./columba --help
```

#### Alignment options

##### Parallelization

Columba can align reads in parallel.
You can choose how many threads are used by using the `-t` or `--threads` option.

```console
-t, --threads           INT   The number of threads to be used. Default is 1.
```

##### Alignment mode

Columba supports two alignment modes: `all` and `best`.

The default mode is `best`.
In this mode, all alignments with the best score (edit or hamming distance) are reported for which the minimum identity between the read and the reference sequence is at least some value.
By default, this value is 95%.
You can set this value by using the `-I` flag.

The other mode is `all` which can be set by adding `-a all` or `--align-mode all` as an option.
In this mode, all alignments up to some edit distance are reported.
By default, this value is 0.
You can set the value by using the `-e` flag.
`all` mode is generally slower than `best` mode.
Note that Columba can currently only handle up to 13 errors per alignment.

Below you can find a summary of the relevant parameters for alignment mode:

```console
  -a, --align-mode        STR   Alignment mode to use. Options are: all, best. Default is best.
  -I, --minIdentity       INT   The minimum identity for alignments in BEST mode. Default is 95.
  -e, --max-distance      STR   The maximum allowed distance (for ALL mode). Default is 0.
```

##### Distance Metric

Columba supports both the Hamming distance and the edit distance.
The Hamming distance will only allow for substitutions, while the edit distance also considers indels.
Note that the edit distance is slower.
You can choose which one to use by using the `-m` or `--metric` option.
By default the edit distance metric is used.

```console
  -m, --metric            STR   Distance metric to use. Options are: edit, hamming. Default is edit.
```

##### Paired end alignment

Columba supports paired-end alignment.
To do paired end alignment run:

```bash
./columba <options> -r base -f <read_file_1> -F <read_file_2>
```

If you have a single interleaved reads file, you can de-interleave them using [this one-line script](https://gist.github.com/nathanhaigh/4544979).

By default, the orientation of the pairs and the fragment size are inferred.
Currently, in a pan-genome context, only pairs of reads that unambiguously to a sequence in the first FASTA file on which the reference is built are considered for inferring the orientation and fragment size.
If this is not the behavior you want, we advise turning inference off using the `-nI` or `--no-inferring`flag
and setting the orientation and minimum and fragment sizes via the `-O` (`--orientation`), `-X` (`--max-insert-size`) and `-N` (`--min-insert-size`) options.

The options for paired-end alignment are listed below:

```console
    -F, --second-reads-file STR     Path to the second reads file (optional). If this is set Columba
                                    will use paired-end alignment.
    -nI, --no-inferring              Do not infer paired-end parameters. Do not infer the paired end
                                    parameters. By default the parameters are inferred. If this option
                                    is set the values provided by -O (default FR), -X (default 500) and 
                                    -N (default 0) are used.
    -O, --orientation       STR     Orientation of the paired end reads. Options are: fr, rf, ff.
                                    Default is fr.
    -X, --max-insert-size   INT     The maximum insert size for paired end reads. Default is 500.
    -N, --min-insert-size   INT     The minimum insert size for paired end reads. Default is 0.
    -D, --discordant        INT     Allow discordant alignments. Optionally you can provide the maximal
                                    number of discordant alignments per pair to allow. Default is 100000.
```

By default, only pairs for which both reads map to the same reference sequence for which the orientation is correct and for which the insert size is between the minimum and maximum insert size (either inferred or given) are reported.
If you are interested in discordant alignment, i.e. alignments with the wrong orientation or insert size or pairs that do not map to the same sequence, you can add the `-D` (`--discordant`) flag.
If no concordant pair can be found the discordant pairs are reported (up to 100 000 per read pair).
Note that this option slows down Columba.

#### Output options

##### SAM format

By default Columba outputs in SAM format to ColumbaOutput.sam.
You can specify another output file by using the `-o` (`--output-file`) option.
By default, each separate alignment of a single read has its own SAM record.
In single-end alignment, you can choose to have all alignments of the same read in one record by using the XA tag.
You can set this by using the `-XA` (`--XA-tag`) flag.

Unmapped reads are given a SAM record with the unmapped flag set.
If you are not interested in unmapped records, you can choose to suppress them with the `-nU` (`--no-unmapped`) flag.
This option can slightly improve Columba's runtime.

In Columba Vanilla, each SAM record has a CIGAR string.
The calculation of this string takes some time.
If you are not interested in the CIGAR string you can suppress its calculation using the `-nC` (`--no-cigar`)flag.
Columba RLC does not output CIGAR strings and instead places an asterisk (*) in the corresponding column.

##### Read Hit Summary (RHS) format

Columba offers a second custom output format called Read Hit Summary (RHS).
This mode can be useful if you are interested in which reference sequences are a hit and not necessarily where in that sequence the match is found.
Currently, this mode is only supported in single-end alignment.

The output file consists of two columns separated by tab characters.
The first column is the read identifier and the second column is a semi-colon-separated list of all the hits.
Each hit is represented as `(refSeq, distance)`, where ``refSeq` is the name of the reference sequence of the hit and distance is the edit or hamming distance to the read.

To use the RHS output format use the `-o` flag to redirect the output to a file with the `rhs` extension.

##### Compression

If you have ZLIB installed on your device, Columba supports the output to gzipped files.
To use this, just use `-o outputFile.sam.gz` or `-o outputFile.rhs.gz`.

##### Logging

Columba logs its progress to the console, you can choose to redirect this to a log file with the `-l` (`--log-file`) option.
The logs contain timestamps and information about Columba's progress.

##### Order of output file

If you run the program with multiple threads, there is no guarantee that the order of reads in your reads file is the same as the order of reads in the output file.
The only guarantee is that hits for the same read(pair) are in a contiguous block.
However, if you want to keep the order, you can use the `-R`flag to force reordering.
Note that this slows down the runtime.

##### Output Options Summary

```console
  -nU, --no-unmapped             Do not output unmapped reads.
  -nC, --no-CIGAR                Do not output CIGAR strings for SAM format. (only in Vanilla)
  -XA, --XA-tag                  Output secondary alignments in XA tag for SAM format.
  -o, --output-file       STR   Path to the output file. Should be .sam or .rhs. 
                                Default is ColumbaOutput.sam.
  -l, --log-file          STR   Path to the log file. Default is stdout.
  -R, --reorder                 Guarantees that output SAM or RHS records are printed in the order
                                corresponding to the order of reads in the original file. Setting
                                this will cause Columba to be somewhat slower and use somewhat more
                                memory.
```

#### Advanced Options

For advanced users who want to experiment with Columba's functionality some other options are available.
Note that the default settings have been carefully selected to bring the best performance.

You can find more info about these settings, as well as on the use of custom search schemes [here](./further-info/advanced_options/README.md).

## Result reproduction

To reproduce the results presented in our papers please refer to [these instructions](./further-info/result_reproduction/README.md).

## Citation

Please cite our work if you find this code useful in your research!

- [Renders, Luca, Kathleen Marchal, and Jan Fostier. "Dynamic partitioning of search patterns for approximate pattern matching using search schemes." Iscience 24.7 (2021).](https://doi.org/10.1016/j.isci.2021.102687)
- [Renders, Luca, Lore Depuydt, and Jan Fostier. "Approximate pattern matching using search schemes and in-text verification." International Work-Conference on Bioinformatics and Biomedical Engineering. Cham: Springer International Publishing, 2022.](https://doi.org/10.1007/978-3-031-07802-6_36)
- [Renders, Luca, Lore Depuydt, Sven Rahmann, and Jan Fostier. "Automated Design of Efficient Search Schemes for Lossless Approximate Pattern Matching." International Conference on Research in Computational Molecular Biology. Cham: Springer Nature Switzerland, 2024.](https://doi.org/10.1007/978-1-0716-3989-4_11)

If you use Columba RLC, please also cite the b-move paper:

- [Lore Depuydt, Luca Renders, Simon Van de Vyver, Lennart Veys, Travis Gagie, and Jan Fostier. "b-move: faster bidirectional character extensions in a run-length compressed index." 24th International Workshop on Algorithms in Bioinformatics (WABI 2024). Leibniz International Proceedings in Informatics (LIPIcs), 2024.](https://doi.org/10.4230/LIPIcs.WABI.2024.10)

## Contact

Questions and suggestions can be directed to:

- <Luca.Renders@UGent.be>
- <Lore.Depuydt@UGent.be> (regarding the RLC version)
- <Jan.Fostier@UGent.be>

## License and Dependencies

See the license file for information about Columba's license.
Columba makes use of the [{fmt} library](https://github.com/fmtlib/fmt) and falls under the exception of its license.
Columba also makes use of the [libsais](https://github.com/IlyaGrebnov/libsais) and [parallel-hashmap](https://github.com/greg7mdp/parallel-hashmap) libraries. Both fall under Apache-2.0 license, which is included in the repository.
Columba RLC with prefix-free parsing makes use of [Big-BWT](https://gitlab.com/manzai/Big-BWT).
