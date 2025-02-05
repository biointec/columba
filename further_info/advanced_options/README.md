# Advanced Options For Alignment With Columba

This README provides detailed instructions for advanced users who want to customize and optimize the functionality of Columba.
It explains how to adjust key parameters such as sparseness factors, k-mer sizes, and in-text verification settings.
The guide also covers how to create and integrate custom search schemes, including file structure requirements, validation steps, and static and dynamic partitioning strategies.
Additionally, it explains how to manage dynamic selection using multiple search schemes, with examples of folder structures and configurations for optimal use of Columba’s search capabilities.
Finally, this guide introduces the optional $\phi$ move tables, a performance optimization available for the RLC flavor, explaining when and how to enable them for improved alignment performance.

## Adjusting Parameters

A special note for those users of Columba Vanilla who have chosen to [build their index with a different sparseness factor](../../README.md#vanilla-option-suffix-array-sparseness).
They must now set that sparseness factor with the `-s` flag.
If you have chosen a particularly large sparseness factor, then it might be a good idea to decrease the in-text verification parameter using the `-i` flag.

If the length of your reads divided by the maximum number of errors +1 is nearly 10 it is a good idea to use a smaller value for the `-K` flag.

```console
  -p, --partitioning          STR   Partitioning strategy to use. Options are: uniform, static,
                                    dynamic. Default is dynamic.
  -K, --kmer-size             INT   The size of k-mers in the hash table (used as seeds during
                                    partitioning). Default is 10.
  -S, --search-scheme         STR   Search scheme to use. Options are: kuch1, kuch2, kianfar,
                                    pigeon, 01*0, custom, naive, multiple, minU, columba. Default is
                                    columba.
  -c, --custom                STR   Path to custom search scheme (overrides default search scheme).
  -nD, --no-dynamic-selection       Do not use dynamic selection with custom search scheme.
  -d, --dynamic-selection     STR   Path to custom search scheme with dynamic selection (overrides
                                    default search scheme).
  -i, --in-text               INT   In-text verification switch point. Should be a positive integer.
                                    Default is 5.
  -s, --sa-sparseness         INT   Sparseness factor to use. Should be an integer in 2^[0, 8].
                                    Default is 4. The suffix array with this sparseness factor must
                                    have been constructed.
```

## Custom search schemes

Columba can handle any valid search scheme.
The following section details how to point Columba to your search scheme.

In the [search_schemes](./search_schemes/) folder, several example search schemes are already present.

To make a custom search scheme you need to create a folder containing at least a file called `name.txt`, which contains the name of your search scheme on the first line.
For every maximum edit/hamming distance a subfolder should be present, which contains at least the file `searches.txt`.
In this file, the searches of your scheme are written line by line.
Each line contains three space-separated arrays: $\pi$, L and U.
Each array is written between curly braces {} and the values are comma-separated.

---
`Note`

The $\pi$ array should be zero-based! The connectivity property should always be satisfied.
The L and U arrays cannot decrease.
All error distributions should be covered!
To check if your search scheme is valid you can use a Python script provided [in this directory](../../validitychecker/).

---

By default dynamic selection between the provided custom search scheme and its symmetric variant (reverse mapping the $\pi$-arrays) is performed. You can turn this off by using the `-nD` flag.

### Static Partitioning

If you want to provide optimal static partitioning you can create a file named `static_partitioning.txt` in the folder of the maximum edit/hamming distance this partitioning is for.
This file should contain one line with percentages (values between 0 and 1) separated by spaces.
The ith percentage corresponds to the starting position (relative to the size of the pattern) of the (i + 1)th part (again this is zero-based).
The starting position of the first part is always zero and should **not** be provided.

### Dynamic Partitioning

Similarly, to provide values for dynamic partitioning you can create a file called `dynamic_partitioning.txt`.
his file should contain two lines.
The first line consists of percentages (again between 0 and 1) that correspond to the seeding positions, relative to the size of the pattern, of all parts, except the first and last part.
The second line should contain space-separated integers corresponding to the weights of each part.

### Example

Consider a search scheme that supports maximal edit/hamming distances 1, 2 and 4.
For distance 1 no static or dynamic partitioning values are known.
For distance 2 only static partitioning values are known and for distance 4 both static and dynamic partitioning values are known.
The folder structure of this search scheme should look like this:

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

### More Examples

In the [search_schemes directory](../../search_schemes/), several search schemes are available as custom search schemes.

## Dynamic Selection - Multiple Search Schemes

To use dynamic selection you can provide a set of search schemes to be used.
For this, you need to create a folder with a file called `name.txt` in it and for all values of `k`, where you want to provide multiple search schemes you need to create a directory.
The files in this subdirectory are named `scheme<x>.txt`, where `<x>` is the number of the search scheme (starting from 1).
Each file contains the searches from one of the search schemes and the file needs to be structured like the `searches.txt` file of a custom search strategy (see [above](#custom-search-schemes)).

### Example Folder Structure

Consider a strategy where multiple search schemes are used that support maximal edit/hamming distances  2, 4 and 6.
For these values of `k`, there exist respectively 2, 3 and 4 search schemes.
The folder structure then should look like this:

``` .
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

### More examples

You can find more examples in [this directory](../../search_schemes/multiple_opt/).

## PHIMOVE Compiler Option (for RLC Flavor Only)

The `PHIMOVE` compiler option enables the use of $\phi$ move tables, a performance optimization for the RLC flavor. These move tables are similar to the LF move tables, which are always used in RLC. While $\phi$ move tables require slightly more memory, they can improve alignment performance, particularly when many occurrences are expected to be found.  

This option is **only applicable for the RLC flavor**.  

By default, `PHIMOVE` is **disabled**. If you would like to enable it for the RLC flavor, simply add the `PHIMOVE` flag as shown below:  

```bash
bash build_script.sh RLC 64 PHIMOVE
```  

or  

```bash
bash build_script.sh RLC 32 PHIMOVE
```  

Alternatively, you can enable $\phi$ move tables when compiling the code directly by adding the `-DPHI_MOVE=ON` flag to CMake:  

```bash
cmake -DRUN_LENGTH_COMPRESSION=ON -DPHI_MOVE=ON ..
```  

If you try to use the `PHIMOVE` flag with the Vanilla flavor, it will be ignored with a warning because $\phi$ move tables are exclusively designed for the RLC flavor.

