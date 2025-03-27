#!/bin/bash

# Script: columba_build_pfp.sh
# Description: This script performs the Columba build process for PFP.
# Author: Lore Depuydt - lore.depuydt@ugent.be

# Capture start time
start_time=$(date +%s)

# We assume that this script is run from the build folder
columba_build_exe="./columba_build"
big_bwt_exe="./../external/Big-BWT/bigbwt"

# Default seed length
seedLength=100

# Optional Big-BWT parameters
ws=0  # Default window size (unset means 0)
mod=0 # Default mod value (unset means 0)

# Array to store fasta files
fasta_files=()

# Function to show usage
showUsage() {
	echo "Usage: $0 [-l <seedLength>] [-w <ws>] [-p <mod>] -r <index_name> [-f <fasta_files>] [-F <fasta_file_list>]"
	echo
	echo "Required arguments:"
	echo "  -r <index_name>   Name/location of the index to be created."
	echo
	echo "Optional arguments:"
	echo "  -f <fasta_files>  Space-separated list of FASTA files."
	echo "  -F <fasta_file_list>  Path to a file containing a list of FASTA files, one per line."
	echo "  -l <seedLength>   Seed length for replacing non-ACGT characters (default: $seedLength). 0 means that no seed is used."
	echo "  -w <ws>           Window size for Big-BWT. If unset, Big-BWT will use its default window size."
	echo "  -p <mod>          Mod value for Big-BWT. If unset, Big-BWT will use its default mod value."
}

# Function to run a command with /usr/bin/time -v and extract time and memory usage
# Usage: runCommandWithTime <command> [<args>...]
runCommandWithTime() {
	local command="$1"
	shift
	# Run the command and capture output, while measuring time and memory usage
	(/usr/bin/time -v "$command" "$@") || {
		local status=$?
		echo "Error: Command '$command $@' failed with exit status $status." >&2
		exit $status
	}
}

# Function to parse command-line options
parseOptions() {
	# Parse command-line options
	while getopts ":l:r:f:F:w:p:" opt; do
		case $opt in
		l)
			seedLength=$OPTARG
			;;
		r)
			index_name=$OPTARG
			;;
		f)
			# Collect all subsequent arguments as the list of FASTA files
			fasta_files+=("$OPTARG")
			while [[ $OPTIND -le $# && ! ${!OPTIND} =~ ^- ]]; do
				fasta_files+=("${!OPTIND}")
				OPTIND=$((OPTIND + 1))
			done
			;;
		F)
			# Read each line in the specified file and add it to the fasta_files array
			if [[ -f $OPTARG ]]; then
				while IFS= read -r line; do
					fasta_files+=("$line")
				done <"$OPTARG"
			else
				echo "Error: File '$OPTARG' not found." >&2
				exit 1
			fi
			;;
		w)
			ws=$OPTARG
			;;
		p)
			mod=$OPTARG
			;;
		\?)
			echo "Invalid option: -$OPTARG" >&2
			showUsage
			exit 1
			;;
		:)
			echo "Option -$OPTARG requires an argument." >&2
			showUsage
			exit 1
			;;
		esac
	done
	# Shift off the options and optional --
	shift $((OPTIND - 1))

	# Ensure required arguments are provided
	if [ -z "$index_name" ] || [ "${#fasta_files[@]}" -eq 0 ]; then
		showUsage
		exit 1
	fi
}

# Main script logic

# Parse command-line options
parseOptions "$@"

echo "Welcome to the Columba build process with prefix-free parsing!"
echo "-------------------------------------------------------------"
echo "Index name: $index_name"
echo "Input FASTA files:"
for file in "${fasta_files[@]}"; do
	echo "  $file"
done
echo "Seed length: $seedLength"
echo "Big-BWT window size: ${ws:-not set}"
echo "Big-BWT mod value: ${mod:-not set}"
echo "-------------------------------------------------------------"

# Start the preprocessing
echo "Start preprocessing the fasta file(s) with Columba..."
runCommandWithTime "$columba_build_exe" --preprocess -l "$seedLength" -r "$index_name" -f "${fasta_files[@]}"
echo "Preprocessing done!"
echo "-------------------------------------------------------------"

base="${index_name}"

# Build Big-BWT command with optional -w and -p arguments
big_bwt_args=("$big_bwt_exe" -e -s -v "$base")
if [ "$ws" -gt 0 ]; then
	big_bwt_args+=(-w "$ws")
fi
if [ "$mod" -gt 0 ]; then
	big_bwt_args+=(-p "$mod")
fi

# Start the prefix-free parsing
echo "Start prefix-free parsing for the original string..."
runCommandWithTime "${big_bwt_args[@]}"
echo "Prefix-free parsing done!"
echo "-------------------------------------------------------------"

# Adjust the arguments for the reverse string
big_bwt_args=("$big_bwt_exe" -e -s -v "${base}.rev")
if [ "$ws" -gt 0 ]; then
	big_bwt_args+=(-w "$ws")
fi
if [ "$mod" -gt 0 ]; then
	big_bwt_args+=(-p "$mod")
fi

echo "Start prefix-free parsing for the reverse string..."
runCommandWithTime "${big_bwt_args[@]}"
echo "Prefix-free parsing done!"
echo "-------------------------------------------------------------"

# Start building the Columba index
echo "Start building the Columba index..."
runCommandWithTime "$columba_build_exe" --pfp -r "$index_name"
echo "Columba index built!"
echo "-------------------------------------------------------------"

# Remove the temporary files
echo "Remove temporary files..."
rm "${base}.bwt"
rm "${base}.rev.bwt"
rm "${base}.ssa"
rm "${base}.rev.ssa"
rm "${base}.esa"
rm "${base}.rev.esa"
rm "${base}.log"
rm "${base}.rev.log"
rm "${base}"
rm "${base}.rev"
echo "Temporary files removed!"
echo "-------------------------------------------------------------"

# Capture end time
end_time=$(date +%s)

# Calculate total elapsed time
total_time=$((end_time - start_time))

echo "The Columba build process with prefix-free parsing is finished!"
echo "Total time elapsed: $total_time seconds."
