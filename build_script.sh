#!/bin/bash

# Check if input is valid
if [ $# -eq 1 ] || [ $# -eq 2 ] || [ $# -eq 3 ]; then
	FLAVOR=$1
	if [ "$FLAVOR" != "Vanilla" ] && [ "$FLAVOR" != "RLC" ]; then
		echo "Error: Invalid flavor specified. Please choose 'Vanilla' or 'RLC' as the first argument."
		exit 1
	fi

	BIT_WIDTH=""
	PHI_MOVE_FLAG=""

	# Parse the remaining arguments flexibly
	for arg in "${@:2}"; do
		if [ "$arg" == "32" ] || [ "$arg" == "64" ]; then
			if [ -n "$BIT_WIDTH" ]; then
				echo "Error: Multiple bit-width values provided. Only one of '32' or '64' is allowed."
				exit 1
			fi
			BIT_WIDTH=$arg
		elif [ "$arg" == "PHIMOVE" ]; then
			if [ -n "$PHI_MOVE_FLAG" ]; then
				echo "Error: 'PHIMOVE' flag specified multiple times."
				exit 1
			fi
			PHI_MOVE_FLAG=$arg
		else
			echo "Error: Invalid argument '$arg'. Expected '32', '64', or 'PHIMOVE'."
			exit 1
		fi
	done

	# Note if PHIMOVE is used with Vanilla
	if [ "$FLAVOR" == "Vanilla" ] && [ -n "$PHI_MOVE_FLAG" ]; then
		echo "Note: Ignoring 'PHIMOVE' flag as it is only supported with the 'RLC' flavor."
		unset PHI_MOVE_FLAG
	fi
else
	echo "Error: Incorrect number of arguments provided."
	echo "Usage: ./build_script.sh <flavor> [bit-width] [PHIMOVE]"
	echo "   <flavor>      : Required. Choose 'Vanilla' or 'RLC'."
	echo "   [bit-width]   : Optional. Specify '32' or '64' (default is 32-bit for Vanilla, 64-bit for RLC)."
	echo "   [PHIMOVE]     : Optional. Specify 'PHIMOVE' as a flag to enable PHI move structures (only supported if flavor is 'RLC')."
	exit 1
fi

# Set options based on provided arguments
option_rlc="-DRUN_LENGTH_COMPRESSION="
if [ "$FLAVOR" == "RLC" ]; then
	option_rlc+="ON"
else
	option_rlc+="OFF"
fi

cmake_command="cmake .. $option_rlc"

# Add bit-width if provided
if [ -n "$BIT_WIDTH" ]; then
	option_32="-DTHIRTY_TWO="
	if [ "$BIT_WIDTH" == "32" ]; then
		option_32+="ON"
	else
		option_32+="OFF"
	fi
	cmake_command="$cmake_command $option_32"
fi

# Add PHIMOVE option if provided and FLAVOR is RLC
if [ "$FLAVOR" == "RLC" ] && [ -n "$PHI_MOVE_FLAG" ]; then
	option_phi="-DPHI_MOVE=ON"
	cmake_command="$cmake_command $option_phi"
fi

# Create build directory
dir="build_$FLAVOR"
if [ -n "$BIT_WIDTH" ]; then
	dir+="_${BIT_WIDTH}"
fi
if [ -n "$PHI_MOVE_FLAG" ]; then
	dir+="_PHIMOVE"
fi

echo "Building Columba $FLAVOR in $dir..."

# Create the directory if it doesn't exist
mkdir -p "$dir"

# Change into the directory
cd "$dir" || {
	echo "Error: Failed to change to directory $dir. Check permissions or directory path."
	exit 1
}

# Run the cmake command, and check for errors
if ! $cmake_command; then
	echo "Error: CMake configuration failed for Columba $FLAVOR in directory $dir. Exiting."
	exit 1
fi

# Run the make command, and check for errors
if ! make; then
	echo "Error: Make process failed for Columba $FLAVOR in directory $dir. Exiting."
	exit 1
fi

# Return to the previous directory
cd .. || {
	echo "Error: Failed to change back to the original directory. Check permissions or directory path."
	exit 1
}

# If everything succeeded
echo "Build complete. Executables for Columba $FLAVOR are in $dir."
