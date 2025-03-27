#!/bin/bash

# Build the Vanilla version
bash build_script.sh Vanilla

# Build the RLC version
bash build_script.sh RLC

echo "Build complete. Executables are in build_Vanilla and build_RLC."
