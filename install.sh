#!/bin/bash

# Stop execution on error
set -e

# Step 1: Clean any previous builds
echo "Cleaning previous builds..."
rm -rf build/ dist/ src.egg-info/

# Step 2: Build the Cython extension
python setup.py build_ext --inplace



# Step 3: Install the package system-wide
echo "Installing package system-wide..."
python -m pip install -e .

# Done
echo "Installation complete!"
