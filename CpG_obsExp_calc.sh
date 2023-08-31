#! /bin/bash
# JRA 2022

# Check the number of command line arguments
if [ $# -ne 1 ]; then
        script_name=$(basename $0)
        echo "Usage: $script_name file.txt"
        exit 1
fi
