#!/bin/sh

notebook_dir="../notebooks"

# iterate through all ipynb files in directory
for notebook in "$notebook_dir"/*.ipynb; do
    # Check if the file exists
    if [ -f "$notebook" ]; then
    
        jupytext --set-formats ipynb,py:percent "$notebook"

        # log the pairing
        echo "Paired: $notebook with Python file"

    else
        "No notebooks found in $notebook_dir"
    
    fi
done