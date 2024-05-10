#!/bin/bash

# Define the input and output directories
input_folder="../input"  # Replace with your actual input folder path
output_folder="../output"  # Replace with your actual output folder path

# Create the output folder if it doesn't exist
mkdir -p "$output_folder"

# Define the output file path
output_file="${output_folder}/merged_results.tsv"

# Initialize a variable to count the number of files processed
file_count=0

# Check if output file already exists and remove it
if [ -f "$output_file" ]; then
    rm "$output_file"
fi

# Find and process each file that matches the pattern
while IFS= read -r file; do
    # Increment file count
    ((file_count++))

    # For the first file, copy as is to keep the header
    if [ $file_count -eq 1 ]; then
        cat "$file" > "$output_file"
    else
        # For subsequent files, skip the header (assumes header is one line)
        tail -n +2 "$file" >> "$output_file"
    fi
done < <(find "$input_folder" -name "*_TSMR_result.tsv")

# Check if any files were found and processed
if [ $file_count -eq 0 ]; then
    echo "No files found to merge."
else
    # Report the number of files processed
    echo "Number of files merged: $file_count"
fi
