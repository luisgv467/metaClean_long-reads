#!/bin/bash

# Function to display usage instructions
usage() {
  echo "Usage: $0 -i <input_file_or_directory> -o <output_file>"
  exit 1
}

# Parse command-line arguments
while getopts ":i:o:" opt; do
  case $opt in
    i) input="$OPTARG" ;;
    o) output_file="$OPTARG" ;;
    *) usage ;;
  esac
done

# Ensure both -i and -o are provided
if [ -z "$input" ] || [ -z "$output_file" ]; then
  usage
fi

# Check if input exists
if [ ! -e "$input" ]; then
  echo "Error: Input '$input' does not exist."
  exit 1
fi

# Create the output file if it doesn't exist
if [ ! -f "$output_file" ]; then
  touch "$output_file"
fi

# Process the input
if [ -d "$input" ]; then
  # Loop through each FASTQ file in the directory
  for file in "$input"/*.fastq.gz; do
    if [ -f "$file" ]; then
      # Count the number of reads
      d=$(zcat "$file" | wc -l | awk '{print $1/4}')
      # Append the results to the output file
      echo "$(basename "$file") $d" >> "$output_file"
    fi
  done
elif [ -f "$input" ]; then
  # Count the number of reads in the single file
  d=$(zcat "$input" | wc -l | awk '{print $1/4}')
  # Append the result to the output file
  echo "$(basename "$input") $d" >> "$output_file"
else
  echo "Error: '$input' is neither a valid file nor directory."
  exit 1
fi

echo "Read counting complete. Results saved in $output_file."
