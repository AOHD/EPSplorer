# Set the directory
directory="./genomes"

# Check if the directory exists
if [ ! -d "$directory" ]; then
  echo "Error: Directory '$directory' not found."
  exit 1
fi

# Create the TSV file
tsv_file="magstats.tsv"
echo -e "bin\tmidas4_tax" > "$tsv_file"

# Iterate over files in the directory
for file in "$directory"/*; do
# Get the file name without extension
filename=$(basename "$file" .fasta)
# Append to the TSV file
echo -e "$filename\t$filename" >> "$tsv_file"
done

echo "TSV file created: $tsv_file"