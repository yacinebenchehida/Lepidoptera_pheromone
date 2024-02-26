#!/usr/bin/env python

import subprocess
import os

# Define the function to process each row
def process_row(species, genome_url, gff3_url, cds_url):
	# Replace spaces in species name with underscores
	species_name = species.replace(' ', '_')

	# Define the output directory where files will be saved
	output_directory = f"/mnt/scratch/projects/biol-specgen-2018/yacine/Pheromones/Blast/Inputs/{species_name}/"

	# Print information for debugging (optional)
	print(f"Processing species: {species_name}")
	print(f"Genome URL: {genome_url}")
	print(f"GFF3 URL: {gff3_url}")
	print(f"CDS URL: {cds_url}")

	# Create the output directory if it doesn't exist
	os.makedirs(output_directory, exist_ok=True)

	# Here you can execute your commands using the extracted information
	# For example:
	# Download genome, GFF3, and CDS files
	subprocess.run(['wget', genome_url, '-O', f"{output_directory}{species_name}-genome.fa.gz"])
	subprocess.run(['wget', gff3_url, '-O', f"{output_directory}{species_name}.gff3.gz"])
	subprocess.run(['wget', cds_url, '-O', f"{output_directory}{species_name}-cds.fa.gz"])

	# Unzip the downloaded files
	subprocess.run(['gunzip', f"{output_directory}{species_name}-genome.fa.gz"])
	subprocess.run(['gunzip', f"{output_directory}{species_name}.gff3.gz"])
	subprocess.run(['gunzip', f"{output_directory}{species_name}-cds.fa.gz"])

# Open the input file and process each line
with open('data.txt', 'r') as file:
	for line in file:
		# Split the line into columns
		columns = line.strip().split('\t')
		# Extract information from each column
		species, genome_url, gff3_url, cds_url = columns
		# Process the row
		process_row(species, genome_url, gff3_url, cds_url)
