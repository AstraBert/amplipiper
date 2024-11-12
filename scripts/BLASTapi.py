
import time
from Bio import SeqIO, Entrez, Blast
from Bio.Blast import NCBIXML
from argparse import ArgumentParser
import gzip
import os
from collections import defaultdict as d


def parse_arguments():
    """Parse command line arguments."""
    parser = ArgumentParser(
        description="BLAST sequence analyzer with taxonomy ID retrieval")
    parser.add_argument(
        "-i", "--infile", help="Path to input file with sequences", required=True)
    parser.add_argument(
        "-o", "--outfile", help="Path to output file for BLAST results", required=True)
    parser.add_argument(
        "-e", "--email", help="Email for Entrez (required for NCBI access)", required=True)
    return parser.parse_args()


def load_data(infile):
    """Load sequences from the input file, unzipping if necessary."""
    if infile.endswith(".gz"):
        with gzip.open(infile, "rt", encoding="latin-1") as file:
            if infile.endswith((".fasta.gz", ".fna.gz", ".fas.gz", ".fa.gz")):
                return {record.id: str(record.seq) for record in SeqIO.parse(file, "fasta")}
            else:
                raise ValueError("Unsupported file format for gzipped file")
    elif infile.endswith((".fasta", ".fna", ".fas", ".fa")):
        with open(infile, "r") as file:
            return {record.id: str(record.seq) for record in SeqIO.parse(file, "fasta")}
    else:
        raise ValueError("Unsupported file format")


def get_taxonomy_id(hit_id):
    """Fetch the genus and species name from NCBI for a given hit ID."""
    try:
        # Extract just the accession part if needed (e.g., remove gi| and any prefixes)
        clean_id = hit_id.split('|')[-2] if '|' in hit_id else hit_id

        handle = Entrez.efetch(db="nucleotide", id=clean_id,
                               rettype="gb", retmode="text")
        for record in SeqIO.parse(handle, "genbank"):
            if record.annotations.get("organism"):
                # Extract full genus and species from the organism field
                organism_name = record.annotations["organism"]
                parts = organism_name.split()
                if len(parts) >= 2:
                    genus = parts[0]
                    species = parts[1]
                    return f"{genus} {species}"
                else:
                    return organism_name  # Return what is available
        handle.close()
    except Exception as e:
        print(f"Failed to retrieve taxonomy ID for {hit_id}: {e}")
        return "N/A"
    finally:
        time.sleep(0.5)  # Introduce a delay of 0.5 seconds between requests


def run_blast_and_save_results(infile, outfile):
    """Run BLAST for sequences in the input file and save results to the output file."""
    print("Loading data...")
    print("Data loaded successfully.")

    result_file_path = os.path.join(outfile)

    print("Starting BLAST for full sequence")
    result_handle = Blast.qblast("blastn", "nt", open(infile).read())

    data = d(lambda: d(list))

    blast_full = NCBIXML.parse(result_handle)
    with open(result_file_path, "a") as result_file:
        print("Parsing BLAST results")
        for blast_record in blast_full:
            C = 1
            for alignment in blast_record.alignments:
                # only use 10 best Hits
                if C > 10:
                    break
                C += 1
                # get the taxonomic information for the hit using NCBI Entrez
                taxonomy_id = get_taxonomy_id(alignment.hit_id)
                HSP = []
                for hsp in alignment.hsps:
                    perc_identity = 100 * \
                        round(hsp.identities / hsp.align_length, 4)
                    HSP.append(perc_identity)
                data[blast_record.query][taxonomy_id].append(max(HSP))

    result_handle.close()
    print("Finished parsing results")
    print("Writing BLAST results")

    result_file = open(result_file_path, "w")
    # Processing data to count occurrences and find the highest percentage for each result
    result_dict = {}
    for test, results in data.items():
        result_dict[test] = {}
        for result, percentages in results.items():
            count = len(percentages)
            highest_percentage = round(max(percentages), 2)
            result_dict[test][result] = {
                'count': count,
                'highest_percentage': highest_percentage
            }

    # Find the maximum number of results across all tests for formatting
    max_results_per_test = max(len(results)
                               for results in result_dict.values())

    # Print header
    header = ["SAMPLE"] + \
        [f"Taxon{i + 1} (sim%)" for i in range(max_results_per_test)]
    result_file.write(",".join(header) + "\n")

    # Print results for each test
    for test, results in result_dict.items():
        sorted_results = sorted(
            results.items(), key=lambda x: x[1]['highest_percentage'], reverse=True)
        row = [test]  # Start row with the test/sample name
        for i, (result, info) in enumerate(sorted_results):
            # Create the formatted string for each result
            entry = f"{result} ({info['highest_percentage']}%); count = {info['count']}"
            row.append(entry)

        # Fill the rest of the row with empty columns if needed
        row += [""] * (max_results_per_test - len(sorted_results))

        # Join the row without adding a comma after the last entry
        formatted_row = ",".join(row)

        # Write the formatted row to the file
        result_file.write(formatted_row + "\n")

    result_file.close()
    print("Finished writing results")


if __name__ == "__main__":
    args = parse_arguments()
    Entrez.email = args.email  # Set the email for Entrez
    run_blast_and_save_results(args.infile, args.outfile+"/final.csv")
