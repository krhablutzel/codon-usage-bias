import os
import requests
from bs4 import BeautifulSoup
import gzip
import directory_scan as ds
import fasta_parser as fp
from Bio import SeqIO
from Bio.Seq import Seq


def count_codons(seq):
    """returns dict w/ counts of each codon"""
    if len(seq) % 3 != 0:
        print(f'Error: seq len {len(seq)}. sequence must have length multiple of 3.')
        return -1

    codons = {}
    for i in range(0, len(seq), 3):
        codon = str(seq[i:i+3])
        if codon in codons:
            codons[codon] += 1
        else:
            codons[codon] = 1

    return codons


def convert_to_fraction(codons):
    """convert each codon count to a fraction of total occurrences for an amino acid"""
    # excludes codons with an unknown base (N)
    uchart = {"Phe": ["UUU", "UUC"],
              "Leu": ["UUA", "UUG", "CUU", "CUC", "CUA", "CUG"],
              "Met": ["AUG"]}

    freqs = {}

    for aa in uchart:
        # count occurrences of aa
        count = 0
        for codon in uchart[aa]:
            if codon in codons:
                count += codons[codon]

        # insert relative codon freq into freq table
        for codon in uchart[aa]:
            if codon in codons:
                freqs[codon] = codons[codon] / count
            else:
                freqs[codon] = 0  # codon never used

    return freqs


def calculate_CUB(data_folder, spec_list):
    for i in range(3):
        species = spec_list[i]
        print(species)

        # combine the fasta into one long sequence
        dna = fp.splice_fasta(f'{data_folder}/{species}.fna')

        # transcribe to RNA
        seq = dna.transcribe()

        # calculate codon usage bias for sequence
        print(str(seq)[:100])
        print(f'len: {len(seq)}')

        codons = count_codons(seq)
        print(codons)

        if codons != -1:
            freqs = convert_to_fraction(codons)
            print(freqs)


def main():
    url = 'https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/'
    data_folder = 'data'
    download = False  # whether to download new data from ncbi
    if download:
        spec_list = [spec[:-1] for spec in ds.collect_data(url, data_folder)]
    else:
        # get every fna filename in data_folder
        # https://python-forum.io/thread-10014.html
        spec_list = [entry.name[:-4] for entry in os.scandir(data_folder) if (entry.is_file() and entry.name[-4:] == ".fna")]
        print(spec_list)
    calculate_CUB(data_folder, spec_list)


main()
