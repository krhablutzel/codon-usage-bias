import os
import requests
from bs4 import BeautifulSoup
import gzip
import directory_scan as ds
import fasta_parser as fp
from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd


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
              "Ile": ["AUU", "AUC", "AUA"],
              "Met": ["AUG"],
              "Val": ["GUU", "GUC", "GUA", "GUG"],
              "Ser": ["UCU", "UCC", "UCA", "UCG", "AGU", "AGC"],
              "Pro": ["CCU", "CCC", "CCA", "CCG"],
              "Thr": ["ACU", "ACC", "ACA", "ACG"],
              "Ala": ["GCU", "GCC", "GCA", "GCG"],
              "Tyr": ["UAU", "UAC"],
              "Stop": ["UAA", "UAG", "UGA"],
              "His": ["CAU", "CAC"],
              "Gln": ["CAA", "CAG"],
              "Asn": ["AAU", "AAC"],
              "Lys": ["AAA", "AAG"],
              "Asp": ["GAU", "GAC"],
              "Glu": ["GAA", "GAG"],
              "Cys": ["UGU", "UGC"],
              "Trp": ["UGG"],
              "Arg": ["CGU", "CGC", "CGA", "CGG", "AGA", "AGG"],
              "Gly": ["GGU", "GGC", "GGA", "GGG"]}

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
                freqs[codon] = [codons[codon] / count]
            else:
                freqs[codon] = [0]  # codon never used

    return freqs


def calculate_CUB(data_folder, spec_list, taxon):
    cub_dfs = []
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
        # print(codons)

        if codons != -1:
            freqs = convert_to_fraction(codons)
            print(freqs)

        data = pd.DataFrame.from_dict(freqs)
        data["Species"] = [species]
        data["Taxon"] = [taxon]
        # data["AccessionNum"] = [] ???
        cub_dfs.append(data)  # store w/ all species dfs
        print(data)

    return pd.concat(cub_dfs)


def main():
    taxa = ["vertebrate_mammalian", "vertebrate_other", "invertebrate", "plant", "viral"]
    all_cub_dfs = []
    for taxon in taxa:
        url = f'https://ftp.ncbi.nlm.nih.gov/genomes/refseq/{taxon}/'
        data_folder = 'data'
        download = True  # whether to download new data from ncbi
        if download:
            spec_list = [spec[:-1] for spec in ds.collect_data(url, data_folder)]
        else:
            # get every fna filename in data_folder
            # https://python-forum.io/thread-10014.html
            spec_list = [entry.name[:-4] for entry in os.scandir(data_folder) if (entry.is_file() and entry.name[-4:] == ".fna")]
            print(spec_list)
        cub_df = calculate_CUB(data_folder, spec_list, taxon)
        all_cub_dfs.append(cub_df)  # store with cub tables for other taxa

    # combine all cub tables into one
    all_cub = pd.concat(all_cub_dfs)

    all_cub.to_csv('cub.csv', index=False)



main()
