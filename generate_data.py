import os
import time
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
    for i in range(0, 300, 3):  # len(seq), 3):
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


def calculate_CUB_stored_file(data_folder, spec_list, taxon):
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
            # print(freqs)

            data = pd.DataFrame.from_dict(freqs)
            data["Species"] = [species]
            data["Taxon"] = [taxon]
            # data["AccessionNum"] = [] ???
            data.to_csv("working.csv", mode='a', index=False)  # append row of data
            print(data)

            return data


def calculate_CUB(seq, species, taxon, accession_num="NA"):
    # calculate codon usage bias for sequence
    print(str(seq)[:100])
    print(f'len: {len(seq)}')

    codons = count_codons(seq)
    # print(codons)

    if codons != -1:
        freqs = convert_to_fraction(codons)
        print(freqs)

        data = pd.DataFrame.from_dict(freqs)
        data["AccessionNum"] = [accession_num]
        data["SeqLen"] = [len(seq)]
        data["Species"] = [species]
        data["Taxon"] = [taxon]
        data.to_csv("working.csv", mode='a', index=False)  # append row of data
        print(data)

        return data


def calculate_CUB_stored_files(data_folder, spec_list, taxon):
    """calculate CUB for multiple locally stored files in a list"""
    taxon_cubs = []
    for i in range(3):
        species = spec_list[i]
        print(species)

        # combine the fasta into one long sequence
        dna = fp.splice_fasta(f'{data_folder}/{species}.fna')

        # transcribe to RNA
        seq = dna.transcribe()

        # calculate CUB
        cub = calculate_CUB(seq, species, taxon)
        taxon_cubs.append(cub)

    return taxon_cubs


def collect_and_calculate_CUB(data_url, taxon):
    # find all species directories on the page
    directories = ds.getDirectories(data_url)

    if directories is None:
        return

    taxon_cubs = []
    # test on three species
    for spec_dir in directories:
        # for skipping to homo sapiens
        short_run = False
        if short_run:
            if spec_dir[0] != "H":
                continue

        # wait a moment before moving to next species
        # to avoid overwhelming the servers and getting blocked
        time.sleep(5)

        species = spec_dir[:-1]

        # download fasta genome from ncbi (w/o saving locally)
        viral = taxon == "viral"  # boolean flag
        result = ds.downloadFasta(data_url, spec_dir, download=False, viral=viral)
        if result is None:
            continue

        accession_num, file_content = result

        # combine the fasta into one long sequence
        dna = fp.splice_fasta_2(file_content)

        # transcribe to RNA
        seq = dna.transcribe()

        # calculate CUB
        cub = calculate_CUB(seq, species, taxon, accession_num)
        taxon_cubs.append(cub)

    return taxon_cubs


def main():
    # taxa = ["vertebrate_mammalian", "vertebrate_other", "invertebrate", "plant", "viral"]
    taxa = ["viral"]
    # todo: viral genomes from latest_assembly_versions/
    all_cubs = []
    for taxon in taxa:
        url = f'https://ftp.ncbi.nlm.nih.gov/genomes/refseq/{taxon}/'
        data_folder = 'data'
        download = False  # whether to download new data from ncbi
        from_memory = False  # whether to calculate based on local fasta files
        if download:
            spec_list = [spec[:-1] for spec in ds.collect_data(url, data_folder)]
            taxon_cubs = calculate_CUB_stored_files(data_folder, spec_list, taxon)
        elif from_memory:
            # get every fna filename in data_folder
            # https://python-forum.io/thread-10014.html
            spec_list = [entry.name[:-4] for entry in os.scandir(data_folder)
                         if (entry.is_file() and entry.name[-4:] == ".fna")]
            print(spec_list)
            taxon_cubs = calculate_CUB_stored_files(data_folder, spec_list, taxon)
        else:  # new best mode - download seq only temporarily to calculate cub
            taxon_cubs = collect_and_calculate_CUB(url, taxon)
        all_taxa = pd.concat(taxon_cubs)
        all_cubs.append(all_taxa)
        all_taxa.to_csv("working_by_taxon.csv", mode='a', index=False)

    # write table of CUBs to out_file
    out_file = "cub.csv"
    out_cub = pd.concat(all_cubs)
    out_cub.to_csv(out_file, index=False)



main()
