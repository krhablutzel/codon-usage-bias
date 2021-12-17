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
import statistics


def calculate_CUB_stored_files(data_folder, spec_list, taxon):
    """calculate CUB for multiple locally stored files in a list"""
    taxon_cubs = []
    for i in range(len(spec_list)):
        species = spec_list[i]
        print(species)

        # combine the fasta into one long sequence
        dna = fp.splice_fasta(f'{data_folder}/{species}.fna')

        # transcribe to RNA
        seq = dna.transcribe()

        # calculate CUB
        cub = fp.calculate_CUB(seq, species, taxon)
        print(cub)
        taxon_cubs.append(cub)

    return taxon_cubs


def collect_and_calculate_CUB(data_url, taxon, exclude_spec=[], include_spec=[]):
    if not include_spec:
        # find all species directories on the page
        include_spec = ds.getDirectories(data_url)

        if include_spec is None:  # blocked from ncbi
            print("found no include species")
            return

        # exclude until working with viruses again
        # with open("data/viral_species.txt", "w") as f:
            # f.write("\n".join(include_spec))  # string of species directories from viral page

    # remove species we've already calculated cub for
    all_species = [species for species in include_spec if species[:-1] not in exclude_spec]

    # print the species we're about to iterate over
    # all_species.sort(reverse=True)  # iterate from Z viruses b/c blocked from A still
    print("iterating on:", all_species[:min(len(all_species), 20)])  # print up to 20 directories

    taxon_cubs = []
    # calculate cub for all species directories
    for spec_dir in all_species:
        # for skipping to homo sapiens
        short_run = False
        if short_run:
            if spec_dir[0] != "H":
                continue

        # wait a moment before moving to next species
        # to avoid overwhelming the servers and getting blocked
        # time.sleep(5)

        species = spec_dir[:-1]

        # download fasta genome from ncbi (w/o saving locally)
        viral = taxon == "viral"  # boolean flag
        result = ds.downloadFasta(data_url, spec_dir, download=False, viral=viral)
        if result is None:
            continue

        accession_num, file_content = result

        # combine the fasta into one long sequence
        dna = fp.splice_fasta_in_memory(file_content)
        dna_biased = fp.splice_fasta_in_memory(file_content, biased=True)

        # transcribe to RNA
        seq = dna.transcribe()  # use biased regions only for now
        seq_biased = dna_biased.transcribe()

        # calculate CUB
        cub = fp.calculate_CUB(seq_biased, species, taxon, accession_num, len(dna), fp.wright_CUB(seq))
        taxon_cubs.append(cub)

    return taxon_cubs


def download_species_by_taxa(data_url, taxon, data_folder):
    taxa_spec = ds.getDirectories(data_url)
    while taxa_spec is None:  # blocked from ncbi
        print("found no include species, trying again")
        time.sleep(40)
        taxa_spec = ds.getDirectories(data_url)

    # save list of species in taxon
    with open(f'{data_folder}/{taxon}.txt', "w") as f:
        f.write("\n".join(taxa_spec)+"\n")  # string of species directories from viral page


def calculate_CUB_queued(data_url, taxon, data_folder, exclude_spec=[]):
    # collect all species from saved file
    with open(f'{data_folder}/{taxon}.txt', "r") as spec_file:
        lines = spec_file.readlines()
        spec_dirs = [line[:-1] for line in lines]  # remove \n

    # remove species we've already calculated cub for
    all_species = [species for species in spec_dirs if species[:-1] not in exclude_spec]

    taxon_cubs = []
    collected_species = []
    # three passes through the list of species to collect data
    # if collection fails, species can go to back of queue to try again
    for i in range(3):
        # don't iterate if collected
        all_species = [species for species in all_species if species not in collected_species]

        # calculate cub for all species directories
        for spec_dir in all_species:
            species = spec_dir[:-1]

            # download fasta genome from ncbi (w/o saving locally)
            viral = taxon == "viral"  # boolean flag
            result = ds.downloadFasta(data_url, spec_dir, download=False, viral=viral)
            if result is None:
                continue

            accession_num, file_content = result

            # combine the fasta into one long sequence
            dna = fp.splice_fasta_in_memory(file_content)
            dna_biased = fp.splice_fasta_in_memory(file_content, biased=True)

            # transcribe to RNA
            seq = dna.transcribe()  # use biased regions only for now
            seq_biased = dna_biased.transcribe()

            # calculated CUB for all regions
            cub = fp.calculate_CUB(seq, species, taxon, accession_num, biased_regions="full")
            taxon_cubs.append(cub)

            # calculate CUB for biased regions only
            cub_biased = fp.calculate_CUB(seq_biased, species, taxon, accession_num, orig_len=len(dna),
                                          orig_nc=fp.wright_CUB(seq), biased_regions="biased")
            taxon_cubs.append(cub_biased)
            collected_species.append(spec_dir)

    return taxon_cubs


def main():
    # taxa = ["vertebrate_mammalian", "vertebrate_other", "invertebrate", "plant", "viral"]
    taxa = ["viral"]
    all_cubs = []
    for taxon in taxa:
        url = f'https://ftp.ncbi.nlm.nih.gov/genomes/refseq/{taxon}/'
        data_folder = 'data'
        download = False  # whether to download new data from ncbi
        from_memory = False  # whether to calculate based on local fasta files
        resume_run = False  # resume from interrupted run
        run_from_species_urls = False  # run from list of species urls (from error log)
        queued_run = True  # run from species queues
        if download:  # download every genome to local storage
            spec_list = [spec[:-1] for spec in ds.collect_data(url, data_folder)]
            taxon_cubs = calculate_CUB_stored_files(data_folder, spec_list, taxon)
        elif from_memory:  # calculate cub from locally stored files
            # get every fna filename in data_folder
            # https://python-forum.io/thread-10014.html
            spec_list = [entry.name[:-4] for entry in os.scandir(data_folder)
                         if (entry.is_file() and entry.name[-4:] == ".fna")]
            print(spec_list)
            taxon_cubs = calculate_CUB_stored_files(data_folder, spec_list, taxon)
        elif resume_run:  # make sure to start from the taxon we stopped at!
            fp.clean_working_csv()
            existing_data = pd.read_csv("data/working_clean.csv")
            exclude_spec = list(existing_data["Species"])
            print("Excluding\n", exclude_spec[:20])
            if run_from_species_urls:
                with open("data/viral_species.txt", "r") as spec_file:
                    lines = spec_file.readlines()
                    include_spec = [line[:-1] for line in lines]  # remove \n
                    print("Including\n", include_spec[:20])
            else:
                include_spec = []
            new_data = collect_and_calculate_CUB(url, taxon, exclude_spec, include_spec)
            taxon_cubs = pd.concat([existing_data, new_data])
        elif run_from_species_urls:  # but not resuming run
            print("running from species urls")
            # include-spec is currently empty???
            with open("data/viral_species.txt", "r") as spec_file:
                lines = spec_file.readlines()
                include_spec = [line[:-1] for line in lines]  # remove \n
                print("Including\n", include_spec[:20])
            taxon_cubs = collect_and_calculate_CUB(url, taxon, include_spec=include_spec)
        elif queued_run:
            # download_species_by_taxa(url, taxon, data_folder + '/species_lists')

            # get species to exclude
            try:
                fp.clean_working_csv()
                existing_data = pd.read_csv(f"{data_folder}/working_clean.csv")
                exclude_spec = list(existing_data["Species"])
            except:
                exclude_spec = []

            # collect cubs
            taxon_cubs = calculate_CUB_queued(url, taxon, data_folder + '/species_lists', exclude_spec)
        else:  # most cohesive mode - download seq only temporarily to calculate cub
            print("~~beginning normal run~~")
            taxon_cubs = collect_and_calculate_CUB(url, taxon)
        all_taxa = pd.concat(taxon_cubs)
        all_cubs.append(all_taxa)
        all_taxa.to_csv("data/working_by_taxon.csv", mode='a', index=False)

    # write table of CUBs to out_file
    out_file = "data/cub.csv"
    out_cub = pd.concat(all_cubs)
    out_cub.to_csv(out_file, index=False)


main()
