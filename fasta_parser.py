# fasta parser

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from io import StringIO

# References:
# https://biopython.org/docs/1.75/api/Bio.Seq.html


def splice_fasta(file_name):
    """Input is .fna file"""
    records = SeqIO.parse(file_name, 'fasta')  # records can only be accessed once
    # print(list(records)[0:2])
    # print(list(records)[0].seq)
    spacer = Seq('')
    # for record in records:
    #     print(record)
    #     print(record.seq)

    # only include proteins within correct reading frame and starting with ATG
    cds_genome = spacer.join([record.seq for record in records
                              if len(record.seq) % 3 == 0 and record.seq[:3] == Seq('ATG')])

    # records = SeqIO.parse(file_name, 'fasta')
    # print(sum([len(record.seq) % 3 for record in records if len(record.seq) % 3 == 0]))

    # print('cds genome:')
    # print(cds_genome)

    # with open('./data/test.fna', "w") as f:
    #     SeqIO.write(cds_genome, f, "fasta")
    return cds_genome


def splice_fasta_2(file_content):
    """input: fasta file content as string"""
    # https://stackoverflow.com/questions/38358191/biopython-parse-from-variable-instead-of-file
    with StringIO(file_content) as f:
        records = SeqIO.parse(f, 'fasta')
        # print(list(records)[0:2])
        # print(list(records)[0].seq)
        spacer = Seq('')
        cds_genome = spacer.join([record.seq for record in records
                              if len(record.seq) % 3 == 0 and record.seq[:3] == Seq('ATG')])
    return cds_genome


def clean_working_csv():
    """Remove duplicate headers from working.csv output file"""
    with open("data/working.csv", "r") as working:
        keep_content = working.readline()  # read header
        new_line = working.readline()  # read next line
        count = 1
        while new_line:  # keep every other next line
            if count % 2 == 1:  # lines 1, 3, etc.
                keep_content += new_line
            new_line = working.readline()
            count += 1

        with open("data/working_clean.csv", "w") as f:
            f.write(keep_content)


def report_error_log():
    """Return lists of species we can't access yet"""
    cannot_access = []
    with open("error_log.txt", "r") as f:
        new_line = f.readline()
        count = 0
        access = False
        while new_line:
            if count % 2 == 0:
                if "access" in new_line:
                    access = True
                else:
                    access = False
            elif access:  # can't access this url
                cannot_access.append(new_line[:-1])  # remove \n at end of url
                access = False
            new_line = f.readline()
            count += 1

    return cannot_access


def combine_datasets():
    # read in datasets
    animal_plant = pd.read_csv("data/working_animals_and_plants.csv")
    virus = pd.read_csv("data/working_virus.csv")
    print(animal_plant.shape)
    print(virus.shape)

    # stack into one
    all_pd = pd.concat([animal_plant, virus])
    print(all_pd.shape)

    # save the combined animal/plants + virus dataset
    all_pd.to_csv("data/all_taxa_12-7.csv", index=False)


# run fasta_parser to clean the working.csv in this directory
clean_working_csv()
combine_datasets()
