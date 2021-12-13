# fasta parser

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from io import StringIO
from statistics import mean
from math import isclose

# References:
# https://biopython.org/docs/1.75/api/Bio.Seq.html


def get_uchart():
    """Return universal codon chart"""
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

    return uchart


def get_n_codons_of_aas(include_aas):
    """Return lists of aas with each n of codons"""
    uchart = get_uchart()

    n_codons = {}

    for aa in include_aas:
        n = len(uchart[aa])
        if n in n_codons:
            n_codons[n].append(aa)
        else:
            n_codons[n] = [aa]

    return n_codons


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
    uchart = get_uchart()
    freqs = {}
    aa_counts = {}

    for aa in uchart:
        # count occurrences of aa
        count = 0
        for codon in uchart[aa]:
            if codon in codons:
                count += codons[codon]

        # store count of aa
        aa_counts[aa] = count

        # insert relative codon freq into freq table
        for codon in uchart[aa]:
            if codon in codons:
                freqs[codon] = [codons[codon] / count]
            else:
                freqs[codon] = [0]  # codon never used

    return freqs, aa_counts


def wright_CUB(seq):
    """Calculate Wright's 'effective number of codons' statistic as measure of synonymous codon usage bias."""
    # count codons in sequence
    codons = count_codons(seq)
    if codons == -1:  # error flag
        return -1

    # get frequencies of codons, total count aa appears
    freqs, aa_counts = convert_to_fraction(codons)

    # get codon chart dict
    uchart = get_uchart()  # which codons go w/ each aa

    # calculate faa for each aa
    faas = {}
    for aa in aa_counts:
        # skip 1-codon aas and stop
        if aa in ["Trp", "Met", "Stop"]:
            continue

        # n is total count of codons for aa in the gene
        n = aa_counts[aa]
        if n <= 1:  # aa doesn't appear frequently enough to calculate faa
            continue

        # p is proportion this codon gets used
        sum_p_2 = 0
        for codon in uchart[aa]:
            p_2 = freqs[codon][0] ** 2  # freq was stored in list
            sum_p_2 += p_2

        # calculate faa
        faa = ((n * sum_p_2) - 1) / (n - 1)
        faas[aa] = faa  # store

    # print(faas)

    # average faa for duet, quartet, and sextet aas (and Ile)
    n_codons = get_n_codons_of_aas(faas.keys())
    if 2 in n_codons:
        f_duet = mean([faas[aa] for aa in n_codons[2]])
    if 3 in n_codons:
        f_trio = faas["Ile"]
    if 4 in n_codons:
        f_quartet = mean([faas[aa] for aa in n_codons[4]])
    if 6 in n_codons:
        f_sextet = mean([faas[aa] for aa in n_codons[6]])

    # error handling
    # okay if we overestimate Nc, we'll just miss some bias bc sequence too short
    if (2 not in n_codons) or (isclose(f_duet, 0)):
        f_duet = 9  # stop including f_duet in bias
    if (3 not in n_codons) or (isclose(f_trio, 0)):
        f_trio = 1  # stop including f_trio in bias
    if (4 not in n_codons) or (isclose(f_quartet, 0)):
        f_quartet = 5  # stop including f_quartet in bias
    if (6 not in n_codons) or (isclose(f_sextet, 0)):
        f_sextet = 3  # stop including f_sextet in bias

    # print(f_duet, f_trio, f_quartet, f_sextet)

    # calculate effective number of codons
    n_c = 2 + (9 / f_duet) + (1 / f_trio) + (5 / f_quartet) + (3 / f_sextet)

    # print("Nc: ", n_c)
    return n_c


def calculate_CUB(seq, species, taxon, accession_num="NA", orig_len="NA", orig_nc="NA"):
    # calculate codon usage bias for sequence
    print(str(seq)[:100])
    print(f'len: {len(seq)}')

    # calculate pct biased regions
    if orig_len == "NA":
        prop_biased_reg = "NA"
    elif len(seq) == 0:
        prop_biased_reg = 0
    else:
        prop_biased_reg = round(len(seq) / orig_len, 4)

    codons = count_codons(seq)
    print(codons)

    if codons != -1:
        freqs, _ = convert_to_fraction(codons)
        # print(freqs)

        data = pd.DataFrame.from_dict(freqs)
        data["AccessionNum"] = [accession_num]
        data["SeqLen"] = [orig_len]
        data["BiasedSeqLen"] = [len(seq)]
        data["PropBiasedRegions"] = [prop_biased_reg]
        data["Nc"] = [orig_nc]
        data["BiasedNc"] = [wright_CUB(seq)]  # Nc for biased region only
        data["Species"] = [species]
        data["Taxon"] = [taxon]
        data.to_csv("data/working.csv", mode='a', index=False)  # append row of data
        print(data)

        return data


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


def splice_fasta_in_memory(file_content, biased=False):
    """input: fasta file content as string"""
    # https://stackoverflow.com/questions/38358191/biopython-parse-from-variable-instead-of-file

    # report Nc
    # with StringIO(file_content) as f:
        # records = SeqIO.parse(f, 'fasta')
        # print([wright_CUB(record.seq) for record in records])

    with StringIO(file_content) as f:
        records = SeqIO.parse(f, 'fasta')
        # print(list(records)[0:2])
        # print(list(records)[0].seq)
        spacer = Seq('')
        if biased: # only include sequence if biased moderately
            cds_genome = spacer.join([record.seq for record in records if len(record.seq) % 3 == 0
                                      and record.seq[:3] == Seq('ATG') and wright_CUB(record.seq.transcribe()) < 50])
        else:
            cds_genome = spacer.join([record.seq for record in records if len(record.seq) % 3 == 0
                                      and record.seq[:3] == Seq('ATG')])

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
    all_pd.to_csv("data/all_taxa_12-11.csv", index=False)


# run fasta_parser to clean the working.csv in this directory
# clean_working_csv()
# combine_datasets()
