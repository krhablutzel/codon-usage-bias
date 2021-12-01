# fasta parser

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
    # https://stackoverflow.com/questions/38358191/biopython-parse-from-variable-instead-of-file
    with StringIO(file_content.decode('utf-8')) as f:
        records = SeqIO.parse(f, 'fasta')
        # print(list(records)[0:2])
        # print(list(records)[0].seq)
        spacer = Seq('')
        cds_genome = spacer.join([record.seq for record in records])
    return cds_genome
