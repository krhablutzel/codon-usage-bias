# fasta parser

from Bio import SeqIO
from Bio.Seq import Seq

records = SeqIO.parse('Ailuropoda_melanoleuca.fna', 'fasta')

# print(list(records)[0:2])

# print(list(records)[0].seq)

spacer = Seq('')
cds_genome = spacer.join([record.seq for record in records])

print(cds_genome)
