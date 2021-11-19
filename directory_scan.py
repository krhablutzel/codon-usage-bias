## test

# import os
import requests
from bs4 import BeautifulSoup
import gzip

##directory = r'../..'
##
## # https://newbedev.com/how-to-iterate-over-files-in-a-given-directory
##for entry in os.scandir(directory):
##    print(entry.path)

def getDirectories(url):
    '''returns the names of all the directories on a given page'''

    # https://stackoverflow.com/questions/11023530/python-to-list-http-files-and-directories/34718858
    page = requests.get(url).text

    # print(page)

    soup = BeautifulSoup(page, 'html.parser')

    # directory text/name ends in /
    directories = [node.getText() for node in soup.find_all('a') if node.getText().endswith('/')]

    return directories

url = 'https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/'
directories = getDirectories(url)

# test on one species

spec_dir = directories[0]



def downloadFasta(spec_dir, out_folder):
    '''download the Fasta for the species refseq CDS'''
    # download from /genomes/refseq/vertebrate_mammilian/[species]/representative/[first directory]
    print(spec_dir)
    
    # find refseq name - the only directory listed on representative page
    ref_seq_parent_url = url + f'{spec_dir}representative/'
    ref_seq = getDirectories(ref_seq_parent_url)[0]
    print(ref_seq)

    # get fasta for ref_seq/ref_seq_cds_from_genomic.fna.gz
    fasta_url = ref_seq_parent_url + ref_seq + ref_seq[:-1] + '_cds_from_genomic.fna.gz'

    # write decompressed fasta
    with open(out_folder + '/' + spec_dir[:-1] + '.fna', 'wb') as out:
        fasta_gz = requests.get(fasta_url)
        file_content = gzip.decompress(fasta_gz.content)
        out.write(file_content)
    
# data = requests.get(fasta_url).text


downloadFasta(spec_dir, '.')

