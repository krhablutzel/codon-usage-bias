# import os
import time
import random
import requests
from bs4 import BeautifulSoup
import gzip
import fasta_parser as fp
from Bio import SeqIO
from Bio.Seq import Seq

##directory = r'../..'
##
## # https://newbedev.com/how-to-iterate-over-files-in-a-given-directory
##for entry in os.scandir(directory):
##    print(entry.path)

VERBOSE = True


def getDirectories(url, level=0):
    """returns the names of all the directories on a given page"""
    if VERBOSE:
        print("getting all directories at this url")
        print(url)

    time.sleep(random.randint(3, 8))  # wait a bit so we don't query too fast (triggers ncbi access block)

    if level >= 3:  # still blocked from NCBI servers
        with open("error_log.txt", "a") as f:
            f.write("cannot access \n" + url + "\n")
            return

    # https://stackoverflow.com/questions/11023530/python-to-list-http-files-and-directories/34718858
    page = requests.get(url).text

    # print(page)

    soup = BeautifulSoup(page, 'html.parser')

    # directory text/name ends in /
    directories = [node.getText() for node in soup.find_all('a') if node.getText().endswith('/')]

    if VERBOSE:
        print(directories)

    if not directories:  # ncbi mistook query for DOS attack and blocked access
        time.sleep(30*(level+1))   # wait a bit before continuing to query
        return getDirectories(url, level=level+1)

    return directories


def downloadFasta(data_url, spec_dir, out_folder=".", download=True, viral=False):
    """download the Fasta for the species refseq CDS"""
    # download from /genomes/refseq/vertebrate_mammilian/[species]/representative/[first directory]
    if VERBOSE:
        print(spec_dir)
    
    # find refseq name - the only directory listed on representative page
    # or latest_assembly_version for viral genomes
    if viral:
        sub_dir = 'latest_assembly_versions/'
    else:
        sub_dir = 'representative/'
    ref_seq_parent_url = data_url + f'{spec_dir}{sub_dir}'
    ref_seq = getDirectories(ref_seq_parent_url)
    if ref_seq is None:
        return
    else:
        ref_seq = ref_seq[0]
    if VERBOSE:
        print(ref_seq)

    # get fasta for ref_seq/ref_seq_cds_from_genomic.fna.gz
    fasta_url = ref_seq_parent_url + ref_seq + ref_seq[:-1] + '_cds_from_genomic.fna.gz'

    if download:
        # write decompressed fasta
        with open(out_folder + '/' + spec_dir[:-1] + '.fna', 'wb') as out:
            fasta_gz = requests.get(fasta_url)
            file_content = gzip.decompress(fasta_gz.content)
            out.write(file_content)
        return ref_seq[:-1]
    else:  # return file content rather than downloading
        fasta_gz = requests.get(fasta_url)
        file_content = gzip.decompress(fasta_gz.content)
        return ref_seq[:-1], file_content.decode("utf-8")


def collect_data(data_url, data_folder):
    # find all species directories on the page
    directories = getDirectories(data_url)

    # test on three species
    for i in range(3):
        # species to download genome of
        spec_dir = directories[i]

        # download the fasta for the species
        downloadFasta(data_url, spec_dir, data_folder)

    return directories


