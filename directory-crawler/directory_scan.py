# directory crawler
# recursively downloads every file in a given web directory

# adapted from BIO334/335 Bioinformatics final project work
# for the Smith CS community mini-hackathon on 19 Nov 2021

import os
import requests
from bs4 import BeautifulSoup
import gzip


def recurse_directories(url, curr_dir="./"):
    if curr_dir == "./":  # store output in folder
        new_dir = curr_dir + "copied_directory/"
        os.mkdir(new_dir)
        curr_dir = new_dir

    # get all directories and files on page
    page = requests.get(url).text
    soup = BeautifulSoup(page, 'html.parser')
    links = [node.get('href') for node in soup.find_all('a')]

    # recurse on each directory
    # and download each file
    for link in links:
        if link[0] == '/':  # skip root references
            continue
        elif link[-1] == '/':  # directory
            new_dir = curr_dir + link  # create directory in file system
            os.mkdir(new_dir)
            recurse_directories(url + link, new_dir)
        else:  # download file
            download(link, url, curr_dir)

    return


def download(name, url, out_folder):
    """download file at the url"""

    # write file
    with open(out_folder + name, 'wb') as out:
        file = requests.get(url + name)
        out.write(file.content)
        return


def download_w_decompress(name, url, out_folder):
    """download file at the url"""

    # decompress and write file
    if url[-3] == '.gz':  # compressed file
        with open(out_folder + name, 'wb') as out:
            gz = requests.get(url + name)
            file_content = gzip.decompress(gz.content)
            out.write(file_content)
            return

    # write file
    with open(out_folder + name, 'wb') as out:
        file = requests.get(url + name)
        out.write(file.content)
        return


ncbi = 'https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Felis_catus/'
recurse_directories(ncbi)
