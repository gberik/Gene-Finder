# -*- coding: utf-8 -*-
"""
Created on Sat Feb  1 22:02:04 2014

@author: pruvolo
"""

from os import path


def load_seq(fasta_file):
    """ Reads a FASTA file and returns the DNA sequence as a string.

    fasta_file: the path to the FASTA file containing the DNA sequence
    returns: the DNA sequence as a string
    """
    retval = ""
    f = open(fasta_file)
    lines = f.readlines()
    for l in lines[1:]:
        retval += l[0:-1]
    f.close()
    return retval


def load_nitrogenase_seq():
    """ This function loads a sequence of DNA that is known to code for
        Nitrogenase.  Nitrogenase is an enzyme that fixes atmospheric
        Nitrogen (N_2)

        returns: the nucleotides in the DNA sequence as a string
    """
    f = open(path.join('.', 'data', 'nitrogenase NifH sequence.txt'), 'r')
    nitrogenase = f.readlines()
    f.close()

    # remove the first line as it is simply a sequence name.
    nitrogenase = nitrogenase[1:]
    for i, line in enumerate(nitrogenase):
        nitrogenase[i] = line[9:].replace(' ', '').replace('\r\n', '')

    nitrogenase = ''.join(nitrogenase).upper()
    return nitrogenase


def extract_next_gene(metagenome_lines, next_line):
    """ A helper function for load_metagenome.  This function
        takes an array of lines from the metagenome file and
        the next_line for processing.

        returns: a tuple consisting of the name of the snippet,
                 the sequence of the snippet, and the line number
                 to process next.
    """
    name = metagenome_lines[next_line].strip()[1:]
    next_line += 1
    start_line = next_line

    while next_line < len(metagenome_lines):
        if metagenome_lines[next_line][0] == '>':
            break
        next_line += 1
    return (name,
            ''.join([l.strip() for l in
                    metagenome_lines[start_line:next_line]]),
            next_line)

def load_contigs():
    """ Loads the DNA contigs for a new bacterial communicty
        returns: a list of DNA snippets consisting of (name, sequence)
                 tuples.  The sequence is represented as an uppercase
                 string of nucleotides
    """
    return load_metagenome_helper('genes segments for software_design_extension.txt')

def load_metagenome_helper(metagenome_file):
    """ Loads the metagenome stored in the specified file.
        returns: a list of DNA snippets consisting of (name, sequence)
                 tuples.  The sequence is represented as an uppercase
                 string of nucleotides
    """
    f = open(path.join('.', 'data', metagenome_file), 'r')

    metagenome_lines = f.readlines()
    f.close()
    next_line = 0
    snippets = []
    while next_line < len(metagenome_lines):
        (label, dna, next_line) = extract_next_gene(metagenome_lines,
                                                    next_line)
        snippets.append((label, dna.upper()))
    return snippets


def load_metagenome():
    """ Loads a metagenome of a bacterial contig.
        returns: a list of DNA snippets consisting of (name, sequence)
                 tuples.  The sequence is represented as an uppercase
                 string of nucleotides
    """
    return load_metagenome_helper('3300000497.a_metagenome_phototrophic community.fna')