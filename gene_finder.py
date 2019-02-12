# -*- coding: utf-8 -*-
"""
GENE_FINDER PROJECT
@author: Galip Sina Berik
"""

import random
from amino_acids import aa, codons, aa_table   # you may find these useful
from load import load_seq
dna = load_seq("./data/X73525.fa")

def shuffle_string(s):
    """Shuffles the characters in the input string
        NOTE: this is a helper function, you do not
        have to modify this in any way """
    return ''.join(random.sample(s, len(s)))

# YOU WILL START YOUR IMPLEMENTATION FROM HERE DOWN ###

def get_complement(nucleotide):
    """ Returns the complementary nucleotide

        nucleotide: a nucleotide (A, C, G, or T) represented as a stri ng
        returns: the complementary nucleotide
    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
    >>> get_complement('T')
    'A'
    >>> get_complement('G')
    'C'
    """
    complement = ''
    for i in range(len(nucleotide)):
        x = nucleotide[i]

        if x == 'A':
            var = 'T'
        elif x =='T':
            var = 'A'
        elif x == 'C':
            var = 'G'
        else:
            var = 'C'

        complement = complement + var
    return complement


def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence

        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    """
    reverse_dna = dna[::-1]
    reverse_complement = get_complement(reverse_dna)
    return reverse_complement



def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start
        codon and returns the sequence up to but not including the
        first in frame stop codon.  If there is no in frame stop codon,
        returns the whole string.
        dna: a DNA sequence
        returns: the open reading frame represented as a string
    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'
    >>> rest_of_ORF("ATGAGATAAG")
    'ATGAGA'
    >>> rest_of_ORF("ATGAGATCGG")
    'ATGAGATCGG'

    """
    ORF = ''
    x = 0
    for x in range(0,len(dna),3):
        i = dna[x:x+3]
        if dna[x:x+3] == 'TGA':
            return ORF
        if dna[x:x+3] == 'TAG':
            return ORF
        if dna[x:x+3] == 'TAA':
            return ORF
        ORF = ORF + i
    return ORF


def find_all_ORFs_oneframe(dna):
    """ Finds all non-nested open reading frames in the given DNA
        sequence and returns them as a list.  This function should
        only find ORFs that are in the default frame of the sequence
        (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    >>> find_all_ORFs_oneframe("ATGCATGATGAATGTAGATAGATGTGCCC")
    ['ATGCATGATGAATGTAGA', 'ATGTGCCC']
    """
    ORF2 = []
    x = 0
    StartCodon = 'ATG'
    while x <= len(dna):  
        i = dna[x:x+3]
        if i == StartCodon:
            ORF = rest_of_ORF(dna[x:])
            ORF2.append(ORF)
            x = x + len(ORF)
        x = x+3
    return ORF2


def find_all_ORFs(dna):
    """
    Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    >>> find_all_ORFs("ATGCATGAATGTAGTAGTGCATAA")
    ['ATGCATGAATGTAGTAGTGCA', 'ATGAATGTAGTAGTGCATAA', 'ATG']
    """

    ORFtotal = []
    ORF1 = find_all_ORFs_oneframe(dna)
    ORF2 = find_all_ORFs_oneframe(dna[1:])
    ORF3 = find_all_ORFs_oneframe(dna[2:])
    ORFtotal.extend(ORF1)
    ORFtotal.extend(ORF2)
    ORFtotal.extend(ORF3)
    return ORFtotal




def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    >>> find_all_ORFs_both_strands("ATGCGATGCTAAAGCTAGCCCGATGAATGTAGCCGATGCAATCA")
    ['ATGCGATGC', 'ATGAATGTAGCCGATGCAATCA', 'ATGCTAAAGCTAGCCCGA', 'ATG', 'ATGCAATCA']
    """
    ORF3 = []
    ORF1 = find_all_ORFs(dna)
    reverse = get_reverse_complement(dna)
    ORF2 = find_all_ORFs(reverse)
    ORF3.extend(ORF1)
    ORF3.extend(ORF2)
    return ORF3

def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    >>> longest_ORF("ATGCGATGCTAAAGCTAGCCCGATGAATGTAGCCGATGCAATCA")
    'ATGAATGTAGCCGATGCAATCA'
    """
    AllORFS = []
    x = find_all_ORFs_both_strands(dna)
    AllORFS.extend(x)
    if not AllORFS == []:
        y = max(AllORFS, key=len)
        return y
    if AllORFS == []:
        return ''

def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """
    i = 0
    dna2 = dna
    PossibleLongest = []
    while i <= num_trials:
        dna2 = shuffle_string(dna)
        longestORF = longest_ORF(dna2)
        PossibleLongest.append(longestORF)
        i = i + 1
    y = max(PossibleLongest, key=len)
    return y

def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).

        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment

        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'
        >>> coding_strand_to_AA("ATGCTAATCT")
        'MLI'
        >>> coding_strand_to_AA("ATGGTCTCT")
        'MVS'

    """
    aminoacidtotal = []
    aminoacid = ''
    x = 0
    dna1 = dna
    y = len(dna1) // 3
    y = y*3
    dna1 = dna1[:y]
    for x in range(0,len(dna1),3):
        i = str(dna1[x:x+3])
        aminoacid = aminoacid + aa_table[i]

    return aminoacid


def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """

    threshold = len(longest_ORF_noncoding(dna, 1500))
    aminoacids = []
    Allorfs = find_all_ORFs(dna)
    x = 0
    for x in range((len(Allorfs))):
        i = Allorfs[x]
        if len(i) > threshold:
            aminoacids.append(coding_strand_to_AA(i))
    return aminoacids

if __name__ == "__main__":
    import doctest
    #doctest.testmod()
    print (gene_finder(dna))
