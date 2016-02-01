# -*- coding: utf-8 -*-
"""
YOUR HEADER COMMENT HERE

@author: YOUR NAME HERE

"""

import random
from amino_acids import aa, codons, aa_table   # you may find these useful
from load import load_seq


def shuffle_string(s):
    """Shuffles the characters in the input string
        NOTE: this is a helper function, you do not
        have to modify this in any way """
    return ''.join(random.sample(s, len(s)))

# YOU WILL START YOUR IMPLEMENTATION FROM HERE DOWN ###


def get_complement(nucleotide):
    """ Returns the complementary nucleotide

        nucleotide: a nucleotide (A, C, G, or T) represented as a string
        returns: the complementary nucleotide

        Also adding all the bases just for complete testing. Will also add extreme cases such as non-strings and non-DNA to test for robustness
    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
    >>> get_complement('T')
    'A'
    >>> get_complement('A')
    'T'
    >>> get_complement(15)
    'This is not a DNA base, good sir'
    >>> get_complement('Z')
    'This is not a DNA base'

    """
    
    if not type(nucleotide) is str:
    	return "This is not a DNA base, good sir"

    normal_list = ['G','C','T','A']	
    complementary_list = ['C','G','A','T']	

    for i in range(0,4):
    	if nucleotide == normal_list[i]:
    		return complementary_list[i]
    return "This is not a DNA base"		


def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence

        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string

        Adding additional cases to prove robustness, as well as additional testing for completeness

    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    >>> get_reverse_complement("CCZZTHSGWNSS")
    'This DNA string is not real'
    >>> get_reverse_complement("TTTTTTTTTTTT")
    'AAAAAAAAAAAA'
    >>> get_reverse_complement("CCGTACATG")
    'CATGTACGG'
    """
    

    
    complementary_reverse_DNA = ""
    reverse_DNA = dna[::-1]

    for base in reverse_DNA:
    	complement = get_complement(base)
    	if len(complement) == 1:
    		complementary_reverse_DNA = complementary_reverse_DNA + complement
    	else:
    		return "This DNA string is not real"	

    return complementary_reverse_DNA	


def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start
        codon and returns the sequence up to but not including the
        first in frame stop codon.  If there is no in frame stop codon,
        returns the whole string.

        dna: a DNA sequence
        returns: the open reading frame represented as a string

        also testing to see if entire string is returned when there is no stop codon
        will also test for outlier cases

    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'
    >>> rest_of_ORF("ATGAGAGA")
    'ATGAGAGA'
    >>> rest_of_ORF(15)
    'Something wrong with this DNA strand'
    """
    
    ORF = dna
    codon = ""

    if not type(dna) is str:
    		return "Something wrong with this DNA strand"

    for index in range(0,len(dna),3):
    	codon = dna[index:index+3]
    	for stopcodon in codons[10]:
    		if(codon == stopcodon):
    			ORF = ORF[0:index]
    			
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

        will also test for robustness of ATG recognization and ORF identification

        also will test when ATG is not at the beginning and also when there are junk codons in the middle to see it can recongize exactly where codons are

    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    >>> find_all_ORFs_oneframe("ATGTAGATGTAGATGTAGATGTAG")
    ['ATG', 'ATG', 'ATG', 'ATG']
    >>> find_all_ORFs_oneframe("CCGCGCATGAGATGATAGCCGGGGATGAGATAG")
    ['ATGAGA', 'ATGAGA']


    """
    

    ORF_list = []

    codon = ""
    index = 0

    while index < len(dna):
    	codon = dna[index:index+3]
    	if(codon == "ATG"):
    		ORF = rest_of_ORF(dna[index:])
    		ORF_list.append(ORF)
    		index = index + len(ORF)
    	index = index +3	

    return ORF_list		

    	




def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs


        testing for outlier cases where there are no good codons

        also will test for completeness by using a lot of start and stop codons

    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    >>> find_all_ORFs("AAAAAAAAAAAAAAAAAAA")
    []
    >>> find_all_ORFs("CCATGAATGTGAGCGATGACTGATTAG")
    ['ATG', 'ATGACTGAT', 'ATGAATGTGAGCGATGAC']
    """
    
    return [i for index in range(0,3) for i in find_all_ORFs_oneframe(dna[index:])]



def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

        I don't think that this method needs further testing because the only difference is that I'm doing find_all_ORFs() except on the opposite side of DNA.
        This method's functionality has already been tested thoroughly in the helper methods it uses, there is no need to rigorously test the accumulation of a group of helper methods
        which have all been well tested and guaranteed functionality

    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """


    return [i for i in find_all_ORFs(dna)] + [i for i in find_all_ORFs(get_reverse_complement(dna))]



def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    # TODO: implement this
    pass


def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """
    # TODO: implement this
    pass


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
    """
    # TODO: implement this
    pass


def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    # TODO: implement this
    pass

if __name__ == "__main__":
    import doctest
    doctest.run_docstring_examples(find_all_ORFs_both_strands, globals(),verbose=True)
