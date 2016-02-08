# -*- coding: utf-8 -*-
"""
This is my gene finder program, which can take a genome and process where certain genes are.

@author: Kevin Zhang

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
    'ValueError: This is not a DNA nucleotide, you're a chicken.'
    >>> get_complement('Z')
    'ValueError: This is not a DNA nucleotide, you're a chicken.'

    """
    
	
    base_list = {'C':'G','G':'C','A':'T','T':'A'}

    if nucleotide in base_list:
    	return base_list[nucleotide]

    else:
    	raise ValueError("This is not a DNA nucleotide, you're a chicken.")	


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
    'ValueError: This is not a DNA nucleotide, you're a chicken.'
    >>> get_reverse_complement("TTTTTTTTTTTT")
    'AAAAAAAAAAAA'
    >>> get_reverse_complement("CCGTACATG")
    'CATGTACGG'
    """
    

    
    complementary_reverse_DNA = ""
    reverse_DNA = dna[::-1]

    for base in reverse_DNA:
    	complementary_reverse_DNA = complementary_reverse_DNA + get_complement(base)


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

        also will test when ATG is not at the beginning and also when there are junk codons in the middle to see it can recognize exactly where codons are

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
    	index = index + 3	

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
    
    return reduce(list.__add__, [find_all_ORFs_oneframe(dna[index:]) for index in range(3)])



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


    return find_all_ORFs(dna) + find_all_ORFs(get_reverse_complement(dna))



def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string

        There is no specification on what to do if there are two ORFs of the same longest length. In that case, my method will simply keep the first.
        There is a test to ensure that is true.

    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    >>> longest_ORF("ATGGGGTGAATGAAATGA")
    'ATGGGG'
    """
    
    ORF_list = find_all_ORFs_both_strands(dna)

    max_ORF = ''

    for ORF in ORF_list:
    	if len(max_ORF) < len(ORF):
    		max_ORF = ORF

    
    return max_ORF		


def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF

        This method can't be tested with doctest because the shuffling is randomized and can't be accurately predicted for a strict test.

        This method can be kind of tested using very strict test cases and very small num_trials, which I do outside of doctest manually. 
        Then I just run it a bajillion times and figure out what the numbers come out to be.

        In this manner, I used 'ATGGGGTGA' on 2 trials, and after testing came out with a lot of 0, 3, 6, and 9. So I know that 0 is when nothing was scrambled right, 
        3 is when it was scrambled 'ATGTGA', 6 is when got scrambled to its original order, and 9 is when it lost the stop codon, so it just counted the whole string.


    """

    max_length = 0
    shuffled_dna = ''

    for i in range(num_trials):
    	shuffled_dna = shuffle_string(dna)
    	max_ORF = longest_ORF(shuffled_dna)
    	if max_length < len(max_ORF):
    		max_length = len(max_ORF)

    return max_length		



def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).

        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment

        Also made some additional tests for some edge cases.         

        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'
        >>> coding_strand_to_AA("AAAAAAAAAA")
        'KKK'
        >>> coding_strand_to_AA("")
        ''
        >>> coding_strand_to_AA("AA")
        ''
    """
    aminoacid = ''
    codon = ''

    if(len(dna) > 2):
	    for index in range(0,len(dna),3):
	    	codon = dna[index:index+3]
	    	if(len(codon) == 3):
	    		aminoacid = aminoacid + aa_table[codon]

    return aminoacid


def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    
    threshold = longest_ORF_noncoding(dna, 10)
    ORF_list = find_all_ORFs_both_strands(dna)

    aa_seq = []

    for ORF in ORF_list:
    	if len(ORF) > threshold:
    		aa_seq.append(coding_strand_to_AA(ORF))

    return aa_seq		

if __name__ == "__main__":
    import doctest
    #doctest.run_docstring_examples(coding_strand_to_AA, globals(),verbose=True)
    #print longest_ORF_noncoding("ATGGGGTGA",2)
    #from load import load_seq
    #dna = load_seq("./data/X73525.fa")
    #print gene_finder(dna)
    from load import load_contigs
    contigs = load_contigs()
    print gene_finder(contigs[7][1])
