# -*- coding: utf-8 -*-
"""
Created on Mon Apr 17 10:55:15 2017

@author: @author: Dimitri Cabaud

The file "Species_on_Anon_planet.txt" contains 100 mini-"genome" sequences from living species found on the planet Anon. 
In this planet, life is functioning with proteins comprising only 8 types of amino acids. 
"Life X" is a new species on that planet and we want to find its known relatives so that we can better understand it biologically. 
"Life X" has a mini-genome that has been sequenced and the genome can be found in the file "Life_X_Query_Seq_2016.txt". 
This program find the Life X's closest relatives among the 100 known species on the planet.
The alignment methods can be chosen (global, local).
"""

from skbio.alignment import local_pairwise_align
from skbio.alignment import global_pairwise_align
from skbio.sequence import GrammaredSequence
from skbio.util import classproperty
import time


def import_data_from_fasta(nameofthefile):
    """
    This function import a fasta file and create two arrays which contains 
    all the names and their sequences associated
    :param string name of the file
    :return list of the sequences   
    :return list of the name  
    """
    species_dna = []
    species_name = []
    line_number = 0
    
    with open(nameofthefile, "r") as f:
        for line in f.readlines():
            if line_number%2 == 0:
                species_name.append(line.strip())
            else:
                species_dna.append(line.strip())
            line_number+=1
    return (species_dna,species_name)

  
class CustomSequence(GrammaredSequence):
    """
    This class define a custom dictionary of amino acids
    """
    @classproperty
    def degenerate_map(cls):
        return {}

    @classproperty
    def definite_chars(cls):
        return set("ABCDEFHG")


    @classproperty
    def default_gap_char(cls):
        return '-'

    @classproperty
    def gap_chars(cls):
        return set('-.')
 

def multi_pairwise_alignment(function_alignment, seq1_name, seq1_dna, array_seq2_name, array_seq2_dna, substitution_matrix, gap_open_penalty = 10, gap_extend_penalty = 2):
    """
    This function makes n pairwise alignment between one sequences and n others 
    :param function method used for the alignment 
    :param string name of the first sequence
    :param string sequence
    :param array of string names
    :param array of string sequences
    :param dico of dico for the subtitution matrix
    :param int gap open default = 0
    :param int gap extend default = 1
    :return tabularmsa alignement   
    :return tabularmsa score
    :return tabularmsa start_end_positions
    """
    alignment, score, start_end_positions = [],[],[]
    seq1 = CustomSequence(seq1_dna[0])
    for dna in range(len(array_seq2_dna)):
        seq2 = CustomSequence(array_seq2_dna[dna])        
        alignment_l, score_l, start_end_positions_l = function_alignment(seq1, seq2, gap_open_penalty, gap_extend_penalty, substitution_matrix)
        alignment.append(alignment_l)
        score.append(score_l)
        start_end_positions.append(start_end_positions_l)
        print("Progress: ",dna,"%")
    return alignment, score, start_end_positions

def closest_relative(number, alignment, score, known_species_dna, known_species_name, x_specie_dna, x_specie_name):
    """
    This function find the n closest alignments and return all their relative 
    informations
    """
    closest_name = []
    closest_dna = []
    closest_score = []
    for i in range(number):
        max_indice = score.index(max(score))
        closest_score.append(max(score))
        closest_name.append(known_species_name[max_indice])
        closest_dna.append(known_species_dna[max_indice])
        del known_species_name[max_indice]
        del known_species_dna[max_indice]
        del score[max_indice]
    return closest_name, closest_dna, closest_score
 
def parameters_definition():
    print("Do you want to define by yourself the algorithm and the parameters of the gap open/extend penalty and the match/mismatch score ?")
    print("Tape 1 for YES or 2 for NO and use the default parameters and press enter")
    print("Default parameters: algo = global_pairwise / gap open penalty = -10 / gap extend penalty = -2 / match score = 2 / mismatch score = -1")
    answer = int(input())
    print('\n')
    match = 2
    mismatch = -1
    function = 2
    gap_open_penalty = -10
    gap_extend_penalty = -2
    if answer == 1:
        print("Please tape 1 for local_pairwise_align / 2 for global_pairwise_align ")
        function = int(input())
        print('\n')
        print("Please tape the value of the gap penalty open and press enter")
        print("For example: -10")
        gap_open_penalty = int(input())
        print('\n')
        print("Please tape the value of the gap penalty extend and press enter")
        print("For example: -2")
        gap_extend_penalty = int(input())
        print('\n')
        print("Now tape the value of the match score and press enter")
        print("For example: 2")
        match = int(input())
        print('\n')
        print("To finish please tape the value of the mismatch score and press enter")
        print("For example: -1")
        mismatch = int(input())
        print('\n')
    print("How many alignments close to the sequence do you want? (max 100)")
    number_of_alignments = int(input())
    print('\n')
    print("Thank you ! Now your request is executing...")
    print('\n')
    return number_of_alignments, function, gap_open_penalty, gap_extend_penalty, match, mismatch

      
"""
Main
"""    
          
number_of_alignments, function, gap_open_penalty, gap_extend_penalty, match, mismatch = parameters_definition() 

substitution_matrix={
'A':{'A':match,'B':mismatch,'C':mismatch,'D':mismatch,'E':mismatch,'F':mismatch,'G':mismatch,'H':mismatch},
'B':{'A':mismatch,'B':match,'C':mismatch,'D':mismatch,'E':mismatch,'F':mismatch,'G':mismatch,'H':mismatch},
'C':{'A':mismatch,'B':mismatch,'C':match,'D':mismatch,'E':mismatch,'F':mismatch,'G':mismatch,'H':mismatch},
'D':{'A':mismatch,'B':mismatch,'C':mismatch,'D':match,'E':mismatch,'F':mismatch,'G':mismatch,'H':mismatch},
'E':{'A':mismatch,'B':mismatch,'C':mismatch,'D':mismatch,'E':match,'F':mismatch,'G':mismatch,'H':mismatch},
'F':{'A':mismatch,'B':mismatch,'C':mismatch,'D':mismatch,'E':mismatch,'F':match,'G':mismatch,'H':mismatch},
'G':{'A':mismatch,'B':mismatch,'C':mismatch,'D':mismatch,'E':mismatch,'F':mismatch,'G':match,'H':mismatch},
'H':{'A':mismatch,'B':mismatch,'C':mismatch,'D':mismatch,'E':mismatch,'F':mismatch,'G':mismatch,'H':match},
} 

start_time = time.time()
  
known_species_dna, known_species_name = import_data_from_fasta("100_known_species_Seq_2016.txt")
x_specie_dna, x_specie_name = import_data_from_fasta("Life_X_Query_Seq_2016.txt")

if function == 1: 
    alignment, score, start_end_positions = multi_pairwise_alignment(local_pairwise_align, x_specie_name, x_specie_dna, known_species_name, known_species_dna, substitution_matrix, -gap_open_penalty, -gap_extend_penalty)
elif function == 2:
    alignment, score, start_end_positions = multi_pairwise_alignment(global_pairwise_align, x_specie_name, x_specie_dna, known_species_name, known_species_dna, substitution_matrix, -gap_open_penalty, -gap_extend_penalty)


closest_name, closest_dna, closest_score = closest_relative(number_of_alignments, alignment, score, known_species_dna, known_species_name, x_specie_dna, x_specie_name) 
print('The ',number_of_alignments,' th closest alignments of Life X are:')
print(closest_name)
print('With a respective score of: ')
print(closest_score)
print("--- %s seconds ---" % (time.time() - start_time)) 

