# -*- coding: utf-8 -*-
"""
DNA alignment with Gap Penalty 
(Neddleman/Wunsch techniques)
Created on Fri Mar 03 15:26:03 2017

@author: Dimitri Cabaud
"""
import numpy as np
import time

def init_matrix(dna_vector1, dna_vector2):
    """
    This function initialize a matrix with a size of the lenght of the 
    two vectors + 1, and fill it with 0.
    :param string sequence of DNA
    :param string sequence of DNA
    :return list of lists the matrix    
    """
    i = len(dna_vector1)+1
    j = len(dna_vector2)+1
    Matrix = np.matrix([[0 for x in range(j)] for y in range(i)])
    return Matrix
    
def score_match(i, j, dna_vector1, dna_vector2, match_score = 2, mismatch_score = -1):
    """
    This function compare a case with coordinate (i,j) of a matrix associated to two 
    DNA sequences with the an other case with coordinate (i-1,j-1)
    :param int abscissa of the matrix (i)
    :param int ordinate of the matrix (j)
    :param string sequence of DNA
    :param string sequence of DNA
    :param int score if there is match by default 2
    :param int score if there is no match by default -1
    :return int score
    """
    if dna_vector1[i-1] is dna_vector2[j-1]:
        return match_score
    else:
        return mismatch_score
    
def max_score(matrix, i, j, score_match, gap_penalty = -2):
    """
    This function send the max between 3 scores
    :param list of lists matrix 
    :param int abscissa of the matrix (i)
    :param int ordinate of the matrix (j)
    :param int score
    :param int gap penalty by default -2
    :return int the max score     
    """
    score_1 = matrix[i-1,j-1] + score_match
    score_2 = matrix[i,j-1] + gap_penalty
    score_3 = matrix[i-1,j] + gap_penalty
    return max(score_1, score_2, score_3)

def create_matrix(dna_vector1, dna_vector2, gap_penalty = -2, match_score = 2, mismatch_score = -1):
    """
    :param string sequence of DNA
    :param string sequence of DNA
    :param int gap penalty by default -2
    :param int score if there is match by default 2
    :param int score if there is no match by default -1
    :return list of lists the matrix
    """
    matrix = init_matrix(dna_vector1, dna_vector2)
    for i in range(1,len(dna_vector1)+1):
        for j in range(1,len(dna_vector2)+1):
            s = score_match(i, j, dna_vector1, dna_vector2, match_score, mismatch_score)
            matrix[i,j] = max_score(matrix, i, j, s, gap_penalty)
    return matrix
    
def comparison_traceback(matrix, i, j, dna_vector1, dna_vector2, working_solution1, working_solution2, 
                         gap_penalty = -2, match_score = 2, mismatch_score = -1):
    """
    This function compare a case with coordinate (i,j) of a matrix associated to two 
    DNA sequences with the an other case with coordinate (i-1,j-1) in order to construct
    the pathway
    :param list of lists matrix 
    :param int abscissa of the matrix (i)
    :param int ordinate of the matrix (j)
    :param string sequence of DNA
    :param string sequence of DNA
    :param string the current amino acid of the first sequence of DNA
    :param string the current amino acid of the second sequence of DNA
    :param int gap penalty by default -2
    :return string the new amino acid of the first sequence of DNA
    :return string the new amino acid of the second sequence of DNA
    :return int abscissa of the matrix (i) associated to the first amino acid
    :return int ordinate of the matrix (j) associated to the second amino acid
    :return list of the amino acid to put in the stack of the first sequence of DNA
    :return list of the amino acid to put in the stack of the second sequence of DNA
    :return list abscissa of the matrix (i) associated to the first amino acid to put in the stack
    :return list ordinate of the matrix (j) associated to the second amino acid to put in the stack
    """
    score_penalty = score_match(i, j, dna_vector1, dna_vector2, match_score, mismatch_score)
    match = 0
    (working_stack_solution1, working_stack_solution2, working_stack_i_to_add, 
    working_stack_j_to_add) = [], [], [], []    
    if matrix[i,j]== matrix[i-1,j-1] + score_penalty:
        working_solution1 = dna_vector1[i-1]
        working_solution2 = dna_vector2[j-1]
        match = 1
        new_i = i-1
        new_j = j-1
    if matrix[i,j]== matrix[i-1,j] + gap_penalty:
        if match == 0:
            working_solution1 = dna_vector1[i-1]
            working_solution2 = '_'
            match = 1
            new_i = i-1
            new_j = j
        else:
            working_stack_solution1 = dna_vector1[i-1]
            working_stack_solution2 = '_'
            working_stack_i_to_add.append(i-1)
            working_stack_j_to_add.append(j)            
    if matrix[i,j]== matrix[i,j-1] + gap_penalty:
        if match == 0:
            working_solution1 = '_'
            working_solution2 = dna_vector2[j-1]
            match = 1
            new_i = i
            new_j = j-1
        else:
            working_stack_solution1 = '_'
            working_stack_solution2 = dna_vector2[j-1]
            working_stack_i_to_add.append(i)
            working_stack_j_to_add.append(j-1)
    return (working_solution1, working_solution2, new_i, new_j, working_stack_solution1, 
    working_stack_solution2, working_stack_i_to_add, working_stack_j_to_add)
    
            
def traceback(matrix, dna_vector1, dna_vector2, i, j, working_solution1, working_solution2, 
              gap_penalty = -2, match_score = 2, mismatch_score = -1):
    """
    This function execute the traceback and return one solution (a pathway) and all the 
    coordinate of the others possibles solutions found
    :param list of lists matrix 
    :param string sequence of DNA
    :param string sequence of DNA
    :param int abscissa of the matrix (i)
    :param int ordinate of the matrix (j)
    :param string the current amino acid of the first sequence of DNA
    :param string the current amino acid of the second sequence of DNA
    :return string one solution for the first sequence of DNA analysed
    :return string one solution for the second sequence of DNA analysed
    :return list the stack of the solutions of the first sequence of DNA
    :return list the stack of the solutions of the second sequence of DNA
    """
    solution1 = working_solution1
    solution2 = working_solution2
    stack_solution1, stack_solution2, stack_i_to_add, stack_j_to_add = [], [], [], []
    while i!=0 and j!=0:
        (working_solution1, working_solution2, i, j, working_stack_solution1, working_stack_solution2, 
        working_stack_i_to_add, working_stack_j_to_add)  = comparison_traceback(matrix, i, j, dna_vector1, 
                                                           dna_vector2, working_solution1, working_solution2,
                                                           gap_penalty, match_score, mismatch_score)
        solution1 += working_solution1
        solution2 += working_solution2
        stack_solution1 += working_stack_solution1
        stack_solution2 += working_stack_solution2
        stack_i_to_add += working_stack_i_to_add
        stack_j_to_add += working_stack_j_to_add
    return solution1, solution2, stack_solution1, stack_solution2, stack_i_to_add, stack_j_to_add

def reverse(list_solution):
    """
    This function reverse all the strings in a list
    :param list of string 
    :return list of string reversed
    """
    return [element[::-1] for element in list_solution]


def complete_sequence(list_solution):
    """
    This function is recursive, it completes all the solutions with the part of precedent string 
    in a list of strings
    :param list of string
    :return list of string completed
    """
    new_list_solution = [list_solution[element-1][:len(list_solution[element-1])-len(list_solution[element])]
                         +list_solution[element] for element in range(1,len(list_solution))]
    new_list_solution.insert(0,list_solution[0])
    return new_list_solution
           
def stack_traceback(matrix, dna_vector1, dna_vector2, gap_penalty = -2, match_score = 2, mismatch_score = -1):
    """
    This function find all the possible alignements of two sequences of DNA using 
    Neddleman/Wunsch techniques
    :param list of lists the matrix
    :param string sequence of DNA
    :param string sequence of DNA
    :return list of string all the possible alignements of the first sequence
    :return list of string all the possible alignements of the second sequence
    """
    (final_solution1, final_solution2, stack_i_to_add, stack_j_to_add, stack_solution1, 
    stack_solution2)  = [], [], [], [], [], []
    index_final = 0
    i = len(dna_vector1)-1
    j = len(dna_vector2)-1
    stack_i = [i]
    stack_j = [j]
    working_solution1 = dna_vector1[-1]
    working_solution2 = dna_vector2[-1]
    while i!=0 and j!=0:
        (solution1, solution2, stack_solution1, stack_solution2, stack_i_to_add, 
        stack_j_to_add) = traceback(matrix, dna_vector1, dna_vector2, i, j, working_solution1, 
                          working_solution2, gap_penalty, match_score, mismatch_score)    
        if not stack_solution1 and not stack_solution2:
            i = 0
            j = 0
            final_solution1.insert(index_final,solution1)
            final_solution2.insert(index_final,solution2)
        else:        
            final_solution1.insert(index_final,solution1)
            final_solution2.insert(index_final,solution2)
            working_solution1 = stack_solution1[-1]
            working_solution2 = stack_solution2[-1]
            stack_i += stack_i_to_add
            stack_j += stack_j_to_add
            i = stack_i[-1]
            j = stack_j[-1]
            index_final+=1
    final_solution1 = complete_sequence(final_solution1)
    final_solution2 = complete_sequence(final_solution2)
    final_solution1 = reverse(final_solution1)
    final_solution2 = reverse(final_solution2)
    return final_solution1, final_solution2
    
def sequences_and_parameters_definition():
    dna_vector1 = 'GGATCGA'
    dna_vector2 = 'GAATTCAGTTA'
    gap_penalty = -2
    match_score = 2
    mismatch_score = -1        
    print("Do you want to input by yourself the DNA sequences to do the alignment job ?")
    print("Type 1 for YES or 2 for NO and use the default sequences and press enter")
    print("Default sequences: GGATCGA / GAATTCAGTTA")
    answer1 = input()  
    if answer1 == 1:
        print("Please type the first sequence with quotations and press enter")
        print("For example: GGATCGA")
        dna_vector1 = input()
        print("Now type the second sequence with quotations and press enter")
        print("For example: GAATTCAGTTA")
        dna_vector2 = input()

    print("Now do you want to define by yourself the parameters of the gap penalty and the match/mismatch score ?")
    print("Type 1 for YES or 2 for NO and use the default parameters and press enter")
    print("Default parameters: gap penalty = -2 / match score = 2 / mismatch score = -1")
    answer2 = input()
    if answer2 == 1:
        print("Please type the value of the gap penalty and press enter")
        print("For example: -2")
        gap_penalty = input()
        print("Now type the value of the match score and press enter")
        print("For example: 2")
        match_score = input()
        print("To finish please type the value of the mismatch score and press enter")
        print("For example: -1")
        mismatch_score = input()
    print("Thank you ! Now your request is executing...")
    return dna_vector1, dna_vector2, gap_penalty, match_score, mismatch_score

    
"""
Main
"""     
dna_vector1, dna_vector2, gap_penalty, match_score, mismatch_score = sequences_and_parameters_definition()   
start_time = time.time()
m = create_matrix(dna_vector1, dna_vector2, gap_penalty, match_score, mismatch_score)
solution1, solution2 = stack_traceback(m, dna_vector1, dna_vector2, gap_penalty, match_score, mismatch_score)
print("The differents results of alignements for the first sequence are: ")
print(solution1)
print("The differents results of alignements for the second sequence are: ")
print(solution2)
print("--- %s seconds ---" % (time.time() - start_time))    






