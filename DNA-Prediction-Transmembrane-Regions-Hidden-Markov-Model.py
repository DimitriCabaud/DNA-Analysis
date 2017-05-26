# -*- coding: utf-8 -*-
"""
Created on Thu May 18 17:17:26 2017

@author: Dimitri Cabaud
"""
from __future__ import division
import collections
import numpy as np
from hmmlearn import hmm
import time


def import_data(nameofthefile):
    
    data = []    
    with open(nameofthefile, "r") as f:
        for i, line in enumerate(f):
            if i%2 != 0:
                data.append(line.strip().split())
    return (data)
  
def clean_data_for_emission_probabilities(data_list):
    data_transmembrane, data_none_transmembrane = [[] for x in range(5000)], [[] for x in range(5000)]
    for i in range(len(data_list)):
        env = 0
        counter = 0
        for j in range(len(data_list[i][0])):
            if data_list[i][0][j] == '[':
                env = 1
                counter = 0
            elif data_list[i][0][j] == '(':
                env = 2
                counter = 0
            elif data_list[i][0][j] == ']' or data_list[i][0][j] == ')':
                env = 0
            elif env == 1:
                data_none_transmembrane[counter].append(data_list[i][0][j])
            elif env == 2:
                data_transmembrane[counter].append(data_list[i][0][j])  
            counter += 1   
    data_transmembrane = [''.join(x) for x in data_transmembrane if x != []]                   
    data_none_transmembrane = [''.join(x) for x in data_none_transmembrane if x != []] 
    return data_transmembrane, data_none_transmembrane
 
def calculate_emission_probabilities(chain):
    list_names, list_probabilities = [], []
    for i in range(len(chain)):
        names, probabilities = [], []
        frequencies = collections.Counter(chain[i])
        for cle,valeur in frequencies.items():
            names.append(cle)
            probabilities.append(valeur/len(chain[i]))
        list_names.append(names)
        list_probabilities.append(probabilities)
    return list_names, list_probabilities
    
def sort_emission_probabilities(list_names, list_probabilities, list_observations):
    list_probabilities_sorted = []
    for i in range(len(list_names)):
        current_list_probabilities = [0] * 20
        for j in range(len(list_names[i])): 
            for k in range(len(list_observations)): 
                if list_observations[k] == list_names[i][j]:
                    current_list_probabilities[k] = list_probabilities[i][j]
        list_probabilities_sorted.append(current_list_probabilities)  
    return list_probabilities_sorted

def calculate_start_probabilities(probabilities_emission_t, probabilities_emission_nt, data_list):
    count_t = 0
    count_nt = 0
    for i in range(len(data_list)):
        if data_list[i][0][0] == '[':
            count_nt += 1 
        elif data_list[i][0][0] == '(':
            count_t += 1
            
    start_probabilities = [0] * (len(probabilities_emission_t)+len(probabilities_emission_nt)+1)
    start_probabilities[0] = count_t/len(data_list)
    start_probabilities[len(probabilities_emission_t)] = count_nt/len(data_list)                      
    return start_probabilities

def calculate_transition_probabilities(data_list, len_probabilities_emission_t, len_probabilities_emission_nt):
    probabilities_transition = [[0] * (len_probabilities_emission_t+len_probabilities_emission_nt+1) for i in range((len_probabilities_emission_t+len_probabilities_emission_nt+1))]
    for i in range(len(data_list)): 
        env = 0
        counter = 0
        var_before = [0]*2
        for j in range(len(data_list[i][0])):
            if data_list[i][0][j] == '[':
                env = 2
                counter = 0
            elif data_list[i][0][j] == '(':
                env = 1
                counter = 0
            elif data_list[i][0][j] == ']':
                env = env
            elif data_list[i][0][j] == ')':
                env = env
            elif env == 1:
                if var_before[0] == 1:           
                    probabilities_transition[var_before[1]][counter]+=1                    
                elif var_before[0] == 2:               
                    probabilities_transition[var_before[1]+len_probabilities_emission_t][counter]+=1
                var_before = [env,counter]
                counter += 1
            elif env == 2:
                if var_before[0] == 1:                   
                    probabilities_transition[var_before[1]][counter+len_probabilities_emission_t]+=1
                elif var_before[0] == 2:                    
                    probabilities_transition[var_before[1]+len_probabilities_emission_t][counter+len_probabilities_emission_t]+=1
                var_before = [env,counter]
                counter += 1
                
    probabilities_transition[len_probabilities_emission_t+len_probabilities_emission_nt] = probabilities_transition[len_probabilities_emission_t+len_probabilities_emission_nt-1]
    probabilities_transition_normalized = [[0] * (len_probabilities_emission_t+len_probabilities_emission_nt+1) for i in range((len_probabilities_emission_t+len_probabilities_emission_nt+1))] 
    for i in range(len(probabilities_transition)):
        for j in range(len(probabilities_transition[i])):    
            if probabilities_transition[i][j]!= 0:
                probabilities_transition_normalized[i][j] = probabilities_transition[i][j] / sum(probabilities_transition[i])
    return probabilities_transition_normalized  
    

def convert_sequence_to_list_position(sequence, list_ref):
    list_position = []
    for i in range(len(sequence[0])):
        for j in range(len(list_ref)):
            if sequence[0][i] == list_ref[j]:
                list_position.append(j)
    list_position_final = [list_position]
    return list_position_final
    
def create_state_list(names_t, names_nt):
    list_state = []
    for i in range(len(names_t)):
        list_state.append('T'+str(i))
    for j in range(len(names_nt)):
        list_state.append('NT'+str(j))
    list_state.append('NT'+str(len(names_t)+len(names_nt)+1))
    return list_state

def convert_data_to_sequence_states(data_list, len_data_transmembrane):
    list_data_sequence_states = []
    for i in range(len(data_list)): 
        env = 0
        counter = 0
        sequence_states = []
        for j in range(len(data_list[i][0])):
            if data_list[i][0][j] == '[':
                env = 2
                counter = 0
            elif data_list[i][0][j] == '(':
                env = 1
                counter = 0
            elif data_list[i][0][j] == ']':
                env = env
            elif data_list[i][0][j] == ')':
                env = env
            elif env == 1:
                sequence_states.append(counter)                    
                counter += 1
            elif env == 2:
                sequence_states.append(counter + len_data_transmembrane)                    
                counter += 1
        list_data_sequence_states.append(sequence_states)
    return list_data_sequence_states
     
def compare_two_sequences_of_states(sequence_pred, sequence_ref, len_data_transmembrane):
    tp = 0
    tn = 0
    fp = 0
    fn = 0
    for i in range(len(sequence_pred)):
        if (int(sequence_pred[i]) < len_data_transmembrane) and (int(sequence_ref[i]) < len_data_transmembrane) : 
            tp+=1
        elif (int(sequence_pred[i]) >= len_data_transmembrane) and (int(sequence_ref[i]) >= len_data_transmembrane) : 
            tn+=1
        elif (int(sequence_pred[i]) < len_data_transmembrane) and (int(sequence_ref[i]) >= len_data_transmembrane) : 
            fp+=1
        elif (int(sequence_pred[i]) >= len_data_transmembrane) and (int(sequence_ref[i]) < len_data_transmembrane) : 
            fn+=1
    return tp, tn, fp, fn

def jackknife_validation(data, observations):
    tp, tn, fp, fn = 0,0,0,0
    len_test = int(len(data)/10)
    train, test = [], []   
    for i in range(10):
        
        print('Progression Jackknife validation : loop number',i,'/10')
        start = len_test*i
        end = len_test*(i+1)
        test = data[start:end]
        train = data[0:start] + data[end:]
    
        data_transmembrane, data_none_transmembrane = clean_data_for_emission_probabilities(train)
    
        names_t, probabilities_emission_t = calculate_emission_probabilities(data_transmembrane)
        names_nt, probabilities_emission_nt = calculate_emission_probabilities(data_none_transmembrane)       
        probabilities_emission_t = sort_emission_probabilities(names_t, probabilities_emission_t, observations)
        probabilities_emission_nt = sort_emission_probabilities(names_nt, probabilities_emission_nt, observations)
        probabilities_emission = probabilities_emission_t + probabilities_emission_nt 
        probabilities_emission.append(probabilities_emission_nt[-1])
        
        start_probabilities = calculate_start_probabilities(probabilities_emission_t, probabilities_emission_nt, train)
        
        probabilities_transition = calculate_transition_probabilities(train, len(probabilities_emission_t), len(probabilities_emission_nt))
        
        states = create_state_list(names_t, names_nt)
        n_states = len(states)
                
        start_probability = np.array(start_probabilities)      
        transition_probability = np.array(probabilities_transition)       
        emission_probability = np.array(probabilities_emission)
        
        model = hmm.MultinomialHMM(n_components=n_states)
        model.startprob_ = start_probability
        model.transmat_ = transition_probability
        model.emissionprob_ = emission_probability
    
        for j in range(len(test)):
            progress = (j/len(test))*100
            print('loading  ',"%.1f" % progress, ' %')
            sequence_test = test[j]
            sequence_test = convert_sequence_to_list_position(sequence_test, observations)
            sequence_test = np.array(sequence_test).T
            
            logprob, transmembrane_predictions = model.decode(sequence_test, algorithm="viterbi")
                    
            test_converted = convert_data_to_sequence_states(test, len(data_transmembrane))
            test_current = test_converted[j]
            tp_current, tn_current, fp_current, fn_current = compare_two_sequences_of_states(transmembrane_predictions, test_current, len(data_transmembrane))       
            tp+=tp_current
            tn+=tn_current
            fp+=fp_current
            fn+=fn_current
    return tp, tn, fp, fn
     
def first_option_user(data, observations, sequence_test):    
    # calculate the emission matrix
    data_transmembrane, data_none_transmembrane = clean_data_for_emission_probabilities(data)
    
    names_t, probabilities_emission_t = calculate_emission_probabilities(data_transmembrane)
    names_nt, probabilities_emission_nt = calculate_emission_probabilities(data_none_transmembrane)
    probabilities_emission_t = sort_emission_probabilities(names_t, probabilities_emission_t, observations)
    probabilities_emission_nt = sort_emission_probabilities(names_nt, probabilities_emission_nt, observations)
    probabilities_emission = probabilities_emission_t + probabilities_emission_nt 
    probabilities_emission.append(probabilities_emission_nt[-1])
    
    # calculate the start probabilities
    start_probabilities = calculate_start_probabilities(probabilities_emission_t, probabilities_emission_nt, data)
    
    # calculate the transition probabilities
    probabilities_transition = calculate_transition_probabilities(data, len(probabilities_emission_t), len(probabilities_emission_nt))
    
    
    # build a HMM model
    states = create_state_list(names_t, names_nt)
    n_states = len(states)
    
    start_probability = np.array(start_probabilities)
    transition_probability = np.array(probabilities_transition)
    emission_probability = np.array(probabilities_emission)
    
    model = hmm.MultinomialHMM(n_components=n_states)
    model.startprob_ = start_probability
    model.transmat_ = transition_probability
    model.emissionprob_ = emission_probability
    
    # predict a sequence of hidden states based on visible states
    # test with protein 1bcc G
    #sequence_test = ['[RQFGHLTRVRHLITYSLSPFEQRPFPHYFSKGVPNVWR](RLRACILRVAPPFLAFYLLYTWG)[TQEFEKSKRKNPAAYVN]']
    sequence_test = convert_sequence_to_list_position(sequence_test, observations)
    sequence_test = np.array(sequence_test).T
    
    logprob, transmembrane_predictions = model.decode(sequence_test, algorithm="viterbi")

    print('\n')    
    print("Sequence test:", ", ".join(map(lambda x: observations[x], sequence_test)))
    print('\n')
    print("Transmembrane Predictions:", ", ".join(map(lambda x: states[x], transmembrane_predictions)))
    print('\n')   
    # validation of this prediction
    data_conversion = convert_data_to_sequence_states(data, len(data_transmembrane))
    current_sequence = data_conversion[16]
    tp, tn, fp, fn = compare_two_sequences_of_states(transmembrane_predictions, current_sequence, len(data_transmembrane))       
    print("The number of true positive TP = ",tp," and true negative TN = ",tn)
    print("The number of false positive FP = ",fp," and false negative FN = ",fn)
    print('\n')
    sen = (tp/(tp+fn))
    spe = (tn/(tn+fp))
    print("The sensitivity (TP/(TP+FN)) = ",sen," and the specificity (TN/(TN+FP)) = ",spe)

def second_option_user(data, observations):
    # cross validation (Jackknife validation)    
    start_time = time.time()
    tp, tn, fp, fn = jackknife_validation(data, observations)
    print("The number of true positive TP = ",tp," and true negative TN = ",tn)
    print("The number of false positive FP = ",fp," and false negative FN = ",fn)
    print('\n')
    sen = (tp/(tp+fn))
    spe = (tn/(tn+fp))
    print("The sensitivity (TP/(TP+FN)) = ",sen," and the specificity (TN/(TN+FP)) = ",spe)
    print("--- %s seconds ---" % (time.time() - start_time)) 
    
def interface():
    observations = ['A','R','N','D','C','E','Q','G','H','I','L','K','M','F','P','S','T','W','Y','V']

    print("Do you want to change the name of the file containing the data ?")
    print("Tape 1 and Enter to change it or 2 to skip (by default: transmembrane_domain_sequence.txt ")
    answer = int(input())
    print('\n')
    if answer == 1:
        print("Please tape the name of the new name and tape Enter to validate")
        nameofthefile = str(input())
        data = import_data(nameofthefile)
        print('\n')
    if answer == 2:
        nameofthefile = 'transmembrane_domain_sequence.txt'
        data = import_data(nameofthefile)
    print("Now do you want to create a Hidden Markov Model with all the data and predict the transmembrane regions of a sequence ? (tape 1)" )
    print("or to check the performance of the model with the Jackknife validation ? (tape 2)" )        
    answer = int(input())
    if answer == 1:
        print("Please tape the sequence you want to analyse for the transmembrane prediction: ")
        print("For example: [RQFGHLTRVRHLITYSLSPFEQRPFPHYFSKGVPNVWR](RLRACILRVAPPFLAFYLLYTWG)[TQEFEKSKRKNPAAYVN]")
		print("where () is a transmembrane region and [] is not")
        sequence_test = [str(input())]
        print('\n')
        first_option_user(data, observations, sequence_test)
    if answer == 2:
        print("The Jackknife validation is quite long, please be patient and take a break ;) " )
        second_option_user(data, observations)

"""
Main
"""
interface()



