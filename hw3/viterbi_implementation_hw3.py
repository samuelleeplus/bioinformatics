#!/usr/bin/env python

# below you can find short explanations of the variables used :
# delta_i is the highest probability along a single path
# obs is observations
# T is length of observations
# states are all your states
# q is state q in all your states
# phi is start/initial probabilies
# trans_p is transition probabilities (A)
# emit_p is emission probabilities (B)


# You will modify viterbi function only!
# Helps visualize the steps of Viterbi. Do not edit print function. 
def print_path_probability_matrix(P,states):
    for q in states: print(q+ ":\t" + "\t".join(("%.2e" % (p)) for p in P[q]))


# Your assignment starts from here. You will fill out the missing lines 
# in the implementation of viterbi algorithm below. 
def viterbi(obs, states, phi, trans_p, emit_p):
    # We will create two variables; P and backpointer. 
    # Both variables (P and backpointer) are dictionaries.
    # You can find additional information on usage of data structures from
    # https://docs.python.org/3/tutorial/datastructures.html

    P = {}
    backpointer = {}
    T=len(obs) 

    # Initialization step
    # In both dictionaries, you will store lists, finally it will become a matrix
    # In P, you will store probabilities at state q in position i (P[q][i])
    # Similarly in backpointer, you will store where you came from to the current 
    # state at the current position
    for q in states:
        P[q] = [phi[q] * emit_p[q][obs[0]]]
        backpointer[q] = ["Start"]

        
    # Recursion step
    #########################################
    # missing lines 1: start recursion step.
    # you will need to calculate probability and state 
    newP=[{}]
    for i in range(1, T):
        newP.append({})
        newPath = {}
        for y in states:
            (prob, state) = max((newP[i-1][y0] * trans_p[y0][y] * emit_p[y][obs[i]], y0) for y0 in states)
            newP[i][y] = prob
            newPath[y] = Path[state] + [y]
    
        backpointer = newPath
        Path = newPath
    P = newP  
    print(P)

    #########################################
    

    # Finalization step
    (delta_i, state)=max((P[q][-1] * trans_p[q]["End"], q) for q in states)
    P["End"]=delta_i
    backpointer["End"]=state

    # Below line is not required. It is just to visualize path probability matrix 
    print_path_probability_matrix(P,states)

    # To retrieve the best path
    #########################################
    # missing lines 2: create a path variable and retrieve the best path. 
    # hint: you need to refer to backpointer 
    
    

    
    
    #########################################

    return ("probability of the best path to observe the sequence: %.2e \nthe best path: %s" % (delta_i, "\t".join(path[::-1])))


# Your coding assignment finishes here. Below is to run your algorithm on a real biological question.
states = ('Exon', 'Splice', 'Intron' )
seq="CTTCATGTGAAAGCAGACGTAAGTCA"
observations=tuple(seq)
start_probability = {'Exon': 1.0, 'Splice': 0.0, 'Intron': 0.0}
transition_probability = {
   'Exon' : {'Exon': 0.9, 'Splice': 0.1, 'Intron': 0.0, 'end': 0.0},
   'Splice' : {'Exon': 0.0, 'Splice': 0.0, 'Intron': 1.0, 'end': 0.0},
   'Intron' : {'Exon': 0.0, 'Splice': 0.0, 'Intron': 0.9, 'end': 0.1}
   }
emission_probability = {
   'Exon' : {'A': 0.25, 'C': 0.25, 'G': 0.25, 'T': 0.25},
   'Splice' : {'A': 0.05, 'C': 0.0, 'G': 0.95, 'T': 0.0},
   'Intron' : {'A': 0.4, 'C': 0.1, 'G': 0.1, 'T': 0.4},
   }

#########################################

print(viterbi(observations, states, start_probability, transition_probability, emission_probability))


# missing lines 3: call viterbi function to calculate probability and the best 
# path for the states observed, and provide results. 
#########################################
