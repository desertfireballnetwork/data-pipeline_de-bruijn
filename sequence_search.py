#!/usr/bin/env python3

# searches a given sequence for a given subsequence using the Python Levenshtein distance module
# Robert Howie, July 2014
# refactored MCT, may 15

#import standard modules
import distance #pip install distance
import argparse
import itertools
#import numpy as np #only for testing

import de_bruijn

def txt_results( all_possible_observations, minimum_values, minima_positions):
    print("\nPossible observations (when taking unknowns into account):\n{0}\n".format(all_possible_observations) )
    print("Best Match distance(s):\n{0}\n".format(minimum_values) )
    print("Positions of best matches:\n{0}\n".format(minima_positions) )

def gui_results( all_scores):
    [plt.plot(i) for i in all_scores]
    plt.show()


def main( arguments):
    """main program, pass an argparse list"""
    if arguments.sequence: #if using arbitrary sequence
        sequence = arguments.sequence        
    else: #if using de bruijn sequence
        #generate De Bruijn sequence and convert to string
        sequence = "".join([str(i) for i in de_bruijn.de_bruijn_sequence(arguments.use_de_bruijn_sequence[0], arguments.use_de_bruijn_sequence[1]).de_bruijn_sequence])
        
    # create sequence set for determining alphabet size
    sequence_set = set(sequence)
    
    # simulates cyclic sequence by doubling string, memory usage could be reduced
    # by only extending string by observation_length-1
    sequence_doubled = sequence * 2
    
    observation = arguments.observation

    #change X's to x's
    observation = "".join("x" if i == "X" else i for i in observation)
    
    if arguments.flip_observation:
        observation = observation[::-1] #flip string
    
    #debug print observation set stripping of any x/X's
    #print(set(i for i in observation if i != "x" and i != "X")) 
    
    if arguments.invert_observation and len(sequence_set) == 2: #if binary observation
        observation = "".join("0" if i == "1" else "1" for i in observation) #invert binary string
    
    #substitution for unknowns (x's)
    #print([i for i,j in enumerate(observation) if j == "x"])
    
    unknown_locations = [i for i,j in enumerate(observation) if j == "x"]
    
    substitutions = itertools.product(*itertools.repeat(list(sequence_set), len(unknown_locations)))
    substitutions_list = sorted([list(i) for i in substitutions])

    #for  in unknown_locations
    #create additional observations for different possibilities
    observation = [list(observation)] * len(substitutions_list)
    
    # nested list comprehension to generate the possible observations from
    # the different substitutions and convert back to strings
    all_possible_observations = ["".join(observation[i][j] if j not in unknown_locations 
                                        else substitutions_list[i][unknown_locations.index(j)]
                                        for j in range(len(observation[i])) )
                                for i in range(len(observation)) ]
    
    #double for fake wrap around
    all_possible_observations_doubled = [i*2 for i in all_possible_observations]
    
    #calculate match ratios
    if arguments.use_levenshtein_distance == False:
        #print("Using Hamming distance")
        all_scores = [[distance.hamming(i, sequence_doubled[j:j+len(i)])
                    for j in range(len(sequence))] for i in all_possible_observations]
    else:
        #print("Using Levenshtein distance")
        all_scores = [[distance.levenshtein(i, sequence_doubled[j:j+len(i)])
                    for j in range(len(sequence))] for i in all_possible_observations]
    
    #np.savetxt("scores.txt",np.array(all_scores)) #save for debugging
    
    #for i in all_scores
    #maximise matches
    minimum_values = [min(i) for i in all_scores]
    minima_positions = [[position for position, value in enumerate(i) if value == min(i)]
                        for i in all_scores]
    
    return all_possible_observations, minimum_values, minima_positions, all_scores

    
if __name__ == "__main__":
    # constants
    default_dbs_alphabet_size = 2
    default_dbs_subsequence_length = 9

    #create argument parser
    parser = argparse.ArgumentParser(description="Searches arbitrary sequences or De Bruijn sequences for observed subsequences (uses Hamming distance by default)")
    #optional argument for de bruijn sequences
    parser.add_argument("-d", "--use-de-bruijn-sequence", nargs=2, type=int, default=[default_dbs_alphabet_size, default_dbs_subsequence_length],
                        help="[-d <alphabet size ({0})> <subsequence length ({1})>]".format(default_dbs_alphabet_size, default_dbs_subsequence_length))
    #optional argument for arbitrary sequences
    parser.add_argument("-s", "--sequence", help="[-s <arbitrary sequence>]", type=str)
    #optional argument for Levenshtein distance instead of Hamming distance
    parser.add_argument("-l", "--use-levenshtein-distance", help="use Levenshtein distance instead of Hamming distance", action="store_true")
    #optional argument for try flipping observed sequence direction
    parser.add_argument("-f", "--flip-observation", help="flip observation direction", action="store_true")
    #optional argument for try inverted observation (binary only)
    parser.add_argument("-i", "--invert-observation", help="invert observation (binary only)", action="store_true")
    #positional argument for observed sequence
    parser.add_argument("observation", help="<observed sequence (use x's for unknowns eg: 0010x01xx)>", type=str)
    arguments = parser.parse_args()
    
    all_possible_observations, minimum_values, minima_positions, all_scores = main( arguments)

    txt_results( all_possible_observations, minimum_values, minima_positions)

    #try:
    import matplotlib.pyplot as plt
    #except ImportError:
    #    pass
    #else:
    gui_results( all_scores)

