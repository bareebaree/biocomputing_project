# -*- coding: utf-8 -*-
"""
Created on Wed Mar 26 21:31:44 2025

@author: james
"""

from Bio import pairwise2
import random

def shuffle_sequence(seq) -> str:
    """Shuffle a sequence while preserving its composition.
    
    PARAMETERS: seq ->(str): Any given sequence from the ID:sequence dictionary
    
    RETURNS:  Shuffled sequence appended to a shuffled sequence list.
    
    """
    seq_list = list(seq)
    random.shuffle(seq_list)
    return "".join(seq_list)

def shuffle_and_align(query_seq: str, best_match_seq: str) -> float:
    """Shuffles the query sequence and computes an alignment score. Calls shuffle_sequence() to do so.
       Scoring scheme = match = +1
       gap/mismatch = 0
       
    PARAMETERS: query_seq -> (str): The query sequence.
    
    RETURNS: The score of each shuffled query 
    
    """
    shuffled_query = shuffle_sequence(query_seq)
    return pairwise2.align.globalxx(shuffled_query, best_match_seq, score_only=True)
