# -*- coding: utf-8 -*-
"""
Created on Tue Mar  4 16:47:01 2025

@author: james
"""

from Bio import pairwise2

best_match = None
best_score = float('-inf')

for name, seq in database.items():
    score = pairwise2.align.globalxx(query_seq, seq, score_only=True)
    
    if score > best_score:
        best_score = score
        best_match = name