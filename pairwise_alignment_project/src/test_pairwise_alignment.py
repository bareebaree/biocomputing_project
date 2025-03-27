
import pytest
from pairwise_functions import shuffle_sequence, shuffle_and_align

# 1. Test shuffle_sequence maintains composition
def test_shuffle_sequence_composition():
    seq = "AAAGTC"
    shuffled = shuffle_sequence(seq)
    assert sorted(seq) == sorted(shuffled)
    assert set(seq) == set(shuffled)
    assert len(seq) == len(shuffled)
    assert seq != shuffled  # It should usually not be the same

# 2. Test shuffle_and_align returns a valid score
def test_shuffle_and_align_score_type():
    query = "ACGT"
    target = "ACGT"
    score = shuffle_and_align(query, target)
    assert isinstance(score, float) or isinstance(score, int)

# 3. Test known score ranges
def test_shuffle_and_align_score_range():
    query = "AAAA"
    target = "TTTT"  # no matches expected
    score = shuffle_and_align(query, target)
    assert score == 0  # No matches, so score should be 0
