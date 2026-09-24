import pytest
from helix.core import analyze_six_frames, clean_sequence, reverse_complement, translate_frame, GENETIC_CODES

def test_clean_fasta_and_rna():
    assert clean_sequence('>x\nAUG GCU 123') == 'ATGGCT'

def test_invalid_character():
    with pytest.raises(ValueError):
        clean_sequence('ATGZAA')

def test_reverse_complement():
    assert reverse_complement('ATGCCG') == 'CGGCAT'

def test_standard_translation():
    table = GENETIC_CODES['1'][1]
    calls = translate_frame('ATGGCTTAA', 0, table)
    assert ''.join(c.aa1 for c in calls) == 'MA*'

def test_six_frames_returned():
    result = analyze_six_frames('ATGGCTTAA', min_orf_length_nt=9)
    assert len(result['frames']) == 6
    assert result['length'] == 9
    assert result['frames'][0]['protein'] == 'MA*'
    assert result['total_orfs'] >= 1

def test_ambiguous_codon_is_unknown():
    result = analyze_six_frames('ATGNNNTAA', min_orf_length_nt=9)
    assert result['frames'][0]['protein'] == 'M?*'
