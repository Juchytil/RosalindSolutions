# My Rosalind Solutions

## Counting Nucleotides

def nuc_count(seq):
    A_count = seq.count('A')
    T_count = seq.count('T')
    G_count = seq.count('G')
    C_count = seq.count('C')
    counts = (A_count,T_count,G_count,C_count)
    return counts

## DNA to RNA

def dna_to_rna(seq):
  rna = seq.replace('T','U')
  return rna

## GC Content Function

def GC_content(sequence):
    seq_len = sequence.count('')
    G_value = sequence.count('G')
    C_value = sequence.count('C')
    GC_value = G_value + C_value
    GC_percent = (GC_value/seq_len)*100
    return GC_percent

## Reverse Compliment Strand

def reverse_compliment_seq(seq):
  compliment_lib = str.maketrans({'A':'T','T':'A','C':'G','G':'C'})
  rev_comp_seq = seq.translate(compliment_lib)[::-1]
  return rev_comp_seq

