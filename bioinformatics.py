# My Rosalind Solutions

sequence = 'ATGCGGGCGAGCGTTTCGGAGGGTATTTATTATCTTTCTATCATTTTTTAGGGGAGGATTTTAGGGGATTATCTCTCGATCGATTATCGATCC'

# Counting Nucleotides (DNA)

def nuc_count(seq):
    A_count = seq.count('A')
    T_count = seq.count('T')
    G_count = seq.count('G')
    C_count = seq.count('C')
    counts = (A_count,T_count,G_count,C_count)
    return counts

## DNA to RNA (RNA)

def dna_to_rna(seq):
  rna = seq.replace('T','U')
  return rna

## Reverse Compliment Strand (REVC)

def reverse_compliment_seq(seq):
  compliment_lib = str.maketrans({'A':'T','T':'A','C':'G','G':'C'})
  rev_comp_seq = seq.translate(compliment_lib)[::-1]
  return rev_comp_seq

## GC Content Function (GC)

def GC_content(sequence):
    seq_len = sequence.count('')
    G_value = sequence.count('G')
    C_value = sequence.count('C')
    GC_value = G_value + C_value
    GC_percent = (GC_value/seq_len)*100
    return GC_percent

## Point mutation count (HAMM)

def mutations_count(seq_a, seq_b):
  mutation = 0
  for i in range(len(seq_a)):
    if seq_a[i] != seq_b[i]:
      mutation += 1
  return mutation

## RNA to protein (PROT)

RNA_codon_table = {
'A': ('GCU', 'GCC', 'GCA', 'GCG'),
'C': ('UGU', 'UGC'),
'D': ('GAU', 'GAC'),
'E': ('GAA', 'GAG'),
'F': ('UUU', 'UUC'),
'G': ('GGU', 'GGC', 'GGA', 'GGG'),
'H': ('CAU', 'CAC'),
'I': ('AUU', 'AUC', 'AUA'),
'K': ('AAA', 'AAG'),
'L': ('UUA', 'UUG', 'CUU', 'CUC', 'CUA', 'CUG'),
'M': ('AUG',),
'N': ('AAU', 'AAC'),
'P': ('CCU', 'CCC', 'CCA', 'CCG'),
'Q': ('CAA', 'CAG'),
'R': ('CGU', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'),
'S': ('UCU', 'UCC', 'UCA', 'UCG', 'AGU', 'AGC'),
'T': ('ACU', 'ACC', 'ACA', 'ACG'),
'V': ('GUU', 'GUC', 'GUA', 'GUG'),
'W': ('UGG',),
'Y': ('UAU', 'UAC'),}

rna_sequence = sequence.replace('T','U')

def rna_to_prot(rna_seq):
  for i in range(0,len(rna_seq),3):
    print(rna_seq[i:i+3])

print(rna_to_prot(rna_sequence))

