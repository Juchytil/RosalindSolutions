# My Rosalind Solutions

### for opening FASTA file if needed
def fasta_to_dict(fasta_file):
  with open(fasta_file, 'r') as file:
    fastas={}
    for line in file.readlines():
      if line [0]=='>':
        header=line[1:].strip()
        fastas[header]=''
      else:
        fastas[header]+=line[:-1]
    return fastas

## Counting Nucleotides (DNA)*

def nuc_count(seq):
    A_count = seq.count('A')
    T_count = seq.count('T')
    G_count = seq.count('G')
    C_count = seq.count('C')
    return A_count,T_count,G_count,C_count

## DNA to RNA (RNA)*

def dna_to_rna(seq):
  rna = seq.replace('T','U')
  return rna

## Reverse Compliment Strand (REVC)*

def reverse_compliment_seq(seq):
  compliment_lib = str.maketrans({'A':'T','T':'A','C':'G','G':'C'})
  rev_comp_seq = seq.translate(compliment_lib)[::-1]
  return rev_comp_seq

## GC Content Function (GC)*

def GC_content(sequence):
    seq_len = len(sequence)
    G_value = sequence.count('G')
    C_value = sequence.count('C')
    GC_value = G_value + C_value
    GC_percent = (GC_value/seq_len)*100
    return GC_percent

## Point mutation count (HAMM)*

def mutations_count(seq_a, seq_b):
  mutation = 0
  for i in range(len(seq_a)):
    if seq_a[i] != seq_b[i]:
      mutation += 1
  return mutation

## RNA to protein (PROT)*

rna_codon_dictionary = {
        "UUU" : "F","CUU" : "L","AUU" : "I","GUU" : "V",
        "UUC" : "F","CUC" : "L","AUC" : "I","GUC" : "V","UUA" : "L","CUA" : "L",
        "AUA" : "I","GUA" : "V","UUG" : "L","CUG" : "L","AUG" : "M","GUG" : "V",
        "UCU" : "S","CCU" : "P","ACU" : "T","GCU" : "A","UCC" : "S","CCC" : "P",
        "ACC" : "T","GCC" : "A","UCA" : "S","CCA" : "P","ACA" : "T","GCA" : "A",
        "UCG" : "S","CCG" : "P","ACG" : "T","GCG" : "A","UAU" : "Y","CAU" : "H",
        "AAU" : "N","GAU" : "D","UAC" : "Y","CAC" : "H","AAC" : "N","GAC" : "D",
        "UAA" : "STOP","CAA" : "Q","AAA" : "K","GAA" : "E","UAG" : "STOP","CAG" : "Q",
        "AAG" : "K","GAG" : "E","UGU" : "C","CGU" : "R","AGU" : "S","GGU" : "G",
        "UGC" : "C","CGC" : "R","AGC" : "S","GGC" : "G","UGA" : "STOP","CGA" : "R",
        "AGA" : "R","GGA" : "G","UGG" : "W","CGG" : "R","AGG" : "R","GGG" : "G"
}

def rna_to_prot(rna_seq):
  protien = ''
  for i in range(0,len(rna_seq)-3,3):
    if rna_codon_dictionary[rna_seq[i:i+3]] == 'STOP':
      break
    protien += rna_codon_dictionary[rna_seq[i:i+3]]
  return protien

### Protein to DNA(For Fun)*

import itertools as it
from itertools import product

prot_to_codon_dictionary = {
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

def prot_to_rnas(prot):
  needed_prot_dictionary = {}
  for letter in prot:
    if letter in prot_to_codon_dictionary:  
        needed_prot_dictionary[letter] = prot_to_codon_dictionary[letter]  
  value_lists = [needed_prot_dictionary[letter] for letter in prot if letter in needed_prot_dictionary]
  combinations = [''.join(combo) for combo in product(*value_lists)]
  return combinations

## Finding a Motif (SUBS)*

def motif_start_loc(seq,motif):
  for i in range(len(seq)):
    if seq[i:i+len(motif)] == motif:
      print(i+1)


## Calculating Protein Mass (PRTM)*

prot_mass_table = {
    'A' : '71.03711',  'C' : '103.00919', 'D' : '115.02694',
    'E' : '129.04259', 'F' : '147.06841', 'G' : '57.02146',
    'H' : '137.05891', 'I' : '113.08406', 'K' : '128.09496',
    'L' : '113.08406', 'M' : '131.04049', 'N' : '114.04293',
    'P' : '97.05276',  'Q' : '128.05858', 'R' : '156.10111',
    'S' : '87.03203',  'T' : '101.04768', 'V' : '99.06841',
    'W' : '186.07931', 'Y' : '163.06333'
}

def prot_mass(seq):
  protien_mass = 0
  for i in seq:
    protien_mass += float(prot_mass_table[i])
  return protien_mass

## RNA Splicing (SPLC)*

def rna_splicing(seq,introns):
  for intron in introns:
    seq = seq.replace(intron,'')
  rna_seq = dna_to_rna(seq)
  protein = rna_to_prot(rna_seq)
  return protein

## Locating Restriction Sites (REVP)*

def restriction_site(seq):
  low_range = 4
  high_range = 12
  for n in range(len(seq)-low_range):
    for i in range(low_range, high_range):
      if seq[n:n+i] == reverse_compliment_seq(seq[n:n+i]) and len(seq[n:n+i]) == i:
        print(n+1,len(seq[n:n+i]))

## Transitions and Transversions (TRAN)*

transition_list = [('A', 'G'), ('T', 'C'), ('C', 'T'), ('G', 'A')]
transversions_list = [('A', 'T'), ('A', 'C'), ('T', 'A'), ('T', 'G'), ('C', 'A'), ('C', 'G'), ('G', 'T'), ('G', 'C')]

def transition_transversion(seq1,seq2):
  transition = 0
  transversion = 0
  for i in range(len(seq1)):
    if (seq1[i], seq2[i]) in transition_list:
        transition += 1
    if (seq1[i], seq2[i]) in transversions_list:
        transversion += 1
  return transition/transversion

## Finding a Shared Motif (LCSM)*

def shared_motif(sequences_dict):
  list_of_seqs = list(sequences_dict.values())
  ref = list_of_seqs[0]
  other_seqs = list_of_seqs[1:]
  motifs = []
  for i in range(len(ref)):
    for n in range(len(ref)):
      if all(ref[i:i+n] in s for s in other_seqs):
        motifs.append(ref[i:i+n])
  return max(motifs, key=len)

## Open Reading Frames (ORF)*

dna_codons = {
	"TTT": 'F', "TTC": 'F', "TTA": 'L',    "TTG": 'L',
	"TCT": 'S', "TCC": 'S', "TCA": 'S',    "TCG": 'S',
	"TAT": 'Y', "TAC": 'Y', "TAA": 'Stop', "TAG": 'Stop',
	"TGT": 'C', "TGC": 'C', "TGA": 'Stop', "TGG": 'W',
	"CTT": 'L', "CTC": 'L', "CTA": 'L',    "CTG": 'L',
	"CCT": 'P', "CCC": 'P', "CCA": 'P',    "CCG": 'P',
	"CAT": 'H', "CAC": 'H', "CAA": 'Q',    "CAG": 'Q',
	"CGT": 'R', "CGC": 'R', "CGA": 'R',    "CGG": 'R',
	"ATT": 'I', "ATC": 'I', "ATA": 'I',    "ATG": 'M',
	"ACT": 'T', "ACC": 'T', "ACA": 'T',    "ACG": 'T',
	"AAT": 'N', "AAC": 'N', "AAA": 'K',    "AAG": 'K',
	"AGT": 'S', "AGC": 'S', "AGA": 'R',    "AGG": 'R',
	"GTT": 'V', "GTC": 'V', "GTA": 'V',    "GTG": 'V',
	"GCT": 'A', "GCC": 'A', "GCA": 'A',    "GCG": 'A',
	"GAT": 'D', "GAC": 'D', "GAA": 'E',    "GAG": 'E',
	"GGT": 'G', "GGC": 'G', "GGA": 'G',    "GGG": 'G'
}

dna_start_codon = "ATG"
dna_stop_codon = [key for key, val in DNA_CODON_MAP.items() if val == "Stop"]
dna_comps = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

rna_codons = {
	"UUU": 'F', "UUC": 'F', "UUA": 'L',    "UUG": 'L',
	"UCU": 'S', "UCC": 'S', "UCA": 'S',    "UCG": 'S',
	"UAU": 'Y', "UAC": 'Y', "UAA": 'Stop', "UAG": 'Stop',
	"UGU": 'C', "UGC": 'C', "UGA": 'Stop', "UGG": 'W',
	"CUU": 'L', "CUC": 'L', "CUA": 'L',    "CUG": 'L',
	"CCU": 'P', "CCC": 'P', "CCA": 'P',    "CCG": 'P',
	"CAU": 'H', "CAC": 'H', "CAA": 'Q',    "CAG": 'Q',
	"CGU": 'R', "CGC": 'R', "CGA": 'R',    "CGG": 'R',
	"AUU": 'I', "AUC": 'I', "AUA": 'I',    "AUG": 'M',
	"ACU": 'T', "ACC": 'T', "ACA": 'T',    "ACG": 'T',
	"AAU": 'N', "AAC": 'N', "AAA": 'K',    "AAG": 'K',
	"AGU": 'S', "AGC": 'S', "AGA": 'R',    "AGG": 'R',
	"GUU": 'V', "GUC": 'V', "GUA": 'V',    "GUG": 'V',
	"GCU": 'A', "GCC": 'A', "GCA": 'A',    "GCG": 'A',
	"GAU": 'D', "GAC": 'D', "GAA": 'E',    "GAG": 'E',
	"GGU": 'G', "GGC": 'G', "GGA": 'G',    "GGG": 'G'
}

rna_start_codon = "AUG"
rna_stop_codon = [key for key, val in RNA_CODON_MAP.items() if val == "Stop"]
rna_comps = {'A': 'U', 'U': 'A', 'C': 'G', 'G': 'C'}

dna = ''.join(dna_list[0].content)
dna_rc = ''.join(dna_comps[c] for c in reversed(dna))

def find_proteins(s):
	starts = []
	for i in range(len(s)-2):
		if s[i:i+3] == dna_start_codon:
			starts.append(i)
	proteins = []
	for start in starts:
		cur_protein = []
		for i in range(start, len(s)-(len(s)-start)%3, 3):
			aa = dna_codons[s[i:i+3]]
			if aa == "Stop":
				proteins.append(''.join(cur_protein))
				break
			cur_protein.append(aa)
	return proteins

# all_orfs(test)

## Rabbits and Recurrence Relations (FIB)*

def rabbits(months,litter):
  month = 0
  pairs = 1
  onemonthpairs = 0
  twomonthpairs = 0
  while month < months:
    pairs += onemonthpairs
    onemonthpairs = 0 + twomonthpairs
    twomonthpairs = 0 + pairs*litter
    month += 1
  return pairs, onemonthpairs, twomonthpairs

## Mortal FIbonaccie Rabbits (FIBD) WORK IN PROGRESS

# def rabbits(months,litter,lifespan):
#   pairs_at_month = []
#   month = 0
#   pairs = 1
#   onemonthpairs = 0
#   twomonthpairs = 0
#   while month < months:
#     pairs += onemonthpairs - pairs_at_month[lifespan]
#     pairs_at_month.append(pairs)
#     onemonthpairs = 0 + twomonthpairs
#     twomonthpairs = 0 + pairs*litter
#     month += 1
#     if len(pairs_at_month) >= lifespan:
#       lifespan += 1
#     print(pairs_at_month)
#   return pairs, onemonthpairs, twomonthpairs

# rabbits(5,1,3)

## Mendels First Law (IPRB)*

def mendelsfirst(dom,het,rec):
  total = dom + het + rec
  dompickone = dom / total
  hetpickone =  het / total
  recpickone = rec / total
  hetdom = hetpickone*(dom/(total - 1))
  hethet = hetpickone*((het - 1)/(total - 1))
  hetrec = hetpickone*(rec/(total - 1))
  recdom = recpickone*(dom/(total - 1))
  rechet = recpickone*(het/(total -  1))
  percentdom = dompickone + (hetdom) + (hethet*0.75) + (hetrec*0.5) + (recdom) + (rechet*0.5)
  return percentdom

## Calculating Expected Offspring (IEV)*

def expected_offspring(genotype1,genotype2,genotype3,genotype4,genotype5,genotype6):
  dom_offspring = 0
  dom_offspring += (genotype1*2) + (genotype2*2) + (genotype3*2) + ((genotype4*2)*0.75) + ((genotype5*2)*0.5) + (genotype6*0)
  return dom_offspring

## Consensus and Profile (CONS)* CORRECT BUT UNABLE TO GET FORMATING FOR ROSALIND WEBSITE

import pandas as pd

def consensus_profile(fastadict):
  row_labels = ['A', 'T', 'C', 'G']
  list_of_seqs = list(fastadict.values())  
  nucleotide_map = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
  result = {i: [0, 0, 0, 0] for i in range(len(list_of_seqs[0]))}
  for seqs in list_of_seqs:
    for i, char in enumerate(seqs):
      if char in nucleotide_map:
          result[i][nucleotide_map[char]] += 1
  profiledf = pd.DataFrame.from_dict(result, orient='index', columns=row_labels).T 
  consensus_sequence = "".join(profiledf.idxmax())
  return profiledf, consensus_sequence
