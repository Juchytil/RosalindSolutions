# My Rosalind Solutions

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

rna_codon_table = {"UUU" : "F", "CUU" : "L", "AUU" : "I", "GUU" : "V",
           "UUC" : "F", "CUC" : "L", "AUC" : "I", "GUC" : "V",
           "UUA" : "L", "CUA" : "L", "AUA" : "I", "GUA" : "V",
           "UUG" : "L", "CUG" : "L", "AUG" : "M", "GUG" : "V",
           "UCU" : "S", "CCU" : "P", "ACU" : "T", "GCU" : "A",
           "UCC" : "S", "CCC" : "P", "ACC" : "T", "GCC" : "A",
           "UCA" : "S", "CCA" : "P", "ACA" : "T", "GCA" : "A",
           "UCG" : "S", "CCG" : "P", "ACG" : "T", "GCG" : "A",
           "UAU" : "Y", "CAU" : "H", "AAU" : "N", "GAU" : "D",
           "UAC" : "Y", "CAC" : "H", "AAC" : "N", "GAC" : "D",
           "UAA" : "STOP", "CAA" : "Q", "AAA" : "K", "GAA" : "E",
           "UAG" : "STOP", "CAG" : "Q", "AAG" : "K", "GAG" : "E",
           "UGU" : "C", "CGU" : "R", "AGU" : "S", "GGU" : "G",
           "UGC" : "C", "CGC" : "R", "AGC" : "S", "GGC" : "G",
           "UGA" : "STOP", "CGA" : "R", "AGA" : "R", "GGA" : "G",
           "UGG" : "W", "CGG" : "R", "AGG" : "R", "GGG" : "G"}

def rna_to_prot(rna_seq):
  protien = ''
  for i in range(0,len(rna_seq),3):
    if rna_codon_table[rna_seq[i:i+3]] == 'STOP':
      break
    protien += rna_codon_table[rna_seq[i:i+3]]
  return protien


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

## Transitions and Transversions (TRAN)

# def transition_transversion_ratio(seqa,seqb):
#   transversion = 0
#   transition = 0
#   for i in range(len(seqa))
#     if seqa[i] and 

## Finding a Shared Motif (LCSM)

## Open Reading Frames (ORF)

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

## Mortal FIbonaccie Rabbits (FIBD)

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