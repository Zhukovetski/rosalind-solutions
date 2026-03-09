import pyperclip as pp

import read_rosalind as rr

seq = rr.read_latest_rosalind_file()
seq = "AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA"
codon_table = { 'UUU': 'F',    'CUU': 'L',    'AUU': 'I',    'GUU': 'V',    
                'UUC': 'F',    'CUC': 'L',    'AUC': 'I',    'GUC': 'V',    
                'UUA': 'L',    'CUA': 'L',    'AUA': 'I',    'GUA': 'V',    
                'UUG': 'L',    'CUG': 'L',    'AUG': 'M',    'GUG': 'V',    
                'UCU': 'S',    'CCU': 'P',    'ACU': 'T',    'GCU': 'A',    
                'UCC': 'S',    'CCC': 'P',    'ACC': 'T',    'GCC': 'A',    
                'UCA': 'S',    'CCA': 'P',    'ACA': 'T',    'GCA': 'A',    
                'UCG': 'S',    'CCG': 'P',    'ACG': 'T',    'GCG': 'A',    
                'UAU': 'Y',    'CAU': 'H',    'AAU': 'N',    'GAU': 'D',    
                'UAC': 'Y',    'CAC': 'H',    'AAC': 'N',    'GAC': 'D',    
                'UAA': 'Stop', 'CAA': 'Q',    'AAA': 'K',    'GAA': 'E',    
                'UAG': 'Stop', 'CAG': 'Q',    'AAG': 'K',    'GAG': 'E',    
                'UGU': 'C',    'CGU': 'R',    'AGU': 'S',    'GGU': 'G',    
                'UGC': 'C',    'CGC': 'R',    'AGC': 'S',    'GGC': 'G',    
                'UGA': 'Stop', 'CGA': 'R',    'AGA': 'R',    'GGA': 'G',    
                'UGG': 'W',    'CGG': 'R',    'AGG': 'R',    'GGG': 'G'}

protein_string = []
for i in range(0, len(seq)-3, 3):
    protein_string.append(codon_table[seq[i:i+3]])
protein_string = "".join(protein_string)
pp.copy(protein_string)
print(protein_string)