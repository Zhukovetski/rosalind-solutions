import read_rosalind as rr
import pyperclip as pp

string = rr.read_latest_rosalind_file()
#string = "GATATATGCATATACTT ATAT"
seq1, seq2 = string.split()

alloc_sub = []
len1, len2 = len(seq1), len(seq2) 
for i in range(0, len1 - len2):
    if seq1[i: i+len2] == seq2:
        alloc_sub.append(str(i + 1))

result = " ".join(alloc_sub)
pp.copy(result)
print(" ".join(alloc_sub))

