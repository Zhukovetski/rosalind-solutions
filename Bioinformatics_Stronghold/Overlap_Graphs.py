import pyperclip as pp

import fasta_reader as fr
import read_rosalind as rr

s = rr.read_latest_rosalind_file()
# s = """>Rosalind_0498
#         AAATAAA
#         >Rosalind_2391
#         AAATTTT
#         >Rosalind_2323
#         TTTTCCC
#         >Rosalind_0442
#         AAATCCC
#         >Rosalind_5013
#         GGGTGGG"""

list_seq = fr.reader(s)
dict_seq = {list_seq[i]: list_seq[i + 1] for i in range(0, len(list_seq) - 1, 2)}
adj_list = []
k = 3

for name1, seq1 in dict_seq.items():
    for name2, seq2 in dict_seq.items():
        if seq1.endswith(seq2[:k]) and name1 != name2:
            adj_list.append(name1[1:] + " " + name2[1:])

result = "\n".join(adj_list)
pp.copy(result)
print(result)
