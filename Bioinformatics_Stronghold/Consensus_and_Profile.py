import read_rosalind as rr
import pyperclip as pp
import fasta_reader as fr

string = rr.read_latest_rosalind_file()
# string = """>Rosalind_1
#             ATCCAGCT
#             >Rosalind_2
#             GGGCAACT
#             >Rosalind_3
#             ATGGATCT
#             >Rosalind_4
#             AAGCAACC
#             >Rosalind_5
#             TTGGAACT
#             >Rosalind_6
#             ATGCCATT
#             >Rosalind_7
#             ATGGCACT"""

list_seq = fr.reader(string)
print(list_seq)
len1 = len(list_seq[1])
len2 = len(list_seq) // 2
matrix = [list(i) for i in list_seq[1::2]]
prof_matrix = []

for i in "ACGT":
    line_amount = []
    line_amount.append(f"{i}:")
    for j in range(len1):
        amount = 0
        for k in range(len2):
            if matrix[k][j] == i:
               amount += 1
        line_amount.append(amount)
    prof_matrix.append(line_amount)

dict_s = {0: "A", 1: "C", 2: "G", 3: "T"}
consen_s = []

for i in range(1, len1+1):
    index = 0
    max_ = 0
    for j in range(4):
        val = prof_matrix[j][i]
        if  val > max_:
            max_ = val
            index = j
    consen_s.append(dict_s[index])
result = "".join(consen_s)

for i in prof_matrix:
    result += "\n" + " ".join(str(j) for j in i) 
print(result)
pp.copy(result)