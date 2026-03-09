import read_rosalind as rr
import pyperclip as pp
import fasta_reader as fr

s = rr.read_latest_rosalind_file()
# s = """>Rosalind_1
#         GATTACA
#         >Rosalind_2
#         TAGACCA
#         >Rosalind_3
#         ATACA"""




list_seq = fr.reader(s)
list_seq = [i for i in list_seq[1::2]]
prev_common = ["A", "T", "G", "C"]
result = []
common_sub = 1
while common_sub:
    common_sub = []
    for s1 in prev_common:
        for s2 in "ATGC":
            s3 = s1 + s2
            for seq in list_seq:
                if seq.find(s3) == -1:
                    break
            else:
                common_sub.append(s3)
    if not common_sub:
        break

    result = common_sub
    prev_common = common_sub
    
pp.copy(result[0])
print(result)

