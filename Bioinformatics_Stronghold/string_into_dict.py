string = """UUU F      CUU L      AUU I      GUU V
UUC F      CUC L      AUC I      GUC V
UUA L      CUA L      AUA I      GUA V
UUG L      CUG L      AUG M      GUG V
UCU S      CCU P      ACU T      GCU A
UCC S      CCC P      ACC T      GCC A
UCA S      CCA P      ACA T      GCA A
UCG S      CCG P      ACG T      GCG A
UAU Y      CAU H      AAU N      GAU D
UAC Y      CAC H      AAC N      GAC D
UAA Stop   CAA Q      AAA K      GAA E
UAG Stop   CAG Q      AAG K      GAG E
UGU C      CGU R      AGU S      GGU G
UGC C      CGC R      AGC S      GGC G
UGA Stop   CGA R      AGA R      GGA G
UGG W      CGG R      AGG R      GGG G """

containt = string.split()
len1, len2 = 0, 0
for s in string[::2]:
    if len(s) > len1:
        len1 = len(s)
for s in string[1::2]:
    if len(s) > len2:
        len2 = len(s)
string = [
    f"'{i}': '{j}',{' ' * (len2 - len(j))}" for i, j in zip(string[::2], string[1::2])
]
for i in range(4, len(string), 4):
    string[i] = f"\n{string[i]}"
string = " ".join(string)
# _dict= {string[i]: string[i+1] for i in range(0, len(string)-1)}
# pp.copy(string)
# print(string)
dict_codon = {}
for i in range(1, len(containt), 2):
    dict_codon[containt[i]] = dict_codon.get(containt[i], 0) + 1
list_codon = [f"'{i}': {j},  " for i, j in dict_codon.items()]
for i in range(3, len(list_codon), 3):
    list_codon[i] = f"\n{list_codon[i]}"
print("".join(list_codon))
