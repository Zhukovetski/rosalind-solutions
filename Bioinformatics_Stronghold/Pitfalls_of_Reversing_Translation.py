import pyperclip as pp

import read_rosalind as rr

contain = rr.read_latest_rosalind_file()
# contain = "MA"

codon_table = {
    "F": 2,
    "L": 6,
    "I": 3,
    "V": 4,
    "M": 1,
    "S": 6,
    "P": 4,
    "T": 4,
    "A": 4,
    "Y": 2,
    "H": 2,
    "N": 2,
    "D": 2,
    "Stop": 3,
    "Q": 2,
    "K": 2,
    "E": 2,
    "C": 2,
    "R": 6,
    "G": 4,
    "W": 1,
}
result = 3
for char in contain:
    result *= codon_table[char]
    result %= 1_000_000
print(result)
pp.copy(result)
