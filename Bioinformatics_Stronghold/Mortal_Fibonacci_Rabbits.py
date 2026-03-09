import pyperclip as pp

import read_rosalind as rr

s = rr.read_latest_rosalind_file()
# s = "6 3"
n, m = map(int, s.split())
ages = [0] * m
ages[0] = 1

for month in range(n - 1):
    newborn = sum(ages[1:])
    for i in range(m - 2, -1, -1):
        ages[i + 1] = ages[i]
    ages[0] = newborn

result = sum(ages)
pp.copy(result)
print(result)
