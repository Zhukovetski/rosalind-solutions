import math

import pyperclip as pp

import read_rosalind as rr

string = rr.read_latest_rosalind_file()
# string = "2 1"
k, N = string.split()
k, N = int(k), int(N)
T = 2**k
total_prob = 0
if N > T / 2:
    for i in range(N, T + 1):
        c = math.comb(T, i)
        total_prob += c * 0.25**i * 0.75 ** (T - i)
else:
    for i in range(N):
        total_prob += math.comb(T, i) * 0.25**i * 0.75 ** (T - i)
    total_prob = round(1 - total_prob, 3)

pp.copy(total_prob)
print(total_prob)
