import pyperclip

import read_rosalind as rr

s = rr.read_latest_rosalind_file()
nums = [int(i) for i in s.split()]
t = sum(nums)
k, m, n = nums

prob_rec = (m * (m - 1) * 0.25 + m * n * 0.5 + n * m * 0.5 + n * (n - 1) * 1) / (
    t * (t - 1)
)


pyperclip.copy(str(1 - prob_rec))
