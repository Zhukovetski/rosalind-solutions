import pyperclip as pp

import read_rosalind as rr

string = rr.read_latest_rosalind_file()
# string = "We tried list and we tried dicts also we tried Zen"

list_str = string.split()
result = {}
for s in list_str:
    result[s] = 1 + result.get(s, 0)
str_result = []
for key, num in result.items():
    str_result.append(f"{key} {num}\n")
str_result = "".join(str_result)
print(str_result)
pp.copy(str_result)
