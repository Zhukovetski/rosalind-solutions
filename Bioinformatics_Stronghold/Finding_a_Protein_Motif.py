import read_rosalind as rr
import pyperclip as pp
import requests as rq
import re
import time

raw_data = rr.read_latest_rosalind_file()
# raw_data = """A2Z669
# B5ZC00
# P07204_TRBM_HUMAN
# P20840_SAG1_YEAST"""

url = "http://www.uniprot.org/uniprot/{}.fasta"
motif = r"(?=(N[^P][ST][^P]))"
list_ID = raw_data.split()
fin_result = []
for id in list_ID:
    i = id.split("_")[0]
    response = rq.get(url.format(i))
    if response.status_code == 200:
        data = "".join(line.strip() for line in response.text.splitlines()[1:])
        result = []
        
        for m in re.finditer(motif, data):
            result.append(str(m.start() + 1))
        if result:
            fin_result.append(id)
            fin_result.append(" ".join(result))
        time.sleep(1)
    else:
        print(f"Ошибка при загрузке {response.status_code}")

fin_result = "\n".join(fin_result)
pp.copy(fin_result)
print(fin_result)