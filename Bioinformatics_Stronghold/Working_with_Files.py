import re

import pyperclip as pp

import read_rosalind as rr

string = rr.read_latest_rosalind_file()
# string = """Bravely bold Sir Robin rode forth from Camelot
# Yes, brave Sir Robin turned about
# He was not afraid to die, O brave Sir Robin
# And gallantly he chickened out
# He was not at all afraid to be killed in nasty ways
# Bravely talking to his feet
# Brave, brave, brave, brave Sir Robin
# He beat a very brave retreat"""

string = re.split(r"(?<=\n)", string)
string = "".join([string[i] for i in range(1, len(string), 2)])
pp.copy(string)
print(string)
