import read_rosalind as rr
import pyperclip as pp

s = rr.read_latest_rosalind_file()
#s = "1 0 0 1 0 1"
nums = [int(i) for i in s.split()]
offset = (nums[0] + nums[1] + nums[2] + nums[3] * 0.75 + nums[4] * 0.5) * 2
pp.copy(offset)
print(offset)

