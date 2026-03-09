
def reader(string):
    list_s = string.split()
    result = []
    result.append(list_s[0])
    seq = []
    for i in list_s[1:]:
        if i[0] == ">":         
            result.append("".join(seq))
            result.append(i)
            seq = []
        else:
            seq.append(i)
    else:
        result.append("".join(seq))
    return result