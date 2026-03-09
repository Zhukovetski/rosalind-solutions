
def solve_task(content):

    A, C, G, T = 0, 0, 0, 0

    for x in content:
        if x == "A":
            A += 1
        elif x == "C":
            C += 1
        elif x == "G":
            G += 1
        else:
            T += 1
    return str(A) + " " + str(C) + " " + str(G) + " " + str(T)