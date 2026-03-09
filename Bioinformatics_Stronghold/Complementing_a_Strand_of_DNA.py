
def solve_task(content):
    return content.translate(str.maketrans("ATCG", "TAGC"))[::-1]