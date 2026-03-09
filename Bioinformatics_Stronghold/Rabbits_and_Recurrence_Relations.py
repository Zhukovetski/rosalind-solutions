def solve_task(s):
    n, k = map(int, s.split())
    prev2, prev1 = 1, 1
    for _ in range(n - 2):
        current = prev1 + (prev2 * k)
        prev2, prev1 = prev1, current
    return prev1
