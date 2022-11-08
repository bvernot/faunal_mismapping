import sys
from collections import Counter

cnt = Counter()

for line in sys.stdin:

    if line.startswith('<'):
        l1 = line.rstrip()
    elif line.startswith('>'):
        l2 = line.rstrip()
        print(l1)
        print(l2)
        print(''.join('.' if l1[i] == l2[i] else 'X' for i in range(len(l1))))
        subs = [l1[i] + l2[i] for i in range(len(l1)) if l1[i] != l2[i]]
        print(subs)
        cnt.update(subs)

        print(cnt)
        print()
        continue
    pass
