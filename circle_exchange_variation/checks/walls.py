"""Exact checks for antidiagonal polarization and graded moment-map walls.

Run with python3 checks/walls.py.
"""
from collections import Counter


def partitions(n, maximum):
    if n == 0:
        yield ()
    for head in range(min(n, maximum), 0, -1):
        for tail in partitions(n-head, head):
            yield (head,) + tail


def reverse(p):
    return Counter({-i: n for i, n in p.items()})


def product(p, r):
    out = Counter()
    for i, a in p.items():
        for j, b in r.items():
            out[i+j] += a*b
    return out


def sub(p, r):
    out = p.copy()
    out.subtract(r)
    return out


def stalk(v, q):
    change = sub(product(q, reverse(v)), product(v, reverse(q)))
    f = sub(q, product({0:1, 1:-1}, change))
    return sub(reverse(f), f)


def check_stalk(v):
    for i in v:
        q = Counter({i:1})
        change = stalk(v,q)
        assert change == reverse(stalk(reverse(v), reverse(q)))
        assert change == sub(Counter(), stalk(reverse(v), reverse(q)))
        if v == reverse(v):
            expected = sub(reverse(q),q)
            expected.subtract(product(expected.copy(), product({0:2,1:-1,-1:-1},v)))
            assert change == expected


def verify_reflections(bits, left, right):
    # Distinguish spectral labels, including equal eigenvalues at the wall.
    labels = list(bits)
    ordered = [j for j in labels if bits[j] and j != right]
    ordered += [left,right]
    ordered += [j for j in labels if not bits[j] and j != left]
    charges = {j: int(j==right)-int(j==left) for j in labels}
    for i, target in enumerate(ordered):
        pos = labels.index(target)
        while pos > i:
            a, b = labels[pos-1:pos+1]
            assert {a,b} != {left,right}
            labels[pos-1:pos+1] = [b,a]
            pos -= 1
    assert [bits[j] for j in labels] == [
        int(i<0 or i==1) for i in bits
    ]
    assert [charges[j] for j in labels] == [
        int(i==1)-int(i==0) for i in bits
    ]


def check_shape(rows):
    boxes = {(i,j) for j,r in enumerate(rows) for i in range(r)}
    v = Counter(j-i for i,j in boxes)
    check_stalk(v)
    indices = range(min(v), max(v)+2)
    bits = {i:int(i<=0)-v[i]+v[i-1] for i in indices}
    assert set(bits.values()) <= {0,1}
    hooks = Counter()
    for x,y in boxes:
        arm = sum(i>x and j==y for i,j in boxes)
        leg = sum(j>y and i==x for i,j in boxes)
        hooks[arm+leg+1,leg] += 1

    walls = Counter()
    for a in indices:
        for d in range(a+1,indices.stop):
            if bits[a] or not bits[d]:
                continue
            w = Counter({i:1 for i in range(a,d)})
            mu = sub(v,w)
            boundary = sub(Counter({0:1}), product({0:2,1:-1,-1:-1},reverse(mu)))
            h = product(w,boundary)
            leg = sum(not bits[i] for i in range(a+1,d))
            assert h[0] == 1
            assert sum(h.values()) == d-a
            assert sum(i*n for i,n in h.items()) == sum(w)
            assert sum(n for i,n in h.items() if i>0) == leg
            verify_reflections(bits,a,d)
            walls[d-a,leg] += 1
    assert walls == hooks
    return sum(walls.values())


if __name__ == "__main__":
    shapes = [s for degree in range(1, 14) for s in partitions(degree, degree)]
    total = sum(check_shape(s) for s in shapes)
    print(f"Verified {total} hook intervals in {len(shapes)} diagrams through degree 13.")
    print("Checked reflections, degree transport, wall determinants, and stalk cancellation.")
