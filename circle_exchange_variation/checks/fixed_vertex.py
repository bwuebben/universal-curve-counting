"""Finite reverse-plane-partition expansion of the polarized vertex."""

from collections import Counter


def multiply_ratio(series, degree, numerator, denominator):
    """Multiply by (1 - numerator*y^degree)/(1 - denominator*y^degree)."""
    result = series.copy()
    for k in range(degree, len(series)):
        result[k] += denominator * result[k - degree] - numerator * series[k - degree]
    return result


def reverse_plane_partitions(boxes, n):
    positions = {box:i for i,box in enumerate(boxes)}
    parents = [[positions[b] for b in ((i-1,j),(i,j-1)) if b in positions] for i,j in boxes]

    def extend(heights, left):
        if len(heights) == len(boxes):
            if not left:
                yield heights
            return
        least = max((heights[i] for i in parents[len(heights)]), default=0)
        for entry in range(least, left+1):
            yield from extend(heights+(entry,), left-entry)

    yield from extend((),n)


def polarization(boxes, heights):
    """Delta H as a Laurent character in y,p,c including its entire y-zero part."""
    result = Counter()

    def add(degree, axial, shift, sign):
        # (p^shift-1)/(1-p), multiplied by s^axial, with s=p/c.
        for e in range(min(0,shift),max(0,shift)):
            result[degree,e+axial,-axial] += sign if shift < 0 else -sign

    for (i,j), d in zip(boxes, heights):
        add(j-i,-j,-d,1)
        for (k,l), e in zip(boxes, heights):
            degree = j-i-l+k
            add(degree,l-j,e-d,-1)
            add(degree+1,l-j,e-d,1)
    assert not result[0,0,0] and not result[0,0,-1], "unexpected invariant Euler"
    return {key:val for key,val in result.items() if val}


def vertex(shape, p, c, bound):
    boxes = [(i,j) for j,width in enumerate(shape) for i in range(width)]
    s = p/c
    empty = [1]+[0]*bound
    weights = [p**0]
    for d in range(1,bound+1):
        weights.append(weights[-1]*(1-c*p**(d-1))/(1-p**d))
    for i,j in boxes:
        ell = sum(k==i and l>j for k,l in boxes)
        h = 1 + ell + sum(k>i and l==j for k,l in boxes)
        empty = [sum(empty[n-h*r]*weights[r]*s**(-ell*r) for r in range(n//h+1))
                 for n in range(bound+1)]
    rows = []
    for n in range(bound+1):
        total = [0]*(bound+1)
        for heights in reverse_plane_partitions(boxes,n):
            factor = c**(-n)
            polynomial = empty.copy()
            for (l,i,j), multiplicity in polarization(boxes,heights).items():
                z = p**i*c**j
                if not l:
                    factor *= ((1-c*z)/(1-z))**multiplicity
                    continue
                if l < 0:
                    factor *= c**multiplicity
                    l, z, multiplicity = -l, 1/(c*z), -multiplicity
                numerator, denominator = (c*z,z) if multiplicity > 0 else (z,c*z)
                for _ in range(abs(multiplicity)):
                    polynomial = multiply_ratio(polynomial,l,numerator,denominator)
            for k, value in enumerate(polynomial):
                total[k] += factor*value
        rows.append(total)
    return rows


def partitions(n, largest):
    if not n:
        yield ()
    else:
        for width in range(min(n,largest),0,-1):
            for rest in partitions(n-width,width):
                yield (width,)+rest
