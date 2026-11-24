"""Enumerate SL(2) orbits of exact prime-power order in the torus fiber."""

from collections import Counter
from itertools import product


def orbit_census(prime, power):
    modulus = prime ** power
    unseen = {
        entry for entry in product(range(modulus), repeat=4)
        if any(value % prime for value in entry)
    }
    census = Counter()
    while unseen:
        start = unseen.pop()
        todo = [start]
        size = 0
        while todo:
            a, b, c, d = todo.pop()
            size += 1
            # These two row additions generate SL(2) over this local ring.
            for other in (
                ((a + c) % modulus, (b + d) % modulus, c, d),
                (a, b, (c + a) % modulus, (d + b) % modulus),
            ):
                if other in unseen:
                    unseen.remove(other)
                    todo.append(other)
        a, b, c, d = start
        determinant = (a*d - b*c) % modulus
        valuation = 0
        while valuation < power and determinant % prime == 0:
            valuation += 1
            determinant //= prime
        census[valuation, size] += 1
    return census


def expected_census(prime, power):
    primitive_columns = prime ** (2 * power) - prime ** (2 * power - 2)
    rows = {}
    for valuation in range(power + 1):
        degree = primitive_columns * prime ** (power - valuation)
        residues = 1 if valuation == 0 else prime ** valuation + prime ** (valuation - 1)
        determinants = 1 if valuation == power else (prime - 1) * prime ** (power - valuation - 1)
        rows[valuation, degree] = residues * determinants
    return rows


if __name__ == "__main__":
    for prime, power in ((2, 1), (2, 2), (2, 3), (3, 1), (3, 2), (5, 1), (7, 1), (11, 1)):
        actual = orbit_census(prime, power)
        assert actual == expected_census(prime, power)
        assert sum(degree * count for (_, degree), count in actual.items()) == (
            prime ** (4 * power) - prime ** (4 * (power - 1))
        )
        print(f"exact order {prime ** power}: {dict(sorted(actual.items()))}")
