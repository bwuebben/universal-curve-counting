"""Check the elementary hook coefficients in the fixed-excess denominator formula."""

from fractions import Fraction
from math import factorial, gcd, isqrt, lcm, prod

from etale_spin_characters import partitions


def todd_denominator(weight):
    return prod(
        p ** (weight // (p - 1))
        for p in range(2, weight + 2)
        if all(p % q for q in range(2, isqrt(p) + 1))
    )


def convolution(left, right, order):
    return [
        sum(left[j] * right[i - j] for j in range(i + 1))
        for i in range(order + 1)
    ]


def connected_rows(genus, max_degree, excess):
    """Ordinary log in degree of the series with its minimum t-power removed."""
    z = [[Fraction(0) for _ in range(excess + 1)] for _ in range(max_degree + 1)]
    connected = [row.copy() for row in z]
    z[0][0] = Fraction(1)
    for degree in range(1, max_degree + 1):
        for shape in partitions(degree):
            block = [Fraction(1)] + [Fraction(0)] * excess
            for i, length in enumerate(shape):
                for j in range(length):
                    hook = length - j + sum(other > j for other in shape[i + 1 :])
                    sine_square = [
                        Fraction(2 * (-1) ** k * hook ** (2 * k + 2), factorial(2 * k + 2))
                        for k in range(excess + 1)
                    ]
                    for _ in range(genus - 1):
                        block = convolution(block, sine_square, excess)
            z[degree] = [a + b for a, b in zip(z[degree], block)]
        # Z times its logarithmic derivative equals the derivative of Z.
        connected[degree] = z[degree].copy()
        for j in range(1, degree):
            product = convolution(connected[j], z[degree - j], excess)
            for k in range(excess + 1):
                connected[degree][k] -= Fraction(j, degree) * product[k]
    return connected


def main():
    excess = 4
    for weight in range(1, 25):
        assert (Fraction(2, factorial(2 * weight + 2)) * todd_denominator(2 * weight)).denominator == 1
    for genus in (2, 3, 4):
        for degree, row in enumerate(connected_rows(genus, 8, excess)[1:], start=1):
            for r in range(1, excess + 1):
                bound = gcd(degree, r) * todd_denominator(2 * r)
                assert bound % row[r].denominator == 0, (genus, degree, r)
            if degree in (2, 3, 5, 7):
                assert row[0].denominator == degree
    example = connected_rows(2, 4, excess)
    witnesses = {
        1: {1: Fraction(-1, 12)},
        2: {1: Fraction(1, 360), 2: Fraction(19, 32)},
        3: {1: Fraction(-1, 20160), 3: Fraction(-8278, 567)},
        4: {1: Fraction(1, 1814400), 2: Fraction(211, 48384), 4: Fraction(233871691, 193536)},
    }
    for r, entries in witnesses.items():
        for degree, expected in entries.items():
            assert example[degree][r] == expected
        common = lcm(*(value.denominator for value in entries.values()))
        assert common == r * todd_denominator(2 * r)
        print(f"Excess {r}: optimal common denominator {common}")
    base = connected_rows(2, 2, 1)[2][1]
    double_core = connected_rows(3, 1, 1)[1][1]
    residual = (base - 8 * double_core) % 1
    assert (base, double_core, residual) == (Fraction(-13, 4), Fraction(-1, 6), Fraction(1, 12))
    print(f"g=2, degree 2, genus 4: equal-spin residue {base % 1}, distinct-spin residue {residual}")
    print("Verified elementary coefficients through degree 8 and excess 4 for genera 2, 3 and 4.")


if __name__ == "__main__":
    main()
