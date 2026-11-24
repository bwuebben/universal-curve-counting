"""Exact second-excess reductions; no computation of the unknown Euler integral.

The geometric component decomposition is in Note 059.  This checks the
coefficient extraction and a compact projective-bundle calculation separately.
The latter is deliberately NOT identified with the stable-map integral.
"""

from fractions import Fraction as Q

from fixed_excess_denominators import connected_rows


def sine_power_row(power):
    """Coefficients through t^4 after removing t^power from (2 sin(t/2))^power."""
    return [Q(1), Q(-power, 24), Q(power * (5 * power - 2), 5760)]


def degree_two_row(genus):
    """Independent closed expansion of 2 s(t)^p s(2t)^p - s(t)^(2p)/2."""
    p = 2 * genus - 2
    s = sine_power_row(p)
    double_s = [2**p * 4**r * s[r] for r in range(3)]
    product = [sum(s[j] * double_s[r - j] for j in range(r + 1)) for r in range(3)]
    square = sine_power_row(2 * p)
    return [2 * product[r] - square[r] / 2 for r in range(3)]


def main():
    for genus in range(2, 11):
        p = 2 * genus - 2
        row = degree_two_row(genus)
        assert row == connected_rows(genus, 2, 2)[2]
        etale_unsigned = Q(2 ** (2 * genus) - 1, 2)
        core_row = sine_power_row(2 * p)
        next_core_row = sine_power_row(2 * p + 2)
        first_exponent_magnitude = (genus - 1) * 2 ** (2 * genus - 3)
        inherited = (
            row[2]
            - etale_unsigned * core_row[2]
            + first_exponent_magnitude * next_core_row[1]
        )
        expected = Q(2) ** (2 * genus - 6) * (genus - 1) * (genus - 2)
        assert inherited == expected
        for epsilon in (-1, 1):
            n_min = Q(-1, 2) - epsilon * etale_unsigned
            n_first = epsilon * first_exponent_magnitude
            inherited_direct = (
                epsilon * row[2]
                + n_min * core_row[2]
                + n_first * next_core_row[1]
            )
            etale = -core_row[2] / 2
            assert inherited_direct - etale == epsilon * expected
        print(f"g={genus}: n_(2,{2 * genus + 1}) = main Euler integral - epsilon*{expected}")

    assert degree_two_row(2) == [Q(15, 2), Q(-13, 4), Q(19, 32)]
    etale_mixed = -sine_power_row(4)[2] / 2
    etale_equal = Q(15, 2) * sine_power_row(4)[2]
    ghost_equal = -Q(1, 6) * Q(-2)
    main_equal = degree_two_row(2)[2] - etale_equal - ghost_equal
    assert (etale_mixed, etale_equal, ghost_equal, main_equal) == (
        Q(-1, 160), Q(3, 32), Q(1, 3), Q(1, 6)
    )
    print("g=2 mixed: L_(2,5) = main Euler integral - 1/160; n_(2,5) = main Euler integral")
    print("g=2 equal: etale=3/32, elliptic-tail component=1/3, main=1/6")

    # B=P(V) over Pic^2(C), with C of genus two (the base is an abelian surface).
    # ch(V)=3-4P; p_*zeta^2=1, p_*zeta^3=4P, p_*zeta^4=8P^2.
    # c2(W tensor T^-1)=P^2/2-P*zeta/2+zeta^2/4, T^2=O(1).
    # Terms with P^3 vanish.  Integral on the root gerbe includes factor 1/2.
    c2 = {(2, 0): Q(1, 2), (1, 1): Q(-1, 2), (0, 2): Q(1, 4)}
    square = {}
    for (p1, z1), a in c2.items():
        for (p2, z2), b in c2.items():
            if p1 + p2 <= 2:
                key = (p1 + p2, z1 + z2)
                square[key] = square.get(key, Q(0)) + a * b
    push = {2: Q(1), 3: Q(4), 4: Q(8)}
    integral = sum(a * push[z] * 2 for (p, z), a in square.items()) / 2
    assert integral == 0
    assert integral != main_equal
    print("Compact projective-bundle extension integral=0; differs from equal-spin main integral=1/6.")
    print("No mixed-spin main Euler integral is evaluated by this check.")


if __name__ == "__main__":
    main()
