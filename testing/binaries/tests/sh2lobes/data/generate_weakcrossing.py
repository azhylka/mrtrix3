#!/usr/bin/env python3
"""Generate weakcrossing_lmax8.mif.

A two-lobe phantom whose second lobe is deliberately faint: it sits above zero, and so is
segmented when the FMLS peak-amplitude threshold is disabled, but falls below the default
threshold of 0.1 and is discarded when the default applies. It therefore distinguishes
`sh2lobes` run with `-fmls_no_thresholds` (16 fixels) from `sh2lobes` run with its defaults
(8 fixels), which the symmetric phantoms cannot do.

The profile is

    S_z + w * S_x        with   S(theta) = cos^8(theta) - cos^8(42 deg)   and   w = 0.16

using the same `alpha = 42 deg`, `lmax = 8` cone as `crossing_lmax8_alpha42.mif`, which is
itself `S_z + S_x`. Rather than re-derive the rotation of S onto the x axis, this script
recovers it from that file as `S_x = crossing - S_z`: only the zonal part `S_z` then has to be
computed here, and the result is guaranteed consistent with the phantom the other tests use.
The recovered `S_x` is checked against `S_z` for equal per-degree norm, which a rotation must
preserve; see README.md for the construction and the coefficient conventions.

Depends only on the standard library. Run from this directory.
"""

import math
import re
import struct
from fractions import Fraction

ALPHA_DEG = 42.0
LMAX = 8
WEIGHT = 0.16
SOURCE = "crossing_lmax8_alpha42.mif"
TARGET = "weakcrossing_lmax8.mif"


def legendre_polynomials(lmax):
    """P_0..P_lmax as exact monomial coefficient lists, lowest power first."""
    polys = [[Fraction(1)], [Fraction(0), Fraction(1)]]
    for n in range(1, lmax):
        # (n+1) P_{n+1} = (2n+1) x P_n - n P_{n-1}
        shifted = [Fraction(0)] + polys[n]
        prev = polys[n - 1] + [Fraction(0)] * (len(shifted) - len(polys[n - 1]))
        polys.append([(Fraction(2 * n + 1) * a - Fraction(n) * b) / (n + 1)
                      for a, b in zip(shifted, prev)])
    return polys[:lmax + 1]


def zonal_sh_coefficients():
    """SH coefficients of  cos^8(theta) - cos^8(alpha)  about +z, in the MRtrix basis."""
    polys = legendre_polynomials(LMAX)

    # a_l = (2l+1)/2 * integral of x^8 P_l(x) over [-1,1]; exact, since the integral of x^k
    #   over [-1,1] is 2/(k+1) for even k and 0 for odd k.
    legendre_coefs = []
    for l, poly in enumerate(polys):
        integral = sum(coef * Fraction(2, k + 8 + 1)
                       for k, coef in enumerate(poly)
                       if (k + 8) % 2 == 0)
        legendre_coefs.append(Fraction(2 * l + 1, 2) * integral)

    # Subtracting the constant cos^8(alpha) alters the l = 0 term alone.
    legendre_coefs[0] -= Fraction(math.cos(math.radians(ALPHA_DEG)) ** 8).limit_denominator(10 ** 12)

    # c_l = a_l * sqrt(4 pi / (2l + 1)), stored at index l(l+1)/2 (the m = 0 term).
    coefs = [0.0] * ((LMAX + 1) * (LMAX + 2) // 2)
    for l in range(0, LMAX + 1, 2):
        coefs[l * (l + 1) // 2] = float(legendre_coefs[l]) * math.sqrt(4.0 * math.pi / (2 * l + 1))
    return coefs


def read_mif(path):
    raw = open(path, "rb").read()
    header = raw[:raw.index(b"END") + 3].decode("latin-1")
    offset = int(re.search(r"file: \. (\d+)", header).group(1))
    dim = [int(x) for x in re.search(r"dim: (\S+)", header).group(1).split(",")]
    count = dim[0] * dim[1] * dim[2] * dim[3]
    return dim, list(struct.unpack("<%df" % count, raw[offset:offset + 4 * count]))


def write_mif(path, dim, values, comment):
    header = ("mrtrix image\n"
              "dim: %s\n"
              "vox: 1,1,1,1\n"
              "layout: +0,+1,+2,+3\n"
              "datatype: Float32LE\n"
              "transform: 1,0,0,0\n"
              "transform: 0,1,0,0\n"
              "transform: 0,0,1,0\n"
              "comment: %s\n"
              "file: . 1024\n"
              "END\n" % (",".join(str(d) for d in dim), comment)).encode("latin-1")
    assert len(header) <= 1024, "header does not fit within the declared 1024-byte offset"
    with open(path, "wb") as f:
        f.write(header)
        f.write(b"\n" * (1024 - len(header)))
        f.write(struct.pack("<%df" % len(values), *values))


def main():
    dim, crossing = read_mif(SOURCE)
    nvox = dim[0] * dim[1] * dim[2]
    ncoef = dim[3]

    s_z = zonal_sh_coefficients()
    assert len(s_z) == ncoef, "coefficient count %d does not match %s (%d)" % (len(s_z), SOURCE, ncoef)

    # Layout +0,+1,+2,+3: axis 0 is fastest and the volume slowest, so coefficient v of voxel i
    #   lives at flat index v * nvox + i. Every voxel of these phantoms is identical; take
    #   voxel 0 and confirm the rest agree.
    crossing_coefs = [crossing[v * nvox] for v in range(ncoef)]
    for v in range(ncoef):
        for i in range(nvox):
            assert abs(crossing[v * nvox + i] - crossing_coefs[v]) < 1e-6, \
                "%s is not uniform across voxels" % SOURCE

    s_x = [c - z for c, z in zip(crossing_coefs, s_z)]

    # A rotation preserves the norm within each degree l; if S_z is right, S_x must match it.
    #   (l, m) sits at index l(l+1)/2 + m with m running from -l to +l, so degree l occupies
    #   the 2l+1 entries starting at l(l-1)/2.
    for l in range(0, LMAX + 1, 2):
        base = l * (l - 1) // 2
        norm_z = math.sqrt(sum(s_z[base + m] ** 2 for m in range(2 * l + 1)))
        norm_x = math.sqrt(sum(s_x[base + m] ** 2 for m in range(2 * l + 1)))
        assert abs(norm_z - norm_x) < 1e-5 * max(1.0, norm_z), \
            "l=%d norm mismatch: S_z %.9f vs recovered S_x %.9f " \
            "(S_z is wrong, or %s is not S_z + S_x)" % (l, norm_z, norm_x, SOURCE)

    weak_coefs = [z + WEIGHT * x for z, x in zip(s_z, s_x)]

    peak = 1.0 - math.cos(math.radians(ALPHA_DEG)) ** 8
    trough = -math.cos(math.radians(ALPHA_DEG)) ** 8
    print("strong lobe peak (+z): %.4f" % (peak + WEIGHT * trough))
    print("weak   lobe peak (+x): %.4f   (default FMLS peak threshold is 0.1)"
          % (trough + WEIGHT * peak))

    values = [weak_coefs[v] for v in range(ncoef) for _ in range(nvox)]
    write_mif(TARGET, dim, values,
              "cos^8(theta) - cos^8(42 deg) along +z plus %g of the same along +x; "
              "the weak +x lobe peaks below the default FMLS peak-amplitude threshold" % WEIGHT)
    print("wrote %s" % TARGET)


if __name__ == "__main__":
    main()
