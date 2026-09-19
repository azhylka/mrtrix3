# sh2lobes test phantoms

Synthetic FOD images whose lobe geometry is known in closed form, used by the `sh2lobes/*` tests
to check the reported lobe width against an analytic value rather than against the output of a
prior software version.

Every image is 2x2x2 voxels on an identity transform with unit voxel size, and every voxel holds
the same spherical harmonic coefficients, so each voxel is an independent repeat of the same test.

## Construction

Each single-fibre phantom is the axially symmetric function

    A(theta) = cos^2n(theta) - cos^2n(alpha)

where `theta` is the angle from the fibre axis. Two properties make this useful:

* `cos^2n(theta)` is a polynomial of degree `2n` in `cos(theta)`, so `A` is **exactly**
  representable in an even spherical harmonic basis of `lmax = 2n` -- there is no truncation
  error to account for.
* `A` is positive for `theta < alpha` and negative beyond it. The FMLS segmenter assigns only
  positive directions to a positive lobe, so the segmented lobe is exactly the cone of
  half-angle `alpha`.

The reported width is not that cone, however. It is the extent of the lobe's *core*: twice the
largest angle between the fixel direction and any direction of the lobe at which the amplitude
reaches a quarter of the lobe's maximum. For these phantoms the peak is `A(0) = 1 - k` with
`k = cos^2n(alpha)`, so the quarter-amplitude contour lies at

    cos^2n(theta) - k = 0.25 * (1 - k)   =>   theta = acos( (0.25 + 0.75k)^(1/2n) )

and the expected width is twice that. It is always narrower than `2 * alpha`, which measured out
to the zero crossing on the shallow tail of the profile.

Subtracting the constant only alters the `l = 0` coefficient. It must stay positive
(`cos^2n(alpha) < 1/(2n+1)`), otherwise the segmenter rejects the voxel outright; this is what
bounds how narrow a lobe can be made at a given `lmax`.

Writing `A` as a Legendre series `sum_l a_l P_l(cos theta)`, the MRtrix zonal coefficients are
`c_l = a_l * sqrt(4 pi / (2l + 1))`, stored at index `l(l+1)/2` (0, 3, 10, 21, 36). To place the
lobe on an axis other than `z`, those coefficients are converted to rotational harmonics
(`Math::SH::SH2RH`) and convolved onto a delta function in the desired direction
(`Math::SH::delta`, `Math::SH::sconv`).

## Files

| file | lmax | axis | expected width |
|---|---|---|---|
| `cone_lmax2_alpha75.mif` | 2 | +z | 113.55 deg (113.01 sampled) |
| `cone_lmax4_alpha60.mif` | 4 | +z | 84.85 deg (84.73 sampled) |
| `cone_lmax6_alpha45.mif` | 6 | +z | 66.36 deg (65.68 sampled) |
| `cone_lmax4_alpha60_oblique.mif` | 4 | (1,2,3)/sqrt(14) | 84.85 deg (84.83 sampled) |
| `positive_definite_lmax2.mif` | 2 | +z | 125.10 deg (125.01 sampled) |
| `crossing_lmax8_alpha42.mif` | 8 | +z and +x | 54.68 deg each (54.00 sampled) |
| `weakcrossing_lmax8.mif` | 8 | +z, and +x at 0.16 weight | not asserted; used for lobe counts |

`positive_definite_lmax2.mif` is `cos^2(theta) + 0.05`, which is positive in every direction. It
is therefore a single lobe spanning the whole sphere: a single-fibre voxel of a real FOD that
never goes negative behaves the same way. Under a measure taken out to the zero crossing this was
the degenerate case, saturating at `180 deg` because there is no zero crossing to find. The
quarter-amplitude contour still resolves it at `125.10 deg`, which is the point of this phantom.

`weakcrossing_lmax8.mif` is `S_z + 0.16 * S_x` for the same `alpha = 42 deg` profile `S`. The
weight is chosen so that the `+x` lobe peaks at `0.052`: above zero, and so segmented when the
FMLS peak-amplitude threshold is disabled, but below the default threshold of `0.1`, and so
discarded when the default applies. It is the phantom that distinguishes the two settings --
every other lobe here peaks above `0.85`, far outside the range where the threshold acts. The
faint lobe is much narrower than a `42 deg` cone, since it only clears zero within `17.6 deg` of
the `x` axis. No width is asserted for this phantom; the tests that use it identify the retained
lobe by peak amplitude instead, `1 - 1.16 * cos^8(42 deg) = 0.8921` against `0.0521` for the
faint one, which the segmenter refines by Newton optimisation and so does not quantise.
It is produced by `generate_weakcrossing.py` (standard library only), which derives the rotated
`S_x` from `crossing_lmax8_alpha42.mif` rather than re-deriving the rotation, and checks the
per-degree norms agree as a rotation requires.

`crossing_lmax8_alpha42.mif` is the sum of two `alpha = 42 deg` profiles along `z` and `x`. The
two cones do not touch (`2 x 36.77 < 90`), so each lobe is bounded by its own zero crossing
rather than by a watershed between the two. This matters: bridging directions between two
*touching* lobes are assigned by `retrospective_assignments` to `adj_lobes.front()`, the
lower-indexed lobe, which would make an otherwise symmetric crossing asymmetric by roughly one
direction-set spacing.

Supporting images:

* `oblique_axis_peaks.mif` -- the reference peak directions `(1,2,3)/sqrt(14)`, for comparing
  against `fixel2peaks` output.
* `halfmask.mif` -- a binary mask selecting the `x = 1` half of the grid (4 of the 8 voxels).
  It has to be *this* half rather than `x = 0`: `FMLS::FODQueueWriter::operator()` skips
  masked-out voxels in a `do`/`while` that exits when the loop is exhausted, and then emits the
  voxel it stopped on regardless, so a mask whose last voxel in scan order is excluded yields one
  spurious segmented voxel outside the mask. Selecting `x = 1` puts `(1,1,1)` inside the mask and
  avoids the issue. The quirk is in the shared segmenter input stage, so `fod2fixel` has it too;
  it simply goes unnoticed on real data, where the voxels surrounding a brain mask hold an empty
  FOD and yield no fixels however they are handled. It cannot be sidestepped by thresholding
  here, because every voxel of these phantoms carries the same non-empty FOD: the spurious voxel
  would contribute its full complement of fixels at any threshold.

## "Sampled" values

The sampled column is what the command actually reports. The width is the largest angle between
the fixel direction and any qualifying direction of the lobe, so it is quantised by the
1281-direction tessellation the segmenter uses (mean nearest-neighbour angle 4.09 degrees) and
always slightly underestimates the analytic value: the outermost direction that clears the
quarter-amplitude threshold generally sits a little inside the true contour. Across the phantoms
above the shortfall is at most 0.7 degrees on the full width. The tests allow 2 degrees, and 4
degrees for the crossing phantom, whose boundary is not a circle of constant `theta` and is
therefore sampled less evenly.

The threshold itself is taken relative to the largest amplitude among the sampled directions of
the lobe, not the Newton-refined peak that `-peak_amp` reports. Both sides of the comparison are
then drawn from the same sampling, and the direction attaining the maximum always qualifies, so
every lobe has a well-defined width.

## Format

Written directly as MRtrix `.mif`: a text header, padding to a 1024-byte offset declared by the
`file: . 1024` entry, then the raw little-endian `float32` payload in `+0,+1,+2,+3` layout
(axis 0 fastest, volume slowest). They are marked `binary` in `.gitattributes`, since
line-ending conversion would corrupt the payload.
