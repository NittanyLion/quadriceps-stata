<p align="center"><img src="logo.svg" alt="Quadriceps logo" width="200"></p>

# quadriceps (Stata)

![authored by: JP](authored_by.svg)

> **Paper:** Joris Pinkse, *Positive weight Hermite and Legendre quadrature rules* (2026) — **[arXiv:2609.26840](https://arxiv.org/abs/2609.26840)**; Zenodo, DOI: **[10.5281/zenodo.22904159](https://doi.org/10.5281/zenodo.22904159)**
>
> **Data deposit:** Zenodo — DOI: **[10.5281/zenodo.22881864](https://doi.org/10.5281/zenodo.22881864)**

Positive-weight cubature rules in several dimensions, for two weights:

| command | weight (default) | one-dimensional cousin |
|---|---|---|
| `ghpos d q` | standard normal density `N(0, I_d)` on `ℝᵈ` | the `q`-node Gauss–Hermite rule |
| `lepos d q` | uniform density on `[0,1]ᵈ` | the `q`-node Gauss–Legendre rule |

A rule of degree `p` is a set of `n` nodes `x_i` in `ℝᵈ` and weights `w_i > 0` with
`sum_i w_i f(x_i) = integral of f(x) ω(x) dx` for every polynomial `f` of total degree `<= p`.
The product of `q`-node one-dimensional Gauss rules does this for `p = 2q - 1` with `qᵈ` nodes.
The rules stored here, the smallest positive-weight rules known to the author, do it with far
fewer; at `d = 5` the saving is more than a factor of ten. They cover `2 <= d <= 5`.
[`RULES.md`](RULES.md) lists every rule with its node count, Möller's lower bound, measured
accuracy and origin.

This is the Stata twin of [Quadriceps.jl](https://github.com/NittanyLion/Quadriceps.jl) (Julia),
[quadriceps-py](https://github.com/NittanyLion/quadriceps-py) (Python) and
[quadriceps-r](https://github.com/NittanyLion/quadriceps-r) (R). The four packages share their
data, their conventions and their function names; the data are refreshed from the Julia package
whenever a smaller rule is found.

## Installation

Stata 16 or later. From Stata:

```stata
net install quadriceps, from(https://raw.githubusercontent.com/NittanyLion/quadriceps-stata/main) replace
```

The package needs nothing beyond Stata and Mata. It installs three commands (`quadriceps`,
`ghpos`, `lepos`), their help files, the Mata source `quadriceps.mata` (compiled on first use),
and its two data files; `help quadriceps` has the details.

## Use

```stata
ghpos 3 4                         // d = 3, q = 4 (degree 7): 27 nodes instead of 64, in frame quadriceps
frame quadriceps: generate double f = w * x1^2 * x2^4
frame quadriceps: summarize f, meanonly
display r(sum)                    // E[Z₁² Z₂⁴] = 3

lepos 2 5, frame(square) replace  // q = 5 (degree 9): 17 nodes on the unit square instead of 25
ghpos 3, p(7)                     // the first rule again, requested by its degree
ghpos 3 4, matrix(R)              // also as the 27 x 4 Stata matrix R (columns x1 x2 x3 w)
ghpos 3 4, mata(X w)              // also as Mata matrices X (27 x 3) and w (27 x 1)
```

The rule goes into a frame, by default the frame `quadriceps` (overwritten on every call), as
variables `x1` … `xd` and `w`, one node per observation. The current data are not touched. The
second argument `q` has its one-dimensional meaning: the number of nodes of the one-dimensional
Gauss rule, which is exact to degree `2q - 1`. `ghpos d q` returns a `d`-dimensional rule of
that same degree `p = 2q - 1`: a replacement for the `qᵈ`-node product grid, and for `d = 1`
the `q`-node Gauss rule itself.

To ask for a degree instead, use the option `p()`: `ghpos d, p(7)`, `lepos d, p(12)`. Any
`p >= 0` is accepted. Rules are stored at odd degrees and a request is served by the smallest
stored rule of degree `>= p`, so an even `p` gets the rule for `p + 1`. Give `q` or `p()`, not
both.

### `nonormalize`

Gauss–Hermite rules integrate against `exp(-x²)`. The rules here are made for the standard
normal density, which is what an expectation needs, so that is the default:

| | default | `nonormalize` (the classical convention) |
|---|---|---|
| `ghpos` | weight `(2π)⁻ᵈᐟ² exp(-‖x‖²/2)`; weights sum to 1 | weight `exp(-‖x‖²)`; weights sum to `πᵈᐟ²` |
| `lepos` | uniform density on `[0,1]ᵈ`; weights sum to 1 | plain integral over `[-1,1]ᵈ`; weights sum to `2ᵈ` |

For `Y ~ N(mu, L L')` use the nodes `mu + L x_i` with the same weights; for a box, rescale the
columns of the Le nodes.

### `pragmatic`

Rules are stored for `2 <= d <= 5`, up to a degree that depends on the family and on `d` (see
[`RULES.md`](RULES.md) or `quadriceps rules`). For any other request:

* without `pragmatic` (the default) the command stops with return code 499 and a message that
  says how far the stored rules go;
* with `pragmatic` it returns the cheapest tensor product of lower-dimensional rules: the split
  of `d` into stored rules and one-dimensional Gauss rules that needs the fewest nodes. The
  result is a valid positive-weight rule of the requested degree. It is not small, but it is
  much smaller than the plain product grid whenever a stored rule can be a factor.

```stata
ghpos 7 5                               // error: no stored rule in seven dimensions
ghpos 7 5, pragmatic                    // (d = 2) x (d = 5): a few thousand nodes; the grid has 78125

quadriceps nnodes gh 10 3, pragmatic    // the node count, without building the rule
quadriceps ruleinfo gh 7 5, pragmatic   // the factors, with their origins
```

With `pragmatic` a request that a stored rule covers returns that stored rule, as without it.

### Other subcommands

* `quadriceps rules [gh|le]` lists the stored rules.
* `quadriceps nnodes family d [q], [p() pragmatic]` gives a node count without building the rule.
* `quadriceps ruleinfo family d [q], [p() pragmatic]` describes the rule and its origin.
* `quadriceps check family d [q], [p() pragmatic]` builds the rule and measures how exact it is.

### Mata

The same rules are available in Mata once the package has been loaded (any call of
`quadriceps`, for instance `quadriceps version`, does that):

```stata
mata: ghpos(3, 4, X=., w=.)                     // fills X (27 x 3) and w (27 x 1)
mata: w' * (X[,1]:^2 :* X[,2]:^4)               // 3
mata: quadriceps_nnodes("le", 6, 3, 1)          // "le", d, p, pragmatic
mata: quadriceps_exactness(X, w, 7, "gh")       // largest relative monomial error
```

`ghpos(d, q, X, w [, normalize, pragmatic])`, `lepos(...)`,
`quadriceps_rule(family, d, p, pragmatic, normalize, X, w)`, `quadriceps_nnodes(...)`,
`quadriceps_ruleinfo(...)` and `quadriceps_exactness(...)`; see `help quadriceps`.

All rules are stored in one binary file, `quadriceps_rules.bin`, with `quadriceps_index.tsv`
next to it as the catalog; [`FORMAT.md`](FORMAT.md) specifies the format, which the Julia,
Python and R packages share byte for byte. Mata reads the file with its buffered binary I/O; no
plugin and no external program is involved.

## Accuracy

Rules are stored in double precision, each the rounding of an 80-digit rule, so what remains of
its error is rounding: the largest relative monomial error over all monomials of degree `<= p`
is below `5.2e-15` for every GH rule (most are `1e-16` to `1e-15`; the worst is `d = 5`,
`p = 21`) and below `4.5e-16` for every Le rule. All weights are positive, and all nodes of the
Le rules lie strictly inside the cube. The gate applied when the data were built, and again by
the tests on every rule, is `1e-11`. The measured value of each rule is in the catalog (`relerr`)
and in [`RULES.md`](RULES.md). One-dimensional Gauss rules are computed by the Golub–Welsch
eigenvalue method in Mata.

## Whose rules these are

<!-- BEGIN GENERATED credits -->
116 of the 142 rules were computed from scratch by the author. A further 11 (Legendre) rules were
obtained by node elimination started from Diallo and Worku's published rules. Finally, 15 are
rules from the literature (copied in, or found again by the author's search and recognized).
<!-- END GENERATED credits -->
`quadriceps ruleinfo` and the `origin` column of [`RULES.md`](RULES.md) say which is which; cite
the source named there when you use such a rule. [`NOTICE.md`](NOTICE.md) has the details and
the license notice that travels with the derived files.

## Tests

```stata
net get quadriceps, from(https://raw.githubusercontent.com/NittanyLion/quadriceps-stata/main) replace
do quadriceps_test.do
```

or, from a checkout, `adopath + "path/to/quadriceps-stata"` and then `do quadriceps_test.do` in
that directory. The script checks every stored rule (size, positive weights, exactness up to its
degree), both conventions, the fallback and the error paths, writes `quadriceps_test.log`, and
ends with a line `N checks, M failed`. It takes a minute or two. (`replace` matters: without it,
`net get` leaves an older copy of the test file in place, silently, with `r(602)`.)

The author does not have Stata. The package was written against the Stata and Mata manuals and
is tested by colleagues who do; if something fails for you, please open an issue with the log.

## License

MIT; see [`LICENSE`](LICENSE). The rules that descend from, or coincide with, published rules
carry their sources' notices in [`NOTICE.md`](NOTICE.md).
