# AGENTS.md

Guidance for coding agents (and people) working in this repository.

## Paper and deposit

The rules are described in Joris Pinkse, *Positive weight Hermite and Legendre quadrature rules* (2026),
on arXiv as **[arXiv:2609.26840](https://arxiv.org/abs/2609.26840)** and deposited on Zenodo:

* paper: **[10.5281/zenodo.22904159](https://doi.org/10.5281/zenodo.22904159)** (concept DOI: it always resolves to the latest version)
* rules, to 80 digits: **[10.5281/zenodo.22881864](https://doi.org/10.5281/zenodo.22881864)** (one record for both weight families)

Cite the paper by its arXiv identifier and its Zenodo DOI together. Keep the block at the top of
`README.md` in step with this one, in all four packages (Quadriceps.jl, quadriceps-py,
quadriceps-r, quadriceps-stata).

## What this is

The Stata twin of Quadriceps.jl: positive-weight cubature rules for the Gaussian weight (`ghpos`)
and the uniform weight on the cube (`lepos`). The Julia package is the master copy. Its
`build/update.sh` refreshes `quadriceps_rules.bin`, `quadriceps_index.tsv`, `RULES.md`,
`FORMAT.md` and `NOTICE.md` here from the Julia package, commits and pushes. The `GENERATED
credits` block of `README.md` is filled in by the same script. **Never edit those by hand** (and
never add per-rule CSV or text files: the rules live in the one binary file), and make behavior
changes in all four packages (Julia, Python, R, Stata) together.

## Written blind

**The author has no Stata license.** Everything here was written against the Stata and Mata
manuals without ever being run, and is tested by colleagues who have Stata: they run
`do quadriceps_test.do` and send back `quadriceps_test.log`. Consequences:

* Do not claim that anything works. A change is "written", not "tested", until a log says so.
* Keep the untested surface small: no feature that the other three packages do not have, no
  clever Stata idioms where a plain one will do, no dependence on community-contributed packages.
* Prefer the documented form of every Mata function over a remembered one; the manuals are at
  `https://www.stata.com/manuals/<chapter>.pdf` (e.g. `m-5bufio.pdf`, `rnet.pdf`).
* When a tester reports a failure, fix the cause, not the check, and ask for a rerun.

Lessons from the testers' logs so far (Stata 19.5 SE, Mac, 2026-09-24):

* A column vector and a row vector are not c-conformable (`(0::r) :+ (0..r)` is an error); a
  column against a matrix with the same number of rows is fine. See [M-2] op_colon.
* Mata functions defined in an ado-file's `mata:` block are private to that ado-file.
* In a one-line `else command`, a `` `=exp' `` macro on the else line is expanded even when the
  branch is not taken.

## Layout

| path | what it holds |
|---|---|
| `quadriceps.ado` | the Stata command `quadriceps` and its subprograms, and `_quadriceps_load`, which compiles `quadriceps.mata` on first use (and again when `quadriceps_mata_version()` disagrees with the ado's version) |
| `quadriceps.mata` | the Mata code: data reader, catalog, planner, Golub–Welsch, tensor product, exactness measure, the glue to frames and matrices. It is a separate file loaded with `run` because Mata functions defined inside an ado-file are private to it ([M-1] Ado); this way `ghpos()` and friends are usable from Mata directly |
| `ghpos.ado`, `lepos.ado` | one-line aliases for `quadriceps gh` and `quadriceps le` |
| `quadriceps.sthlp` | the help file (SMCL); `ghpos.sthlp` and `lepos.sthlp` include it |
| `quadriceps_rules.bin`, `quadriceps_index.tsv` | the data (format `QUADRICEPS1`, see `FORMAT.md`) — generated; byte-identical to `data/rules.bin` and `data/index.tsv` of Quadriceps.jl, renamed because `net install` puts every file into a flat letter directory of the ado-path |
| `quadriceps_test.do` | the tests: every stored rule, both conventions, the fallback, the error paths; ancillary file, fetched with `net get` |
| `stata.toc`, `quadriceps.pkg` | the `net install` manifest; on every release bump `Distribution-Date` in `quadriceps.pkg`, the `*! version` lines, and the version string in both `_quadriceps_load` and `quadriceps_mata_version()` (they must agree, or the loader recompiles on every call) |
| `RULES.md`, `FORMAT.md`, `NOTICE.md` | generated copies from Quadriceps.jl |

## Conventions that must hold

* **`q` and `p`.** `q` is the number of nodes of the one-dimensional Gauss rule; the rule
  returned has degree `p = 2q - 1`, and for `d = 1` it is the `q`-node Gauss rule. The degree
  can be given instead as `p()`. Exactly one of the two. Internally everything works in `p`.
* **Normalized frame inside.** Files, tensor products and `quadriceps_exactness()` use
  `N(0, I_d)` for GH and the uniform density on `[0,1]ᵈ` for Le, weights summing to 1.
  `nonormalize` is applied once, at the end, in `quadriceps_rule()`.
* **Normalized is the default**, unlike the classical Gauss rules. This is deliberate.
* **`pragmatic`.** Without it: return code 499 when no stored rule covers the request. With it:
  the cheapest tensor product of stored rules and Gauss rules; a stored rule wins ties.
* **Positive weights only**, relative monomial error below `1e-11`; the tests enforce both.
* **Output**: a frame with variables `x1..xd`, `w` (default frame `quadriceps`, always
  overwritten; any other frame needs `replace`), optionally a Stata matrix `n x (d+1)` and Mata
  matrices. Never Stata matrices alone: they are capped at `c(max_matsize)` rows, which the
  largest stored rules exceed.
* Stata 16 or later (frames). Mata only; no plugin, no external program. American spelling.
* Return codes: 198 for argument errors, 499 for "no stored rule", 498 for anything else the
  package refuses, 110 for an existing frame.

## Repository

`authored_by.svg` is the author's shield (the same file as in the other packages); the README
shows it after the badges. Do not replace it with a generated shields.io badge. The logo is
`logo.svg`, a copy of the one in Quadriceps.jl.

There is no CI: GitHub Actions has no Stata. Public, `github.com/NittanyLion/quadriceps-stata`,
branch `main`. MIT license (`LICENSE`); `NOTICE.md` carries the notices of the rules that descend
from published ones.

SSC: to submit, email Kit Baum (baum@bc.edu) the files listed in `quadriceps.pkg`; there is no
web form. `net install` from GitHub works without SSC.
