# Erdős–Straus Conjecture ($n \equiv 1 \pmod{24}$): Sieve Dimension, Search Depth, and Composite-Witness Sieves

[![Verification CI](https://github.com/Benjamin-Dahan/Conjecture-Erdos-Straus/actions/workflows/ci.yml/badge.svg)](https://github.com/Benjamin-Dahan/Conjecture-Erdos-Straus/actions/workflows/ci.yml)
![Language](https://img.shields.io/badge/C-C11%20%2B%20128--bit%20Safe-00599C?logo=c&logoColor=white)
![Python](https://img.shields.io/badge/Python-3.12%20Exact%20Arithmetic-3776AB?logo=python&logoColor=white)
[![arXiv 2608.24035v2](https://img.shields.io/badge/arXiv-2608.24035v2-b31b1b.svg)](https://arxiv.org/abs/2608.24035)
[![arXiv 2611.08717v1](https://img.shields.io/badge/arXiv-2611.08717v1-b31b1b.svg)](https://arxiv.org/abs/2611.08717)

Research papers and exact integer verification code studying the Erdős–Straus equation

$$\frac{4}{n} = \frac{1}{x} + \frac{1}{y} + \frac{1}{z} \qquad (x, y, z \in \mathbb{N}^*)$$

for primes $n \equiv 1 \pmod{24}$—the residual congruence class not covered by elementary modular filters—through explicit, **factorization-free** divisor witnesses and analytic sieve theory.

---

## Papers Included in [`papers/`](papers/)

1. **Sieve dimension and search depth for the Erdős–Straus conjecture, $n \equiv 1 \pmod{24}$** — *Benjamin Dahan* (submitted / under review at *The Ramanujan Journal*)
   - [`papers/2608.24035v2_sieve_dimension_depth.pdf`](papers/2608.24035v2_sieve_dimension_depth.pdf) · [`LaTeX source`](papers/2608.24035v2_sieve_dimension_depth.tex)
   - **Key results:**
     - **Two-parameter Type II & Type I normal forms:** $4/n$ admits a Type II decomposition iff some $an+u$ ($\gcd(u,a)=1$) has a divisor $s \equiv -1 \pmod{4ua}$ (`Theorem 3.9`), and a Type I decomposition iff some $an+u$ has a divisor $s \equiv -n \pmod{4ua}$ (`Theorem 3.14`), enabling simultaneous factorization-free search in $\mathcal{O}(J^{3/2}\sqrt{n} + J^2 \log n)$ arithmetic operations.
     - **Characterization of the hyperbolic slice ($a=1$):** Proves the hyperbolic family $n = 4\alpha u v - u - v$ is strictly proper (`Theorem 3.6` / `Corollary 3.7`), identifying the exact **14 exceptional primes** $\mathcal{H}(5\times 10^9)$ outside $a=1$ (`Proposition 3.8`, extending Xu's 9 wild primes below $7.2\times 10^5$) and resolving all 14 with parameters $u,a \le 5$ (`Proposition 3.13`).
     - **Sieve dimension & large sieve constants:** Unconditional fixed-depth dimension $|E_{\mathcal{P}}(N)| \ll_{\mathcal{P}} N(\log N)^{-1-|\mathcal{P}|/2}$, exact two-variable Dirichlet evaluation yielding $\kappa_1^* = \frac{9}{2\pi^2} = 0.455945\dots$ ($3\times$ improvement in exponent), sub-quadratic $o(N^2)$ batch interval algorithm, and proof that single-prime witnesses hit a universal $\tau_3$ ceiling at $(\log N)^{2/3}$.

2. **Composite witness sieves, Kneser–Plancherel containers, and the exponent $5/7$ for the Erdős–Straus conjecture** — *Benjamin Dahan* (preprint)
   - [`papers/2611.08717v1_composite_witness_sieves.pdf`](papers/2611.08717v1_composite_witness_sieves.pdf) · [`LaTeX source`](papers/2611.08717v1_composite_witness_sieves.tex)
   - **Key results:**
     - Breaks Vaughan's 1970 barrier $|E(N)| \ll N\exp(-c(\log N)^{2/3})$ unconditionally on the full exceptional set via semi-composite witnesses $s = pq \equiv -1 \pmod{12ua}$:
       $$|E(N)| \ll N \exp\!\left(-c' \frac{(\log N)^{5/7}}{(\log\log N)^{2/7}}\right)$$
     - Reaches exponent $1 - o(1)$ on the structured small-sheet component $E_2^{\mathrm{small}}(N; J(N))$ via Kneser multiplicative stabilizer kernels, quotient containers, Pontryagin–Plancherel character containers, and a Krawtchouk–Montgomery generating-polynomial large sieve.

---

## Verified C11 & Python Tooling

This repository replaces an earlier 2025 computational prototype with overflow-safe C11 (`uint64_t` + `__uint128_t` identity verification) and an exact Python `fractions.Fraction` reproducibility suite aligned directly with **Paper 1, Section 7**.

### Project Structure

```
.
├── papers/
│   ├── 2608.24035v2_sieve_dimension_depth.{tex,pdf}
│   └── 2611.08717v1_composite_witness_sieves.{tex,pdf}
├── src/
│   ├── c/
│   │   └── erdos_straus_sieve.c   # C11 Algorithm 1, Segmented Type I/II Sieve, & H(N) Finder
│   └── python/
│       └── verify_paper1.py       # Section 7 exact integer verification suite
├── tests/
│   └── test_c_cli.py              # CLI & edge-case tests
├── Makefile
└── .github/workflows/ci.yml
```

### Build & Run Verification

```bash
make CFLAGS="-Wall -Wextra -Werror -std=c11 -O3"
make test
```

### Measured Local Benchmarks (`gcc -O3`, `x86_64`)

- **Exhaustive hyperbolic exceptional set $\mathcal{H}(10^7)$ (`--find-H 10000000`):**
  Completed in **`0.31 s`**, checking all **82,887** primes $n \le 10^7$ with $n \equiv 1 \pmod{24}$ and reproducing the exact **14 exceptional primes**:
  `409, 577, 5569, 9601, 23929, 83449, 102001, 329617, 712321, 1134241, 1724209, 1726201, 5212561, 8813281`.
- **Full verification up to $10^6$ (`--verify-range 1000000 --depth 64`):**
  Completed in **`9.79 s`**, verifying all **9,732** primes $n \le 10^6$ ($n \equiv 1 \pmod{24}$): **9,732 / 9,732** resolved by both the two-parameter sieve ($J=64$) and Algorithm 1 (maximum shift $c = 63$).
- **Exact Python reproduction (`src/python/verify_paper1.py`):**
  All Section 7 checks pass in exact rational/integer arithmetic, including Example 3.4 ($n=20353$), first 600 primes shift counts, $n=100129$ rejection iteration counts (`1670`, `6439`, `8697`), all **4,661** two-parameter witness reconstructions, $\kappa_0$ / $\kappa_1^*$ tables, and $c=7$ blind-set counts.

### CLI Usage Examples

```bash
# Compare Algorithm 1, Hyperbolic Sieve (a=1), and Simultaneous Type I/II Sieve on a prime:
./build/erdos_straus_sieve --check 409 --depth 5

# Force Type II only in the two-parameter sieve (exhibits modulator a=2 resolution of 409):
./build/erdos_straus_sieve --check 409 --depth 5 --type2-only

# Compute the exact exceptional set H(10^7) outside the a=1 hyperbolic slice:
./build/erdos_straus_sieve --find-H 10000000

# Verify every prime n <= 1,000,000 (n ≡ 1 mod 24):
./build/erdos_straus_sieve --verify-range 1000000 --depth 64
```

---

## Author

**Benjamin Dahan** — Université Paris Dauphine – PSL (`benjamin.dahan@dauphine.eu`)
