#include <ctype.h>
#include <inttypes.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*
 * Factorization-free Erdős–Straus (4/n = 1/x + 1/y + 1/z) solvers and
 * hyperbolic exceptional-set sieve accompanying:
 *   [1] B. Dahan, "Sieve dimension and search depth for the Erdős–Straus conjecture,
 *       n ≡ 1 (mod 24)", arXiv:2608.24035v2 (2026).
 *   [2] B. Dahan, "Composite witness sieves, Kneser–Plancherel containers, and the
 *       exponent 5/7 for the Erdős–Straus conjecture", arXiv:2611.08717v1 (2026).
 */

typedef struct {
    bool found;
    uint64_t x;
    uint64_t y;
    uint64_t z;
    uint64_t shift_c;
    uint64_t mod_a;
    uint64_t u;
    uint64_t v;
    uint64_t witness_s;
    uint64_t iterations;
    char branch[16]; /* "TypeII", "TypeI", or "Alg1" */
} EsSolution;

static uint64_t gcd_u64(uint64_t a, uint64_t b) {
    while (b != 0) {
        uint64_t t = a % b;
        a = b;
        b = t;
    }
    return a;
}

static uint64_t isqrt_u64(uint64_t n) {
    if (n < 2) {
        return n;
    }
    uint64_t x = n;
    uint64_t y = (x + 1) >> 1;
    while (y < x) {
        x = y;
        y = (x + n / x) >> 1;
    }
    return x;
}

/* Overflow-safe exact verification: 4/n == 1/x + 1/y + 1/z. */
static bool verify_identity(uint64_t n, uint64_t x, uint64_t y, uint64_t z) {
    if (n == 0 || x == 0 || y == 0 || z == 0) {
        return false;
    }
    /* Scale numerator and denominator by n*x*y*z/g? Better: __uint128_t if products fit,
       else compare 4/n - 1/x = (4x-n)/(nx) == (y+z)/(yz) after reducing each side! */
    if (4ULL * x <= n) {
        return false;
    }
    uint64_t lhs_num = 4ULL * x - n;
    uint64_t lhs_den = n * x;
    uint64_t g1 = gcd_u64(lhs_num, lhs_den);
    lhs_num /= g1;
    lhs_den /= g1;

    /* Reduce y, z separately so (y+z)/g_yz and (y/g_yz)*z stay within 128-bit range. */
    uint64_t g_yz = gcd_u64(y, z);
    uint64_t y_red = y / g_yz;
    uint64_t z_red = z / g_yz;
    uint64_t rhs_num = y_red + z_red;
    uint64_t g2 = gcd_u64(rhs_num, y_red);
    rhs_num /= g2;
    uint64_t y_part = y_red / g2;
    uint64_t g3 = gcd_u64(rhs_num, z_red);
    rhs_num /= g3;
    uint64_t z_part = z_red / g3;
    uint64_t g4 = gcd_u64(rhs_num, g_yz);
    rhs_num /= g4;
    uint64_t g_part = g_yz / g4;

#if defined(__SIZEOF_INT128__)
    unsigned __int128 den_r = (unsigned __int128)y_part * z_part * g_part;
    return (rhs_num == lhs_num) && (den_r == (unsigned __int128)lhs_den);
#else
    if (rhs_num != lhs_num || lhs_den % y_part != 0) {
        return false;
    }
    uint64_t rem = lhs_den / y_part;
    if (rem % z_part != 0) {
        return false;
    }
    return (rem / z_part) == g_part;
#endif
}

static bool is_prime_u64(uint64_t n) {
    if (n < 2) {
        return false;
    }
    if (n % 2 == 0) {
        return n == 2;
    }
    if (n % 3 == 0) {
        return n == 3;
    }
    for (uint64_t d = 5; d <= n / d; d += 6) {
        if (n % d == 0 || n % (d + 2) == 0) {
            return false;
        }
    }
    return true;
}

/*
 * Algorithm 1 (§3.8 of arXiv:2608.24035v2):
 * Divisor-based factorization-free search for primes n ≡ 1 (mod 4) [specifically mod 24].
 * Scans shifts c = 3, 7, 11, ... < n. Setting exact_mode=true disables the
 * modulator heuristic a > (n+3)/10 and keeps only the unconditional V > K^2 bound.
 */
static EsSolution solve_alg1(uint64_t n, bool exact_mode, uint64_t single_c) {
    EsSolution sol = {0};
    snprintf(sol.branch, sizeof(sol.branch), "%s", "Alg1");

    uint64_t c_start = (single_c > 0) ? single_c : 3;
    uint64_t c_end = (single_c > 0) ? single_c : (n - 1);

    for (uint64_t c = c_start; c <= c_end; c += 4) {
        uint64_t q = n / c;
        uint64_t r = n % c;
        int64_t b_base = (int64_t)(q % 4);
        uint64_t K = (n + c) / 4;

        for (uint64_t delta = 0;; delta++) {
            sol.iterations++;
            int64_t b = b_base - 4 * (int64_t)delta;
            int64_t denom = 3 - b;
            if (denom <= 0) {
                continue;
            }

            /* V = (c * denom - r) / 4 >= 0 since c * denom ≡ r (mod 4) */
            uint64_t c_denom = c * (uint64_t)denom;
            if (c_denom <= r) {
                continue;
            }
            uint64_t V = (c_denom - r) / 4;
            if (V == 0) {
                continue;
            }

            uint64_t g = gcd_u64(K, V);
            uint64_t d2 = V / g;
            if (K % d2 == 0) {
                uint64_t d1 = K / g;
                uint64_t sigma = (q + 4 - (uint64_t)b) / 4;
                uint64_t a = sigma / g;

                sol.found = true;
                sol.shift_c = c;
                sol.mod_a = a;
                sol.u = d1;
                sol.v = d2;
                sol.x = K;
                sol.y = n * sigma;
                sol.z = n * a * (K / d2);
                return sol;
            }

#if defined(__SIZEOF_INT128__)
            bool v_exceeds_k2 = ((unsigned __int128)V > (unsigned __int128)K * K);
#else
            bool v_exceeds_k2 = (V / K > K) || (V / K == K && V % K > 0);
#endif
            if (v_exceeds_k2) {
                break;
            }
            if (!exact_mode) {
                /* Guard: (q + 4 - b) > 4 * g * (n + 3) / 10 */
#if defined(__SIZEOF_INT128__)
                unsigned __int128 lhs = (unsigned __int128)(q + 4 - (uint64_t)b) * 10ULL;
                unsigned __int128 rhs = (unsigned __int128)4ULL * g * (n + 3);
                if (lhs > rhs) {
                    break;
                }
#else
                if ((q + 4 - (uint64_t)b) > (4ULL * g * (n + 3)) / 10ULL) {
                    break;
                }
#endif
            }
        }
    }
    return sol;
}

/* Reconstruct Type II decomposition from Theorem thm:gen_criterion:
 * s | (a*n + u), s ≡ -1 (mod 4*u*a).
 */
static void reconstruct_type2(uint64_t n, uint64_t a, uint64_t u, uint64_t s, EsSolution *sol) {
    uint64_t target = a * n + u;
    uint64_t v = target / s;
    uint64_t c = (u + v) / a;
    uint64_t K = (n + c) / 4;
    sol->found = true;
    sol->shift_c = c;
    sol->mod_a = a;
    sol->u = u;
    sol->v = v;
    sol->witness_s = s;
    sol->x = K;
    sol->y = n * a * (K / u);
    sol->z = n * a * (K / v);
    snprintf(sol->branch, sizeof(sol->branch), "%s", "TypeII");
}

/* Reconstruct Type I decomposition from Theorem thm:type_I_two_parameter:
 * s | (a*n + u), gcd(u,a) = 1, s ≡ -n (mod 4*u*a).
 */
static void reconstruct_type1(uint64_t n, uint64_t a, uint64_t u, uint64_t s, EsSolution *sol) {
    uint64_t target = a * n + u;
    uint64_t v = target / s;
    uint64_t alpha = (s + n) / (4ULL * u * a);
    sol->found = true;
    sol->shift_c = 0;
    sol->mod_a = a;
    sol->u = u;
    sol->v = v;
    sol->witness_s = s;
    sol->x = alpha * u * v;
    sol->y = a * alpha * u;
    sol->z = a * alpha * v * n;
    snprintf(sol->branch, sizeof(sol->branch), "%s", "TypeI");
}

/*
 * Two-phase geometric segmented sieve (Algorithm 2 / Corollary cor:simultaneous_sieve):
 * Scans progression [a*n + 1, a*n + J] without factorization.
 * If allow_type1 == false and a_max == 1, executes the exact hyperbolic sieve (a=1, Type II).
 */
static EsSolution solve_segmented_sieve(
    uint64_t n, uint64_t a_max, uint64_t u_max, bool allow_type1) {
    EsSolution sol = {0};

    for (uint64_t a = 1; a <= a_max; a++) {
        uint64_t base = a * n;
        uint64_t hi = base + u_max;
        uint64_t R = isqrt_u64(hi);

        /* Phase A: small divisor s <= R */
        uint64_t s_start = allow_type1 ? 3ULL : 3ULL;
        uint64_t s_step = allow_type1 ? 2ULL : 4ULL; /* Type II needs s ≡ 3 (mod 4); Type I allows any odd s >= 3 */
        for (uint64_t s = s_start; s <= R; s += s_step) {
            uint64_t M = ((base + s) / s) * s; /* smallest multiple of s >= base + 1 */
            for (; M <= hi; M += s) {
                sol.iterations++;
                uint64_t u = M - base;
                if (u > u_max) {
                    break;
                }
                if (a > 1 && gcd_u64(u, a) != 1) {
                    continue;
                }
                uint64_t mod = 4ULL * u * a;
                if ((s + 1) % mod == 0) {
                    reconstruct_type2(n, a, u, s, &sol);
                    return sol;
                }
                if (allow_type1 && ((s + (n % mod)) % mod == 0)) {
                    reconstruct_type1(n, a, u, s, &sol);
                    return sol;
                }
            }
        }

        /* Phase B: small cofactor v <= R */
        for (uint64_t v = 1; v <= R; v++) {
            uint64_t M = ((base + v) / v) * v;
            for (; M <= hi; M += v) {
                sol.iterations++;
                uint64_t u = M - base;
                uint64_t s = M / v;
                if (a > 1 && gcd_u64(u, a) != 1) {
                    continue;
                }
                uint64_t mod = 4ULL * u * a;
                if ((s + 1) % mod == 0) {
                    reconstruct_type2(n, a, u, s, &sol);
                    return sol;
                }
                if (allow_type1 && ((s + (n % mod)) % mod == 0)) {
                    reconstruct_type1(n, a, u, s, &sol);
                    return sol;
                }
            }
        }
    }
    return sol;
}

/*
 * Exhaustive hyperbolic exceptional set H(N) up to bound N (Theorem thm:hyper_char & Prop prop:nine_primes):
 * Computes every prime n <= N, n ≡ 1 (mod 24), NOT representable as 4*alpha*u*v - u - v
 * for any positive integers alpha, u, v (no truncation on u).
 *
 * Uses a compact bitset indexed ONLY over integers k where n = 24*k + 1 <= N,
 * reducing memory to N/192 bytes and visiting only multiples matching n ≡ 1 (mod 24).
 */
static int run_find_H(uint64_t N, bool quiet) {
    if (N < 73) {
        if (!quiet) {
            printf("H(%" PRIu64 ") = {} (count: 0)\n", N);
        }
        return 0;
    }

    uint64_t max_k = (N - 1) / 24;
    size_t nbytes = (size_t)((max_k + 8) / 8);
    uint8_t *marked = (uint8_t *)calloc(nbytes, 1);
    uint8_t *is_comp = (uint8_t *)calloc(nbytes, 1);
    if (!marked || !is_comp) {
        fprintf(stderr, "Error: memory allocation failed for N=%" PRIu64 "\n", N);
        free(marked);
        free(is_comp);
        return 1;
    }

    /* Simple sieve on arithmetic progression n = 24*k + 1 */
    is_comp[0] |= 1U; /* n = 1 is not prime */
    uint64_t limit_p = isqrt_u64(N);
    uint8_t *small_comp = (uint8_t *)calloc((size_t)limit_p + 2, 1);
    if (!small_comp) {
        fprintf(stderr, "Error: memory allocation failed\n");
        free(marked);
        free(is_comp);
        return 1;
    }
    for (uint64_t p = 2; p <= limit_p; p++) {
        if (small_comp[p]) {
            continue;
        }
        for (uint64_t m = p * p; m <= limit_p; m += p) {
            small_comp[m] = 1;
        }
        if (p <= 3) {
            continue;
        }
        /* Find first n = 24*k + 1 >= p^2 divisible by p */
        uint64_t p_mod24 = p % 24;
        /* p * p_inv ≡ 1 (mod 24) because for gcd(p,24)==1, p^2 ≡ 1 (mod 24)! */
        uint64_t inv24 = p_mod24;
        uint64_t k = (p * inv24 - 1) / 24; /* Gives multiple of p ≡ 1 mod 24 equal to p*inv24 or p^2 */
        uint64_t n0 = 24 * k + 1;
        if (n0 < p * p) {
            uint64_t diff = (p * p - n0 + (24 * p - 1)) / (24 * p);
            k += diff * p;
        } else {
            /* If inv24 > p, step back */
            while (k >= p && 24 * (k - p) + 1 >= p * p) {
                k -= p;
            }
        }
        for (; k <= max_k; k += p) {
            is_comp[k >> 3] |= (uint8_t)(1U << (k & 7));
        }
    }
    free(small_comp);

    /* Mark all n = 24*k + 1 = v*s - u with s = 4*alpha*u - 1, u <= ceil((N+1)/3). */
    uint64_t u_bound = (N + 1 + 2) / 3;
    for (uint64_t u = 1; u <= u_bound; u++) {
        uint64_t u_plus_1_mod24 = (u + 1) % 24;
        uint64_t u_plus_1_mod3 = u_plus_1_mod24 % 3;
        uint64_t s_max = N + u;
        uint64_t step_s = 4ULL * u;

        for (uint64_t s = step_s - 1; s <= s_max; s += step_s) {
            uint64_t s_mod3 = s % 3;
            uint64_t v0;
            uint64_t v_step;

            /* Every odd s satisfies s^{-1} ≡ s (mod 8). */
            uint64_t v_mod8 = (u_plus_1_mod24 * (s & 7ULL)) & 7ULL;
            if (s_mod3 == 0) {
                if (u_plus_1_mod3 != 0) {
                    continue;
                }
                v0 = (v_mod8 == 0) ? 8ULL : v_mod8;
                v_step = 8ULL;
            } else {
                uint64_t v_mod3 = (u_plus_1_mod3 * s_mod3) % 3ULL;
                /* CRT combine mod 8 and mod 3 */
                v0 = v_mod8;
                while ((v0 % 3ULL) != v_mod3) {
                    v0 += 8ULL;
                }
                if (v0 == 0) {
                    v0 = 24ULL;
                }
                v_step = 24ULL;
            }

            uint64_t n = v0 * s - u;
            if (n > N) {
                continue;
            }
            uint64_t k = (n - 1) / 24;
            uint64_t dk = (v_step * s) / 24;
            for (; k <= max_k; k += dk) {
                marked[k >> 3] |= (uint8_t)(1U << (k & 7));
            }
        }
    }

    uint64_t total_primes = 0;
    uint64_t h_count = 0;
    if (!quiet) {
        printf("Exceptional set H(%" PRIu64 ") [primes n <= %" PRIu64 ", n ≡ 1 (mod 24), outside a=1 hyperbolic family]:\n", N, N);
    }
    for (uint64_t k = 0; k <= max_k; k++) {
        if ((is_comp[k >> 3] & (1U << (k & 7))) == 0) {
            total_primes++;
            if ((marked[k >> 3] & (1U << (k & 7))) == 0) {
                uint64_t n = 24 * k + 1;
                h_count++;
                if (!quiet) {
                    printf("  [%" PRIu64 "] n = %" PRIu64 "\n", h_count, n);
                }
            }
        }
    }
    if (!quiet) {
        printf("Summary: %" PRIu64 " exceptional prime(s) out of %" PRIu64 " primes n ≡ 1 (mod 24) in [1, %" PRIu64 "].\n",
               h_count, total_primes, N);
    }

    free(marked);
    free(is_comp);
    return 0;
}

static void print_solution(uint64_t n, const EsSolution *s) {
    if (!s->found) {
        printf("n=%" PRIu64 ": no solution found (iterations=%" PRIu64 ")\n", n, s->iterations);
        return;
    }
    bool ok = verify_identity(n, s->x, s->y, s->z);
    printf("n=%" PRIu64 " [%s] 4/%" PRIu64 " = 1/%" PRIu64 " + 1/%" PRIu64 " + 1/%" PRIu64
           " (a=%" PRIu64 ", u=%" PRIu64 ", v=%" PRIu64 ", c=%" PRIu64 ", iter=%" PRIu64 ") verified=%s\n",
           n, s->branch, n, s->x, s->y, s->z,
           s->mod_a, s->u, s->v, s->shift_c, s->iterations, ok ? "YES" : "NO");
}

static void usage(const char *prog) {
    fprintf(stderr,
            "Usage:\n"
            "  %s --check <n> [--depth J] [--exact-alg1] [--type2-only] [--single-c C]\n"
            "      Verify 4/n via Algorithm 1, Hyperbolic Sieve (a=1), and Two-Parameter Sieve.\n"
            "  %s --verify-range <N> [--depth J]\n"
            "      Verify every prime n <= N, n ≡ 1 (mod 24).\n"
            "  %s --find-H <N> [--quiet]\n"
            "      Compute exact exceptional set H(N) outside the a=1 hyperbolic family (Theorem 3.6 & Proposition 3.8).\n",
            prog, prog, prog);
}

int main(int argc, char **argv) {
    if (argc < 3) {
        usage(argv[0]);
        return 1;
    }

    const char *mode = argv[1];
    uint64_t param = strtoull(argv[2], NULL, 10);
    uint64_t depth = 64;
    bool exact_alg1 = false;
    bool allow_type1 = true;
    uint64_t single_c = 0;
    bool quiet = false;

    for (int i = 3; i < argc; i++) {
        if (strcmp(argv[i], "--depth") == 0 && i + 1 < argc) {
            depth = strtoull(argv[++i], NULL, 10);
        } else if (strcmp(argv[i], "--exact-alg1") == 0) {
            exact_alg1 = true;
        } else if (strcmp(argv[i], "--type2-only") == 0) {
            allow_type1 = false;
        } else if (strcmp(argv[i], "--single-c") == 0 && i + 1 < argc) {
            single_c = strtoull(argv[++i], NULL, 10);
        } else if (strcmp(argv[i], "--quiet") == 0) {
            quiet = true;
        } else {
            usage(argv[0]);
            return 1;
        }
    }

    if (strcmp(mode, "--check") == 0) {
        if (!is_prime_u64(param) || (param % 24) != 1) {
            fprintf(stderr, "Error: n=%" PRIu64 " must be a prime ≡ 1 (mod 24).\n", param);
            return 1;
        }
        EsSolution s_alg1 = solve_alg1(param, exact_alg1, single_c);
        EsSolution s_hyp = solve_segmented_sieve(param, 1, (param + 3) / 3, false);
        EsSolution s_two = solve_segmented_sieve(param, depth, depth, allow_type1);

        printf("=== Verification for n = %" PRIu64 " ===\n", param);
        printf("Algorithm 1 (%s):       ", exact_alg1 ? "unaccelerated" : "accelerated");
        print_solution(param, &s_alg1);
        printf("Hyperbolic Sieve (a=1): ");
        print_solution(param, &s_hyp);
        printf("Two-Param Sieve (J=%" PRIu64 "): ", depth);
        print_solution(param, &s_two);
        if (single_c > 0) {
            return 0;
        }
        return (s_alg1.found && verify_identity(param, s_alg1.x, s_alg1.y, s_alg1.z)) ? 0 : 2;
    }

    if (strcmp(mode, "--verify-range") == 0) {
        uint64_t tested = 0;
        uint64_t two_param_resolved = 0;
        uint64_t alg1_resolved = 0;
        uint64_t hyp_exceptions = 0;
        uint64_t max_shift = 0;

        for (uint64_t n = 73; n <= param; n += 24) {
            if (!is_prime_u64(n)) {
                continue;
            }
            tested++;
            EsSolution s_two = solve_segmented_sieve(n, depth, depth, true);
            if (s_two.found && verify_identity(n, s_two.x, s_two.y, s_two.z)) {
                two_param_resolved++;
            }
            EsSolution s_alg1 = solve_alg1(n, exact_alg1, 0);
            if (s_alg1.found && verify_identity(n, s_alg1.x, s_alg1.y, s_alg1.z)) {
                alg1_resolved++;
                if (s_alg1.shift_c > max_shift) {
                    max_shift = s_alg1.shift_c;
                }
            } else {
                fprintf(stderr, "FAIL: Algorithm 1 failed for n=%" PRIu64 "\n", n);
                return 2;
            }
            if (!s_two.found) {
                fprintf(stderr, "NOTE: two-parameter sieve depth J=%" PRIu64 " missed n=%" PRIu64 "\n", depth, n);
            }
            /* Check hyperbolic membership if small enough or fast check */
            EsSolution s_hyp = solve_segmented_sieve(n, 1, (n <= 1000000ULL ? (n + 1) / 3 : 200), false);
            if (!s_hyp.found) {
                hyp_exceptions++;
            }
        }
        printf("Verified %" PRIu64 " primes n <= %" PRIu64 ", n ≡ 1 (mod 24):\n", tested, param);
        printf("  Two-param sieve (J=%" PRIu64 "): %" PRIu64 "/%" PRIu64 " resolved\n", depth, two_param_resolved, tested);
        printf("  Algorithm 1:             %" PRIu64 "/%" PRIu64 " resolved (max shift c=%" PRIu64 ")\n", alg1_resolved, tested, max_shift);
        printf("  Outside hyperbolic cap:  %" PRIu64 " (when exhaustively scanned for n<=1e6 / cap=200 otherwise)\n", hyp_exceptions);
        return (alg1_resolved == tested) ? 0 : 2;
    }

    if (strcmp(mode, "--find-H") == 0) {
        return run_find_H(param, quiet);
    }

    usage(argv[0]);
    return 1;
}
