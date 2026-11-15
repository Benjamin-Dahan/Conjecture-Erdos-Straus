#!/usr/bin/env python3
import os
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BIN_CANDIDATES = [
    os.path.join(ROOT, "build", "erdos_straus_sieve"),
    os.path.join(ROOT, "build", "erdos_straus_sieve.exe"),
]
BIN = next((p for p in BIN_CANDIDATES if os.path.exists(p)), None)
if BIN is None:
    print("FAIL: binary not built.")
    sys.exit(1)


def main():
    # 1. Check standard prime (73), Example 3.4 (20353), and wild prime outside a=1 (409, 102001)
    for n in (73, 409, 20353, 102001, 8813281):
        out = subprocess.check_output([BIN, "--check", str(n), "--depth", "10"], text=True)
        assert "verified=YES" in out, out
    # Confirm 409 shows no solution in Hyperbolic Sieve (a=1), Type I in simultaneous sieve, and Type II at a=2 with --type2-only
    out_409 = subprocess.check_output([BIN, "--check", "409", "--depth", "5"], text=True)
    assert "Hyperbolic Sieve (a=1): n=409: no solution found" in out_409, out_409
    assert "Two-Param Sieve (J=5): n=409 [TypeI]" in out_409 and "verified=YES" in out_409, out_409
    out_409_t2 = subprocess.check_output([BIN, "--check", "409", "--depth", "5", "--type2-only"], text=True)
    assert "Two-Param Sieve (J=5): n=409 [TypeII]" in out_409_t2 and "verified=YES" in out_409_t2, out_409_t2

    # 2. Check n=100129 single-c rejection counts match Theorem thm:reject_cost measurements
    for c, expected_it in ((3, 1670), (7, 6439), (19, 8697)):
        out = subprocess.check_output([BIN, "--check", "100129", "--single-c", str(c)], text=True)
        assert f"iterations={expected_it})" in out, out

    # 3. Verify all primes <= 100,000 (n ≡ 1 mod 24)
    out_range = subprocess.check_output([BIN, "--verify-range", "100000", "--depth", "64"], text=True)
    assert "100000" in out_range and "409" not in out_range.splitlines()[0], out_range

    # 4. Find H(10^6) and match Xu's 9 wild primes
    out_h = subprocess.check_output([BIN, "--find-H", "1000000"], text=True)
    assert "9 exceptional prime(s) out of 9732 primes" in out_h, out_h

    print("ALL C CLI TESTS PASSED.")


if __name__ == "__main__":
    main()
