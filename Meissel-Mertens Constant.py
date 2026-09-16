#!/usr/bin/env python3
"""
Meissel-Mertens Constant Calculator (HPC OEIS Edition)
======================================================
Calculates Meissel-Mertens constant (M_1) to exactly [N] significant digits using 
Möbius inversion series over logarithmic Zeta functions, 12-core parallel chunking, 
C-accelerated gmpy2 math, and strict OEIS truncation formatting.
"""
from __future__ import annotations

import argparse
import gc
import math
import multiprocessing as mp
import os
import sys
import time

import mpmath



os.environ['MPMATH_GMPY2'] = '1'

sys.set_int_max_str_digits(0)

NUM_WORKERS = 12


def mobius(n) -> int:
    """Mobius.
    
    Args:
        n:
    
    Returns:
        int: Result of type int
    
    """
    if n == 1:
        return 1
    p = 0
    d = 2
    temp = n
    while d * d <= temp:
        if temp % d == 0:
            p += 1
            temp //= d
            if temp % d == 0:
                return 0
        d += 1
    if temp > 1:
        p += 1
    return -1 if p % 2 != 0 else 1


def worker_mobius_chunk(args) -> Any:
    """Worker function for mobius chunk.
    
    Args:
        args:
    
    Returns:
        Any: The computed result
    
    """
    start, end, dps = args
    mpmath.mp.dps = dps
    ctx = mpmath.mp

    partial_sum = ctx.mpf(0)
    for k in range(start, end):
        mu = mobius(k)
        if mu != 0:
            term = (ctx.mpf(mu) / ctx.mpf(k)) * ctx.ln(ctx.zeta(k))
            partial_sum += term
    return partial_sum


def save_oeis_files(constant_name, digits_str, target_digits):
    """Save oeis files to file.
    
    Args:
        constant_name:
        digits_str:
        target_digits:
    
    """
    clean_digits = digits_str.replace(".", "")[:target_digits]
    
    raw_filename = f"{constant_name}_{target_digits}_digits.txt"
    with open(raw_filename, "w", encoding="utf-8") as f:
        f.write(clean_digits)
    print(f"Saved raw digit output to {raw_filename}")

    b_filename = f"b_file_{constant_name}_{target_digits}.txt"
    with open(b_filename, "w", encoding="utf-8") as f:
        for idx, digit in enumerate(clean_digits, start=1):
            f.write(f"{idx} {digit}\n")
    print(f"Saved OEIS b-file output to {b_filename}")


def compute_meissel_mertens_hpc(target_digits) -> Any:
    """Compute meissel mertens hpc using optimized algorithms.
    
    Args:
        target_digits:
    
    Returns:
        Any: The computed result
    
    """
    dps_working = target_digits + 50
    mpmath.mp.dps = dps_working
    ctx = mpmath.mp

    res = ctx.euler
    terms = int(target_digits * 3.5) + 20
    chunk_size = math.ceil((terms - 2) / NUM_WORKERS)

    chunks = []
    for i in range(NUM_WORKERS):
        start = 2 + i * chunk_size
        end = min(terms, 2 + (i + 1) * chunk_size)
        if start < terms:
            chunks.append((start, end, dps_working))

    with mp.Pool(processes=NUM_WORKERS) as pool:
        results = pool.map(worker_mobius_chunk, chunks)

    for chunk_sum in results:
        res += chunk_sum

    del results
    gc.collect()

    meissel_str = ctx.nstr(res, dps_working)
    clean_digits = meissel_str.replace(".", "")[:target_digits]

    del res
    gc.collect()

    save_oeis_files("Meissel_Mertens", clean_digits, target_digits)
    return clean_digits


def main():
    """Entry point — parse arguments and run the main computation.
    
    """
    parser = argparse.ArgumentParser(description="HPC Meissel-Mertens OEIS Calculator")
    parser.add_argument("-n", "--digits", type=int, default=3000, help="Target digits (default: 1000)")
    args = parser.parse_args()

    t0 = time.time()
    digits = compute_meissel_mertens_hpc(args.digits)
    t1 = time.time()

    print(f"Execution finished in {t1 - t0:.4f} seconds using {NUM_WORKERS} cores.")

if __name__ == "__main__":
    main()
