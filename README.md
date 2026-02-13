# Erdős-Straus Conjecture: HPC & Theoretical Approach

![Language](https://img.shields.io/badge/Language-C_%7C_Python-blue.svg) ![Focus](https://img.shields.io/badge/Focus-High_Performance_Computing-green.svg) ![Status](https://img.shields.io/badge/Status-Research_Complete-orange.svg)

## 📄 Abstract
This research project investigates the **Erdős-Straus conjecture** ($4/n = 1/x + 1/y + 1/z$) through a novel theoretical approach based on **modular arithmetic** and high-performance computational verification.

## ⚡ Performance & Optimization
The project compares a baseline Python implementation versus an optimized C implementation using parallel processing.

* **Python (NumPy/SymPy):** Used for prototyping and understanding the algorithmic logic.
* **C + OpenMP:** Rewritten algorithms utilizing low-level memory management and parallel processing.
* **Result:** Achieved a **95% execution speedup**, allowing verification of hypotheses up to $N = 10^{10}$.

## 🧮 Theoretical Insight (The "3361 Case")
A key finding of this research (detailed in the attached PDF) is the handling of edge cases where standard modular congruences fail. notably **$n = 3361$**, which is the only known exception to the standard algorithmic coverage in the first $10^{10}$ integers, requiring a specific solving path.

## 📂 Project Structure
* `Research_Paper_Erdos_Straus.pdf`: Full mathematical paper explaining the modular arithmetic approach and results.
* `erdos_straus_optimized.c`: High-performance C implementation for large numbers ($N \ge 10^{10}$).
* `erdos_straus_basic.py`: Python implementation for educational purpose and logic verification.

## 🚀 How to Run

### 1. C Implementation (Optimized)
Compile the C code (ensure you have GCC installed):
```bash
gcc erdos_straus_optimized.c -o erdos_optimizer
