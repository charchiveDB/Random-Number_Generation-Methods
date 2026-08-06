# Random-Number_Generation-Methods

# The Ziggurat Method: Random Number Generation

An R implementation of the **Ziggurat algorithm** for efficiently generating random variables from the **Normal** and **Exponential** distributions.

## Overview

This project implements the Ziggurat method described by Marsaglia and Tsang (2000), one of the fastest algorithms for sampling from unimodal probability distributions. The algorithm partitions a probability density into equal-area rectangular regions ("layers"), allowing over 99% of generated samples to be accepted without evaluating the probability density function.

Our implementation recreates the algorithm in **R**, adapting the original C implementation while preserving the core ideas behind the method.

## Features

- Implementation of `rnormzig()` for Normal random variable generation
- Implementation of `rexpzig()` for Exponential random variable generation
- Automatic construction of Ziggurat tables using `zigtable()`
- Distribution-specific tail sampling
- Chi-square goodness-of-fit validation
- Comparison against R's built-in `rnorm()` and `rexp()`

## Algorithm

The implementation consists of three primary components:

1. **Construct the Ziggurat**
   - Divide the density into 256 equal-area regions.
   - Compute rectangle boundaries iteratively using the inverse PDF.
   - Determine the optimal rectangle area using R's `uniroot()`.

2. **Generate Samples**
   - Randomly choose one rectangle.
   - Draw a uniform sample inside the rectangle.
   - Immediately accept samples in the fast-accept region.
   - Otherwise perform an acceptance-rejection test or invoke the distribution-specific tail generator.

3. **Distribution-Specific Adjustments**
   - **Normal:** Randomly assign a positive or negative sign using symmetry.
   - **Exponential:** Sample directly since the distribution is only defined on positive values.

## Validation

The implementation was validated using a **chi-square goodness-of-fit test**.

### Normal Distribution

- 1,000,000 generated samples
- 200 bins over [-7, 7]

| Generator | χ² Statistic |
|-----------|-------------:|
| `rnormzig()` | 174.20 |
| R `rnorm()` | 184.42 |

### Exponential Distribution

- 1,000,000 generated samples
- 200 bins over [0, 7]

| Generator | χ² Statistic |
|-----------|-------------:|
| `rexpzig()` | 167.75 |
| R `rexp()` | 199.76 |

The 95% critical value (199 degrees of freedom) is **232.912**.

All generated samples produced χ² values below the critical threshold, demonstrating that both implementations accurately reproduce the intended probability distributions.

## Challenges

The original Ziggurat implementation was written in **C**, making direct translation into R nontrivial.

Major challenges included:

- Converting zero-based indexing to R's one-based indexing.
- Replacing low-level bitwise optimizations unavailable in R.
- Adapting C-specific sampling optimizations while maintaining readability.
- Managing the performance limitations of interpreted R code.

Although the resulting implementation is slower than optimized C versions, it clearly demonstrates the theory behind the Ziggurat algorithm and serves as an educational implementation.

## Repository Structure

```
.
├── rnormzig.R          # Normal Ziggurat generator
├── rexpzig.R           # Exponential Ziggurat generator
├── zigtable.R          # Computes rectangle boundaries
├── chisquare_test.R    # Goodness-of-fit validation
├── examples.R          # Example usage
└── README.md
```

## References

1. Marsaglia, G., & Tsang, W. W. (2000). *The Ziggurat Method for Generating Random Variables*. Journal of Statistical Software, 5(8).

2. Leong, P. H. W., et al. (2005). *A Comment on the Implementation of the Ziggurat Method*. Journal of Statistical Software, 12(7).

3. NumPy Developers. *Random Number Generator Documentation*.

## Authors

Charlotte Huang and collaborators

