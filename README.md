[![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=dnafinder/quade&file=quadetest.m)

# quadetest

Quade test for nonparametric two-way ANOVA in complete block designs, implemented in MATLAB.

The Quade test is an alternative to the Friedman test. It removes block effects by ranking within blocks, but also weights blocks according to the magnitude of their within-block ranges. This makes it often more powerful than the Friedman test. Conceptually, Quade is to Friedman as Wilcoxon signed-rank is to the sign test, and it reduces to Wilcoxon when there are only two treatments.

## Syntax

- [stats, mc] = quadetest(X)
- [stats, mc] = quadetest(X, 'Alpha', ALPHA, 'PostHoc', POSTHOC, 'Display', DISPLAY)

## Inputs

### Required

- X  
  Numeric matrix of size (Blocks × Treatments), real, finite, non-NaN, non-empty.  
  Rows represent blocks, columns represent treatments.

### Name–Value options

- 'Alpha'  
  Significance level (scalar in (0,1), default 0.05).

- 'PostHoc'  
  Logical-like flag to enable/disable post-hoc multiple comparisons (default true).  
  When true and the global Quade test rejects the null hypothesis, pairwise comparisons between treatments are performed.

- 'Display'  
  Logical-like flag to control command-window output (default true).  
  If false, no text is printed; only stats and mc are returned.

For 'PostHoc' and 'Display', the following are accepted as true/false:  
true/false, 1/0, "on"/"off", "yes"/"no", "true"/"false" (case-insensitive).

## Outputs

### STATS structure

- stats.nObs       – number of observations (blocks × treatments)  
- stats.blocks     – number of blocks (rows)  
- stats.treatments – number of treatments (columns)  
- stats.W          – Quade statistic  
- stats.F_df_num   – numerator degrees of freedom (treatments − 1)  
- stats.F_df_denom – denominator degrees of freedom (blocks − 1) × (treatments − 1)  
- stats.F_p        – two-tailed p-value from the F approximation  
- stats.alpha      – significance level used  
- stats.rejectNull – true if the treatments do not have identical effects (F_p < alpha)

If no output argument is specified, the function prints the same information in the Command Window when 'Display' is true.

### MC structure (post-hoc)

Returned when 'PostHoc' is true and stats.rejectNull is true; otherwise empty.

- mc.method      – 'Quade-Conover-type LSD'  
- mc.Rdiff       – Treatments × Treatments matrix of absolute differences among treatment scores (Ti)  
- mc.cv          – critical value for the differences  
- mc.pvalue      – empty (this procedure is expressed in terms of a critical difference, not explicit p-values)  
- mc.significant – logical matrix (lower triangle), true where absolute differences exceed the critical value

## Method (breve)

- For each block (row) of X, the observations are ranked using tiedrank, yielding a rank matrix R.  
- The within-block range is computed for each block and ranked, giving a weight vector Q.  
- A weighted, centred rank matrix rij is constructed as (R − (c + 1)/2) .* Q, where c is the number of treatments.  
- Treatment scores Ti are obtained by summing rij across blocks, and a Quade statistic W is computed, which is approximated by an F distribution with (c − 1, (r − 1)(c − 1)) degrees of freedom.  
- If the F-based p-value is less than Alpha, the null hypothesis of identical treatment effects is rejected.  
- When the global test is significant and 'PostHoc' is true, a Quade-Conover-type LSD procedure is applied: absolute differences among treatment scores are compared to a critical value based on the t distribution and the Quade residual term; differences exceeding this value are flagged as significant.

## Example

    x = [115 142  36  91  28;
          28  31   7  21   6;
         220 311 108  51 117;
          82  56  24  46  33;
         256 298 124  46  84;
         294 322 176  54  86;
          98  87  55  84  25];

    [stats, mc] = quadetest(x);

## Citation

If you use this function in a scientific publication, please cite it as:

Cardillo G. (2009). QUADETEST: Quade test for non parametric two way ANalysis Of VAriance. Available on GitHub: https://github.com/dnafinder/quade
