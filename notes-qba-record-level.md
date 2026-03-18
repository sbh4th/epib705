# Notes on Record-Level QBA and the Multiple Imputation Connection

*From a conversation with Claude, March 2026*

---

## Background

Fox et al. (2021) describe a **record-level** approach to quantitative bias analysis (QBA) as an
alternative to summary-level (contingency table) adjustment. The claimed advantage is that
record-level analysis retains information on other covariates and supports multiple simultaneous
bias adjustments. The broad steps are:

1. Identify the source of bias
2. Select bias parameters
3. Assign probability distributions to each bias parameter
4. Sample from those distributions, apply simple bias-analysis formulas, and sample the
   bias-adjusted effect estimate
5. Save the bias-adjusted estimate and repeat

The computational cost is steep: a study of N records run for 10,000 iterations produces a working
dataset of N × 10,000 rows.

---

## The Core Tension: Collapsing to a Contingency Table

The Fox et al. record-level procedure still requires collapsing the data to a contingency table at
the bias-correction step. The reason is mechanical: the correction formulas for misclassification
(PPV/NPV), unmeasured confounding, and selection bias were all derived for 2×2 tables and take cell
counts as inputs. There is no closed-form version of these corrections that operates directly on a
regression coefficient from a multivariable model.

**The problem:** if you are estimating the effect of binary X on binary Y while controlling for W
and Z, collapsing X and Y to a marginal 2×2 table (ignoring W and Z) and applying the bias
correction there means you have corrected the *marginal* X–Y association, not the *conditional*
association you actually care about. These are different quantities whenever W or Z confound or
modify the effect.

The Fox et al. record-level approach partially addresses this by generating corrected probabilities
at the record level (e.g. PPV and NPV per record) and merging them back before the final regression.
But the bias-correction step upstream still used the marginal table. This is the unresolved tension.

**When does it matter most?** When the source of bias is differential across levels of W or Z — e.g.
misclassification of X varies by W, or the unmeasured confounder's prevalence varies by Z. In those
cases, using a marginal 2×2 correction introduces its own error.

**Proper fix:** Stratify the bias correction by the relevant covariates — collapse to a 2×2 table
*within strata* of W and Z, apply the formulas separately in each stratum, generate stratum-specific
corrected probabilities, merge back, and then run the multivariable regression. This is what Fox et
al. mean by "collapsing with reweighting is possible."

---

## A Cleaner Approach: Direct Record-Level Simulation

For unmeasured confounding specifically, you can bypass the contingency table machinery entirely and
simulate U directly at the record level:

1. Specify a distribution for U (e.g. binary with some prevalence)
2. Specify the relationship between U and X (e.g. as an odds ratio)
3. Specify the relationship between U and Y (e.g. as an odds ratio)
4. For each iteration, generate a value of U for every record consistent with those parameters
5. Fit the regression of Y on X + W + Z + U
6. The coefficient on X is now adjusted for U

No contingency table, no PPV/NPV formulas, no merging. This is simply **multiple imputation of an
unobserved variable** under an assumed data-generating model.

The same logic applies to misclassification. Rather than computing PPV/NPV from a marginal table,
you compute the posterior probability that the true exposure is 1 given the observed exposure and
the assumed sensitivity/specificity — directly via Bayes' theorem applied record-by-record — and
sample from that posterior. The Fox et al. procedure does the same thing, but routes through
PPV/NPV as an intermediate step; that is one way to derive the posterior, not the only way.

**Why does Fox et al. take the harder route?**
- The framework was developed for summary-level data first; the record-level chapter extends that
  machinery rather than rethinking it from scratch.
- Contingency table parameters (sensitivity, specificity, confounder prevalence among exposed vs.
  unexposed) are quantities epidemiologists reason about directly.
- It maintains formal compatibility with summary-level results.

But the contingency table route is not *necessary* — it is largely a historical artifact.

---

## The Multiple Imputation Connection

The direct simulation approach is essentially multiple imputation of missing data (the unobserved
U, or the true exposure). The bias parameters map directly onto imputation model parameters.
Framing QBA this way allows you to handle multivariable models, effect modification, differential
misclassification, and complex covariate structures naturally — tools like `mice` in R are directly
applicable.

### Key citations

- **Beesley et al.** — showed that multiple imputation removes bias due to outcome misclassification
  across a range of sensitivity/specificity scenarios; highlighted advantages of the missing-data
  framing for epidemiologists already familiar with imputation methods.

- **Shepherd et al.** — demonstrated that multiple imputation combined with marginal structural
  models produces consistent estimates in the presence of both confounding and measurement error,
  and showed feasibility in large observational databases.

- **Blum et al.** — proposed a multi-bias approach involving simultaneous adjustment of unmeasured
  confounding, misclassification, and selection bias via imputation and regression weighting, where
  the imputed value corresponds to the probability of the missing data.

- **`unmconf` R package (Bayesian)** — implements Bayesian regression directly modeling one or two
  unmeasured confounders for normal, binary, Poisson, and gamma outcomes; credible intervals achieve
  near-nominal coverage when the unmeasured confounder is explicitly modeled.

- **Simulation-based QBA for study design** — a recent proposal for individual-level simulation-based
  QBA at the design stage, noted for its flexibility and relative ease of implementation for
  nonrandomized studies.

---

## Comparison Code

A companion R script (`qba-record-level-comparison.R`, to be created) illustrates both approaches
on a synthetic dataset where the truth is known (binary Y, binary X, measured covariates W and Z,
unmeasured binary confounder U):

- **Method 1 (Fox et al.):** collapses to the marginal X–Y 2×2 table → applies confounder bias
  correction formulas → merges per-cell probabilities of U back to record level → samples U → fits
  Y ~ X + W + Z + U.
- **Method 2 (direct):** never touches a contingency table → computes P(U=1 | X_i, Y_i) per record
  via Bayes' theorem → samples U from that posterior → fits the same final regression.

Both should recover a bias-adjusted interval that includes the true OR and sits between the naive
estimate and the oracle. The direct approach is more easily extended to differential bias or
effect-modified corrections by modifying only the per-record likelihood step.

---

## What "Oracle" and "True OR" Mean in the Simulation

These are two distinct things:

- **True OR** is a fixed constant — the exact value hardcoded into the data-generating model
  (e.g. `log(1.5)` in the `plogis()` call). It never varies; it is the parameter we are trying
  to recover.

- **Oracle OR** is a finite-sample *estimate* — the result of fitting the regression on the one
  realised dataset with U observed (cheating). Because the study has finite N, this estimate will
  not land exactly on the true OR; it will differ by sampling variability. Re-running with a
  different random seed would give a slightly different oracle each time. With N → ∞ it converges
  to the true OR.

The oracle answers the question: *"How well could the bias correction possibly do on this particular
dataset?"* — not the truth itself, but the best achievable estimate given the data at hand. A
well-functioning QBA method should produce a distribution of bias-adjusted ORs centred near the
oracle, with additional spread reflecting uncertainty about U.

---

## When Stratified Correction Matters (and When It Doesn't)

Running `qba-record-level-comparison.R` with modest differential bias (P(U=1|X=1) ranging from
0.50 to 0.70 across W strata) produces nearly identical results for all three methods. The reason
is that **W is included in the final regression**. Even when Methods 1 and 2 use the wrong marginal
P(U=1|X), the regression on Y ~ X + W + Z + U_sim partially absorbs the imputation errors via the
W coefficient, because W predicts both the true U (through the differential prevalence) and Y
directly. The stratified correction (Method 2b) produces slightly better imputations, but the
regression compensates enough that the downstream ORs converge.

### The condition under which stratification clearly helps

Stratified correction gives a visibly better answer when:

1. **The differential bias is extreme** — P(U=1|X=1) swings from near-0 to near-1 across W strata,
   so the marginal prior is badly wrong for both subgroups simultaneously.
2. **The unmeasured confounder is strong** — a large OR(U→Y) means that even modest mis-imputation
   of U produces meaningful residual confounding.
3. **Both conditions hold together** — the regression then cannot absorb the systematic imputation
   errors, even with W in the model.

The script `qba-differential-bias-demo.R` constructs exactly this scenario:

| Cell | True P(U=1\|X,W) | Marginal prior used | Error |
|------|-----------------|---------------------|-------|
| X=1, W=1 | 0.92 | 0.50 | −0.42 |
| X=1, W=0 | 0.08 | 0.50 | +0.42 |
| X=0, W=1 | 0.50 | 0.27 | −0.23 |
| X=0, W=0 | 0.04 | 0.27 | +0.23 |

With OR(U→Y) = 3 and N = 5,000, this mis-imputation is large enough that the regression cannot
compensate, and Method 2b (stratified) should sit noticeably closer to the oracle than Methods 1
and 2. A companion diagnostic plot shows the *imputed* P(U=1) distributions by W stratum directly,
making the mis-imputation visible before looking at the final ORs.

### The broader principle

When the variable that drives differential bias (W here) is in the final regression model, the
marginal correction is surprisingly robust — the regression does part of the work. Stratified
correction matters most when:

- The bias-driving variable is **not** in the final model (unmeasured, or deliberately excluded),
  so the regression cannot compensate at all; or
- The differential bias is so extreme that the mis-imputed U values are essentially noise for one
  subgroup, attenuating the recovered U effect and leaving residual confounding that the W
  coefficient alone cannot undo.
