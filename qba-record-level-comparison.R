# =============================================================================
# qba-record-level-comparison.R
#
# Compares two approaches to probabilistic bias analysis (PBA) for an
# unmeasured binary confounder U, using record-level data:
#
#   Method 1 — Fox et al. (2021): collapse to a marginal 2×2 table of
#              (X, Y), compute P(U=1 | X-cell, Y-cell) from that table,
#              merge back to records, sample U, run multivariable regression.
#
#   Method 2 — Direct simulation: compute P(U=1 | X_i, Y_i) per record
#              via Bayes' theorem, sample U, run the same regression.
#              No contingency table required.
#
#   Method 2b — Direct simulation, conditioned on W: extends Method 2 by
#              computing P(U=1 | X_i, Y_i, W_i) using stratum-specific
#              observed risks. Shows how the direct approach naturally
#              handles differential bias without any extra machinery.
#
# Data-generating model (truth is known):
#   W ~ Bernoulli(0.5)           — measured covariate / confounder
#   Z ~ Bernoulli(0.4)           — measured covariate
#   X | W ~ Bernoulli(plogis(-0.5 + 0.8*W))
#   U | X ~ Bernoulli(p1 if X=1, p0 if X=0)  — UNMEASURED
#     where p1 = P(U=1|X=1) varies by W (differential bias):
#       P(U=1 | X=1, W=1) = 0.70,  P(U=1 | X=1, W=0) = 0.50
#       P(U=1 | X=0, W=1) = 0.35,  P(U=1 | X=0, W=0) = 0.20
#   Y | X,W,Z,U ~ Bernoulli(plogis(-2 + log(1.5)*X + 0.5*W + 0.3*Z + log(2)*U))
#
# Because P(U=1|X) varies by W, the marginal 2×2 correction (Methods 1 & 2)
# will be a W-averaged approximation, while the W-stratified correction
# (Method 2b) recovers the truth more faithfully.
#
# Bias parameters supplied to each method:
#   Marginal (Methods 1 & 2): p1 = mean P(U=1|X=1) across W, likewise p0
#   Stratified (Method 2b):   separate (p1_W1, p0_W1) and (p1_W0, p0_W0)
#   OR_UY = 2.0 for all methods (assumed known and correctly specified)
#
# =============================================================================

library(tidyverse)

set.seed(123)

# ── 0. Constants ──────────────────────────────────────────────────────────────

N       <- 2000    # study size
n_iter  <- 1000    # Monte Carlo iterations

true_or_xy <- 1.5  # true OR for X → Y (conditional on W, Z, U)
true_or_uy <- 2.0  # true OR for U → Y

# True P(U=1 | X, W) — differential by W
true_p_u <- tribble(
  ~X, ~W, ~p_u,
   1,  1,  0.70,
   1,  0,  0.50,
   0,  1,  0.35,
   0,  0,  0.20
)

# Marginal P(U=1|X) averaged over W=Bernoulli(0.5)
# (what Methods 1 & 2 will use — an approximation)
marginal_p1 <- 0.5 * 0.70 + 0.5 * 0.50   # P(U=1|X=1), averaged over W
marginal_p0 <- 0.5 * 0.35 + 0.5 * 0.20   # P(U=1|X=0), averaged over W

cat(sprintf("Marginal bias parameters supplied to Methods 1 & 2:\n"))
cat(sprintf("  P(U=1|X=1) = %.2f\n", marginal_p1))
cat(sprintf("  P(U=1|X=0) = %.2f\n", marginal_p0))
cat(sprintf("  OR(U→Y)    = %.2f\n\n", true_or_uy))


# ── 1. Generate synthetic data ────────────────────────────────────────────────

dat <- tibble(id = 1:N,
              W  = rbinom(N, 1, 0.5),
              Z  = rbinom(N, 1, 0.4)) |>
  mutate(
    X = rbinom(N, 1, plogis(-0.5 + 0.8 * W))
  ) |>
  left_join(true_p_u, by = c("X", "W")) |>
  mutate(
    U = rbinom(N, 1, p_u),
    Y = rbinom(N, 1, plogis(-2 + log(true_or_xy)*X + 0.5*W + 0.3*Z + log(true_or_uy)*U))
  ) |>
  select(id, X, Y, W, Z, U)


# ── 2. Reference models ───────────────────────────────────────────────────────

fit_naive  <- glm(Y ~ X + W + Z,     data = dat, family = binomial())
fit_oracle <- glm(Y ~ X + W + Z + U, data = dat, family = binomial())

or_naive  <- exp(coef(fit_naive)["X"])
or_oracle <- exp(coef(fit_oracle)["X"])

cat(sprintf("Naive OR  (ignores U):    %.3f\n", or_naive))
cat(sprintf("Oracle OR (observes U):   %.3f\n", or_oracle))
cat(sprintf("True OR:                  %.3f\n\n", true_or_xy))


# ── 3. Core helper: P(U=1 | X, Y) via Bayes' theorem ─────────────────────────
#
# Given P(U=1|X) = p_x, OR(U→Y) = or_uy, and the observed marginal risk
# r_x = P(Y=1|X=x), we first solve for the baseline risk
# q = P(Y=1|U=0, X=x) via the quadratic:
#
#   r_x = p_x * [or_uy*q / (1 + (or_uy-1)*q)] + (1-p_x)*q
#
# Then apply Bayes:
#   P(U=1|X=x,Y=y) = P(Y=y|U=1,X=x)*p_x /
#                    [P(Y=y|U=1,X=x)*p_x + P(Y=y|U=0,X=x)*(1-p_x)]
#
# Vectorised: x_val, y_val, p_x, r_x may all be vectors of the same length.

posterior_u1 <- function(y_val, p_x, or_uy, r_x) {
  # Quadratic coefficients for q = P(Y=1|U=0,X)
  a    <- (1 - p_x) * (or_uy - 1)
  b    <- p_x * or_uy + (1 - p_x) - r_x * (or_uy - 1)
  cc   <- -r_x
  disc <- pmax(b^2 - 4 * a * cc, 0)

  q <- if_else(
    abs(a) < 1e-10,
    r_x,                                   # or_uy == 1 edge case
    (-b + sqrt(disc)) / (2 * a)
  ) |> pmax(1e-6) |> pmin(1 - 1e-6)

  # Conditional risks
  py_u1 <- (or_uy * q / (1 + (or_uy - 1) * q)) |> pmax(1e-6) |> pmin(1 - 1e-6)
  py_u0 <- q

  # Likelihoods for observed Y
  lik_u1 <- if_else(y_val == 1, py_u1, 1 - py_u1)
  lik_u0 <- if_else(y_val == 1, py_u0, 1 - py_u0)

  # Posterior
  (lik_u1 * p_x) / (lik_u1 * p_x + lik_u0 * (1 - p_x))
}


# ── 4. Method 1: Fox et al. (marginal 2×2 table → record level) ──────────────

run_fox <- function(data, bias_p1, bias_p0, bias_or_uy) {
  # Marginal observed risks P(Y=1|X) from the full dataset
  r1 <- mean(data$Y[data$X == 1])
  r0 <- mean(data$Y[data$X == 0])

  # Compute P(U=1|X,Y) for the four cells of the 2×2 table
  cell_probs <- tribble(
    ~X, ~Y, ~p_u_cell,
     1,  1,  posterior_u1(1, bias_p1, bias_or_uy, r1),
     1,  0,  posterior_u1(0, bias_p1, bias_or_uy, r1),
     0,  1,  posterior_u1(1, bias_p0, bias_or_uy, r0),
     0,  0,  posterior_u1(0, bias_p0, bias_or_uy, r0)
  )

  data |>
    left_join(cell_probs, by = c("X", "Y")) |>
    mutate(U_sim = rbinom(n(), 1, p_u_cell)) |>
    (\(d) glm(Y ~ X + W + Z + U_sim, data = d, family = binomial()))() |>
    coef() |> (\(b) exp(b["X"]))()
}


# ── 5. Method 2: Direct per-record simulation (marginal) ──────────────────────
#
# Identical information to Method 1 — uses the same marginal r_x and the
# same bias parameters — but skips the table entirely and computes the
# Bayesian posterior record-by-record.

run_direct <- function(data, bias_p1, bias_p0, bias_or_uy) {
  r1 <- mean(data$Y[data$X == 1])
  r0 <- mean(data$Y[data$X == 0])

  data |>
    mutate(
      p_x    = if_else(X == 1, bias_p1, bias_p0),
      r_x    = if_else(X == 1, r1,      r0),
      p_u    = posterior_u1(Y, p_x, bias_or_uy, r_x),
      U_sim  = rbinom(n(), 1, p_u)
    ) |>
    (\(d) glm(Y ~ X + W + Z + U_sim, data = d, family = binomial()))() |>
    coef() |> (\(b) exp(b["X"]))()
}


# ── 6. Method 2b: Direct simulation conditioned on W ─────────────────────────
#
# Extension of Method 2 that uses W-stratified observed risks r_x_w and
# W-stratified bias parameters. Because the bias is actually differential
# by W in this DGP, this method should outperform Methods 1 & 2.
#
# In practice you would elicit separate bias parameters for each W stratum;
# here we supply the true stratum-specific values to give the method its
# best-case advantage.

run_direct_stratified <- function(data, bias_p_u_df, bias_or_uy) {
  # Stratum-specific observed risks P(Y=1 | X, W)
  r_xw <- data |>
    group_by(X, W) |>
    summarise(r_xw = mean(Y), .groups = "drop")

  data |>
    left_join(bias_p_u_df, by = c("X", "W")) |>   # p_u = P(U=1|X,W)
    left_join(r_xw,         by = c("X", "W")) |>   # r_xw = P(Y=1|X,W)
    mutate(
      p_u_post = posterior_u1(Y, p_u, bias_or_uy, r_xw),
      U_sim    = rbinom(n(), 1, p_u_post)
    ) |>
    (\(d) glm(Y ~ X + W + Z + U_sim, data = d, family = binomial()))() |>
    coef() |> (\(b) exp(b["X"]))()
}


# ── 7. Run all three methods ──────────────────────────────────────────────────

cat("Running", n_iter, "iterations across three methods...\n")

results <- map_dfr(seq_len(n_iter), \(i) {
  tibble(
    iteration  = i,
    fox        = run_fox(
                   dat, marginal_p1, marginal_p0, true_or_uy),
    direct     = run_direct(
                   dat, marginal_p1, marginal_p0, true_or_uy),
    direct_strat = run_direct_stratified(
                   dat, true_p_u, true_or_uy)
  )
}, .progress = TRUE)


# ── 8. Summarise ─────────────────────────────────────────────────────────────

summary_tbl <- results |>
  pivot_longer(c(fox, direct, direct_strat),
               names_to = "method", values_to = "or") |>
  group_by(method) |>
  summarise(
    median = median(or),
    mean   = mean(or),
    p2.5   = quantile(or, 0.025),
    p97.5  = quantile(or, 0.975),
    .groups = "drop"
  ) |>
  mutate(method = recode(method,
    fox          = "Method 1: Fox et al. (marginal 2×2)",
    direct       = "Method 2: Direct (marginal)",
    direct_strat = "Method 2b: Direct (stratified on W)"
  ))

cat("\n── Summary ───────────────────────────────────────────────\n")
cat(sprintf("  True OR:                %.3f\n", true_or_xy))
cat(sprintf("  Naive OR  (ignores U):  %.3f\n", or_naive))
cat(sprintf("  Oracle OR (observes U): %.3f\n\n", or_oracle))
print(summary_tbl, n = Inf)


# ── 9. Plot ───────────────────────────────────────────────────────────────────

plot_dat <- results |>
  pivot_longer(c(fox, direct, direct_strat),
               names_to = "method", values_to = "or") |>
  mutate(method = recode(method,
    fox          = "Method 1: Fox et al.\n(marginal 2×2 → record level)",
    direct       = "Method 2: Direct simulation\n(marginal, per-record Bayes)",
    direct_strat = "Method 2b: Direct simulation\n(stratified on W)"
  ))

ref_lines <- tribble(
  ~or,        ~label,    ~colour,
  or_naive,   "Naive",   "grey40",
  or_oracle,  "Oracle",  "steelblue",
  true_or_xy, "True OR", "black"
)

ggplot(plot_dat, aes(x = or, fill = method, colour = method)) +
  geom_density(alpha = 0.35, linewidth = 0.6) +
  geom_vline(data = ref_lines,
             aes(xintercept = or, colour = NULL, linetype = label),
             colour = c("grey40", "steelblue", "black"),
             linewidth = 0.8) +
  scale_linetype_manual(values = c("Naive" = "dashed",
                                   "Oracle" = "dotdash",
                                   "True OR" = "solid")) +
  scale_fill_manual(values  = c("#E69F00", "#56B4E9", "#009E73")) +
  scale_colour_manual(values = c("#E69F00", "#56B4E9", "#009E73")) +
  labs(
    title    = "Bias-adjusted OR distributions: Fox et al. vs. direct simulation",
    subtitle = sprintf(
      "N = %d, %d iterations  |  True OR = %.1f  |  Naive = %.3f  |  Oracle = %.3f",
      N, n_iter, true_or_xy, or_naive, or_oracle
    ),
    x        = "Bias-adjusted OR for X (adjusted for W, Z, U)",
    y        = "Density",
    fill     = "Method",
    colour   = "Method",
    linetype = "Reference"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom",
        legend.box      = "vertical")

ggsave("qba-record-level-comparison.png", width = 10, height = 6, dpi = 150)
cat("\nPlot saved to qba-record-level-comparison.png\n")
