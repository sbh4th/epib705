# =============================================================================
# qba-differential-bias-demo.R
#
# Demonstrates the scenario where stratified QBA (Method 2b) clearly
# outperforms the standard marginal correction (Methods 1 & 2).
#
# WHY THE PREVIOUS SCRIPT DIDN'T SHOW A DIFFERENCE:
#   In qba-record-level-comparison.R the differential bias by W was modest
#   (P(U=1|X=1) ranging from 0.50 to 0.70 across W strata). When the bias
#   is that mild and W is in the final regression, the model partially absorbs
#   the imputation errors via the W coefficient. Methods converge.
#
# THE CONDITION THAT MAKES STRATIFICATION MATTER:
#   P(U=1|X) must vary SO dramatically across W that the marginal prior is
#   badly wrong for BOTH subgroups at once:
#
#     P(U=1|X=1, W=1) = 0.92   ← U almost certain among exposed, W=1
#     P(U=1|X=1, W=0) = 0.08   ← U rare among exposed, W=0
#     Marginal P(U=1|X=1)       = 0.50  ← wrong by ~0.40 for BOTH groups
#
#   Even with W in the final regression, such severe mis-imputation creates
#   a spurious U_sim–W correlation that the regression cannot fully untangle.
#
# SETUP:
#   W ~ Bernoulli(0.5)         — measured, included in final model
#   X | W ~ Bernoulli(plogis(-0.5 + 0.8*W))
#   U | X, W:  (extreme differential prevalence)
#     P(U=1|X=1,W=1) = 0.92,  P(U=1|X=1,W=0) = 0.08
#     P(U=1|X=0,W=1) = 0.50,  P(U=1|X=0,W=0) = 0.04
#   Y | X,W,U ~ Bernoulli(plogis(-2 + log(1.5)*X + 0.8*W + log(3)*U))
#     True OR(X→Y) = 1.5,  OR(U→Y) = 3  (strong unmeasured confounder)
#
# WHAT WE EXPECT:
#   Methods 1 & 2 (marginal) will use P(U=1|X=1)=0.50 for all X=1 subjects,
#   massively over-imputing U=1 for W=0 subjects and under-imputing for W=1.
#   This attenuates the recovered U effect and leaves residual confounding by U.
#   Method 2b (stratified) uses the correct W-specific priors and should
#   centre near the oracle.
# =============================================================================

library(tidyverse)

set.seed(2026)

# ── 0. Constants ──────────────────────────────────────────────────────────────

N       <- 5000
n_iter  <- 1000

true_or_xy <- 1.5
true_or_uy <- 3.0    # strong confounder — makes mis-imputation costly

true_p_u <- tribble(
  ~X, ~W,  ~p_u,
   1,  1,   0.92,
   1,  0,   0.08,
   0,  1,   0.50,
   0,  0,   0.04
)

# Marginal P(U=1|X), averaging over W ~ Bernoulli(0.5)
marginal_p1 <- 0.5 * 0.92 + 0.5 * 0.08   # = 0.50
marginal_p0 <- 0.5 * 0.50 + 0.5 * 0.04   # = 0.27

cat("── Bias parameters ───────────────────────────────────────\n")
cat(sprintf("  True P(U=1|X=1,W=1) = %.2f\n", 0.92))
cat(sprintf("  True P(U=1|X=1,W=0) = %.2f\n", 0.08))
cat(sprintf("  Marginal P(U=1|X=1) = %.2f  ← wrong by ±0.42 for each W group\n\n",
            marginal_p1))


# ── 1. Generate data ──────────────────────────────────────────────────────────

dat <- tibble(
  id = 1:N,
  W  = rbinom(N, 1, 0.5),
  Z  = rbinom(N, 1, 0.4)
) |>
  mutate(X = rbinom(N, 1, plogis(-0.5 + 0.8 * W))) |>
  left_join(true_p_u, by = c("X", "W")) |>
  mutate(
    U = rbinom(N, 1, p_u),
    Y = rbinom(N, 1, plogis(-2 + log(true_or_xy)*X + 0.8*W + log(true_or_uy)*U))
  ) |>
  select(id, X, Y, W, Z, U)


# ── 2. Reference models ───────────────────────────────────────────────────────

or_naive  <- glm(Y ~ X + W + Z,     data = dat, family = binomial()) |>
  coef() |> (\(b) exp(b["X"]))()

or_oracle <- glm(Y ~ X + W + Z + U, data = dat, family = binomial()) |>
  coef() |> (\(b) exp(b["X"]))()

cat(sprintf("Naive OR  (ignores U):    %.3f\n", or_naive))
cat(sprintf("Oracle OR (observes U):   %.3f\n", or_oracle))
cat(sprintf("True OR:                  %.3f\n\n", true_or_xy))


# ── 3. Posterior P(U=1|X,Y) helper ───────────────────────────────────────────
# (same function as qba-record-level-comparison.R)

posterior_u1 <- function(y_val, p_x, or_uy, r_x) {
  a    <- (1 - p_x) * (or_uy - 1)
  b    <- p_x * or_uy + (1 - p_x) - r_x * (or_uy - 1)
  cc   <- -r_x
  disc <- pmax(b^2 - 4 * a * cc, 0)

  q <- if_else(
    abs(a) < 1e-10,
    r_x,
    (-b + sqrt(disc)) / (2 * a)
  ) |> pmax(1e-6) |> pmin(1 - 1e-6)

  py_u1 <- (or_uy * q / (1 + (or_uy - 1) * q)) |> pmax(1e-6) |> pmin(1 - 1e-6)
  py_u0 <- q

  lik_u1 <- if_else(y_val == 1, py_u1, 1 - py_u1)
  lik_u0 <- if_else(y_val == 1, py_u0, 1 - py_u0)

  (lik_u1 * p_x) / (lik_u1 * p_x + lik_u0 * (1 - p_x))
}


# ── 4. Methods ────────────────────────────────────────────────────────────────

run_fox <- function(data, bias_p1, bias_p0, bias_or_uy) {
  r1 <- mean(data$Y[data$X == 1])
  r0 <- mean(data$Y[data$X == 0])

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

run_direct <- function(data, bias_p1, bias_p0, bias_or_uy) {
  r1 <- mean(data$Y[data$X == 1])
  r0 <- mean(data$Y[data$X == 0])

  data |>
    mutate(
      p_x   = if_else(X == 1, bias_p1, bias_p0),
      r_x   = if_else(X == 1, r1,      r0),
      p_u   = posterior_u1(Y, p_x, bias_or_uy, r_x),
      U_sim = rbinom(n(), 1, p_u)
    ) |>
    (\(d) glm(Y ~ X + W + Z + U_sim, data = d, family = binomial()))() |>
    coef() |> (\(b) exp(b["X"]))()
}

run_direct_stratified <- function(data, bias_p_u_df, bias_or_uy) {
  r_xw <- data |>
    group_by(X, W) |>
    summarise(r_xw = mean(Y), .groups = "drop")

  data |>
    left_join(bias_p_u_df, by = c("X", "W")) |>
    left_join(r_xw,         by = c("X", "W")) |>
    mutate(
      p_u_post = posterior_u1(Y, p_u, bias_or_uy, r_xw),
      U_sim    = rbinom(n(), 1, p_u_post)
    ) |>
    (\(d) glm(Y ~ X + W + Z + U_sim, data = d, family = binomial()))() |>
    coef() |> (\(b) exp(b["X"]))()
}


# ── 5. Run iterations ─────────────────────────────────────────────────────────

cat("Running", n_iter, "iterations...\n")

results <- map_dfr(seq_len(n_iter), \(i) {
  tibble(
    iteration    = i,
    fox          = run_fox(dat, marginal_p1, marginal_p0, true_or_uy),
    direct       = run_direct(dat, marginal_p1, marginal_p0, true_or_uy),
    direct_strat = run_direct_stratified(dat, true_p_u, true_or_uy)
  )
}, .progress = TRUE)


# ── 6. Summary ────────────────────────────────────────────────────────────────

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

cat("\n── Results ───────────────────────────────────────────────\n")
cat(sprintf("  True OR:                %.3f\n", true_or_xy))
cat(sprintf("  Naive OR  (ignores U):  %.3f\n", or_naive))
cat(sprintf("  Oracle OR (observes U): %.3f\n\n", or_oracle))
print(summary_tbl, n = Inf)


# ── 7. Diagnostic: show WHY the marginal correction fails ─────────────────────
#
# Plot the imputed P(U=1) values from one iteration, split by W stratum.
# The marginal correction assigns the same posterior to all subjects with
# the same (X, Y) cell, regardless of W. The stratified correction assigns
# different posteriors to W=0 and W=1 subjects in the same (X, Y) cell.

r1 <- mean(dat$Y[dat$X == 1])
r0 <- mean(dat$Y[dat$X == 0])
r_xw <- dat |> group_by(X, W) |> summarise(r_xw = mean(Y), .groups = "drop")

diag_dat <- dat |>
  # Marginal posteriors
  mutate(
    p_x_marg = if_else(X == 1, marginal_p1, marginal_p0),
    r_x_marg = if_else(X == 1, r1, r0),
    p_u_marg = posterior_u1(Y, p_x_marg, true_or_uy, r_x_marg)
  ) |>
  # Stratified posteriors
  left_join(true_p_u |> rename(p_x_strat = p_u), by = c("X", "W")) |>
  left_join(r_xw, by = c("X", "W")) |>
  mutate(p_u_strat = posterior_u1(Y, p_x_strat, true_or_uy, r_xw)) |>
  pivot_longer(c(p_u_marg, p_u_strat),
               names_to = "correction", values_to = "p_u_imputed") |>
  mutate(
    correction = recode(correction,
      p_u_marg  = "Marginal correction",
      p_u_strat = "Stratified correction"
    ),
    W_label = if_else(W == 1, "W = 1  (true P(U=1|X=1) = 0.92)",
                               "W = 0  (true P(U=1|X=1) = 0.08)"),
    X_label = if_else(X == 1, "Exposed (X=1)", "Unexposed (X=0)")
  )

ggplot(diag_dat, aes(x = p_u_imputed, fill = correction, colour = correction)) +
  geom_density(alpha = 0.35, linewidth = 0.6) +
  facet_grid(W_label ~ X_label) +
  scale_fill_manual(values  = c("Marginal correction" = "#E69F00",
                                "Stratified correction" = "#009E73")) +
  scale_colour_manual(values = c("Marginal correction" = "#E69F00",
                                 "Stratified correction" = "#009E73")) +
  labs(
    title    = "Diagnostic: imputed P(U=1) by W stratum and correction method",
    subtitle = paste0(
      "Marginal correction uses P(U=1|X=1)=0.50 for all X=1 subjects regardless of W.\n",
      "Stratified correction uses W-specific priors (0.92 for W=1, 0.08 for W=0)."
    ),
    x      = "Imputed P(U = 1)",
    y      = "Density",
    fill   = NULL,
    colour = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom",
        strip.text = element_text(face = "bold"))

ggsave("qba-differential-bias-diagnostic.png", width = 9, height = 7, dpi = 150)
cat("\nDiagnostic plot saved to qba-differential-bias-diagnostic.png\n")


# ── 8. Main results plot ──────────────────────────────────────────────────────

plot_dat <- results |>
  pivot_longer(c(fox, direct, direct_strat),
               names_to = "method", values_to = "or") |>
  mutate(method = recode(method,
    fox          = "Method 1: Fox et al.\n(marginal 2×2 → record level)",
    direct       = "Method 2: Direct simulation\n(marginal, per-record Bayes)",
    direct_strat = "Method 2b: Direct simulation\n(stratified on W)"
  ))

ggplot(plot_dat, aes(x = or, fill = method, colour = method)) +
  geom_density(alpha = 0.35, linewidth = 0.6) +
  geom_vline(xintercept = or_naive,    linetype = "dashed",  colour = "grey40",   linewidth = 0.8) +
  geom_vline(xintercept = or_oracle,   linetype = "dotdash", colour = "steelblue",linewidth = 0.8) +
  geom_vline(xintercept = true_or_xy,  linetype = "solid",   colour = "black",    linewidth = 0.8) +
  annotate("text", x = or_naive   + 0.01, y = Inf, vjust = 2,   hjust = 0,
           label = "Naive",    colour = "grey40",    size = 3.5) +
  annotate("text", x = or_oracle  + 0.01, y = Inf, vjust = 3.5, hjust = 0,
           label = "Oracle",   colour = "steelblue", size = 3.5) +
  annotate("text", x = true_or_xy + 0.01, y = Inf, vjust = 1,   hjust = 0,
           label = "True OR",  colour = "black",     size = 3.5) +
  scale_fill_manual(values  = c("#E69F00", "#56B4E9", "#009E73")) +
  scale_colour_manual(values = c("#E69F00", "#56B4E9", "#009E73")) +
  labs(
    title    = "Bias-adjusted OR: marginal vs. stratified QBA under extreme differential bias",
    subtitle = sprintf(
      "N = %d, %d iterations  |  True OR = %.1f  |  P(U=1|X=1) = 0.92 (W=1) vs 0.08 (W=0)",
      N, n_iter, true_or_xy
    ),
    x      = "Bias-adjusted OR for X (adjusted for W, Z, U)",
    y      = "Density",
    fill   = "Method",
    colour = "Method"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom",
        legend.box      = "vertical")

ggsave("qba-differential-bias-results.png", width = 10, height = 6, dpi = 150)
cat("Results plot saved to qba-differential-bias-results.png\n")
