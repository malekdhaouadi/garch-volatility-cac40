# GARCH / EGARCH volatility modelling: CAC 40 and EURO STOXX 50

GARCH(1,1) and EGARCH(1,1) with Student-t innovations (`rugarch`, R): in-sample fit and diagnostics, rolling out-of-sample evaluation (QLIKE, Diebold-Mariano), a long-sample comparison of 21-day volatility forecasts with VSTOXX implied volatility (2008-2026), and a 1-day VaR / Expected Shortfall backtest with tail-model experiments.

![Conditional volatility](figures/cond_vol_garch_egarch.png)

## Summary

- In sample, EGARCH is preferred (AIC -6.504 vs -6.458); the leverage term is negative and significant (α₁ = -0.169) and residual diagnostics are clean.
- 1-day out-of-sample (500 days): QLIKE 1.427 (EGARCH) vs 1.438 (GARCH), Diebold-Mariano p = 0.67, so no significant difference.
- Against VSTOXX (21-day horizon, 2008-2026): variance risk premium of +3.2 (EURO STOXX 50, from 2011) to +3.9 (CAC 40 proxy, from 2008) vol points. VSTOXX tracks realized volatility better than either model (p <= 0.026); GARCH and EGARCH are indistinguishable once mean levels are removed.
- VaR 2008-2026: Student-t GARCH and EGARCH under-cover at 1% (exceedance rates 1.48% and 1.81%, Kupiec p <= 0.002) but beat 250-day historical simulation on quantile loss.
- Skewed-t and empirical-tail variants bring the 1% exceedance rate to 1.04-1.15% on a hold-out (2017-2026) that was not used for selection (baseline: 1.55%). Exceedances of the selected GARCH remain predictable (DQ p = 0.014).

## Data

- CAC 40 (`^FCHI`) adjusted close, Yahoo Finance (`quantmod`). In-sample and 1-day evaluation: 2015-01-05 to 2026-10-26 (N = 3,023). Long-sample VSTOXX comparison: from 2000-01-04 (N = 6,853). Missing prices dropped, not interpolated.
- EURO STOXX 50 (`^STOXX50E`), Yahoo, from 2007-04-02 (N = 4,906).
- VSTOXX (V2TX), daily closes 1999-2026, downloaded manually from stoxx.com to `data/v2tx_history.csv`.

CAC 40 returns (2015-2026): skewness -0.83, excess kurtosis 11.19, Jarque-Bera p < 2.2e-16. Squared returns: Ljung-Box Q(10) = 1144.9, ARCH-LM(5) = 478.5, both p < 2.2e-16.

## Models

Mean: constant. Innovations: Student-t (the VaR experiments below also use a skewed Student-t and zero-mean variants). Estimation: ML, robust standard errors.

- GARCH: σ²ₜ = ω + α·ε²ₜ₋₁ + β·σ²ₜ₋₁
- EGARCH: ln σ²ₜ = ω + β·ln σ²ₜ₋₁ + α₁·zₜ₋₁ + γ₁·(|zₜ₋₁| − E|z|)

`rugarch` convention: `alpha1` is the sign (leverage) term, `gamma1` the magnitude term.

## Estimates (CAC 40, 2015-2026)

| | GARCH | EGARCH |
|---|---|---|
| ω | ≈0 (p = 0.52) | -0.2927 (0.0041) |
| α / α₁ | 0.1410 (0.0200) | -0.1686 (0.0119) |
| β | 0.8316 (0.0540) | 0.9692 |
| γ₁ | n/a | 0.1322 (0.0091) |
| ν (Student-t) | 5.19 (0.70) | 5.67 (0.54) |
| Persistence / half-life | 0.9725 / 24.9 d | 0.9692 / 22.2 d |
| Log-likelihood | 9766.1 | 9836.2 |
| AIC / BIC | -6.4579 / -6.4479 | -6.5036 / -6.4916 |

EGARCH is preferred on AIC, BIC and likelihood. α₁ < 0 and significant: leverage effect (news impact: z = -1 raises ln σ² by about 0.30, z = +1 lowers it by about 0.04). The reported EGARCH standard error on β is not meaningful (ω and β nearly collinear).

![News impact curve](figures/news_impact_egarch.png)

## Residual diagnostics (p-values)

| Test | GARCH | EGARCH |
|---|---|---|
| Ljung-Box z, lag 10 / 20 | 0.985 / 0.876 | 0.957 / 0.786 |
| Ljung-Box z², lag 10 / 20 | 0.515 / 0.847 | 0.916 / 0.993 |
| ARCH-LM, q = 5 / 10 | 0.255 / 0.577 | 0.653 / 0.931 |

EGARCH: sign bias joint p = 0.46; Nyblom joint statistic 1.25 (10% critical value 1.49). Adjusted Pearson goodness-of-fit rejects the Student-t (p < 0.01). `arch_lm_test` in `R/utils.R` is validated against `FinTS::ArchTest`.

## Out-of-sample, 1-day ahead (CAC 40, last 500 days)

1,000-day moving window, refit every 50 days. Proxy: squared demeaned return.

| Loss | GARCH | EGARCH | DM p-value |
|---|---|---|---|
| QLIKE | 1.4380 | 1.4274 | 0.67 |
| MSE | 4.252e-08 | 4.020e-08 | 0.016 |

EGARCH is not significantly better under QLIKE; the MSE difference is driven by a few extreme days.

## Implied volatility: 21-day forecasts vs VSTOXX

VSTOXX is a 30-calendar-day implied volatility, so model forecasts are averaged over 21 trading days. At each origin, models are re-estimated every 50 days on a 1,000-day moving window. Realized volatility = realized over the next 21 trading days. Vol points, annualized.

Two runs:

- **EURO STOXX 50 vs VSTOXX** (like for like): first origin 2011-03-23 (Yahoo history starts 2007-04, and the window needs 1,000 days), 3,886 origins.
- **CAC 40 vs VSTOXX** (proxy, different index): first origin 2008-01-02, 4,769 origins; the only run covering 2008-09. Daily returns of the CAC 40 and EURO STOXX 50 correlate at 0.969 over 2007-2026 (0.976 in 2008-09, 0.987 in 2020).

| | EURO STOXX 50 | | | CAC 40 (proxy) | | |
|---|---|---|---|---|---|---|
| | GARCH | EGARCH | VSTOXX | GARCH | EGARCH | VSTOXX |
| Correlation with VSTOXX | 0.852 | 0.795 | n/a | 0.879 | 0.795 | n/a |
| Mean error vs realized | +2.47 | +0.05 | +3.23 | +2.01 | -0.77 | +3.93 |
| RMSE vs realized | 7.95 | 7.42 | 7.66 | 8.63 | 8.65 | 8.70 |
| Correlation with realized | 0.552 | 0.559 | 0.646 | 0.617 | 0.570 | 0.689 |
| RMSE, mean removed | 7.55 | 7.42 | 6.94 | 8.39 | 8.62 | 7.76 |

**Variance risk premium** (VSTOXX − subsequent realized): +3.23 (Newey-West s.e. 0.40, lag 20; positive on 78% of origins) on the EURO STOXX 50 run, +3.93 (s.e. 0.42; 81%) on the CAC 40 run.

**Diebold-Mariano on squared error vs realized, HAC lag 20** (stat > 0: second model better):

| Comparison | EURO STOXX 50 | CAC 40 |
|---|---|---|
| GARCH vs EGARCH, raw | +2.26 (p = 0.024) | -0.07 (p = 0.945) |
| GARCH vs EGARCH, mean removed | +0.75 (p = 0.454) | -0.75 (p = 0.456) |
| EGARCH vs VSTOXX, mean removed | +2.23 (p = 0.026) | +2.95 (p = 0.003) |
| GARCH vs VSTOXX, mean removed | +3.03 (p = 0.002) | +3.32 (p = 0.001) |

- VSTOXX tracks changes in realized volatility better than either model; GARCH and EGARCH are indistinguishable once mean levels are removed.
- EGARCH's lower raw RMSE on the EURO STOXX 50 comes from lower bias (mean error +0.05 vs +2.47), not better tracking; on the CAC 40 the raw difference is not distinguishable from zero (p = 0.945).

**By regime** (CAC 40 vs VSTOXX; vol points):

| Period | n | Premium | RMSE GARCH | RMSE EGARCH | RMSE VSTOXX |
|---|---|---|---|---|---|
| 2008-09 (GFC) | 508 | +3.5 | 12.92 | 15.28 | 13.26 |
| 2010-12 | 766 | +4.4 | 8.39 | 7.68 | 8.39 |
| 2013-19 | 1,777 | +3.3 | 6.30 | 5.71 | 6.13 |
| 2020 (COVID) | 255 | +0.9 | 17.13 | 18.45 | 16.54 |
| 2021-26 | 1,463 | +5.1 | 6.96 | 5.76 | 7.45 |

EGARCH has the lowest RMSE in three of the five regimes and the highest in the two crises. The premium is smallest in 2020. The EURO STOXX 50 run shows the same pattern from 2011 onward.

![Forecasts vs VSTOXX, CAC 40](figures/ext_iv_vs_garch_CAC40.png)
![Forecasts vs VSTOXX, EURO STOXX 50](figures/ext_iv_vs_garch_STOXX50E.png)

## VaR and Expected Shortfall backtest (CAC 40)

1-day-ahead VaR and ES, 4,812 out-of-sample days (2008-01-02 to 2026-10-26). GARCH and EGARCH (Student-t, 1,000-day moving window, refit every 50 days) against a 250-day historical simulation (HS) benchmark. VaR/ES for the Student-t models are computed analytically from the forecast mean, volatility and estimated ν. Tests: Kupiec (unconditional coverage), Christoffersen (independence, conditional coverage), McNeil-Frey ES test (mean of standardized exceedance residuals; one-sided p for ES too optimistic), tick (quantile) loss with Diebold-Mariano.

**1% level** (expected exceedances: 48.1)

| | GARCH | EGARCH | HS250 |
|---|---|---|---|
| Exceedances (rate) | 71 (1.48%) | 87 (1.81%) | 73 (1.52%) |
| Kupiec p | 0.002 | 4e-7 | 8e-4 |
| Independence p | 0.40 | 0.097 | 8e-4 |
| Conditional coverage p | 0.006 | 7e-7 | 1.4e-5 |
| ES test: mean residual (one-sided p) | -0.011 (0.46) | -0.146 (0.087) | -0.456 (0.002) |
| Tick loss (x1e4) | 4.26 | 4.44 | 5.17 |

**5% level** (expected exceedances: 240.6)

| | GARCH | EGARCH | HS250 |
|---|---|---|---|
| Exceedances (rate) | 317 (6.59%) | 289 (6.01%) | 259 (5.38%) |
| Kupiec p | 1.4e-6 | 0.002 | 0.23 |
| Independence p | 0.002 | 0.17 | 2.5e-9 |
| Conditional coverage p | 8e-8 | 0.003 | 9e-9 |
| ES test: mean residual (one-sided p) | -0.060 (0.071) | -0.190 (1e-4) | -0.188 (0.003) |
| Tick loss (x1e4) | 14.6 | 14.5 | 16.1 |

- All three methods produce too many exceedances at 1% (Kupiec rejected); EGARCH is worst (1.81%).
- Both GARCH models have lower tick loss than HS at both levels (Diebold-Mariano p <= 0.001). GARCH vs EGARCH is not significant (p = 0.125 at 1%, 0.452 at 5%).
- Exceedance clustering: independence is rejected for HS at both levels, for GARCH at 5% only, and for EGARCH at neither.
- ES: the GARCH ES is not rejected at 1% and only marginally at 5% (p = 0.07); EGARCH ES is rejected at 5%; HS underestimates tail losses at both levels.
- By year (1% VaR): 2008 and 2020 account for 26 of EGARCH's 87 exceedances (15 and 11, Basel red zone is 10+ per 250 days); GARCH has 5 and 8. Outside these two years the rate is still above 1% (GARCH 1.35%, EGARCH 1.42%).

![1% VaR exceedances](figures/var_exceedances_CAC40.png)
![Cumulative exceedances](figures/var_cum_exceedances_CAC40.png)

## VaR tail experiments (CAC 40)

The baseline above under-covers: 1.55% of days in 2017-2026 breach the 1% VaR. Sixteen specifications were tried: GARCH/EGARCH x constant/zero mean x Student-t/skewed Student-t x 1,000/2,000-day window, each with two tail estimates: parametric (quantile of the fitted distribution) and filtered historical simulation (FHS: empirical quantile of the window's standardized residuals, rescaled by the forecast volatility). That gives 32 candidates plus the 250-day historical simulation (HS250), all on the same 4,812 days.

Design fixed before running: candidates are ranked on 2008-2016 only (mean |log(exceedance rate / alpha)| at 1% and 5%); results are reported on 2017-2026 (2,512 days), which is not used for selection. A dynamic quantile (DQ) test is added, and a 2.5% level. The two top candidates tied on the selection score (0.0574); the tie was broken by listing order, so the selected model is the GARCH below, and the tied EGARCH is reported alongside it. The baseline ranks 30th of 33 on the selection period.

Hold-out 2017-2026, values at the 1% / 2.5% / 5% levels:

| | Rate (%) | Kupiec p | DQ p | ES test p | Tick loss (x1e4) |
|---|---|---|---|---|---|
| Baseline: GARCH, const. mean, Student-t, 1,000 d | 1.55 / 3.90 / 6.53 | 0.010 / 3e-5 / 8e-4 | 5e-4 / 2e-6 / 7e-9 | 0.55 / 0.65 / 0.13 | 3.86 / 7.68 / 12.5 |
| Selected: GARCH, zero mean, skewed-t, 1,000 d | 1.04 / 2.55 / 5.10 | 0.86 / 0.88 / 0.83 | 0.014 / 0.034 / 8e-5 | 0.70 / 0.56 / 0.53 | 3.86 / 7.63 / 12.5 |
| EGARCH, const. mean, skewed-t, 2,000 d | 1.15 / 2.39 / 4.82 | 0.45 / 0.72 / 0.67 | 0.41 / 0.79 / 0.028 | 0.49 / 0.21 / 0.28 | 3.71 / 7.42 / 12.1 |
| EGARCH, zero mean, skewed-t, 2,000 d | 1.11 / 2.35 / 4.62 | 0.57 / 0.62 / 0.37 | 0.44 / 0.81 / 0.061 | 0.49 / 0.25 / 0.24 | 3.71 / 7.42 / 12.1 |
| HS250 | 1.59 / 2.83 / 5.45 | 0.006 / 0.30 / 0.30 | <1e-14 | 0.042 / 0.011 / 0.035 | 4.64 / 8.48 / 13.6 |

Expected exceedances on the hold-out: 25.1 / 62.8 / 125.6.

- The skewed-t specifications bring exceedance rates close to nominal (1.04-1.15% at the 1% level, against 1.55% for the baseline). Decomposition for GARCH (1,000-day window, 1% level, hold-out): the baseline has 39 exceedances; zero mean alone gives 34, skewed-t alone 30, FHS (empirical residual tail) alone 26, zero mean and skewed-t together 26. Each step changes only a few exceedances (sampling noise about +/-5), so only the baseline against the combined fixes is a visible gap. None of these GARCH variants passes the DQ test at 1% (p <= 0.014), whatever the mean, distribution or tail method; tick loss barely differs (3.84 to 3.91). All 33 candidates are in `outputs/var_experiments_all.csv`.
- Coverage is not independence: the selected GARCH matches the average rate, but its exceedances remain predictable (independence p = 0.029 at 1%, DQ p = 0.014 and 8e-5 at 1% and 5%). The EGARCH 2,000-day candidates pass the DQ test at 1% and 2.5% and are marginal at 5%. Among constant-mean skewed-t candidates at 1% (hold-out), GARCH with a 1,000-day window fails the DQ test (p = 6.6e-4) while EGARCH 1,000 d (p = 0.27), GARCH 2,000 d (p = 0.45) and EGARCH 2,000 d (p = 0.41) pass, with 26 to 31 exceedances each: either the asymmetry or the longer window is enough, and the two effects cannot be separated here.
- All skewed-t candidates pass the ES test; HS250 is rejected at all three levels.
- Tick loss at 1% is 3.71 for the EGARCH 2,000-day candidates against 3.86 for the baseline and the selected model (not tested for significance).
- Exceedances of the 1% VaR in 2008 and 2020 for the EGARCH skewed-t, 2,000-day candidate: 2 and 7 (the Student-t EGARCH of the previous section: 15 and 11).

![Cumulative exceedances, selected vs baseline](figures/var_experiments_cum_exceedances.png)

## Limitations

- 21-day windows overlap, so the effective sample is about n/21 (roughly 190 to 230 independent observations); HAC variances account for this approximately. No multiple-testing adjustment is applied to any reported p-value.
- Mean-removed losses use full-sample means: a diagnostic of tracking, not a tradable forecast.
- The 2008-09 result relies on the CAC 40 as a proxy for the EURO STOXX 50; each regime contains one or few episodes.
- Implied volatility is risk-neutral, GARCH is physical-measure.
- Squared returns are a noisy variance proxy (no intraday realized variance). The 1-day evaluation covers only 500 days.
- Student-t rejected by goodness-of-fit; the skewed-t is used only in the VaR experiments, not in the in-sample fit or the volatility comparison.
- VaR/ES: about 48 expected exceedances at 1%, so backtests have limited power; ES values assume the fitted distribution.
- VaR experiments: 33 candidates were evaluated; the selection tie was broken by listing order; the hold-out has about 25 expected exceedances at 1% (26 observed is consistent with true rates of roughly 0.7% to 1.5%).

## Reproduce

```r
renv::restore()
# place data/v2tx_history.csv (stoxx.com, V2TX "History EUR Price")
# config.yml: data_start: "2000-01-01", eval_start: "2008-01-01"
for (f in c("01_data", "02_diagnostics", "03_fits", "04_forecast", "06_extended_eval", "07_var_backtest", "08_var_experiments"))
  source(sprintf("R/%s.R", f))
```

Parameters (window, refit frequency, seed, distribution) are in `config.yml`.