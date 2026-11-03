# GARCH / EGARCH volatility modelling: CAC 40 and EURO STOXX 50

GARCH(1,1) and EGARCH(1,1) with Student-t innovations (`rugarch`, R): in-sample fit and diagnostics, rolling out-of-sample evaluation (QLIKE, Diebold-Mariano), and a long-sample comparison of 21-day volatility forecasts with VSTOXX implied volatility (2008-2026).

![Conditional volatility](figures/cond_vol_garch_egarch.png)

## Data

- CAC 40 (`^FCHI`) adjusted close, Yahoo Finance (`quantmod`). In-sample and 1-day evaluation: 2015-01-05 to 2026-10-26 (N = 3,023). Long-sample VSTOXX comparison: from 2000-01-04 (N = 6,853). Missing prices dropped, not interpolated.
- EURO STOXX 50 (`^STOXX50E`), Yahoo, from 2007-04-02 (N = 4,906).
- VSTOXX (V2TX), daily closes 1999-2026, downloaded manually from stoxx.com to `data/v2tx_history.csv`.

CAC 40 returns (2015-2026): skewness -0.83, excess kurtosis 11.19, Jarque-Bera p < 2.2e-16. Squared returns: Ljung-Box Q(10) = 1144.9, ARCH-LM(5) = 478.5, both p < 2.2e-16.

## Models

Mean: constant. Innovations: Student-t. Estimation: ML, robust standard errors.

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

## Limitations

- 21-day windows overlap, so the effective sample is about n/21 (roughly 190 to 230 independent observations); HAC variances account for this approximately. About eight p-values are reported without multiple-testing adjustment.
- Mean-removed losses use full-sample means: a diagnostic of tracking, not a tradable forecast.
- The 2008-09 result relies on the CAC 40 as a proxy for the EURO STOXX 50; each regime contains one or few episodes.
- Implied volatility is risk-neutral, GARCH is physical-measure.
- Squared returns are a noisy variance proxy (no intraday realized variance). The 1-day evaluation covers only 500 days.
- Student-t rejected by goodness-of-fit; skewed-t not tested.

## Reproduce

```r
renv::restore()
# place data/v2tx_history.csv (stoxx.com, V2TX "History EUR Price")
# config.yml: data_start: "2000-01-01", eval_start: "2008-01-01"
for (f in c("01_data", "02_diagnostics", "03_fits", "04_forecast", "06_extended_eval"))
  source(sprintf("R/%s.R", f))
```

Parameters (window, refit frequency, seed, distribution) are in `config.yml`.