\# GARCH volatility modelling: CAC 40



Volatility forecasting on daily CAC 40 log returns (Jan 2015 to Oct 2026, 3,023 observations)

using GARCH(1,1) and EGARCH(1,1) with Student-t innovations, implemented in R with `rugarch`.



\## Key findings

\- Returns show strong volatility clustering (Ljung-Box on squared returns, p < 2.2e-16)

&#x20; and fat tails (excess kurtosis 11.2, skewness -0.83).

\- EGARCH is preferred in-sample (AIC -6.504 vs -6.458; log-likelihood 9836 vs 9766).

&#x20; The leverage term is negative and significant: falls raise volatility more than rises.

\- Standardized residuals show no remaining autocorrelation, ARCH effects or sign bias.

\- Out-of-sample (500 rolling one-day forecasts), EGARCH has a slightly lower QLIKE

&#x20; (1.427 vs 1.438), but the difference is not statistically significant

&#x20; (Diebold-Mariano p = 0.67).



\## Pipeline

| Script | Purpose |

|---|---|

| `R/01\_data.R` | Download prices, compute log returns |

| `R/02\_diagnostics.R` | Summary stats, Jarque-Bera, Ljung-Box, ARCH-LM |

| `R/03\_fits.R` | Fit models, parameter tables, residual diagnostics |

| `R/04\_forecast.R` | Rolling out-of-sample evaluation, 30-day forecast |



Parameters (window, refit frequency, horizon) are set in `config.yml`.

Dependencies are locked with `renv`.



\## Reproduce

&#x20;   renv::restore()

&#x20;   source("R/01\_data.R"); source("R/02\_diagnostics.R")

&#x20;   source("R/03\_fits.R"); source("R/04\_forecast.R")



\## Status

Comparison with implied volatility (VCAC / VSTOXX) is in progress.

