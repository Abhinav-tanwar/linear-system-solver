# Linear System Solver

Tolerance-aware row reduction, inverse validation, and linear-system residuals.

## Run

Requires R 4.4+ and packages MASS and jsonlite.

```r
install.packages(c("MASS", "jsonlite"))
```

```sh
Rscript scripts/run.R
Rscript tests/checks.R
```

Generated metrics and figures are saved in `results/`. Read [input provenance](data/README.md) and [method changes](docs/provenance.md).

## Executed results

Numerical checks verify residuals rather than relying on exact floating-point equality. Full measured diagnostics are in [metrics.json](results/metrics.json).

## Reproducibility and scope

Seeds are fixed in source. Validation uses numerical invariants or input/split checks. This is a maintained reconstruction of the supplied work; methodological corrections and newly generated data mean some results intentionally differ from the original reports. Original PDFs are omitted to avoid publishing student identifiers.

Additional source coverage includes matrix structure and transpose identities. Spectral diagnostics, where applicable, distinguish image/null bases, rectangular rank–nullity, outer-product rank, and symmetric orthogonality.

## Reports

`report.html` is a ready-to-open summary. To render the R Markdown version after running the analysis, install `rmarkdown` and `knitr`, install Pandoc, and run `rmarkdown::render("report.Rmd")` from the repository root. The optional R Markdown rendering was not executed in the local environment because Pandoc is unavailable.
