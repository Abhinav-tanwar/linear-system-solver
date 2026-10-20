source("R/analysis.R")

out <- tempfile()

dir.create(out)

r <- run_analysis(normalizePath("."), out)

stopifnot(length(r) > 0, file.exists(out))

stopifnot(r$solve_residual < 1e-08, r$inverse_residual < 1e-08)

stopifnot(max(abs(rref(matrix(c(1, 2, 2, 4), 2)) - matrix(c(1, 0, 2, 0), 2))) < 1e-08)

cat("Validation passed\n")

stopifnot(all(is.finite(unlist(r))))

