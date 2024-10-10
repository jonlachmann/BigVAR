test_that("multiplication works", {
  data(Y)
  structs <- c("Basic", "BasicEN", "Lag", "SparseLag", "OwnOther", "SparseOO", "HLAGC", "HLAGOO", "HLAGELEM", "Tapered", "BGR", "MCP", "SCAD")
  structs_x <- c("Basic", "BasicEN", "Lag", "SparseLag", "OwnOther", "SparseOO", "EFX", "HLAGC", "HLAGOO", "HLAGELEM", "Tapered","BGR", "MCP", "SCAD")
  for (struct in structs) {
    test_that(paste0("BigVAR.fit with struct ", struct, " produces expected results."), {
      fit <- BigVAR.fit(Y, struct = struct, p = 2, lambda = 1)
      expect_snapshot_value(fit, style = "deparse")
    })
  }

  for (struct in structs_x) {
    test_that(paste0("BigVAR.fit with struct ", struct, " produces expected results with exogenous data."), {
      VARX <- list(k = 2, s = 2)
      fit  <- BigVAR.fit(Y, p,"Basic", lambda = 1e-2, VARX = VARX, intercept = TRUE)
      expect_snapshot_value(fit, style = "deparse")
    })
  }
})



