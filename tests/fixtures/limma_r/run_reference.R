# Sorties de référence R pour tests/test_limma_reference.py (voir make_fixtures.py).
# (fonction sourcée depuis le code DEqMS 1.31.0, sans installer le paquet).
suppressMessages(library(limma))
deqms_src <- Sys.getenv("DEQMS_R")
env <- new.env(); sys.source(deqms_src, envir = env)
spectraCounteBayes <- env$spectraCounteBayes
for (d in commandArgs(TRUE)) {
  E <- as.matrix(read.csv(gzfile(file.path(d, "expr.csv.gz")), header = FALSE))
  X <- as.matrix(read.csv(file.path(d, "design.csv"), header = FALSE))
  C <- as.matrix(read.csv(file.path(d, "contrasts.csv"), header = FALSE))
  colnames(X) <- paste0("c", seq_len(ncol(X)))
  rownames(C) <- colnames(X); colnames(C) <- paste0("k", seq_len(ncol(C)))
  rownames(E) <- paste0("p", seq_len(nrow(E)))
  fit <- lmFit(E, X)
  fit2 <- contrasts.fit(fit, C)
  fit3 <- suppressWarnings(eBayes(fit2))
  o <- file.path(d, "R")
  dir.create(o, showWarnings = FALSE)
  w <- function(x, f) {
    con <- gzfile(file.path(o, paste0(f, ".gz")), "w")
    write.table(x, con, sep = ",", row.names = FALSE, col.names = FALSE)
    close(con)
  }
  w(fit3$coefficients, "coef.csv"); w(fit3$sigma, "sigma.csv")
  w(fit3$t, "t.csv"); w(fit3$p.value, "p.csv")
  w(apply(fit3$p.value, 2, p.adjust, method = "BH"), "padj.csv")
  w(c(fit3$df.prior[1], fit3$s2.prior[1], fit3$df.total[1], fit3$df.residual[1]), "hyper.csv")
  w(p.adjust(as.vector(fit3$p.value), method = "BH"), "padj_global.csv")
  if (file.exists(file.path(d, "counts.csv"))) {
    cnt <- scan(file.path(d, "counts.csv"), quiet = TRUE)
    fit3$count <- cnt
    sca <- tryCatch(spectraCounteBayes(fit3, coef_col = seq_len(ncol(C))),
                    error = function(e) { cat(basename(d), "DEqMS ERROR:", conditionMessage(e), "\n"); NULL })
    if (!is.null(sca)) {
    w(sca$sca.t, "sca_t.csv"); w(sca$sca.p, "sca_p.csv")
    w(apply(sca$sca.p, 2, p.adjust, method = "BH"), "sca_padj.csv")
    w(sca$sca.dfprior, "sca_dfprior.csv"); w(fitted(sca$model), "sca_loess.csv")
    }
  }
  cat(basename(d), "OK df.prior =", fit3$df.prior[1], "\n")
}
