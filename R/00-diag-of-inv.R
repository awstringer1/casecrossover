# Get the diagonal of Q^-1 via inverting sparse cholesky decomp
diagOfInv = function(x, verbose=FALSE,constrA = NULL,i = NULL) {
  # x: matrix to compute the diagonal of the inverse of
  # verbose: print progress?
  # constrA: optional n x k matrix of linear constraints to correct for. n = dim(x), k = number of constraints
  # i: optional vector of indices of x for which to return the diagonal. only used if correcting
  # for linear constraints. The whole matrix x is needed even if you only want a small subset i
  # of marginal variances. However, correcting for linear constraints needs only the subset i,
  # which if small can drastically decrease run time.

  if (is.null(i)) i <- 1:dim(x)[1]

  if(verbose) {
    cat("cholesky\n")
  }
  cholHere = Matrix::expand(Matrix::Cholesky(x,
                                             LDL=FALSE,
                                             super=FALSE))
  if(verbose) {
    cat("solve\n")
  }
  cholHere$Linv = Matrix::solve(cholHere$L)

  if(verbose) {
    cat("multiply\n")
  }
  # sum the columns of Linv^2
  cholHere$LinvDf = data.table::data.table(
    col = rep(1:nrow(cholHere$Linv), diff(cholHere$Linv@p)),
    x = cholHere$Linv@x^2
  )

  # varDiag = cholHere$LinvDf[, .(sum = sum(x)), by = col]
  # Change the use of "." because R package doesn't like nonstandard eval:
  varDiag = cholHere$LinvDf[, list(sum = sum(x)), by = col]


  if(verbose) {
    cat("permute\n")
  }

  # do the permutation transform
  varDiagMat = Diagonal(nrow(varDiag), varDiag$sum)
  varDiagMatP = crossprod(cholHere$P, varDiagMat) %*% cholHere$P

  if (is.null(constrA)) {
    vars <- varDiagMatP@x
  } else {
    # Correct for linear constraints
    WW <- Matrix::solve(x,constrA)
    # VV <- Matrix::solve(t(constrA) %*% WW,WW)
    VV <- t(Matrix::solve(t(WW) %*% constrA,t(WW)))
    # Correct for the ones that are actually being returned
    correction <- rep(0,length(varDiagMatP@x))
    for (j in i) correction[j] <- VV[j, ] %*% WW[j, ]

    vars <- varDiagMatP@x - correction
  }
  vars[i]
}
