# pmm2_common.R - Spilni utylity dlia vsikh PMM2 modelei


#' Universalnyi alhorytm PMM2 dlia vsikh typiv modelei
#'
#' @param b_init Pochatkovi otsinky parametriv
#' @param X Matrytsia dyzainu
#' @param y Vektor vidhuku
#' @param m2,m3,m4 Tsentralni momenty
#' @param max_iter Maksymalna kilkist iteratsii
#' @param tol Dopusk dlia zbizhnosti
#' @param regularize Chy dodavaty rehuliaryzatsiiu
#' @param reg_lambda Parametr rehuliaryzatsii
#' @param verbose Chy vyvodyty informatsiiu pro prohres
#' @param poly_terms Poperedno obchysleni koefitsiienty polinoma (spysok z elementamy \code{A}, \code{B}, \code{C});
#'   dozvoliaie peredavaty vlasni znachennia dlia spetsialnykh stsenariiv, inakshe vony obchysliuiutsia z momentiv
#' 
#' @return Spysok z rezultatamy otsiniuvannia
#' @keywords internal
pmm2_algorithm <- function(b_init, X, y, m2, m3, m4,
                           max_iter = 50, tol = 1e-6,
                           regularize = TRUE, reg_lambda = 1e-8,
                           verbose = FALSE,
                           poly_terms = NULL) {
  # Potochni otsinky parametriv
  b_cur <- b_init
  converged <- FALSE
  iterations <- 0

  # Obchyslyty koefitsiienty PMM2 polinoma
  # Vyvedennia vidpovidaie rivnianniu (10) zi statti:
  # A = c3 * sigma^3, B = (c4 + 3) * sigma^4 - sigma^4 - 2 c3 * sigma^3 * y,
  # C = c3 * y^2 * sigma^3 - ((c4 + 3) * sigma^4 - sigma^4) * y - sigma^2 * c3 * sigma^3.
  # Pidstanovka sigma^2 = m2, m3 = c3 sigma^3, m4 = (c4 + 3) sigma^4 daie navedeni nyzhche formuly.
  if (is.null(poly_terms)) {
    A <- m3
    B <- m4 - m2^2 - 2 * m3 * y
    C <- m3 * y^2 - (m4 - m2^2) * y - m2 * m3
  } else {
    A <- poly_terms$A
    B <- poly_terms$B
    C <- poly_terms$C
    if (length(B) != nrow(X) || length(C) != nrow(X)) {
      stop("poly_terms$B ta poly_terms$C povynni maty dovzhynu, shcho dorivniuie kilkosti riadkiv X")
    }
  }

  # Vidstezhuvaty istoriiu zbizhnosti, iakshcho verbose
  if (verbose) {
    conv_history <- numeric(max_iter)
  }

  # Osnovnyi tsykl iteratsii
  for (iter in seq_len(max_iter)) {
    iterations <- iter

    # Obchyslyty prohnozovani znachennia
    y_pred <- as.vector(X %*% b_cur)

    # Obchyslyty Z1 = A*y_pred^2 + B*y_pred + C
    Z1 <- A*(y_pred^2) + B*y_pred + C

    # Sformuvaty vektor Z dlia kozhnoho parametra
    p <- length(b_cur)
    Z <- numeric(p)
    for (r in 1:p) {
      Z[r] <- sum(Z1 * X[, r])
    }

    # Obchyslyty pokhidnu JZ11 = 2*A*y_pred + B
    JZ11 <- 2 * A * y_pred + B

    # Sformuvaty matrytsiu Yakobiana
    JZs <- matrix(0, p, p)
    for (i in 1:p) {
      for (j in 1:p) {
        JZs[i, j] <- sum(JZ11 * X[, i] * X[, j])
      }
    }

    # Dodaty rehuliaryzatsiiu, iakshcho potribno
    if (regularize) {
      diag(JZs) <- diag(JZs) + reg_lambda
    }

    # Rozv'iazaty systemu JZs * delta = Z
    delta <- tryCatch({
      solve(JZs, Z)
    }, error = function(e) {
      if (verbose) {
        cat("Pomylka pry rozv'iazanni liniinoi systemy:", conditionMessage(e), "\n")
        cat("Dodaiemo sylnishu rehuliaryzatsiiu\n")
      }
      diag(JZs) <- diag(JZs) + 1e-4
      solve(JZs, Z)
    })

    # Onovyty parametry
    b_new <- b_cur - delta
    diff_par <- sqrt(sum((b_new - b_cur)^2))

    # Zberehty istoriiu zbizhnosti, iakshcho verbose
    if (verbose) {
      conv_history[iter] <- diff_par
      if (iter %% 5 == 0 || iter == 1) {
        cat("Iteratsiia", iter, ": Zmina parametriv =",
            formatC(diff_par, digits = 8), "\n")
      }
    }

    b_cur <- b_new

    # Pereviryty zbizhnist
    if (diff_par < tol) {
      converged <- TRUE
      if (verbose) cat("Zbizhnist dosiahnuta pislia", iter, "iteratsii\n")
      break
    }
  }

  # Poperedzhennia, iakshcho dosiahnuto maksymalnu kilkist iteratsii bez zbizhnosti
  if (!converged && verbose) {
    cat("Poperedzhennia: Alhorytm ne zbihsia pislia", max_iter, "iteratsii\n")
  }

  # Obchyslyty kintsevi zalyshky
  final_res <- as.numeric(y - X %*% b_cur)

  # Pobuduvaty istoriiu zbizhnosti, iakshcho verbose
  if (verbose && iterations > 1) {
    if (requireNamespace("graphics", quietly = TRUE)) {
      graphics::plot(1:iterations, conv_history[1:iterations], type = "b",
                     xlab = "Iteratsiia", ylab = "Zmina parametriv",
                     main = "Istoriia zbizhnosti")
      graphics::abline(h = tol, col = "red", lty = 2)
    }
  }

  # Povernuty rezultaty
  list(
    b = as.numeric(b_cur),
    convergence = converged,
    iterations = iterations,
    residuals = final_res
  )
}
