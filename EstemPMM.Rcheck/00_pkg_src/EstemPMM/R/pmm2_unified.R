#' Unified PMM2 Estimator for Nonlinear Regression Models
#'
#' Цей файл містить уніфіковану реалізацію методу PMM2, придатну для будь-яких
#' нелінійних регресійних моделей (включаючи SARIMAX), для яких можна обчислити
#' залишки та їх похідні (Якобіан) по параметрах.
#'
#' Реалізовано два підходи:
#' 1. Iterative: Повна ітеративна процедура (Nonlinear PMM2).
#' 2. One-step (Global): Однокрокова корекція після класичного оцінювача (наприклад, MLE).


#' Обчислення чисельного Якобіану функції залишків
#'
#' Використовує numDeriv::jacobian для обчислення матриці похідних.
#'
#' @param fn_residuals Функція function(theta), що повертає вектор залишків
#' @param theta Поточні значення параметрів
#' @param method Метод чисельного диференціювання ("Richardson", "simple")
#' @return Матриця Якобіану (n x p)
#' @keywords internal
compute_numerical_jacobian <- function(fn_residuals, theta, method = "Richardson") {
    if (!requireNamespace("numDeriv", quietly = TRUE)) {
        stop("Package 'numDeriv' is required for numerical Jacobian. Please install it.", call. = FALSE)
    }
    
    # Обчислюємо Якобіан: d(residuals)/d(theta)
    # Для регресії e = y - f(theta), тому d(e)/d(theta) = -d(f)/d(theta)
    # numDeriv::jacobian обчислює d(fn)/d(theta)
    J_residuals <- numDeriv::jacobian(fn_residuals, theta, method = method)
    
    # Для PMM2 solver нам потрібно J = d(f)/d(theta) = -d(e)/d(theta)
    # Тому змінюємо знак
    J <- -J_residuals
    
    return(J)
}


#' Обчислення ваг та компонентів PMM2
#'
#' @param residuals Вектор залишків
#' @return Список з моментами та параметрами PMM2 (gamma3, gamma4, weights)
compute_pmm2_components <- function(residuals) {
    n <- length(residuals)

    # Моменти
    m1 <- mean(residuals)
    m2 <- mean((residuals - m1)^2)
    m3 <- mean((residuals - m1)^3)
    m4 <- mean((residuals - m1)^4)

    if (m2 < 1e-10) {
        return(NULL)
    } # Вироджений випадок

    # Стандартизовані моменти
    gamma3 <- m3 / m2^(1.5)
    gamma4 <- m4 / m2^2 - 3

    # Перевірка знаменника PMM2
    denom <- 2 + gamma4
    if (abs(denom) < 1e-6) denom <- 1e-6

    list(
        m2 = m2,
        gamma3 = gamma3,
        gamma4 = gamma4,
        denom = denom
    )
}

#' Розв'язувач кроку PMM2 (PMM2 Step Solver)
#'
#' Розв'язує лінеаризовану систему для знаходження оновлення параметрів.
#' Базується на розкладанні Тейлора: e(theta) ~ e(theta_k) - J * delta
#'
#' @param residuals Поточні залишки
#' @param J Матриця Якобіана (n x p), де J[i, j] = - d(e_i)/d(theta_j)
#'          УВАГА: Знак мінус важливий. Якщо J це d(y_hat)/d(theta), то це ок.
#'          Якщо J це d(e)/d(theta), то у формулі має бути мінус.
#'          Тут припускаємо стандартне визначення регресії: y = f(theta) + e
#'          Тоді e = y - f(theta). d(e)/d(theta) = - d(f)/d(theta).
#'          Ми очікуємо J = d(f)/d(theta) (градієнт функції регресії).
#' @param pmm_stats Статистики з compute_pmm2_components
#'
#' @return Вектор оновлення delta
solve_pmm2_step <- function(residuals, J, pmm_stats) {
    # У класичному PMM2 для лінійної регресії y = Xb + e:
    # b_pmm = (X'X)^(-1) X' (y - c * s * (gamma3/denom)) ... це спрощено
    #
    # Більш загально, PMM2 максимізує (gamma3)^2 / (2 + gamma4).
    # Градієнт цієї функції по theta веде до системи рівнянь.
    #
    # Для реалізації ми використовуємо підхід "Iterative Reweighted Least Squares"
    # або еквівалентний метод корекції зсуву.
    #
    # Формула оновлення (спрощена для PMM2):
    # delta = (J'J)^(-1) J' * (residuals + correction)
    # Де correction залежить від асиметрії розподілу.

    # Але чекайте, PMM2 це не просто зсув. Це метод, що використовує моменти вищих порядків.
    # Основна ідея PMM2: мінімізувати дисперсію оцінки, використовуючи інформацію про розподіл.
    # Для асиметричних розподілів це означає врахування gamma3.

    # Згідно з теорією PMM (Polynomially Modified Maximum Likelihood аналог):
    # Оптимальне рівняння оцінювання: sum( w(e_t) * grad(f_t) ) = 0
    # Де w(e) - поліноміальна вагова функція. Для PMM2 (поліном 2-го порядку):
    # w(z) = z + a * (z^2 - 1)
    #
    # Коефіцієнт 'a' визначається як:
    # a = - gamma3 / (2 + gamma4)  (або з іншим знаком залежно від визначення)
    #
    # Перевіримо знак. Якщо gamma3 > 0 (правий хвіст), ми хочемо зменшити вагу великих позитивних помилок?
    # Ні, PMM2 намагається наблизити score function log-likelihood'у.
    # Score function g(z) = -f'(z)/f(z).
    # Ми апроксимуємо g(z) поліномом.

    # Використовуємо формулу з наших попередніх реалізацій (EstemPMM):
    # correction_factor = - (gamma3 / (2 + gamma4)) * (residuals^2 / sqrt(m2) - sqrt(m2))
    # Або простіше, працюємо з стандартизованими залишками z = e / sigma

    sigma <- sqrt(pmm_stats$m2)
    z <- residuals / sigma

    # Коефіцієнт при квадратичному члені
    lambda <- -pmm_stats$gamma3 / pmm_stats$denom

    # "Модифіковані" залишки, які ми хочемо зробити ортогональними до градієнта
    # pseudo_residuals = z + lambda * (z^2 - 1)
    # Ми хочемо J' * pseudo_residuals = 0
    #
    # Лінеаризація:
    # z(new) ~ z(old) - (J/sigma) * delta
    # Підставляємо в рівняння:
    # J' * [ (z - J/sigma * delta) + lambda * ((z - J/sigma * delta)^2 - 1) ] = 0
    #
    # Це дає квадратичне рівняння відносно delta, але ми можемо спростити,
    # ігноруючи квадратичні члени delta (Newton-Raphson крок).
    #
    # J' * [ z - J/sigma * delta + lambda * (z^2 - 1 - 2*z*(J/sigma)*delta) ] = 0
    # J' * [ z + lambda*(z^2 - 1) ] = J' * [ J/sigma * delta + lambda * 2 * z * J/sigma * delta ]
    # J' * [ z + lambda*(z^2 - 1) ] = J' * diag(1 + 2*lambda*z) * (J/sigma) * delta
    #
    # Ліва частина (LHS) - це градієнт цільової функції PMM2.
    # Права частина - це наближений Гессіан.

    # Ваги для Гессіана
    W_diag <- 1 + 2 * lambda * z

    # Щоб уникнути нестабільності, іноді W_diag замінюють на його очікування (=1),
    # що перетворює метод на Modified Newton або Scoring.
    # Спробуємо спрощений варіант (Scoring), де E[1 + 2*lambda*z] = 1 + 0 = 1.
    # Тоді Hessian ~ J'J / sigma.
    #
    # LHS = J' * (z + lambda * (z^2 - 1))
    # (J'J / sigma) * delta = J' * (z + lambda * (z^2 - 1))
    # delta = sigma * (J'J)^(-1) J' * (z + lambda * (z^2 - 1))
    # delta = (J'J)^(-1) J' * (sigma * z + sigma * lambda * (z^2 - 1))
    # delta = (J'J)^(-1) J' * (residuals + sigma * lambda * (z^2 - 1))

    # Вектор корекції
    correction <- sigma * lambda * (z^2 - 1)

    # Ефективні "спостереження" для кроку МНК
    y_effective <- residuals + correction

    # Крок МНК: delta = (J'J)^(-1) J' y_effective
    # Використовуємо qr.solve для стабільності
    delta <- tryCatch(
        {
            qr.solve(J, y_effective)
        },
        error = function(e) {
            # Fallback if singular
            solve(t(J) %*% J + diag(1e-6, ncol(J))) %*% t(J) %*% y_effective
        }
    )

    as.vector(delta)
}


#' Універсальний PMM2 оцінювач (Iterative)
#'
#' @param theta_init Початкові значення параметрів
#' @param fn_residuals Функція function(theta), що повертає вектор залишків
#' @param fn_jacobian Функція function(theta), що повертає матрицю Якобіана (n x p).
#'                    J[i,j] = d(y_hat_i)/d(theta_j) = -d(epsilon_i)/d(theta_j)
#'                    Якщо NULL, використовується чисельний Якобіан через numDeriv
#' @param max_iter Максимальна кількість ітерацій
#' @param tol Точність збіжності
#' @param verbose Вивід прогресу
#'
#' @return Список з результатами (theta, residuals, convergence, etc.)
#' @export
pmm2_nonlinear_iterative <- function(theta_init, fn_residuals, fn_jacobian = NULL,
                                     max_iter = 100, tol = 1e-6, verbose = FALSE) {
    theta <- theta_init
    p <- length(theta)

    if (verbose) cat("Starting Iterative PMM2...\n")

    for (iter in 1:max_iter) {
        # 1. Обчислення залишків та Якобіана в поточній точці
        res <- fn_residuals(theta)
        
        if (is.null(fn_jacobian)) {
            # Використовуємо чисельний Якобіан
            J <- compute_numerical_jacobian(fn_residuals, theta)
        } else {
            J <- fn_jacobian(theta)
        }

        # Перевірка розмірностей
        if (length(res) != nrow(J)) stop("Mismatch between residuals length and Jacobian rows")
        if (length(theta) != ncol(J)) stop("Mismatch between theta length and Jacobian cols")

        # 2. Статистики PMM2
        stats <- compute_pmm2_components(res)
        if (is.null(stats)) {
            warning("PMM2 stats computation failed (degenerate residuals).")
            break
        }

        # 3. Обчислення кроку оновлення
        delta <- solve_pmm2_step(res, J, stats)

        # 4. Оновлення параметрів
        theta_new <- theta + delta

        # 5. Перевірка збіжності
        max_change <- max(abs(delta))
        if (verbose) {
            cat(sprintf(
                "Iter %d: Max Delta = %.6f, Objective ~ %.6f\n",
                iter, max_change, stats$gamma3^2 / stats$denom
            ))
        }

        if (max_change < tol) {
            theta <- theta_new
            if (verbose) cat("Converged.\n")
            break
        }

        theta <- theta_new
    }

    # Фінальний перерахунок
    final_res <- fn_residuals(theta)
    final_stats <- compute_pmm2_components(final_res)

    list(
        coefficients = theta,
        residuals = final_res,
        iterations = iter,
        converged = (iter < max_iter),
        pmm_stats = final_stats,
        method = "Iterative PMM2"
    )
}


#' Універсальний PMM2 оцінювач (One-step Global)
#'
#' Застосовує одноразову корекцію PMM2 до результатів класичного оцінювання.
#'
#' @param theta_classical Оцінки параметрів, отримані класичним методом (наприклад, MLE)
#' @param fn_residuals Функція function(theta), що повертає вектор залишків
#' @param fn_jacobian Функція function(theta), що повертає матрицю Якобіана.
#'                    Якщо NULL, використовується чисельний Якобіан через numDeriv
#' @param verbose Вивід прогресу
#'
#' @return Список з результатами (theta, residuals, etc.)
#' @export
pmm2_nonlinear_onestep <- function(theta_classical, fn_residuals, fn_jacobian = NULL, verbose = FALSE) {
    if (verbose) cat("Starting One-step PMM2 correction...\n")

    # 1. Обчислення залишків та Якобіана в точці класичної оцінки
    # Це "Global" частина - ми оцінюємо структуру задачі один раз
    res <- fn_residuals(theta_classical)
    
    if (is.null(fn_jacobian)) {
        # Використовуємо чисельний Якобіан
        J <- compute_numerical_jacobian(fn_residuals, theta_classical)
    } else {
        J <- fn_jacobian(theta_classical)
    }

    # 2. Статистики PMM2
    stats <- compute_pmm2_components(res)
    if (is.null(stats)) {
        warning("PMM2 stats computation failed.")
        return(list(coefficients = theta_classical, method = "One-step PMM2 (Failed)"))
    }

    # 3. Обчислення єдиного кроку оновлення
    delta <- solve_pmm2_step(res, J, stats)

    # 4. Оновлення
    theta_new <- theta_classical + delta

    if (verbose) {
        cat("Classical Theta:", paste(round(theta_classical, 4), collapse = ", "), "\n")
        cat("PMM2 Correction:", paste(round(delta, 4), collapse = ", "), "\n")
        cat("Final Theta:    ", paste(round(theta_new, 4), collapse = ", "), "\n")
    }

    # Фінальні залишки (опціонально, для звіту)
    final_res <- fn_residuals(theta_new)

    list(
        coefficients = theta_new,
        residuals = final_res,
        correction = delta,
        pmm_stats = stats, # Статистики на початковому етапі (важливо для аналізу)
        method = "One-step PMM2"
    )
}
