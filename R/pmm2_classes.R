# pmm2_classes.R - Ієрархія класів для моделей PMM2

#' Базовий клас S4 для зберігання результатів PMM2 моделей
#'
#' @slot coefficients числовий вектор оцінених параметрів
#' @slot residuals числовий вектор кінцевих залишків
#' @slot m2 числовий другий центральний момент початкових залишків
#' @slot m3 числовий третій центральний момент початкових залишків
#' @slot m4 числовий четвертий центральний момент початкових залишків
#' @slot convergence логічний або цілий код, що вказує, чи алгоритм збігся
#' @slot iterations числова кількість виконаних ітерацій
#' @slot call оригінальний виклик функції
#'
#' @exportClass BasePMM2
setClass("BasePMM2",
         slots = c(coefficients = "numeric",
                   residuals    = "numeric",
                   m2           = "numeric",
                   m3           = "numeric",
                   m4           = "numeric",
                   convergence  = "logical",
                   iterations   = "numeric",
                   call         = "call"))

#' Клас S4 для зберігання результатів PMM2 регресійної моделі
#'
#' @slot coefficients числовий вектор оцінених параметрів
#' @slot residuals числовий вектор кінцевих залишків
#' @slot m2 числовий другий центральний момент початкових залишків
#' @slot m3 числовий третій центральний момент початкових залишків
#' @slot m4 числовий четвертий центральний момент початкових залишків
#' @slot convergence логічний або цілий код, що вказує, чи алгоритм збігся
#' @slot iterations числова кількість виконаних ітерацій
#' @slot call оригінальний виклик функції
#'
#' @exportClass PMM2fit
setClass("PMM2fit",
         contains = "BasePMM2")

#' Базовий клас S4 для зберігання результатів PMM2 моделей часових рядів
#'
#' @slot coefficients числовий вектор оцінених параметрів
#' @slot residuals числовий вектор кінцевих залишків
#' @slot m2 числовий другий центральний момент початкових залишків
#' @slot m3 числовий третій центральний момент початкових залишків
#' @slot m4 числовий четвертий центральний момент початкових залишків
#' @slot convergence логічний або цілий код, що вказує, чи алгоритм збігся
#' @slot iterations числова кількість виконаних ітерацій
#' @slot call оригінальний виклик функції
#' @slot model_type символьний рядок, що вказує тип моделі
#' @slot intercept числове значення перехоплення
#' @slot original_series числовий вектор оригінального часового ряду
#' @slot order список параметрів порядку
#'
#' @exportClass TS2fit
setClass("TS2fit",
         contains = "BasePMM2",
         slots = c(model_type      = "character",
                   intercept       = "numeric",
                   original_series = "numeric",
                   order           = "list"))

#' Клас S4 для зберігання результатів PMM2 моделі AR
#'
#' @exportClass ARPMM2
setClass("ARPMM2", contains = "TS2fit")

#' Клас S4 для зберігання результатів PMM2 моделі MA
#'
#' @exportClass MAPMM2
setClass("MAPMM2", contains = "TS2fit")

#' Клас S4 для зберігання результатів PMM2 моделі ARMA
#'
#' @exportClass ARMAPMM2
setClass("ARMAPMM2", contains = "TS2fit")

#' Клас S4 для зберігання результатів PMM2 моделі ARIMA
#'
#' @exportClass ARIMAPMM2
setClass("ARIMAPMM2", contains = "TS2fit")

#' Узагальнений метод summary для об'єктів PMM2fit
#'
#' @param object об'єкт класу "PMM2fit"
#' @param formula (опціонально) формула, використана для моделі
#' @param data (опціонально) використані дані
#' @param B кількість бутстреп-реплікацій для статистичного висновку
#' @param ... додаткові аргументи (не використовуються)
#'
#' @return Виводить резюме на консоль; повертає об'єкт (невидимо).
#'
#' @export
setMethod("summary", "PMM2fit",
          function(object, formula=NULL, data=NULL, B=100, ...) {
            cat("Результати оцінювання методом PMM2\n")
            if(!is.null(object@call)) {
              cat("Виклик:\n")
              print(object@call)
              cat("\n")
            }

            cat("Коефіцієнти:\n")
            print(object@coefficients)

            cat("\nЦентральні моменти початкових залишків:\n")
            cat("  m2 =", object@m2, "\n")
            cat("  m3 =", object@m3, "\n")
            cat("  m4 =", object@m4, "\n\n")

            vf <- pmm2_variance_factor(object@m2, object@m3, object@m4)
            if(!is.na(vf$g)) {
              cat("Теоретичні характеристики PMM2 (S = 2):\n")
              cat("  c3 =", vf$c3, "\n")
              cat("  c4 =", vf$c4, "\n")
              cat("  g  =", vf$g, " (очікуване відношення Var[PMM2]/Var[OLS])\n\n")
            }

            cat("Інформація про алгоритм:\n")
            cat("  Статус збіжності:", object@convergence, "\n")
            if("iterations" %in% slotNames(object)) {
              cat("  Ітерацій:", object@iterations, "\n\n")
            } else {
              cat("\n")
            }

            # Якщо користувач хоче побачити p-значення, викликаємо pmm2_inference:
            if(!is.null(formula) && !is.null(data)) {
              cat("Наближений статистичний висновок через бутстреп (B=", B, "):\n", sep="")
              inf_tab <- pmm2_inference(object, formula, data, B=B)
              print(inf_tab)
            } else {
              cat("Для перегляду p-значень, передайте formula= та data=\n")
            }
            invisible(object)
          }
)

#' Узагальнений метод summary для об'єктів TS2fit
#'
#' @param object об'єкт класу "TS2fit" або підкласу
#' @param ... додаткові аргументи (не використовуються)
#'
#' @return Виводить резюме на консоль; повертає об'єкт (невидимо).
#'
#' @export
setMethod("summary", "TS2fit",
          function(object, ...) {
            model_type <- object@model_type
            ar_order <- object@order$ar
            ma_order <- object@order$ma
            d <- object@order$d

            cat("Результати оцінювання часового ряду методом PMM2\n")

            # Вивести тип моделі
            cat("Тип моделі: ")
            if (model_type == "ar") {
              cat("AR(", ar_order, ")\n", sep = "")
            } else if (model_type == "ma") {
              cat("MA(", ma_order, ")\n", sep = "")
            } else if (model_type == "arma") {
              cat("ARMA(", ar_order, ",", ma_order, ")\n", sep = "")
            } else if (model_type == "arima") {
              cat("ARIMA(", ar_order, ",", d, ",", ma_order, ")\n", sep = "")
            }

            if(!is.null(object@call)) {
              cat("Виклик:\n")
              print(object@call)
              cat("\n")
            }

            # Вивести коефіцієнти з розділенням на AR та MA частини
            cat("Коефіцієнти:\n")
            if (ar_order > 0) {
              ar_coefs <- object@coefficients[1:ar_order]
              names(ar_coefs) <- paste0("ar", 1:ar_order)
              cat("AR: ")
              print(ar_coefs)
            }

            if (ma_order > 0) {
              ma_coefs <- object@coefficients[(ar_order+1):(ar_order+ma_order)]
              names(ma_coefs) <- paste0("ma", 1:ma_order)
              cat("MA: ")
              print(ma_coefs)
            }

            # Вивести перехоплення, якщо є
            if (object@intercept != 0) {
              cat("Перехоплення: ", object@intercept, "\n")
            }

            cat("\nЦентральні моменти початкових залишків:\n")
            cat("  m2 =", object@m2, "\n")
            cat("  m3 =", object@m3, "\n")
            cat("  m4 =", object@m4, "\n\n")

            vf <- pmm2_variance_factor(object@m2, object@m3, object@m4)
            if(!is.na(vf$g)) {
              cat("Теоретичні характеристики PMM2 (S = 2):\n")
              cat("  c3 =", vf$c3, "\n")
              cat("  c4 =", vf$c4, "\n")
              cat("  g  =", vf$g, " (очікуване відношення Var[PMM2]/Var[OLS])\n\n")
            }

            cat("Інформація про алгоритм:\n")
            cat("  Статус збіжності:", object@convergence, "\n")
            if("iterations" %in% slotNames(object)) {
              cat("  Ітерацій:", object@iterations, "\n\n")
            } else {
              cat("\n")
            }

            invisible(object)
          }
)
