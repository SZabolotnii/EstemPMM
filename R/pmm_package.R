# pmm_package.R - Головний файл пакету з залежностями та імпортами

#' @importFrom methods is new slotNames
#' @importFrom graphics abline hist legend lines par
#' @importFrom stats acf arima cov dnorm lm lm.fit lowess model.frame model.matrix model.response na.fail na.pass pnorm qqline qqnorm quantile arima.sim diffinv rnorm sd
#' @importFrom utils tail
NULL

#' PMM2: Поліноміальний метод максимізації для робастної регресії та часових рядів
#'
#' Пакет PMM2 пропонує робастні методи для оцінювання параметрів лінійних
#' моделей та моделей часових рядів, робастних до негаусівських помилок.
#'
#' @section Функції лінійної регресії:
#'
#' \code{\link{lm_pmm2}} - Підгонка лінійної моделі за допомогою PMM2
#'
#' \code{\link{compare_with_ols}} - Порівняння PMM2 з OLS
#'
#' @section Функції часових рядів:
#'
#' \code{\link{ts_pmm2}} - Загальна функція для підгонки моделей часових рядів за допомогою PMM2
#'
#' \code{\link{ar_pmm2}} - Підгонка AR моделей
#'
#' \code{\link{ma_pmm2}} - Підгонка MA моделей
#'
#' \code{\link{arma_pmm2}} - Підгонка ARMA моделей
#'
#' \code{\link{arima_pmm2}} - Підгонка ARIMA моделей
#'
#' \code{\link{compare_ts_methods}} - Порівняння PMM2 з класичними методами
#'
#' @section Статистичний висновок:
#'
#' \code{\link{pmm2_inference}} - Бутстреп-висновок для лінійних моделей
#'
#' \code{\link{ts_pmm2_inference}} - Бутстреп-висновок для моделей часових рядів
#'
#' @section Утиліти:
#'
#' \code{\link{pmm_skewness}} - Обчислення асиметрії
#'
#' \code{\link{pmm_kurtosis}} - Обчислення ексцесу
#'
#' \code{\link{compute_moments}} - Обчислення моментів та кумулянтів
#'
#' @docType package
#' @name pmm2
NULL

#' Клас S4 PMM2fit
#'
#' Клас для зберігання результатів оцінювання лінійних моделей за допомогою PMM2
#'
#' @section Слоти:
#' \describe{
#'   \item{coefficients}{Оцінені коефіцієнти}
#'   \item{residuals}{Кінцеві залишки}
#'   \item{m2}{Другий центральний момент}
#'   \item{m3}{Третій центральний момент}
#'   \item{m4}{Четвертий центральний момент}
#'   \item{convergence}{Статус збіжності}
#'   \item{iterations}{Кількість виконаних ітерацій}
#'   \item{call}{Оригінальний виклик}
#' }
#'
#' @docType class
#' @name PMM2fit-class
NULL

#' Клас S4 TS2fit
#'
#' Базовий клас для зберігання результатів оцінювання моделей часових рядів за допомогою PMM2
#'
#' @section Слоти:
#' \describe{
#'   \item{coefficients}{Оцінені коефіцієнти}
#'   \item{residuals}{Кінцеві залишки}
#'   \item{m2}{Другий центральний момент}
#'   \item{m3}{Третій центральний момент}
#'   \item{m4}{Четвертий центральний момент}
#'   \item{convergence}{Статус збіжності}
#'   \item{iterations}{Кількість виконаних ітерацій}
#'   \item{call}{Оригінальний виклик}
#'   \item{model_type}{Тип моделі}
#'   \item{intercept}{Перехоплення}
#'   \item{original_series}{Оригінальний часовий ряд}
#'   \item{order}{Порядки моделі}
#' }
#'
#' @docType class
#' @name TS2fit-class
NULL
