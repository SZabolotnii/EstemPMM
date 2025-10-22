# pmm2_package.R - Holovnyi fail paketu z zalezhnostiamy ta importamy

#' @importFrom methods is new slotNames
#' @importFrom graphics abline hist legend lines par plot
#' @importFrom stats acf aggregate arima as.formula cov dnorm lm lm.fit lowess model.frame model.matrix model.response na.fail na.pass pnorm qqline qqnorm quantile arima.sim diffinv rnorm sd
#' @importFrom utils head tail
NULL

#' PMM2: Polinomialnyi metod maksymizatsii dlia robastnoi rehresii ta chasovykh riadiv
#'
#' Paket PMM2 proponuie robastni metody dlia otsiniuvannia parametriv liniinykh
#' modelei ta modelei chasovykh riadiv, robastnykh do nehausivskykh pomylok.
#'
#' @section Funktsii liniinoi rehresii:
#'
#' \code{\link{lm_pmm2}} - Pidhonka liniinoi modeli za dopomohoiu PMM2
#'
#' \code{\link{compare_with_ols}} - Porivniannia PMM2 z OLS
#'
#' @section Funktsii chasovykh riadiv:
#'
#' \code{\link{ts_pmm2}} - Zahalna funktsiia dlia pidhonky modelei chasovykh riadiv za dopomohoiu PMM2
#'
#' \code{\link{ar_pmm2}} - Pidhonka AR modelei
#'
#' \code{\link{ma_pmm2}} - Pidhonka MA modelei
#'
#' \code{\link{arma_pmm2}} - Pidhonka ARMA modelei
#'
#' \code{\link{arima_pmm2}} - Pidhonka ARIMA modelei
#'
#' \code{\link{compare_ts_methods}} - Porivniannia PMM2 z klasychnymy metodamy
#'
#' @section Statystychnyi vysnovok:
#'
#' \code{\link{pmm2_inference}} - Butstrep-vysnovok dlia liniinykh modelei
#'
#' \code{\link{ts_pmm2_inference}} - Butstrep-vysnovok dlia modelei chasovykh riadiv
#'
#' @section Utylity:
#'
#' \code{\link{pmm_skewness}} - Obchyslennia asymetrii
#'
#' \code{\link{pmm_kurtosis}} - Obchyslennia ekstsesu
#'
#' \code{\link{compute_moments}} - Obchyslennia momentiv ta kumuliantiv
#' @keywords internal
"_PACKAGE"

#' Klas S4 PMM2fit
#'
#' Klas dlia zberihannia rezultativ otsiniuvannia liniinykh modelei za dopomohoiu PMM2
#'
#' @section Sloty:
#' \describe{
#'   \item{coefficients}{Otsineni koefitsiienty}
#'   \item{residuals}{Kintsevi zalyshky}
#'   \item{m2}{Druhyi tsentralnyi moment}
#'   \item{m3}{Tretii tsentralnyi moment}
#'   \item{m4}{Chetvertyi tsentralnyi moment}
#'   \item{convergence}{Status zbizhnosti}
#'   \item{iterations}{Kilkist vykonanykh iteratsii}
#'   \item{call}{Oryhinalnyi vyklyk}
#' }
#'
#' @docType class
#' @name PMM2fit-class
NULL

#' Klas S4 TS2fit
#'
#' Bazovyi klas dlia zberihannia rezultativ otsiniuvannia modelei chasovykh riadiv za dopomohoiu PMM2
#'
#' @section Sloty:
#' \describe{
#'   \item{coefficients}{Otsineni koefitsiienty}
#'   \item{residuals}{Kintsevi zalyshky}
#'   \item{m2}{Druhyi tsentralnyi moment}
#'   \item{m3}{Tretii tsentralnyi moment}
#'   \item{m4}{Chetvertyi tsentralnyi moment}
#'   \item{convergence}{Status zbizhnosti}
#'   \item{iterations}{Kilkist vykonanykh iteratsii}
#'   \item{call}{Oryhinalnyi vyklyk}
#'   \item{model_type}{Typ modeli}
#'   \item{intercept}{Perekhoplennia}
#'   \item{original_series}{Oryhinalnyi chasovyi riad}
#'   \item{order}{Poriadky modeli}
#' }
#'
#' @docType class
#' @name TS2fit-class
NULL
