% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pmm_main.R
\name{lm_pmm2}
\alias{lm_pmm2}
\title{pmm2: Головна функція для PMM2 (S=2)}
\usage{
lm_pmm2(
  formula,
  data,
  max_iter = 50,
  tol = 1e-06,
  regularize = TRUE,
  reg_lambda = 1e-08,
  na.action = na.fail,
  weights = NULL,
  verbose = FALSE
)
}
\arguments{
\item{formula}{Формула R для моделі}

\item{data}{data.frame, що містить змінні у формулі}

\item{max_iter}{ціле: максимальна кількість ітерацій для алгоритму}

\item{tol}{числове: допуск для збіжності}

\item{regularize}{логічне: додати мале значення до діагоналі для числової стабільності}

\item{reg_lambda}{числове: параметр регуляризації (якщо regularize=TRUE)}

\item{na.action}{функція для обробки відсутніх значень, за замовчуванням - na.fail}

\item{weights}{опціональний вектор ваг (поки не реалізовано)}

\item{verbose}{логічне: чи виводити інформацію про прогрес}
}
\value{
Об'єкт S4 \code{PMM2fit}
}
\description{
Підганяє лінійну модель за допомогою методу максимізації поліномів (порядок 2),
який є робастним щодо негаусівських помилок.
}
\details{
Алгоритм PMM2 працює наступним чином:
\enumerate{
\item Підганяє звичайну регресію найменших квадратів (OLS) для отримання початкових оцінок
\item Обчислює центральні моменти (m2, m3, m4) із залишків OLS
\item Ітеративно покращує оцінки параметрів за допомогою підходу на основі градієнта
}

PMM2 особливо корисний, коли терми помилок не є гаусівськими.
}
\examples{
\dontrun{
# Генерувати дані вибірки з t-розподіленими помилками
n <- 100
x <- rnorm(n)
y <- 2 + 3*x + rt(n, df=3)
dat <- data.frame(y=y, x=x)

# Підігнати модель за допомогою PMM2
fit <- lm_pmm2(y ~ x, data=dat)

# Резюме та статистичний висновок
summary(fit, formula=y~x, data=dat)
}
}
