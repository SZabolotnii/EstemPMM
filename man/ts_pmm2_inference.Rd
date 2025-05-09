% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pmm_inference.R
\name{ts_pmm2_inference}
\alias{ts_pmm2_inference}
\title{Бутстреп-висновок для моделей часових рядів PMM2}
\usage{
ts_pmm2_inference(
  object,
  x = NULL,
  B = 200,
  seed = NULL,
  block_length = NULL,
  method = c("residual", "block"),
  parallel = FALSE,
  cores = NULL,
  debug = FALSE
)
}
\arguments{
\item{object}{об'єкт класу TS2fit}

\item{x}{(опціонально) оригінальний часовий ряд; якщо NULL, використовує object@original_series}

\item{B}{кількість бутстреп-реплікацій}

\item{seed}{(опціонально) для відтворюваності}

\item{block_length}{довжина блоку для блокового бутстрепу; якщо NULL, використовує евристичне значення}

\item{method}{тип бутстрепу: "residual" або "block"}

\item{parallel}{логічне, чи використовувати паралельні обчислення}

\item{cores}{кількість ядер для паралельних обчислень}

\item{debug}{логічне, чи виводити додаткову діагностичну інформацію}
}
\value{
data.frame з стовпцями: Estimate, Std.Error, t.value, p.value
}
\description{
Бутстреп-висновок для моделей часових рядів PMM2
}
