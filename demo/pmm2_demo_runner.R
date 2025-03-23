# Завантаження пакету для розробки
library(devtools)
setwd("~/R/EstemPMM")  # Перейти до кореневої папки пакету
load_all()  # Завантажить всі функції з R/

# Завантаження необхідних пакетів
library(ggplot2)
library(gridExtra)
library(dplyr)
library(parallel)
library(reshape2)  # Додаємо для melt функції

# Основна функція для запуску демонстрації
run_demo <- function(quick = TRUE) {
  if (quick) {
    quick_test()
  } else {
    cat("Виберіть демонстрацію для запуску:\n")
    cat("1. Симуляції Монте-Карло (займає багато часу)\n")
    cat("2. Аналіз реальних даних Auto MPG\n")
    cat("3. Порівняння передбачувальної точності з крос-валідацією\n")
    cat("4. Швидке тестування\n")
    cat("0. Вихід\n\n")

    choice <- as.integer(readline("Ваш вибір (0-4): "))

    switch(choice + 1,
           cat("Вихід з демонстрації\n"),
           {
             cat("Запуск симуляцій Монте-Карло...\n")
             all_results <- run_all_simulations()
             cat("Симуляції завершено!\n")
           },
           {
             cat("Запуск аналізу реальних даних...\n")
             results <- apply_to_mpg_data()
             cat("Аналіз завершено!\n")
           },
           {
             cat("Запуск аналізу передбачувальної точності...\n")
             prediction_results <- evaluate_prediction_accuracy()
             cat("Аналіз завершено!\n")
           },
           {
             quick_test()
           },
           cat("Некоректний вибір\n")
    )
  }
}

# Запуск демонстрації
# За замовчуванням - швидке тестування
# 11run_demo(quick = TRUE)

# Для інтерактивного вибору демонстрації, розкоментуйте:
run_demo(quick = FALSE)
