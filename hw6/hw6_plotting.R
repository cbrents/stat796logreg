library(readr)
library(splines)
library(ggplot2)

evans <- read_csv("Data/evans.csv")

glm_chd0 <- glm(chd~smoked + age, data=evans, family=binomial)

# Coefficients for intercept, 'smoked', and age
coef(glm_chd0)

# Plot relationship between log-odds and age
ggplot() + theme_bw() +
  geom_line(aes(x=age,
                y=coef(glm_chd0)[1] + coef(glm_chd0)[3]*age),
            data=evans) +
  xlab("Age (in years)") + ylab("Log odds of CHD")


glm_chd <- glm(chd~smoked + ns(age, 3), data=evans, family=binomial)

# Coefficients for intercept, 'smoked', and the 3 splines
coef(glm_chd)
# Using predict() with 'type="terms"' gives the value of beta*x for each variable in the model
# It combines the spline terms for you
head(predict(glm_chd, type="terms"))

# Plot relationship between log-odds and age
ggplot() + theme_bw() +
  geom_line(aes(x=age,
                y=predict(glm_chd, type="terms")[,"ns(age, 3)"]),
            data=evans) +
  xlab("Age (in years)") + ylab("Log odds of CHD")
