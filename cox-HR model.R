library(MASS)
library(dplyr)
library(ggplot2)
library(broom)
library(survival)
library(survminer)
# Assume 'mydata' is your dataset
# 'outcome' is the binary outcome variable
# 'exposure' is the main exposure of interest
# 'covariate1', 'covariate2' are the confounding variables
data<-read.csv("r.txt", header = TRUE,sep="\t")

Q1 <- quantile(data$CR, 0.25)
Q3 <- quantile(data$CR, 0.75)

data <- data %>%
  mutate(CRq = case_when(
    CR < Q1 ~ "Low",
    CR > Q3 ~ "High",
    TRUE ~ "Low"
  ))

data <- data %>%
  mutate(ascites2 = case_when(
    ascites =="L" ~ "H",
    TRUE ~ ascites
  ))

Q1 <- quantile(data$Na, 0.25)
Q3 <- quantile(data$Na, 0.75)
data <- data %>%
  mutate(Naq = case_when(
    Na < Q1 ~ "Low",
    Na > Q3 ~ "High",
    TRUE ~ "middle"
  ))

Q1 <- quantile(data$PEA, 0.25)
Q2 <- quantile(data$PEA, 0.50)
Q3 <- quantile(data$PEA, 0.75)
data <- data %>%
  mutate(PEAq = case_when(
    PEA < Q2~ "Low",
    PEA > Q2~ "High",
    TRUE ~ "middle"
  ))

data$PEAq <- factor(data$PEAq)

data$PEAq <- relevel(data$PEAq, ref = "Low")
# Cox model
cox_model <- coxph(Surv(Time, Group) ~ PEA + Age + Sex + Amonia + Na + CR + ascites2 + meld, data = data)

# Summary of the model
summary(cox_model)
tidy_model <- tidy(cox_model, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95)
write.csv (tidy_model ,"COX-HR-tidy_model.csv")
