################ PER IN SITU ANALYSIS 2023-2024 ################

rm(list=ls())

# Data loading
library(readxl)

datos <- read_excel("PER in situ.xlsx")


colnames(datos)[colnames(datos) == "#hive"] <- "hive"
colnames(datos)[colnames(datos) == "Sucrose concentration (%w/v)"] <- "sucrose"

# Factorise variables and set category levels 
datos$category <- factor(datos$category, levels = c("NF", "PF", "C"))
datos$sucrose <- factor(datos$sucrose, levels = c("0", "0.1", "0.3", "1", "3", "10", "30"))
datos$Season <- factor(datos$Season)
datos$hive <- factor(datos$hive)

### Descriptive statistics

library(pastecs) 
library(ggplot2)

str(datos)
summary(datos)

table(datos$category, datos$PER)

############################ Models ############################

# Response Variable: PER (Proboscis Extension Response) -> 1/0 -> Bernoulli distribution

# Model without interaction
m1 <- glm(PER ~ category + sucrose , data = datos, family = binomial)

# Note: The default link function for the binomial family is 'logit'

# Model with interaction
m1int <- glm(PER ~ category*sucrose , data = datos, family = binomial)

# Model with sucrose concentration (sucrose) as a continuous variable
datosnum <- datos
# Convert 'sucrose' to character first, then to numeric to preserve values
datosnum$sucrose <- as.character(datosnum$sucrose) 
datosnum$sucrose <- as.numeric(datosnum$sucrose)   
str(datosnum) 

m2 <- glm(PER ~ category + sucrose , data = datosnum, family = binomial)

# Generalized Linear Mixed Model (GLMM) with 'hive' (colony) as a random effect
library(lme4)

m3 <- glmer(PER ~ category + sucrose + (1 | hive), data = datos, family = binomial)

# Compare models using Akaike Information Criterion (AIC)
AIC(m1, m1int, m2, m3)
# Note: The mixed model (m3) is selected due to the lowest AIC value.

# Goodness-of-fit evaluation (Hosmer-Lemeshow test)
library(ResourceSelection)
hl <- hoslem.test(datos$PER, fitted(m3), g=10)
hl 

########################## Assumptions #########################

# Note: Since the response variable is dichotomous (Bernoulli), overdispersion 
# or underdispersion is not applicable, and standard residuals vs. fitted 
# plots are generally uninformative.

########################### ANALYSIS ###########################

summary(m3)
coefficients(m3)

# Note: Odds are calculated as exp(estimate) and represent the ratio of the 
# probability of success to the probability of failure.

# Odds Ratios (OR) represent the magnitude of the effect: 
# the effect of x+1 is OR times the effect of x.

library(sjPlot)
tab_model(m3)
drop1(m3, test = "Chisq")

# The OR for a given category indicates how the odds of proboscis extension 
# change compared to the reference category.

library(car)
Anova(m3, test = "Chisq")
summary(m3)

## Coefficients
# Linear Predictor scale (log-odds)
coef(m3)
confint(m3)

# Odds scale
exp(coef(m3))
exp(confint(m3))

# Set the control category (C) as the reference level
levels(datos$category)
datos$category <- relevel(datos$category, ref = "C")

tab_model(m3) 
# Interpretation note: Review the odds ratio output to determine the percentage 
# increase/decrease in the probability of PER for nectar followers (NF) and 
# pollen followers (PF) relative to the non-follower control group (C). 
# Evaluate the intraclass correlation coefficient (ICC) to report the variance 
# explained by the colony random effect.

##################### Post-hoc Comparisons #####################

library(emmeans)
options(emmeans= list(emmeans = list(infer = c(TRUE, FALSE)),contrast = list(infer = c(TRUE, TRUE))))

# Comparisons on the linear predictor scale (Tukey adjustment by default)
comp <- emmeans(m3, pairwise ~ category, adjust = "tukey") 
summary(comp)
plot(comp, comparisons = TRUE)

# Comparisons on the response variable scale (probabilities)
compVR <- emmeans(m3, pairwise ~ category, type = "response") 
summary(compVR)
plot(compVR, comparisons = TRUE)

### Simple main effects:
CompVR1 <- emmeans(m3, pairwise ~ category | sucrose, type = "response")
CompVR1 

# Biological conclusion: Nectar followers show a significantly higher PER 
# probability across sucrose concentrations compared to the control group (higher sensitivity). 
# Significant differences are also observed between PF and NF, and between PF and C.

############################# Plots ############################

library(ggplot2)
comp1 <- as.data.frame(CompVR1$emmeans)

# Initial exploratory plot
Graf1 <- ggplot(comp1, aes(x=sucrose, y=prob, fill= category)) +
  geom_point(aes(colour= category)) + geom_line(aes(colour= category)) +
  labs(x="[Sucrose]") + labs(y="PER probability") + geom_smooth(se= FALSE, aes(x=sucrose, y=prob, group=category, colour=category)) + geom_errorbar(aes(ymin=prob-SE, ymax=prob+SE, colour= category), width=0.2)

Graf1

### PUBLICATION-QUALITY PLOT ###

library(ggplot2)
library(dplyr)

# 1. Format dataframe: reorder levels 
# (Los datos crudos ya vienen con NF, PF, C, así que solo seteamos el orden)
comp1 <- comp1 %>%
  mutate(
    # Reorder levels for legend consistency
    category = factor(category, levels = c("C", "NF", "PF"))
  )

# 2. Define color palette (consistent with qPCR figures)
colores_qPCR <- c(
  "C"  = "#619CFF", # Blue
  "NF" = "#F8766D", # Red
  "PF" = "#00BA38"  # Green
)

# 3. Build the final base plot
Graf_PER <- ggplot(comp1, aes(x = sucrose, y = prob, color = category, group = category)) +
  geom_line(linetype = "dotted", linewidth = 0.5, alpha = 0.6) +
  geom_point(size = 1.75) +
  geom_errorbar(aes(ymin = prob - SE, ymax = prob + SE), 
                width = 0.05, size = 0.4) +
  scale_color_manual(values = colores_qPCR) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    axis.text.x = element_text(face = "bold"),
    axis.text.y = element_text(),
    axis.title.y = element_text(margin = margin(r = 15), size = 12),
    axis.title.x = element_text(margin = margin(t = 10), size = 12),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.margin = margin(10, 30, 10, 10) # Added right margin space for significance brackets
  ) +
  labs(
    x = "Sucrose concentration (%)",
    y = "PER probability",
    title = NULL
  )

Graf_PER
