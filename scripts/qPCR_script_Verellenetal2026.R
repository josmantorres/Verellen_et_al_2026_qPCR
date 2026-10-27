###############################################################################
# IMPORTANT DATA ACCESSIBILITY NOTE:
# The dataset "qPCR data.xlsx" provided in this repository contains the RAW 
# Ct values extracted directly from the thermocycler. 
# However, the variable named 'ct' in the following script represents the 
# relative gene expression (Fold Change) calculated via the 2^-ΔΔCt method, 
# as detailed in the methodology section of the manuscript.
# To fully reproduce these statistical models and figures, users must first 
# transform the raw Ct values from the dataset into relative expression values 
# using the reference genes provided.
###############################################################################

###############################################################################
############################### ANTENNAL LOBE #################################
###############################################################################

rm(list=ls())

library(readxl)
# Ensure "qPCR data.xlsx" is in your current working directory
DatosqPCR <- read_excel("qPCR data.xlsx", sheet = "Lobulo antenal")
# View(DatosqPCR)

datos <- DatosqPCR
# View(datos)

### Exploratory Data Analysis (EDA)
summary(datos)
class(datos$ct)

datos$categoria <- factor(datos$categoria, levels = c("SN", "SP", "NS"))
datos$sac <- NULL
datos$gen <- factor(datos$gen)
datos$muestra <- factor(datos$muestra)
datos$ct <- as.numeric(datos$ct)

library(ggplot2)
# install.packages("ggpubr") 
library(ggpubr)

medias.datos <- aggregate(ct ~ categoria + gen, datos, mean)
medias.datos

g1 <- ggplot(medias.datos, aes(x=categoria, y=ct, colour=gen, group=gen))
g2 <- g1 + geom_line(aes(linetype=gen), size=.6) + geom_point(aes(shape=gen), size=3) 
g2 # Note: Check interaction between Pollen Followers (PF/SP) and sNPF.

### MODELING

library(glmmTMB)

m1 <- glmmTMB(ct ~  categoria*gen + (1|muestra), data = datos, family = gaussian)

##### m1 Assumptions #####
# Normality of model residuals
library(car)
rp1 <- residuals(m1, type="pearson")
qqPlot(rp1) # Checking fit with qqplot 

shapiro.test(rp1)
# Note: Residuals do not fit a normal distribution (expected). Checking histogram.
hist(datos$ct)

## Testing Gamma distribution:
library(glmmTMB)
m2 <- glmmTMB(ct ~ categoria*gen + (1|muestra), data = datos, family = Gamma(link = "log"))

## Assumptions for Gamma model
library(DHARMa)
supmgamma <- simulateResiduals(m2)
plot(supmgamma, asFactor = FALSE) 
# Note: Significant combined quantile test indicates potential non-linearity or 
# that the link function does not perfectly capture the data structure. Testing alternatives.

# Testing alternative links/families: 
m3 <- glmmTMB(ct ~ categoria*gen + (1|muestra/gen), data = datos, family = Gamma(link = "identity"))
supmgamma2 <- simulateResiduals(m3)
plot(supmgamma2, asFactor = FALSE) 

m4 <- glmmTMB(ct ~ categoria*gen + (1|muestra/gen), data = datos, family = gaussian(link = "log"))
supmgamma3 <- simulateResiduals(m4)
plot(supmgamma3, asFactor = FALSE) 

# Note: Alternative link functions do not improve residual vs predicted plots.
# Checking best fit via AIC.
AIC(m2, m3, m4)
# Selected m2 based on AIC values and theoretical fit for this data type.

### Analysis
library(car)
Anova(m2, type = 3)
# Note: Significant interaction implies that gene expression varies according to behavioral category.

### Post-hoc comparisons
library(emmeans)
compLA <- emmeans(m2, pairwise ~ categoria | gen)
compLA

library(emmeans)
library(ggplot2)
library(dplyr)

emm <- emmeans(m2, ~ categoria * gen) %>% as.data.frame()
ggplot(emm, aes(x = gen, y = emmean, color = categoria, group = categoria)) +
  geom_point(position = position_dodge(0.2)) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), 
                width = 0.2, position = position_dodge(0.2)) +
  labs(y = "Relative expression level", x = "Gene")

### Publication-quality plot with significance brackets

library(dplyr)

# Convert contrasts to data frame
compLA <- as.data.frame(compLA$contrasts)

# Split contrast into group1 and group2
compLA <- compLA %>%
  mutate(
    group1 = sapply(strsplit(as.character(contrast), " - "), `[`, 1),
    group2 = sapply(strsplit(as.character(contrast), " - "), `[`, 2),
    gen = as.factor(gen),
    p.adj = p.value 
  )

df_sigLA <- compLA %>%
  filter(p.adj < 0.05) %>%
  group_by(gen) %>%
  mutate(
    # Bracket height within each gene
    y.position = max(datos$ct[datos$gen == gen[1]]) + row_number() * 0.9
  ) %>%
  ungroup()

df_sigLA <- df_sigLA %>%
  mutate(label = case_when(
    p.adj < 0.001 ~ "***",
    p.adj < 0.01  ~ "**",
    p.adj < 0.05  ~ "*"
  ))

ggplot(datos, aes(x = categoria, y = ct, fill = categoria)) + 
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.7) +
  facet_wrap(~ gen, scales = "free_y") +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold"),       
    axis.text.y = element_text(),                    
    axis.title.y = element_text(margin = margin(r = 15), size = 12), 
    axis.title.x = element_text(margin = margin(t = 10), size = 12)  
  ) +
  labs(
    x = "Follower category",
    y = "Relative expression (Fold Change vs. Control)",
    title = "Antennal lobe"
  ) +
  stat_pvalue_manual(df_sigLA, label = "label", tip.length = 0.03,
                     bracket.size = 0.5, size = 6, inherit.aes = FALSE) + 
  scale_x_discrete(limits = c("NS","SN", "SP"),labels = c("SN" = "NF", "SP" = "PF", "NS" = "C")) 


###############################################################################
################################### BRAIN #####################################
###############################################################################

library(readxl)
DatosqPCR <- read_excel("qPCR data.xlsx", sheet = "Cerebro")
datosC <- DatosqPCR

### Exploratory Data Analysis (EDA)
summary(datosC)
class(datosC$ct)

datosC$categoria <- factor(datosC$categoria, levels = c("SN", "SP", "NS"))
datosC$gen <- factor(datosC$gen)
datosC$muestra <- factor(datosC$muestra)
datosC$ct <- as.numeric(datosC$ct)

library(ggplot2)
medias.datosC <- aggregate(ct ~ categoria + gen, datosC, mean)
g1C <- ggplot(medias.datosC, aes(x=categoria, y=ct, colour=gen, group=gen))
g2C <- g1C + geom_line(aes(linetype=gen), size=.6) + geom_point(aes(shape=gen), size=3) 
g2C # Note: Complex interaction pattern observed.

### MODELING
library(glmmTMB)
m1C <- glmmTMB(ct ~  categoria*gen + (1|muestra), data = datosC, family = gaussian)

##### m1C Assumptions #####
library(car)
rp1C <- residuals(m1C, type="pearson")
qqPlot(rp1C) 
shapiro.test(rp1C)
# Note: Residuals do not fit a normal distribution. Checking histogram.
hist(datosC$ct)

## Testing Gamma distribution:
m2C <- glmmTMB(ct ~ categoria*gen + (1|muestra), data = datosC, family = Gamma(link = "log"))

## Assumptions for Gamma model
library(DHARMa)
supmgammaC <- simulateResiduals(m2C)
plot(supmgammaC, asFactor = FALSE)       
# Note: Assumptions are reasonably met.

### Analysis
library(car)
Anova(m2C, type = 3)

# Note: Global interaction is not significant. However, specific a priori hypotheses 
# require evaluating differences between categories within each gene.
# Given 3 genes and 3 categories (9 comparisons), Benjamini-Hochberg (BH) 
# correction is applied to control the False Discovery Rate (FDR), adhering to 
# standard practices in transcriptomics (e.g., DESeq2, edgeR).

### Post-hoc comparisons
library(emmeans)

# Comparisons by gene
comp <- emmeans(m2C, pairwise ~ categoria | gen)

# Extract adjusted p-values (BH)
summary(comp$contrasts, infer = TRUE, adjust = "BH")

### Plotting
library(ggplot2)

G1C <- ggplot(datosC, aes(x = categoria, y = ct, fill = categoria)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +  
  geom_jitter(width = 0.1, alpha = 0.7) +          
  facet_wrap(~ gen, scales = "free_y") +           
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold")
  ) +
  labs(
    x = "Category",
    y = "Expression level",
    title = "Expression by category and gene"
  )
G1C

### Publication-quality plot with significance brackets

library(dplyr)
library(ggpubr)

df_comp <- as.data.frame(comp$contrasts)

df_comp <- df_comp %>%
  mutate(
    group1 = sapply(strsplit(as.character(contrast), " - "), `[`, 1),
    group2 = sapply(strsplit(as.character(contrast), " - "), `[`, 2),
    gen = as.factor(gen),
    p.adj = p.value 
  )

df_sig <- df_comp %>%
  filter(p.adj < 0.05) %>%
  group_by(gen) %>%
  mutate(
    y.position = max(datosC$ct[datosC$gen == gen[1]]) + row_number() * 0.5
  ) %>%
  ungroup()

df_sig <- df_sig %>%
  mutate(label = case_when(
    p.adj < 0.001 ~ "***",
    p.adj < 0.01  ~ "**",
    p.adj < 0.05  ~ "*"
  ))

ggplot(datosC, aes(x = categoria, y = ct, fill = categoria)) + 
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.7) +
  facet_wrap(~ gen, scales = "free_y") +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold"),       
    axis.text.y = element_text(),                    
    axis.title.y = element_text(margin = margin(r = 15), size = 12), 
    axis.title.x = element_text(margin = margin(t = 10), size = 12)  
  ) +
  labs(
    x = "Follower category",
    y = "Relative expression (Fold Change vs. Control)",
    title = "Brain"
  ) +
  stat_pvalue_manual(df_sig, label = "label", tip.length = 0.03,
                     bracket.size = 0.5, size = 6, inherit.aes = FALSE) + 
  scale_x_discrete(limits = c("NS","SN", "SP"),labels = c("SN" = "NF", "SP" = "PF", "NS" = "C")) 


###############################################################################
################################## ANTENNAE ###################################
###############################################################################

library(readxl)
DatosqPCR <- read_excel("qPCR data.xlsx", sheet = "Antena")
datosA <- DatosqPCR

### Exploratory Data Analysis (EDA)
summary(datosA)
class(datosA$ct)

datosA$categoria <- factor(datosA$categoria, levels = c("SN", "SP", "NS"))
datosA$gen <- factor(datosA$gen)
datosA$muestra <- factor(datosA$muestra)
datosA$ct <- as.numeric(datosA$ct)

library(ggplot2)
medias.datosA <- aggregate(ct ~ categoria + gen, datosA, mean)

g1A <- ggplot(medias.datosA, aes(x=categoria, y=ct, colour=gen, group=gen))
g2A <- g1A + geom_line(aes(linetype=gen), size=.6) + geom_point(aes(shape=gen), size=3) 
g2A 
# Note: Excluded outlier samples (abnormally high Ct for Rp18s internal control) 
# in the Control-NPF group to maintain data integrity.
datosA <- datosA[-c(24, 25, 26, 78, 79, 80), ]

### MODELING
library(glmmTMB)
m1A <- glmmTMB(ct ~  categoria*gen + (1|muestra), data = datosA, family = gaussian)

##### m1A Assumptions #####
library(car)
rp1A <- residuals(m1A, type="pearson")
qqPlot(rp1A) 
shapiro.test(rp1A)
# Note: Residuals do not fit a normal distribution.
hist(datosA$ct)
hist(log(datosA$ct), breaks = 20)

## Testing Gamma distribution:
m2A <- glmmTMB(ct ~ categoria*gen + (1|muestra), data = datosA, family = Gamma(link = "log"))

## Assumptions for Gamma model
library(DHARMa)
supmgammaA <- simulateResiduals(m2A)
plot(supmgammaA, asFactor = FALSE) 
# Note: Assumptions are not met. Testing alternative distributions.

###### Testing log-normal distribution 
m_lognormA <- glmmTMB(log(ct) ~ categoria*gen + (1|muestra),
                      data = datosA,
                      family = gaussian(link="identity"))

supmlognormA <- simulateResiduals(m_lognormA)
plot(supmlognormA , asFactor = FALSE)              
# Note: Assumptions are still violated.

### Testing Tweedie distribution
# install.packages(c("glmmTMB", "statmod", "tweedie"))
library(glmmTMB)
library(statmod)
library(tweedie)

m_twd <- glmmTMB(
  ct ~ categoria*gen + (1|muestra),
  data = datosA,
  family = glmmTMB::tweedie(link = "log")
)

supmtweedie <- simulateResiduals(m_twd)
plot(supmtweedie , asFactor = FALSE)              
# Note: QQplot is acceptable, but residuals vs predicted plot indicates issues. 
# This may be due to the Tweedie p-parameter estimating how variance scales with the mean.

## Modeling variance heterogeneity: 
m_twd_disp1 <- glmmTMB(
  ct ~ categoria*gen + (1|muestra),
  dispformula = ~ 1 + categoria*gen,
  data = datosA,
  family = glmmTMB::tweedie(link = "log")
)

# Checking assumptions for dispersion model
supmtweedie2 <- simulateResiduals(m_twd_disp1)
plot(supmtweedie2 , asFactor = FALSE) 
# Note: Assumptions are met. Proceeding with this model.

### Analysis
library(car)
Anova(m_twd_disp1, type = 3)

### Post-hoc comparisons
library(emmeans)
compA <- emmeans(m_twd_disp1, pairwise ~ categoria | gen)

# Extract adjusted p-values (BH)
summary(compA$contrasts, infer = TRUE, adjust = "BH")

### Plotting
library(ggplot2)

G1A <- ggplot(datosA, aes(x = categoria, y = ct, fill = categoria)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +  
  geom_jitter(width = 0.1, alpha = 0.7) +          
  facet_wrap(~ gen, scales = "free_y") +           
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold")
  ) +
  labs(
    x = "Category",
    y = "Expression level",
    title = "Expression by category and gene"
  )
G1A

### Publication-quality plot with significance brackets

library(dplyr)
library(ggpubr)

df_compA <- as.data.frame(compA$contrasts)

df_compA <- df_compA %>%
  mutate(
    group1 = sapply(strsplit(as.character(contrast), " - "), `[`, 1),
    group2 = sapply(strsplit(as.character(contrast), " - "), `[`, 2),
    gen = as.factor(gen),
    p.adj = p.value 
  )

df_sig <- df_compA %>%
  filter(p.adj < 0.05) %>%      
  mutate(
    gen = as.character(gen),
    group1 = trimws(group1),
    group2 = trimws(group2)
  )

df_sig <- df_sig %>% mutate(gen = as.character(gen))
datosA <- datosA %>% mutate(gen = as.character(gen))

df_sig <- df_sig %>%
  left_join(
    datosA %>%
      group_by(gen) %>%
      summarise(max_ct = max(ct, na.rm = TRUE)),
    by = "gen"
  ) %>%
  group_by(gen) %>%
  mutate(
    y.position = max_ct * 1.15 + row_number() * (max_ct * 0.10)
  ) %>%
  ungroup()

df_sig$y.position <- df_sig$y.position + seq(0, 0.10, length.out = nrow(df_sig))

df_sig <- df_sig %>%
  mutate(label = case_when(
    p.adj < 0.001 ~ "***",
    p.adj < 0.01  ~ "**",
    p.adj < 0.05  ~ "*"
  ))

ggplot(datosA, aes(x = categoria, y = ct, fill = categoria)) + 
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.7) +
  facet_wrap(~ gen, scales = "free_y") +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold"),       
    axis.text.y = element_text(),                    
    axis.title.y = element_text(margin = margin(r = 15), size = 12), 
    axis.title.x = element_text(margin = margin(t = 10), size = 12)  
  ) +
  labs(
    x = "Follower category",
    y = "Relative expression (Fold Change vs. Control)",
    title = "Antennae"
  ) +
  stat_pvalue_manual(df_sig, label = "label", tip.length = 0.03,
                     bracket.size = 0.5, size = 6, inherit.aes = FALSE) + 
  scale_x_discrete(limits = c("NS","SN", "SP"), labels = c("SN" = "NF", "SP" = "PF", "NS" = "C")) +
  scale_fill_manual(values = c("NS" = "#619CFF",  # Control in Blue
                               "SN" = "#F8766D",  # NF in Red
                               "SP" = "#00BA38")) # PF in Green
