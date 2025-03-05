# PLSSEM-RCodes

```{r}
#install.packages("readxl")
#install.packages("mvnfast")
#install.packages("parameters")
#install.packages("MVN")
#install.packages("lavaan")
#install.packages("sem")
#install.packages("psych")
#install.packages("semTools")
#install.packages("semPlot")
library(readxl)
library(psych)
library(ggplot2)
library(car)
library(lavaan)
library(semPlot)
library(scales)
library(dplyr)
library(gridExtra)
library(semTools)
library(mvnfast)
library(parameters)
library(MVN)
library(plspm)
library(gmnl)
```

**MVN library to test normality and kurtosis**

```{r}
RHCSurvey <- read_excel("D:/Users/70896713H/OneDrive - Universidad de Alcala/CSIC/Papers/Morales, Ovando & Mayans 2025/CombinedDataset_complete_.xlsx",sheet = "RHC", skip = 2)
norm <- mvn((RHCSurvey[6:22]), mvnTest = "mardia")
norm.multivar <- as.data.frame(norm$multivariateNormality)
norm <- as.data.frame(norm$Descriptives)
norm
```

**Qqplot and histogram**

```{r}
qqplot <- mvn(data = (RHCSurvey[4:31]), mvnTest = "mardia", univariatePlot="qqplot")
histogram <- mvn(data = (RHCSurvey[4:31]), mvnTest = "mardia", univariatePlot="histogram")
```

***Internal consistency: Cronbach's Alpha***

```{r}
psych::alpha((RHCSurvey[4:31]))
```

**Confirmatory Factor Analysis**

*KMO and Fligner-Killeen test*

```{r}
KMO(RHCSurvey[,c(4:31)])
fligner.test(RHCSurvey[,c(4:31)])
```
*Factor loadings'*

```{r}
scree(RHCSurvey[,c(4:31)], main="Sediment graph")
#Factor loadings.31
Behaviour <- psych::fa(RHCSurvey[,c(4, 5, 7, 8, 10, 11, 13:16, 26:31)], nfactor = 5,rotate = "varimax", fm="pa") %>%
model_parameters(sort = TRUE, threshold = "max")
Behaviour

#Adter running the CFA, we select the cross loadings higher than .70
fac.SEM <- psych::fa(RHCSurvey[,c(4, 5, 7, 8, 10, 11, 13:16, 26:31)], nfactor = 1,rotate = "varimax", fm="pa") %>% model_parameters(sort = TRUE, threshold = "max")
fac.SEM
```

*Structural model*


```{r}
path_matrix <- rbind(
  c(0, 0, 0, 0, 0, 0), # Risk
  c(1, 0, 0, 0, 0, 0), # Trust
  c(0, 0, 0, 0, 0, 0), # PerceivedBehavC
  c(0, 0, 0, 0, 0, 0), # SubjectiveNorms
  c(0, 1, 0, 0, 0, 0), # Attitudes
  c(1, 1, 1, 1, 1, 0)  # Behaviour
)
colnames(path_matrix) <- rownames(path_matrix) <- c("Risk", "Trust", "PerceivedBehavC", "SubjectiveNorms", "Attitudes", "Behaviour")
# Measure model
blocks <- list(
  c("RSK2", "RSK3"),                  # Risk
  c("TR1", "TR2" ),                   # Trust
  c("PBC7", "PBC8", "PBC9"),          # PerceivedBehavC
  c("SN1", "SN2", "SN3"),             # SubjectiveNorms
  c("ATT2", "ATT3","ATT4", "ATT5"),   # Attitudes
  c("PastBehav", "OutOffset")         # Behaviour
)
# Define all variables as reflexive
modes <- c("A", "A", "A", "A", "A", "A")
#PLS-PMPBC9
pls_model <- plspm(RHCSurvey, path_matrix, blocks, modes)
summary(pls_model)
```

**VIF**

```{r}
calc_vif <- function(pls_model) {
    # Extraer las comunalidades y los nombres de las variables
    communalities <- pls_model$outer_model$communality
    variable_names <- pls_model$outer_model$name
    
    # Calcular VIF
    vif_values <- 1 / (1 - communalities)
    
    # Crear un dataframe con los nombres de las variables y los valores de VIF
    vif_df <- data.frame(Variable = variable_names, VIF = vif_values)
    
    return(vif_df)
}

# Calcular VIF y mostrarlo
vif_values <- calc_vif(pls_model)
print(vif_values)

```
**Global fit**

```{r}
pls_model$gof
```
**AVE and composite reliability**

```{r}
# AVE
ave_values <- pls_model$outer_model %>%
  group_by(block) %>%
  summarize(AVE = sum((loading^2)) / (sum((loading^2)) + sum((1 - loading^2))))
# Composite Reliability
composite_reliability <- pls_model$outer_model %>%
  group_by(block) %>%
  summarize(CR = sum((loading)^2) / (sum((loading)^2) + sum((1 - loading)^2)))
print(ave_values)
print(composite_reliability)
```
**Fornell - Larcker criterion**

```{r}
#AVE
ave_values <- pls_model$outer_model %>%
  group_by(block) %>%
  summarize(AVE = sum(loading^2) / length(loading))

#Fornell-Larcker Matrix
fornell_larcker <- matrix(NA, ncol=length(blocks), nrow=length(blocks))
rownames(fornell_larcker) <- colnames(fornell_larcker) <- c("Risk", "Trust", "Attitudes", "SubjectiveNorms", "PerceivedBehavC",    "Behaviour")

#Fill the diagonal with AVE values
diag(fornell_larcker) <- sqrt(ave_values$AVE)
latent_correlations <- cor(pls_model$scores)
# Latent correlations
for (i in 1:length(blocks)) {
  for (j in 1:length(blocks)) {
    if (i != j) {
      fornell_larcker[i, j] <- latent_correlations[i, j]
    }
  }
}
print(fornell_larcker)
```
**Heterotrait - monotrait Matrix**

```{r}
calculate_htmt <- function(pls_model, RHCSurvey, blocks) {
  # latent constructs
  scores <- pls_model$scores
  # HTMT matrix
  n_blocks <- length(blocks)
  htmt_matrix <- matrix(NA, nrow = n_blocks, ncol = n_blocks)
  rownames(htmt_matrix) <- colnames(htmt_matrix) <- names(blocks)
  for (i in 1:(n_blocks - 1)) {
    for (j in (i + 1):n_blocks) {
      indicators_i <- blocks[[i]]
      indicators_j <- blocks[[j]]
      # HTMT correlations
      heterotrait_cor <- cor(RHCSurvey[, indicators_i], RHCSurvey[, indicators_j])
      monotrait_cor_i <- cor(RHCSurvey[, indicators_i])
      monotrait_cor_j <- cor(RHCSurvey[, indicators_j])
      # Calcular HTMT
      htmt_value <- mean(abs(heterotrait_cor)) / (sqrt(mean(abs(monotrait_cor_i[lower.tri(monotrait_cor_i)]))) * sqrt(mean(abs(monotrait_cor_j[lower.tri(monotrait_cor_j)]))))
            htmt_matrix[i, j] <- htmt_value
      htmt_matrix[j, i] <- htmt_value
    }
  }
  return(htmt_matrix)
}
htmt_matrix <- calculate_htmt(pls_model, RHCSurvey, blocks)
print(htmt_matrix)
```
**Boot-strapping**

```{r}
set.seed(123)
pls_model_boot <- plspm(RHCSurvey, path_matrix, blocks, modes, boot.val = TRUE, br = 10000)
boot_results <- pls_model_boot$boot
print(boot_results)
```

***Latent scores' extraction for the logistic model***

```{r}
latent_scores <- pls_model$scores
latent_scores_df <- as.data.frame(latent_scores)
SemConstructs <- cbind(RHCSurvey, latent_scores_df)
file_path <- "SemConstructs.xlsx"
write_xlsx(SemConstructs, path = file_path)
```

**WTP Calculations by segment and type of project**
```{r}
#Mean WTP 
wtp_mean_A <- WTPAP %>%
  group_by(ProjectType, SpeciesType, Location) %>%
  summarise(PriceProject = mean(UnitPrice))
# Mean WTP
wtp_promedio_B <- WTPBP %>%
  group_by(ProjectType, SpeciesType, Location) %>%
  summarise(PriceProject = mean(UnitPrice))

print(wtp_mean_A)
print(wtp_mean_B)
```
