## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- eval=FALSE--------------------------------------------------------------
#  install.packages("mvDFA")

## ---- eval=FALSE--------------------------------------------------------------
#  if(!"devtools" %in% installed.packages()) install.packages("devtools")
#  devtools::install_github("jpirmer/mvDFA")

## ----setup--------------------------------------------------------------------
library(mvDFA)

## -----------------------------------------------------------------------------
Sigma <- Sigma <- matrix(.5, 4, 4)
diag(Sigma) <- 1
Sigma[3,4] <- Sigma[4,3] <- -.3
Sigma

## -----------------------------------------------------------------------------
set.seed(2023)
X <- simulate_cMTS(N = 10^3, H = c(.2, .5, .7, .9), 
                   Sigma = Sigma, simulation_process = "FGN0", 
                   decomposition = "chol", cor_increments = FALSE)
head(X)

## ---- fig.align='center', fig.height=3, fig.width=4---------------------------
x1 <- simulate_cMTS(N = 3*10^2, H = c(.5), Sigma = as.matrix(1), 
                    simulation_process = "FGN0",
                    cor_increments = FALSE)
plot(x1$X1, main = "H = 0.5 and FGN0", type = "l")
x2 <- simulate_cMTS(N = 3*10^2, H = c(.5), Sigma = as.matrix(1), 
                    simulation_process = "FGN.fft",
                    cor_increments = FALSE)
plot(x2$X1, main = "H = 0.5 and FGN.fft", type = "l")

## ---- message=FALSE, warning=FALSE--------------------------------------------
mvDFA_result <- mvDFA(X = X, steps = 50, cores = 1, degree = 1)
mvDFA_result

## -----------------------------------------------------------------------------
mvDFA_result[6:10]

## -----------------------------------------------------------------------------
mvDFA_result$S
mvDFA_result$RMS_tot
mvDFA_result$RMS_gen
head(mvDFA_result$CovRMS_s)

## ---- fig.align='center', fig.height=5, fig.width=5---------------------------
# total
df_tot <- data.frame(mvDFA_result[c("S", "RMS_tot")])
reg_tot <- lm(I(log10(RMS_tot)) ~ 1 + I(log10(S)), data = df_tot)
coef(reg_tot)[2]
mvDFA_result$Ltot

plot(log10(df_tot))
abline(reg_tot)

## ---- fig.align='center', fig.height=5, fig.width=5---------------------------
# generalized
df_gen <- data.frame(mvDFA_result[c("S", "RMS_gen")])
reg_gen <- lm(I(log10(RMS_gen)) ~ 1 + I(log10(S)), data = df_gen)
coef(reg_gen)[2]/4
mvDFA_result$Lgen

plot(log10(df_gen))
abline(reg_gen)

