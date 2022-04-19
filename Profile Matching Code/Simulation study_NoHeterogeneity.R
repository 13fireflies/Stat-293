## Packages
library(designmatch)
library(MASS)
source("oneprob_profmatch.R")
source("problemparameters_profmatch.R")
source("profmatch.R")

## New functions
expit = function(x){return(exp(x)/(1+exp(x)))}

## Simulations setup
#set.seed(1234)
n_cohort <- 1500
Balvars.List <- c("X_1", "X_2", "X_3", "X_4", "X_5", "X_6")
Covars.List <- list(c("X_1", "X_2", "X_3", "X_4", "X_5", "X_6"),
                    c("X_1.2", "X_2.2", "X_3", "X_4.2", "X_5.2", "X_6"),
                    c("X_1.X_3", "X_2.2", "X_4", "X_5", "X_6"))

sigmalist <- c("030", "100")
clist <- c("1", "2", "3")

argslist <- vector("character", 2^3*3)
for (i in 1:length(sigmalist)){
  for (j in 1:length(clist)){
   	ind <- 3 * (i - 1) + j
   	argslist[ind] <- paste0("00010500", sigmalist[i], clist[j])
   }
}

nreps <- 10000
for (a in 1:1){
  args <- argslist[a]
  nboots <- as.numeric(substr(args, 1, 4))
  sigma <- as.numeric(substr(args, 5, 7))
  c <- as.numeric(substr(args, 8, 8))
  
  sigma <- sqrt(sigma)
  
  ## Simulation function
  onerun <- function(r){
    ## Generating the covariates
    X_123.mat <- mvrnorm(n = n_cohort,
                         mu = c(0, 0, 0),
                         Sigma = rbind(c(2, 1, -1),
                                       c(1, 1, -0.5),
                                       c(-1, -0.5, 1)))
    X_1 <- X_123.mat[,1]
    X_2 <- X_123.mat[,2]
    X_3 <- X_123.mat[,3]
    X_4 <- runif(n_cohort, -3, 3)
    X_5 <- rchisq(n_cohort, 1)
    X_6 <- rbinom(n_cohort, 1, 0.5)
    X_7 <- rbinom(n_cohort, 1, 0.5)
    
    ## Creating interaction/2nd moment variables
    X_1.2 <- X_1^2
    X_2.2 <- X_2^2
    X_3.2 <- X_3^2
    X_4.2 <- X_4^2
    X_5.2 <- X_5^2
    X_6.2 <- X_6^2
    X_1.X_2 <- X_1 * X_2
    X_1.X_3 <- X_1 * X_3
    X_1.X_4 <- X_1 * X_4
    X_1.X_5 <- X_1 * X_5
    X_1.X_6 <- X_1 * X_6
    X_2.X_3 <- X_2 * X_3
    X_2.X_4 <- X_2 * X_4
    X_2.X_5 <- X_2 * X_5
    X_2.X_6 <- X_2 * X_6
    X_3.X_4 <- X_3 * X_4
    X_3.X_5 <- X_3 * X_5
    X_3.X_6 <- X_3 * X_6
    X_4.X_5 <- X_4 * X_5
    X_4.X_6 <- X_4 * X_6
    X_5.X_6 <- X_5 * X_6
    
    ## Selection
    S1.var <- X_1 + 2*X_2 - 2 * X_3 - X_4 - 0.5 * X_5 + X_6 + rnorm(n_cohort, 0, sigma)
    S <- ifelse(S1.var > 0, 1, 0)
    
    ## Treatment assignment
    A <- rbinom(n_cohort, 1, 0.5)
    #prob.A <- expit(-0.5 * X_1 + -0.5 * X_2 + 1 * X_3 + 1) * S
    #A <- rbinom(n_cohort, 1, prob.A)
    
    ## Outcome models (continuous)
    Y.1.0 <- X_1 + X_2 + X_3 - X_4 + X_5 + X_6 + rnorm(n_cohort, 0, 1)
    Y.1.1 <- X_1 + X_2 + X_3 - X_4 + X_5 + X_6  + rnorm(n_cohort, 0, 1)
    Y.2.0 <- X_1 + X_2 + 0.2 * X_3 * X_4 - sqrt(X_5) + rnorm(n_cohort, 0, 1)
    Y.2.1 <- X_1 + X_2 + 0.2 * X_3 * X_4 - sqrt(X_5) + rnorm(n_cohort, 0, 1)
    Y.3.0 <- (X_1 + X_2 + X_5)^2 + rnorm(n_cohort, 0, 1)
    Y.3.1 <- (X_1 + X_2 + X_5)^2 + rnorm(n_cohort, 0, 1)
    
    Y.1 <- ifelse(A == 1, Y.1.1, Y.1.0)
    Y.2 <- ifelse(A == 1, Y.2.1, Y.2.0)
    Y.3 <- ifelse(A == 1, Y.3.1, Y.3.0)
    
    ## Combining
    Data.Mat <- cbind(X_1, X_2, X_3, X_4, X_5, X_6,
                      S,
                      A,
                      Y.1, Y.2, Y.3,
                      X_1.2, X_2.2, X_3.2, X_4.2, X_5.2, X_6.2,
                      X_1.X_2, X_1.X_3, X_1.X_4, X_1.X_5, X_1.X_6,
                      X_2.X_3, X_2.X_4, X_2.X_5, X_2.X_6,
                      X_3.X_4, X_3.X_5, X_3.X_6,
                      X_4.X_5, X_4.X_6,
                      X_5.X_6)
    Data.Mat <- as.data.frame(Data.Mat)
    
    ## Overwriting treatment to 0 if not selected
    Data.Mat$A <- ifelse(Data.Mat$S == 0, 0, Data.Mat$A)
    
    ## Subsetted data sets for later use
    Data.Mat.S0 <- subset(Data.Mat, Data.Mat$S == 0)
    Data.Mat.S1 <- subset(Data.Mat, Data.Mat$S == 1)
    Data.Mat.S1.A0 <- subset(Data.Mat.S1, Data.Mat.S1$A == 0)
    Data.Mat.S1.A1 <- subset(Data.Mat.S1, Data.Mat.S1$A == 1)
    
    TSATE.1 <- mean(Y.1.1[S == 0] - Y.1.0[S == 0])
    TSATE.2 <- mean(Y.2.1[S == 0] - Y.2.0[S == 0])
    TSATE.3 <- mean(Y.3.1[S == 0] - Y.3.0[S == 0])
    
    ## Now running the profile matching approach
    Covars <- Covars.List[[c]]
    PMatch.Results <- matrix(data=0, nrow = 1, ncol=31)
    mom_tols <- vector(mode="numeric", length=length(Covars))
    
    for (i in 1:length(Covars)){
      Covar <- Covars[i]
      mom_tols[i] <- 0.05 * sd(Data.Mat.S0[Covar][,1]) ## Defining the tolerance as 0.05 times the pooled standard deviation
    }
    mom_targets <- colMeans(Data.Mat.S0[Covars])
    t_max <- 60*30
    solver <- "gurobi"
    approximate <- 0
    solver <- list(name = solver, t_max = t_max, approximate = approximate, round_cplex = 0, trace = 0)
    t_ind <- Data.Mat.S1$A
    mom_covs <- as.matrix(Data.Mat.S1[Covars])
    mom <- list(covs = mom_covs, tols = mom_tols, targets = mom_targets)
    df.match.results <- profmatch(t_ind, mom = mom, solver = solver)
    id <- df.match.results$id
    rm(df.match.results)
    Data.Mat.S1.pmatched <- Data.Mat.S1[id,]  
    Data.Mat.S1.pmatched.A0 <- subset(Data.Mat.S1.pmatched, Data.Mat.S1.pmatched$A == 0)
    Data.Mat.S1.pmatched.A1 <- subset(Data.Mat.S1.pmatched, Data.Mat.S1.pmatched$A == 1)
    ATE.Y1.PMatch <- mean(Data.Mat.S1.pmatched.A1$Y.1) - mean(Data.Mat.S1.pmatched.A0$Y.1)
    ATE.Y2.PMatch <- mean(Data.Mat.S1.pmatched.A1$Y.2) - mean(Data.Mat.S1.pmatched.A0$Y.2)
    ATE.Y3.PMatch <- mean(Data.Mat.S1.pmatched.A1$Y.3) - mean(Data.Mat.S1.pmatched.A0$Y.3)
    ESS.PMatch <- nrow(Data.Mat.S1.pmatched)
    
    ## Doubly-robust profile matching approach
    covars.formula <- paste0(Covars, sep = "+", collapse="")
    covars.formula <- substr(covars.formula, 1, nchar(covars.formula) - 1)
    if (nrow(subset(Data.Mat.S1.pmatched, Data.Mat.S1.pmatched$A == 1 & Data.Mat.S1.pmatched$S == 1)) > 0 &
        nrow(subset(Data.Mat.S1.pmatched, Data.Mat.S1.pmatched$A == 0 & Data.Mat.S1.pmatched$S == 1)) > 0){
      omod.1.Y1.A1 <- lm(as.formula(paste0("Y.1 ~", covars.formula)), 
                         data = subset(Data.Mat.S1.pmatched, Data.Mat.S1.pmatched$A == 1 & Data.Mat.S1.pmatched$S == 1))
      omod.1.Y2.A1 <- lm(as.formula(paste0("Y.2 ~", covars.formula)), 
                         data = subset(Data.Mat.S1.pmatched, Data.Mat.S1.pmatched$A == 1 & Data.Mat.S1.pmatched$S == 1))
      omod.1.Y3.A1 <- lm(as.formula(paste0("Y.3 ~", covars.formula)), 
                         data = subset(Data.Mat.S1.pmatched, Data.Mat.S1.pmatched$A == 1 & Data.Mat.S1.pmatched$S == 1))
      omod.1.Y1.A0 <- lm(as.formula(paste0("Y.1 ~", covars.formula)), 
                         data = subset(Data.Mat.S1.pmatched, Data.Mat.S1.pmatched$A == 0 & Data.Mat.S1.pmatched$S == 1))
      omod.1.Y2.A0 <- lm(as.formula(paste0("Y.2 ~", covars.formula)), 
                         data = subset(Data.Mat.S1.pmatched, Data.Mat.S1.pmatched$A == 0 & Data.Mat.S1.pmatched$S == 1))
      omod.1.Y3.A0 <- lm(as.formula(paste0("Y.3 ~", covars.formula)), 
                         data = subset(Data.Mat.S1.pmatched, Data.Mat.S1.pmatched$A == 0 & Data.Mat.S1.pmatched$S == 1))
      op.1.Y1.0 <- predict(omod.1.Y1.A0, Data.Mat.S1.pmatched, type="response")
      op.1.Y2.0 <- predict(omod.1.Y2.A0, Data.Mat.S1.pmatched, type="response")
      op.1.Y3.0 <- predict(omod.1.Y3.A0, Data.Mat.S1.pmatched, type="response")
      op.1.Y1.1 <- predict(omod.1.Y1.A1, Data.Mat.S1.pmatched, type="response")
      op.1.Y2.1 <- predict(omod.1.Y2.A1, Data.Mat.S1.pmatched, type="response")
      op.1.Y3.1 <- predict(omod.1.Y3.A1, Data.Mat.S1.pmatched, type="response")  
      
      DR.ATE.Y1 <- mean(op.1.Y1.1) - mean(op.1.Y1.0)
      DR.ATE.Y2 <- mean(op.1.Y2.1) - mean(op.1.Y2.0)
      DR.ATE.Y3 <- mean(op.1.Y3.1) - mean(op.1.Y3.0)
      
      k <- 1
      Bal.PMatch <- vector("numeric", length(Balvars.List))
      for (j in 1:length(Balvars.List)){
        Var <- Balvars.List[j]
        Bal.PMatch[k] <- abs(mean(Data.Mat.S1.pmatched.A0[Var][,1]) - mean(Data.Mat.S0[Var][,1])) / sd(Data.Mat.S0[Var][,1])
        Bal.PMatch[(k+1)] <- abs(mean(Data.Mat.S1.pmatched.A1[Var][,1]) - mean(Data.Mat.S0[Var][,1])) / sd(Data.Mat.S0[Var][,1])
        k <- k + 2
      }
    
      Results.ThisPMatch <- c(ATE.Y1.PMatch,
                              ATE.Y2.PMatch,
                              ATE.Y3.PMatch,
                              ESS.PMatch,
                              Bal.PMatch,
                              NA,
                              NA,
                              NA,
                              NA,
                              NA,
                              NA,
                              DR.ATE.Y1,
                              DR.ATE.Y2,
                              DR.ATE.Y3,
                              NA,
                              NA,
                              NA,
                              NA,
                              NA,
                              NA)
      PMatch.Results <- Results.ThisPMatch
    }
    else{
      PMatch.Results <- rep(NA, 31)
    }
    
    ## IOW
    Covars <- Covars.List[[c]]
    covars.formula <- paste0(Covars, sep = "+", collapse="")
    covars.formula <- substr(covars.formula, 1, nchar(covars.formula) - 1)
    DR <- function(covars){
      omod.1.Y1.A1 <- lm(as.formula(paste0("Y.1 ~", covars.formula)), 
                         data = subset(Data.Mat, Data.Mat$A == 1 & Data.Mat$S == 1))
      omod.1.Y2.A1 <- lm(as.formula(paste0("Y.2 ~", covars.formula)), 
                         data = subset(Data.Mat, Data.Mat$A == 1 & Data.Mat$S == 1))
      omod.1.Y3.A1 <- lm(as.formula(paste0("Y.3 ~", covars.formula)), 
                         data = subset(Data.Mat, Data.Mat$A == 1 & Data.Mat$S == 1))
      omod.1.Y1.A0 <- lm(as.formula(paste0("Y.1 ~", covars.formula)), 
                         data = subset(Data.Mat, Data.Mat$A == 0 & Data.Mat$S == 1))
      omod.1.Y2.A0 <- lm(as.formula(paste0("Y.2 ~", covars.formula)), 
                         data = subset(Data.Mat, Data.Mat$A == 0 & Data.Mat$S == 1))
      omod.1.Y3.A0 <- lm(as.formula(paste0("Y.3 ~", covars.formula)), 
                         data = subset(Data.Mat, Data.Mat$A == 0 & Data.Mat$S == 1))
      op.1.Y1.0 <- predict(omod.1.Y1.A0, Data.Mat, type="response")
      op.1.Y2.0 <- predict(omod.1.Y2.A0, Data.Mat, type="response")
      op.1.Y3.0 <- predict(omod.1.Y3.A0, Data.Mat, type="response")
      op.1.Y1.1 <- predict(omod.1.Y1.A1, Data.Mat, type="response")
      op.1.Y2.1 <- predict(omod.1.Y2.A1, Data.Mat, type="response")
      op.1.Y3.1 <- predict(omod.1.Y3.A1, Data.Mat, type="response")  
      
      DR.mu1 <- sum(Data.Mat.IOW$w.a1 * (Data.Mat.IOW$Y.1 - op.1.Y1.1) + (1 - Data.Mat.IOW$S) * op.1.Y1.1)/sum(1 - Data.Mat.IOW$S)
      DR.mu0 <- sum(Data.Mat.IOW$w.a0 * (Data.Mat.IOW$Y.1 - op.1.Y1.0) + (1 - Data.Mat.IOW$S) * op.1.Y1.0)/sum(1 - Data.Mat.IOW$S)
      DR.ATE.Y1 <- DR.mu1 - DR.mu0
      
      DR.mu1 <- sum(Data.Mat.IOW$w.a1 * (Data.Mat.IOW$Y.2 - op.1.Y2.1) + (1 - Data.Mat.IOW$S) * op.1.Y2.1)/sum(1 - Data.Mat.IOW$S)
      DR.mu0 <- sum(Data.Mat.IOW$w.a0 * (Data.Mat.IOW$Y.2 - op.1.Y2.0) + (1 - Data.Mat.IOW$S) * op.1.Y2.0)/sum(1 - Data.Mat.IOW$S)
      DR.ATE.Y2 <- DR.mu1 - DR.mu0
      
      DR.mu1 <- sum(Data.Mat.IOW$w.a1 * (Data.Mat.IOW$Y.3 - op.1.Y3.1) + (1 - Data.Mat.IOW$S) * op.1.Y3.1)/sum(1 - Data.Mat.IOW$S)
      DR.mu0 <- sum(Data.Mat.IOW$w.a0 * (Data.Mat.IOW$Y.3 - op.1.Y3.0) + (1 - Data.Mat.IOW$S) * op.1.Y3.0)/sum(1 - Data.Mat.IOW$S)
      DR.ATE.Y3 <- DR.mu1 - DR.mu0
      
      return(c(DR.ATE.Y1, DR.ATE.Y2, DR.ATE.Y3))
    }
    pmod <- glm(as.formula(paste0("S ~ ", covars.formula)), 
                data = Data.Mat, 
                family = binomial(link = "probit"))
    ps <- predict(pmod, Data.Mat, type="response")
    pa <- 0.5
    
    ## Computing the IO weights
    Data.Mat.IOW <- Data.Mat
    Data.Mat.IOW$w.a1 <-  (1 - ps)/(ps * pa) * Data.Mat.IOW$S * Data.Mat.IOW$A
    Data.Mat.IOW$w.a0 <-  (1 - ps)/(ps * (1 - pa)) * Data.Mat.IOW$S * (1 - Data.Mat.IOW$A)
    
    ## Computing the ATEs for the IO weighted estimators
    ATE.Y1.IOW <- (sum(Data.Mat.IOW$w.a1))^(-1) * sum(Data.Mat.IOW$w.a1 *  Data.Mat.IOW$Y.1) - 
      (sum(Data.Mat.IOW$w.a0))^(-1) * sum(Data.Mat.IOW$w.a0 *  Data.Mat.IOW$Y.1)
    ATE.Y2.IOW <- (sum(Data.Mat.IOW$w.a1))^(-1) * sum(Data.Mat.IOW$w.a1 *  Data.Mat.IOW$Y.2) - 
      (sum(Data.Mat.IOW$w.a0))^(-1) * sum(Data.Mat.IOW$w.a0 *  Data.Mat.IOW$Y.2)
    ATE.Y3.IOW <- (sum(Data.Mat.IOW$w.a1))^(-1) * sum(Data.Mat.IOW$w.a1 *  Data.Mat.IOW$Y.3) - 
      (sum(Data.Mat.IOW$w.a0))^(-1) * sum(Data.Mat.IOW$w.a0 *  Data.Mat.IOW$Y.3)
    
    ## Effective sample size for the IO weighted estimators
    Data.Mat.IOW$w <- Data.Mat.IOW$A * Data.Mat.IOW$S * Data.Mat.IOW$w.a1 + (1 - Data.Mat.IOW$A) * Data.Mat.IOW$S * Data.Mat.IOW$w.a0
    ESS.IOW <- sum(Data.Mat.IOW$w)^2 / (sum(Data.Mat.IOW$w^2))
    
    ## Balance for the IO weighted estimators
    k <- 1
    Bal.IOW <- vector(mode="numeric", length=length(Balvars.List))
    for (j in 1:length(Balvars.List)){
      Var <- Balvars.List[j]
      Bal.IOW[k] <- abs(
        weighted.mean(subset(Data.Mat.IOW, Data.Mat.IOW$S == 1 & Data.Mat.IOW$A == 0)[Var][,1], 
                      subset(Data.Mat.IOW, Data.Mat.IOW$S == 1 & Data.Mat.IOW$A == 0)$w) -
          mean(Data.Mat.S0[Var][,1])) / sd(Data.Mat.S0[Var][,1])
      Bal.IOW[(k+1)] <- abs(
        weighted.mean(subset(Data.Mat.IOW, Data.Mat.IOW$S == 1 & Data.Mat.IOW$A == 1)[Var][,1], 
                      subset(Data.Mat.IOW, Data.Mat.IOW$S == 1 & Data.Mat.IOW$A == 1)$w) -
          mean(Data.Mat.S0[Var][,1])) / sd(Data.Mat.S0[Var][,1])
      k <- k + 2
    }
    
    ## Combining IOW results
    Results.IOW <- c(ATE.Y1.IOW, ATE.Y2.IOW, ATE.Y3.IOW,
                     ESS.IOW, Bal.IOW)
    
    ## Doubly-robust estimators
    ## Outcome models
    DR1 <- DR("Covars")
    ATE.Y1.DR <- DR1[1]
    ATE.Y2.DR <- DR1[2]
    ATE.Y3.DR <- DR1[3]
    DR.Results <- c(ATE.Y1.DR,
                    ATE.Y2.DR,
                    ATE.Y3.DR)
    
    ## Confidence intervals for profile matching
    ATE.Y1.PMatch.boot <- vector("numeric", nboots)
    ATE.Y2.PMatch.boot <- vector("numeric", nboots)
    ATE.Y3.PMatch.boot <- vector("numeric", nboots)
    ATE.Y1.PMatchDR.boot <- vector("numeric", nboots)
    ATE.Y2.PMatchDR.boot <- vector("numeric", nboots)
    ATE.Y3.PMatchDR.boot <- vector("numeric", nboots)
    
    ## Bootstrapped SEs
    ATE.Y1.IOW.boot <- vector("numeric", nboots)
    ATE.Y2.IOW.boot <- vector("numeric", nboots)
    ATE.Y3.IOW.boot <- vector("numeric", nboots)
    ATE.Y1.DR.boot <- vector("numeric", nboots)
    ATE.Y2.DR.boot <- vector("numeric", nboots)
    ATE.Y3.DR.boot <- vector("numeric", nboots)
    #for (b in 1:nboots){
    oneboot <- function(b){
     Data.Mat.S0.boot <- Data.Mat.S0[sample(nrow(Data.Mat.S0), nrow(Data.Mat.S0), replace = TRUE), ]
     mom_tols <- vector(mode="numeric", length=length(Covars))
      for (i in 1:length(Covars)){
          Covar <- Covars[i]
          mom_tols[i] <- 0.05 * sd(Data.Mat.S0.boot[Covar][,1]) ## Defining the tolerance as 0.05 times the pooled standard deviation
        }
      mom_targets <- colMeans(Data.Mat.S0.boot[Covars])
      Data.Mat.S1.boot <- Data.Mat.S1[sample(nrow(Data.Mat.S1), nrow(Data.Mat.S1), replace = TRUE), ]
      t_ind <- Data.Mat.S1.boot$A
      mom_covs <- as.matrix(Data.Mat.S1.boot[Covars])
      mom <- list(covs = mom_covs, tols = mom_tols, targets = mom_targets)
      df.match.results <- profmatch(t_ind, mom = mom, solver = solver)
      id <- df.match.results$id
      rm(df.match.results)
      Data.Mat.S1.pmatched <- Data.Mat.S1.boot[id,]
      Data.Mat.S1.pmatched.A0 <- subset(Data.Mat.S1.pmatched, Data.Mat.S1.pmatched$A == 0)
      Data.Mat.S1.pmatched.A1 <- subset(Data.Mat.S1.pmatched, Data.Mat.S1.pmatched$A == 1)
      #ATE.Y1.PMatch.boot[b] <- mean(Data.Mat.S1.pmatched.A1$Y.1) - mean(Data.Mat.S1.pmatched.A0$Y.1)
      #ATE.Y2.PMatch.boot[b] <- mean(Data.Mat.S1.pmatched.A1$Y.2) - mean(Data.Mat.S1.pmatched.A0$Y.2)
      #ATE.Y3.PMatch.boot[b] <- mean(Data.Mat.S1.pmatched.A1$Y.3) - mean(Data.Mat.S1.pmatched.A0$Y.3) 
      ATE.Y1.PMatch.boot <- mean(Data.Mat.S1.pmatched.A1$Y.1) - mean(Data.Mat.S1.pmatched.A0$Y.1)
      ATE.Y2.PMatch.boot <- mean(Data.Mat.S1.pmatched.A1$Y.2) - mean(Data.Mat.S1.pmatched.A0$Y.2)
      ATE.Y3.PMatch.boot <- mean(Data.Mat.S1.pmatched.A1$Y.3) - mean(Data.Mat.S1.pmatched.A0$Y.3) 
      
      if (nrow(subset(Data.Mat.S1.pmatched, Data.Mat.S1.pmatched$A == 1 & Data.Mat.S1.pmatched$S == 1)) > 0 &
          nrow(subset(Data.Mat.S1.pmatched, Data.Mat.S1.pmatched$A == 0 & Data.Mat.S1.pmatched$S == 1)) > 0){
        omod.1.Y1.A1 <- lm(as.formula(paste0("Y.1 ~", covars.formula)), 
                           data = subset(Data.Mat.S1.pmatched, Data.Mat.S1.pmatched$A == 1 & Data.Mat.S1.pmatched$S == 1))
        omod.1.Y2.A1 <- lm(as.formula(paste0("Y.2 ~", covars.formula)), 
                           data = subset(Data.Mat.S1.pmatched, Data.Mat.S1.pmatched$A == 1 & Data.Mat.S1.pmatched$S == 1))
        omod.1.Y3.A1 <- lm(as.formula(paste0("Y.3 ~", covars.formula)), 
                           data = subset(Data.Mat.S1.pmatched, Data.Mat.S1.pmatched$A == 1 & Data.Mat.S1.pmatched$S == 1))
        omod.1.Y1.A0 <- lm(as.formula(paste0("Y.1 ~", covars.formula)), 
                           data = subset(Data.Mat.S1.pmatched, Data.Mat.S1.pmatched$A == 0 & Data.Mat.S1.pmatched$S == 1))
        omod.1.Y2.A0 <- lm(as.formula(paste0("Y.2 ~", covars.formula)), 
                           data = subset(Data.Mat.S1.pmatched, Data.Mat.S1.pmatched$A == 0 & Data.Mat.S1.pmatched$S == 1))
        omod.1.Y3.A0 <- lm(as.formula(paste0("Y.3 ~", covars.formula)), 
                           data = subset(Data.Mat.S1.pmatched, Data.Mat.S1.pmatched$A == 0 & Data.Mat.S1.pmatched$S == 1))
        op.1.Y1.0 <- predict(omod.1.Y1.A0, Data.Mat.S1.pmatched, type="response")
        op.1.Y2.0 <- predict(omod.1.Y2.A0, Data.Mat.S1.pmatched, type="response")
        op.1.Y3.0 <- predict(omod.1.Y3.A0, Data.Mat.S1.pmatched, type="response")
        op.1.Y1.1 <- predict(omod.1.Y1.A1, Data.Mat.S1.pmatched, type="response")
        op.1.Y2.1 <- predict(omod.1.Y2.A1, Data.Mat.S1.pmatched, type="response")
        op.1.Y3.1 <- predict(omod.1.Y3.A1, Data.Mat.S1.pmatched, type="response")  
        
        #ATE.Y1.PMatchDR.boot[b] <- mean(op.1.Y1.1) - mean(op.1.Y1.0)
        #ATE.Y2.PMatchDR.boot[b] <- mean(op.1.Y2.1) - mean(op.1.Y2.0)
        #ATE.Y3.PMatchDR.boot[b] <- mean(op.1.Y3.1) - mean(op.1.Y3.0)
        ATE.Y1.PMatchDR.boot <- mean(op.1.Y1.1) - mean(op.1.Y1.0)
        ATE.Y2.PMatchDR.boot <- mean(op.1.Y2.1) - mean(op.1.Y2.0)
        ATE.Y3.PMatchDR.boot <- mean(op.1.Y3.1) - mean(op.1.Y3.0)
      }
      else{
        #ATE.Y1.PMatchDR.boot[b] <- NA
        #ATE.Y2.PMatchDR.boot[b] <- NA
        #ATE.Y3.PMatchDR.boot[b] <- NA
        ATE.Y1.PMatchDR.boot <- NA
        ATE.Y2.PMatchDR.boot <- NA
        ATE.Y3.PMatchDR.boot <- NA
      }
      Data.Mat.boot <- rbind(Data.Mat.S1.boot, Data.Mat.S0.boot)
      
      pmod.boot <- glm(as.formula(paste0("S ~ ", covars.formula)), 
                       data = Data.Mat.boot, 
                       family = binomial(link = "probit"))
      ps.boot <- predict(pmod.boot, Data.Mat.boot, type="response")
      pa.boot <- 0.5
      
      ## Computing the IO weights
      Data.Mat.IOW.boot <- Data.Mat.boot
      Data.Mat.IOW.boot$w1.a1 <-  (1 - ps.boot)/(ps.boot * pa.boot) * Data.Mat.IOW.boot$S * Data.Mat.IOW.boot$A
      Data.Mat.IOW.boot$w1.a0 <-  (1 - ps.boot)/(ps.boot * (1 - pa.boot)) * Data.Mat.IOW.boot$S * (1 - Data.Mat.IOW.boot$A)
      DR.boot <- function(covars){
        covars.formula <- paste0(covars, sep = "+", collapse="")
        covars.formula <- substr(covars.formula, 1, nchar(covars.formula) - 1)
        omod.1.Y1.A1 <- lm(as.formula(paste0("Y.1 ~", covars.formula)), 
                           data = subset(Data.Mat.boot, Data.Mat.boot$A == 1 & Data.Mat.boot$S == 1))
        omod.1.Y2.A1 <- lm(as.formula(paste0("Y.2 ~", covars.formula)), 
                           data = subset(Data.Mat.boot, Data.Mat.boot$A == 1 & Data.Mat.boot$S == 1))
        omod.1.Y3.A1 <- lm(as.formula(paste0("Y.3 ~", covars.formula)), 
                           data = subset(Data.Mat.boot, Data.Mat.boot$A == 1 & Data.Mat.boot$S == 1))
        omod.1.Y1.A0 <- lm(as.formula(paste0("Y.1 ~", covars.formula)), 
                           data = subset(Data.Mat.boot, Data.Mat.boot$A == 0 & Data.Mat.boot$S == 1))
        omod.1.Y2.A0 <- lm(as.formula(paste0("Y.2 ~", covars.formula)), 
                           data = subset(Data.Mat.boot, Data.Mat.boot$A == 0 & Data.Mat.boot$S == 1))
        omod.1.Y3.A0 <- lm(as.formula(paste0("Y.3 ~", covars.formula)), 
                           data = subset(Data.Mat.boot, Data.Mat.boot$A == 0 & Data.Mat.boot$S == 1))
        op.1.Y1.0 <- predict(omod.1.Y1.A0, Data.Mat.boot, type="response")
        op.1.Y2.0 <- predict(omod.1.Y2.A0, Data.Mat.boot, type="response")
        op.1.Y3.0 <- predict(omod.1.Y3.A0, Data.Mat.boot, type="response")
        op.1.Y1.1 <- predict(omod.1.Y1.A1, Data.Mat.boot, type="response")
        op.1.Y2.1 <- predict(omod.1.Y2.A1, Data.Mat.boot, type="response")
        op.1.Y3.1 <- predict(omod.1.Y3.A1, Data.Mat.boot, type="response")  
        
        DR.mu1 <- sum(Data.Mat.IOW.boot$w1.a1 * (Data.Mat.IOW.boot$Y.1 - op.1.Y1.1) + (1 - Data.Mat.IOW.boot$S) * op.1.Y1.1)/sum(1 - Data.Mat.IOW.boot$S)
        DR.mu0 <- sum(Data.Mat.IOW.boot$w1.a0 * (Data.Mat.IOW.boot$Y.1 - op.1.Y1.0) + (1 - Data.Mat.IOW.boot$S) * op.1.Y1.0)/sum(1 - Data.Mat.IOW.boot$S)
        DR.ATE.Y1 <- DR.mu1 - DR.mu0
        
        DR.mu1 <- sum(Data.Mat.IOW.boot$w1.a1 * (Data.Mat.IOW.boot$Y.2 - op.1.Y2.1) + (1 - Data.Mat.IOW.boot$S) * op.1.Y2.1)/sum(1 - Data.Mat.IOW.boot$S)
        DR.mu0 <- sum(Data.Mat.IOW.boot$w1.a0 * (Data.Mat.IOW.boot$Y.2 - op.1.Y2.0) + (1 - Data.Mat.IOW.boot$S) * op.1.Y2.0)/sum(1 - Data.Mat.IOW.boot$S)
        DR.ATE.Y2 <- DR.mu1 - DR.mu0
        
        DR.mu1 <- sum(Data.Mat.IOW.boot$w1.a1 * (Data.Mat.IOW.boot$Y.3 - op.1.Y3.1) + (1 - Data.Mat.IOW.boot$S) * op.1.Y3.1)/sum(1 - Data.Mat.IOW.boot$S)
        DR.mu0 <- sum(Data.Mat.IOW.boot$w1.a0 * (Data.Mat.IOW.boot$Y.3 - op.1.Y3.0) + (1 - Data.Mat.IOW.boot$S) * op.1.Y3.0)/sum(1 - Data.Mat.IOW.boot$S)
        DR.ATE.Y3 <- DR.mu1 - DR.mu0
        
        return(c(DR.ATE.Y1, DR.ATE.Y2, DR.ATE.Y3))
      }
      ## Computing the ATEs for the IO weighted estimators
      #ATE.Y1.IOW.boot[b] <- (sum(Data.Mat.IOW.boot$w1.a1))^(-1) * sum(Data.Mat.IOW.boot$w1.a1 *  Data.Mat.IOW.boot$Y.1) - (sum(Data.Mat.IOW.boot$w1.a0))^(-1) * sum(Data.Mat.IOW.boot$w1.a0 *  Data.Mat.IOW.boot$Y.1)
      #ATE.Y2.IOW.boot[b] <- (sum(Data.Mat.IOW.boot$w1.a1))^(-1) * sum(Data.Mat.IOW.boot$w1.a1 *  Data.Mat.IOW.boot$Y.2) - (sum(Data.Mat.IOW.boot$w1.a0))^(-1) * sum(Data.Mat.IOW.boot$w1.a0 *  Data.Mat.IOW.boot$Y.2)
      #ATE.Y3.IOW.boot[b] <- (sum(Data.Mat.IOW.boot$w1.a1))^(-1) * sum(Data.Mat.IOW.boot$w1.a1 *  Data.Mat.IOW.boot$Y.3) - (sum(Data.Mat.IOW.boot$w1.a0))^(-1) * sum(Data.Mat.IOW.boot$w1.a0 *  Data.Mat.IOW.boot$Y.3)
      ATE.Y1.IOW.boot <- (sum(Data.Mat.IOW.boot$w1.a1))^(-1) * sum(Data.Mat.IOW.boot$w1.a1 *  Data.Mat.IOW.boot$Y.1) - (sum(Data.Mat.IOW.boot$w1.a0))^(-1) * sum(Data.Mat.IOW.boot$w1.a0 *  Data.Mat.IOW.boot$Y.1)
      ATE.Y2.IOW.boot <- (sum(Data.Mat.IOW.boot$w1.a1))^(-1) * sum(Data.Mat.IOW.boot$w1.a1 *  Data.Mat.IOW.boot$Y.2) - (sum(Data.Mat.IOW.boot$w1.a0))^(-1) * sum(Data.Mat.IOW.boot$w1.a0 *  Data.Mat.IOW.boot$Y.2)
      ATE.Y3.IOW.boot <- (sum(Data.Mat.IOW.boot$w1.a1))^(-1) * sum(Data.Mat.IOW.boot$w1.a1 *  Data.Mat.IOW.boot$Y.3) - (sum(Data.Mat.IOW.boot$w1.a0))^(-1) * sum(Data.Mat.IOW.boot$w1.a0 *  Data.Mat.IOW.boot$Y.3)
      
      DR.boot1 <- DR.boot(Covars)
      #ATE.Y1.DR.boot[b] <- DR.boot1[1]
      #ATE.Y2.DR.boot[b] <- DR.boot1[2]
      #ATE.Y3.DR.boot[b] <- DR.boot1[3]
      ATE.Y1.DR.boot <- DR.boot1[1]
      ATE.Y2.DR.boot <- DR.boot1[2]
      ATE.Y3.DR.boot <- DR.boot1[3]
      print(paste0("Replicate number = ", r, "; boot number = ", b))
      return(c(ATE.Y1.PMatch.boot, ATE.Y2.PMatch.boot, ATE.Y3.PMatch.boot,
             ATE.Y1.PMatchDR.boot, ATE.Y2.PMatchDR.boot, ATE.Y3.PMatchDR.boot,
             ATE.Y1.IOW.boot, ATE.Y2.IOW.boot, ATE.Y3.IOW.boot, 
             ATE.Y1.DR.boot, ATE.Y2.DR.boot, ATE.Y3.DR.boot))
    }
    boot.out <- matrix(0, nrow=nboots, ncol=12)
    for (b in 1:nboots){
    	boot.out[b, ] <- oneboot(b)
    }
    
    ATE.Y1.PMatch.boot <- boot.out[,1]
    ATE.Y2.PMatch.boot <- boot.out[,2]
    ATE.Y3.PMatch.boot <- boot.out[,3]
    ATE.Y1.PMatchDR.boot <- boot.out[,4]
    ATE.Y2.PMatchDR.boot <- boot.out[,5]
    ATE.Y3.PMatchDR.boot <- boot.out[,6]
    ATE.Y1.IOW.boot <- boot.out[,7]
    ATE.Y2.IOW.boot <- boot.out[,8]
    ATE.Y3.IOW.boot <- boot.out[,9]
    ATE.Y1.DR.boot <- boot.out[,10]
    ATE.Y2.DR.boot <- boot.out[,11]
    ATE.Y3.DR.boot <- boot.out[,12]
      
    ATE.Y1.PMatch.lower <- ATE.Y1.PMatch - 1.96 * sd(ATE.Y1.PMatch.boot, na.rm = TRUE)
    ATE.Y1.PMatch.upper <- ATE.Y1.PMatch + 1.96 * sd(ATE.Y1.PMatch.boot, na.rm = TRUE)
    ATE.Y2.PMatch.lower <- ATE.Y2.PMatch - 1.96 * sd(ATE.Y2.PMatch.boot, na.rm = TRUE)
    ATE.Y2.PMatch.upper <- ATE.Y2.PMatch + 1.96 * sd(ATE.Y2.PMatch.boot, na.rm = TRUE)
    ATE.Y3.PMatch.lower <- ATE.Y3.PMatch - 1.96 * sd(ATE.Y3.PMatch.boot, na.rm = TRUE)
    ATE.Y3.PMatch.upper <- ATE.Y3.PMatch + 1.96 * sd(ATE.Y3.PMatch.boot, na.rm = TRUE)
    ATE.Y1.PMatchDR.lower <- DR.ATE.Y1 - 1.96 * sd(ATE.Y1.PMatchDR.boot, na.rm = TRUE)
    ATE.Y1.PMatchDR.upper <- DR.ATE.Y1 + 1.96 * sd(ATE.Y1.PMatchDR.boot, na.rm = TRUE)
    ATE.Y2.PMatchDR.lower <- DR.ATE.Y2 - 1.96 * sd(ATE.Y2.PMatchDR.boot, na.rm = TRUE)
    ATE.Y2.PMatchDR.upper <- DR.ATE.Y2 + 1.96 * sd(ATE.Y2.PMatchDR.boot, na.rm = TRUE)
    ATE.Y3.PMatchDR.lower <- DR.ATE.Y3 - 1.96 * sd(ATE.Y3.PMatchDR.boot, na.rm = TRUE)
    ATE.Y3.PMatchDR.upper <- DR.ATE.Y3 + 1.96 * sd(ATE.Y3.PMatchDR.boot, na.rm = TRUE)
    PMatch.Results[c(17:22, 26:31)] <- c(ATE.Y1.PMatch.lower, ATE.Y1.PMatch.upper,
                                         ATE.Y2.PMatch.lower, ATE.Y2.PMatch.upper,
                                         ATE.Y3.PMatch.lower, ATE.Y3.PMatch.upper,
                                         ATE.Y1.PMatchDR.lower, ATE.Y1.PMatchDR.upper,
                                         ATE.Y2.PMatchDR.lower, ATE.Y2.PMatchDR.upper,
                                         ATE.Y3.PMatchDR.lower, ATE.Y3.PMatchDR.upper)
    
    ATE.Y1.IOW.CI <- c(ATE.Y1.IOW - 1.96 * sd(ATE.Y1.IOW.boot, na.rm = TRUE),
                       ATE.Y1.IOW + 1.96 * sd(ATE.Y1.IOW.boot, na.rm = TRUE))
    ATE.Y2.IOW.CI <- c(ATE.Y2.IOW - 1.96 * sd(ATE.Y2.IOW.boot, na.rm = TRUE),
                       ATE.Y2.IOW + 1.96 * sd(ATE.Y2.IOW.boot, na.rm = TRUE))
    ATE.Y3.IOW.CI <- c(ATE.Y3.IOW - 1.96 * sd(ATE.Y3.IOW.boot, na.rm = TRUE),
                       ATE.Y3.IOW + 1.96 * sd(ATE.Y3.IOW.boot, na.rm = TRUE))
    ATE.Y1.DR.CI <- c(ATE.Y1.DR - 1.96 * sd(ATE.Y1.DR.boot, na.rm = TRUE),
                      ATE.Y1.DR + 1.96 * sd(ATE.Y1.DR.boot, na.rm = TRUE))
    ATE.Y2.DR.CI <- c(ATE.Y2.DR - 1.96 * sd(ATE.Y2.DR.boot, na.rm = TRUE),
                      ATE.Y2.DR + 1.96 * sd(ATE.Y2.DR.boot, na.rm = TRUE))
    ATE.Y3.DR.CI <- c(ATE.Y3.DR - 1.96 * sd(ATE.Y3.DR.boot, na.rm = TRUE),
                      ATE.Y3.DR + 1.96 * sd(ATE.Y3.DR.boot, na.rm = TRUE))
    IOW.DR.CI.Results <- c(ATE.Y1.IOW.CI,
                           ATE.Y2.IOW.CI,
                           ATE.Y3.IOW.CI,
                           ATE.Y1.DR.CI,
                           ATE.Y2.DR.CI,
                           ATE.Y3.DR.CI)
    IOW.results <- c(Results.IOW, DR.Results, IOW.DR.CI.Results)
    return(c(TSATE.1, TSATE.2, TSATE.3, PMatch.Results, IOW.results))
  }
  
  sim.out <- matrix(0, nrow=nreps, ncol=65)
  for (r in 1:nreps){
    sim.out[r, ] <- onerun(r)
  }
  write.csv(sim.out, 
            paste0("Results_NoHet_Sigma", sigma^2, "_BalDsgn_", c, ".csv"))
}
