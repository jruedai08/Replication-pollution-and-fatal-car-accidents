#MidTerm - Jrue Dai - March 6th

# Set up the packages and directory
library(pacman)
p_load(dplyr,haven, readr, knitr, readxl, psych, ggplot2, stats4)
p_load(magrittr, qwraps2, car, lmtest, stargazer,sandwich) 
p_load(tidyr, margins)
setwd('D:\\Berkeley\\Lectures\\ARE 212\\Midterm')


#Exercise 1


#1


# load the data
mydata <- read_dta('midterm2024.dta')


# filter out the missing variables
mydatan <- mydata %>% drop_na(pm25,tlml,eastwind,numberAccidents)
mydatan <- mydatan %>% drop_na(aqi,inversion,meanrainfall)


#2


# extract data from 2000
mydatan1 <- filter(mydatan,yearnum==2000)


#3


# plot pm25 against inversion
(scatter1_1 <- ggplot(mydatan1, aes(x = inversion, y = pm25)) +
    geom_point() + # plot the scatter
    labs(title = "pm25 against inversion for year 2000",
         x = "inversion",
         y = "pm25"))

# plot pm25 against east wind
(scatter1_2 <- ggplot(mydatan1, aes(x = eastwind, y = pm25)) +
    geom_point() + # plot the scatter
    labs(title = "pm25 against eastwind for year 2000",
         x = "eastwind",
         y = "pm25"))


#Exercise 2


#1


X2 <- cbind(1,mydatan1$inversion,mydatan1$eastwind)

# ols
beta_ols2 <- solve(t(X2) %*% X2) %*% t(X2) %*% mydatan1$pm25

# beta1 hat
(beta_ols2[2])
# Since estimated beta1 is 0.1055207,
# It means that pm25 would increase 0.1055207 unit everytime inversion increase 1 Kelvin.


#2


n2 <- dim(X2)[1]
k2 <- dim(X2)[2]

# homoskedasticity
e2 <- mydatan1$pm25 - X2 %*% beta_ols2
Vb_homo2 <- solve(t(X2) %*% X2) * (as.numeric(t(e2) %*% e2)/(n2-k2))
se_homo2 <- sqrt(diag(Vb_homo2))

# heteroskedasticity
Xe2 <- cbind(e2,mydatan1$inversion*e2,mydatan1$eastwind*e2)
Vb_hete2 <- solve(t(X2) %*% X2) %*% t(Xe2) %*% Xe2 %*% solve(t(X2) %*% X2)
se_hete2 <- sqrt(diag(Vb_hete2))

# t statistics under homoskedasticity
t_homo2 <- beta_ols2[2]/se_homo2[2]
# t statistics assuming homoskedasticity is 2.2051

# t statistics under heteroskedasticity
t_hete2 <- beta_ols2[2]/se_hete2[2]
# t statistics assuming heteroskedasticity is 2.6885

# Compare these two t statistics, the one under heteroskedasticity is larger
# Which means it has smaller p value


#3


# According to the table, the p-value of t statistics under homoskedasticity is between 0.05 and 0.01
# And the p-value of t statistics under heteroskedasticity is larger then the one under homoskedasticity


#4


# Covariance between the coefficient of constant and eastwind
(Vb_hete2[1,3])


#Exercise 3


#1


X3 <- cbind(1,mydatan1$inversion,mydatan1$eastwind,mydatan1$tlml,mydatan1$meanrainfall,mydatan1$aqi)

#ols
beta_ols3 <- solve(t(X3) %*% X3) %*% t(X3) %*% mydatan1$pm25

# beta 1 hat
(beta_ols3[2])
# Since estimated beta1 is 0.08596511,
# It means that pm25 would increase 0.08596511 unit everytime inversion increase 1 Kelvin.


#2


n3 <- dim(X3)[1]
k3 <- dim(X3)[2]

# standard error under homoskedasticity
e3 <- mydatan1$pm25 - X3 %*% beta_ols3
vb_homo3 <- solve(t(X3) %*% X3) * (as.numeric(t(e3) %*% e3)/(n3-k3))
se_homo3 <- sqrt(diag(vb_homo3))

# t statistics of inversion under homoskedasticity
t_homo3 <- beta_ols3[2]/se_homo3[2]


#3


# The marginal effect would change
# Since rainfall is negative related to pm25,
# And if there is higher inversion temperature, there are less likely to have rainfall,
# Which means there is a negative relation between rainfall and inversion,
# Since then there will be a postive bias,
# the margianl effect of inversion would be larger if there's no rainfall in the model

# Back up the correlation between rainfall and inversion
(corr3 <- solve(t(mydatan1$inversion) %*% mydatan1$inversion) 
  %*% t(mydatan1$inversion) %*% mydatan1$meanrainfall)
# The correlation between inversion and rainfall is negative


#4


# the residual maker of control variables
X3_1 <- cbind(1,mydatan1$eastwind,mydatan1$tlml,mydatan1$meanrainfall,mydatan1$aqi)
M3_1 <- diag(n3) - X3_1 %*% solve(t(X3_1) %*% X3_1) %*% t(X3_1)

# the residual of regression pm25 on control variables
e3_pm25 <- M3_1 %*% mydatan1$pm25

# the residual of regression inversion on control variables
e3_inversion <- M3_1 %*% mydatan1$inversion
mydatan1 <- mutate(mydatan1,
                   residual_pm25=e3_pm25,
                   residual_inversion=e3_inversion)

# plot
(scatter3_1 <- ggplot(mydatan1,aes(x=residual_pm25,y=residual_inversion)) +
    geom_point() + # plot the scatter
    labs(title = "pm25 against inversion temperature for year 2000",
         x = "temperature",
         y = "pm25"))

# Based on the Frisch-Waugh-Lovell Theorem,
# to identify beta1
# we could first regress pm25 on all the variables except inversion and get the residual
# then we could regress inversion on other variables and get the residual
# and finally regress the residual from pm25 on residual from inversion,
# the coefficient of residual from inversion is the estimated beta1


#5


# null hypothesis: temperature, rainfall and aqi not matter
# I'll use fitted F test here

# unrestricted
e3_u <- mydatan1$pm25 - X3 %*% beta_ols3
SSR3_u <- t(e3_u) %*% e3_u

# restricted (equal to equation 1 in exercise 2)
e3_r <- mydatan1$pm25 - X2 %*% beta_ols2
SSR3_r <- t(e3_r) %*% e3_r

(F3 <- ((SSR3_r-SSR3_u)/(k3-k2))/(SSR3_u/(n3-k3)))
# Since F statistics is 12.55094 and follows F(3,2035)
# From the table we know it could reject the null hypothesis at 1 percent significance level,
# which means these three variables matter at 1 percent significance level


#Exercise 4


#1


X4 <- cbind(1,mydatan1$pm25,mydatan1$tlml,mydatan1$meanrainfall,mydatan1$aqi)

# ols
alpha_ols4 <- solve(t(X4) %*% X4) %*% t(X4) %*% mydatan1$numberAccidents

# alpha1 hat
(alpha_ols4[2])
# Since estimated alpha1 is 0.009059152,
# It means that the number of fatal accidents would increase 0.009059152 unit
# everytime pm25 increase one unit


#2


# Since there is a potential ommited variable here,
# Assumption 3 (strict exogenity) is not present, 
# which asks the mean of disturbance equals to 0 conditional on all the observations


#3


# first stage
# We've already estimated the coefficient in exercise 3, question 1
pm25_hat <- X3 %*% beta_ols3

# second stage
X4_hat <- cbind(1,pm25_hat,mydatan1$tlml,mydatan1$meanrainfall,mydatan1$aqi)
alpha_2sls4 <- solve(t(X4_hat) %*% X4_hat) %*% t(X4_hat) %*% mydatan1$numberAccidents

# the effect of pm25 on the number of accidents
(alpha_2sls4[2])


#4


# white robust heteroskedasticity consistent standard error (X4 is true X here)
e4_2sls <- mydatan1$numberAccidents - X4 %*% alpha_2sls4
Xe4_2sls <- cbind(e4_2sls,
                  pm25_hat*e4_2sls,
                  mydatan1$tlml*e4_2sls,
                  mydatan1$meanrainfall*e4_2sls,
                  mydatan1$aqi*e4_2sls)
Vb_hete4 <- solve(t(X4_hat) %*% X4_hat) %*% t(Xe4_2sls) %*% Xe4_2sls %*% solve(t(X4_hat) %*% X4_hat)
(se_hete4 <- sqrt(diag(Vb_hete4)))


#5


# t statistics under heteroskedasticity
(t4_hete <- alpha_2sls4[2]/se_hete4[2])
# t statistics is -2.002356,
# and from the table we know we could reject the null hypothesis at 5 percent significance level


#exercise 5


# unrestricted
e5_u <- mydatan1$pm25 - X3 %*% beta_ols3
SSR5_u <- t(e5_u) %*% e5_u

# restricted
Z5 <- cbind(1,mydatan1$tlml,mydatan1$meanrainfall,mydatan1$aqi)
beta_1stage5 <- solve(t(Z5) %*% Z5) %*% t(Z5) %*% mydatan1$pm25
e5_r <- mydatan1$pm25 - Z5 %*% beta_1stage5
SSR5_r <- t(e5_r) %*% e5_r

# degree of freedom
n5 <- dim(X3)[1]
k5_1 <- dim(X3)[2]
k5_2 <- dim(Z5)[2]

# F statistics
(F5 <- ((SSR5_r-SSR5_u)/(k5_1-k5_2))/(SSR5_u/(n5-k5_1)))
# Since F statistics is 20.55099, larger than 10
# There isn't a weak instrument problem


#exercise 6


#1


# first stage
# We've already estimated the coefficient in exercise 3, question 1 (equation 2)
v6 <- mydatan1$pm25 - X3 %*% beta_ols3

# second stage (control function approach)
X6_hat <- cbind(1,mydatan1$pm25,mydatan1$tlml,mydatan1$meanrainfall,mydatan1$aqi,v6)
beta6_control <- solve(t(X6_hat) %*% X6_hat) %*% t(X6_hat) %*% mydatan1$accidents

# white robust heteroskedasticity consistent standard error
e6_control <- mydatan1$accidents - X6_hat %*% beta6_control
Xe6_control <- cbind(e6_control,
                     mydatan1$pm25*e6_control,
                     mydatan1$tlml*e6_control,
                     mydatan1$meanrainfall*e6_control,
                     mydatan1$aqi*e6_control,
                     v6*e6_control)
Vb_hete6 <- solve(t(X6_hat) %*% X6_hat) %*% t(Xe6_control) %*% Xe6_control %*% solve(t(X6_hat) %*% X6_hat)
(se_hete6 <- sqrt(diag(Vb_hete6)))

# t statistics
(t6_hete <- beta6_control[2]/se_hete6[2])
# Since t statistics is -0.8227688, we couldn't reject the null hypothesis,
# which means the effect of pm25 on accident doesn't occur


#2


# get and add predictions
mydatan1 <- mutate(mydatan1,accidents_fit=X6_hat %*% beta6_control)

# plot
(scatter6_1 <- ggplot(mydatan1, aes(x=pm25,y=accidents_fit)) +
  geom_point() + # plot the scatter
  geom_hline(yintercept=0,color='red') + # plot y=0
  geom_hline(yintercept=1,color='red') + # plot y=1
  labs(title = "Predicted Probability of accidents happening against pm25 level",
       x = "pm25",
       y = "predicted probability of accidents happening"))

# The problem is that our estimated probability might get out of the interval between 0 and 1,
# which is unrealistic.


#3


# logit using control function
logit6_control <- glm(mydatan1$accidents~mydatan1$pm25+mydatan1$tlml+mydatan1$meanrainfall+mydatan1$aqi+v6,
                      family = binomial(link = "logit"))
summary(logit6_control)

# likelihood
loglik6_control <- logLik(logit6_control)

#4


# estimated probability (X6 don't include residual from first stage)
mydatan1 <- mutate(mydatan1,accidents_fit_logit=logit6_control$fitted.values)

# plot
(scatter6_2 <- ggplot(mydatan1, aes(x=pm25,y=accidents_fit_logit)) +
    geom_point() + # plot the scatter
    geom_hline(yintercept=0,color='red') + # plot y=0
    geom_hline(yintercept=1,color='red') + # plot y=1
    labs(title = "Predicted Probability of accidents happening against pm25 level",
         x = "pm25",
         y = "predicted probability of accidents happening"))

# We could see all the estimated probability of accidents is between 0 and 1 now.


#5


# step 1
# Null hypothesis: the coefficient of temperature, rainfall and aqi are all zero
# Alternative hypothesis: not all the coefficient of temperature, rainfall and aqi are zero

# Step 2
# unrestricted
# We've done it in exercise 6, question 3
loglike6_u <- loglik6_control

# Step 3
# restricted

# first stage is equal to equation 1 in exercise 2
V6_r <- mydatan1$pm25 - X2 %*% beta_ols2

# logit model
logit6_r <- glm(mydatan1$accidents~mydatan1$pm25+V6_r,
                family = binomial(link = 'logit'))
summary(logit6_r)
loglike6_r <- logLik(logit6_r)

# Step 4
# LR test
(LR <- 2*(as.numeric(loglike6_u)-as.numeric(loglike6_r)))

# Step 5
# Since the degree of freedom is 3 and the LR is 45.45594
# From the table we know it's larger than 11.34,
# which means we could reject the null hypothesis at 1% significance level


#exercise 7


#a


Na <- 100
set.seed(12345)

# data generation for one sample
Ta <- 1:Na
e1a <- rnorm(n=Na/2,mean=0,sd=2)
e2a <- rnorm(n=Na/2,mean=0,sd=8)
ea <- c(e1a,e2a)
truebeta_a <- c(20,1.5)
Ya <- truebeta_a[1]+truebeta_a[2]*Ta+ea

# ols
Xa <- cbind(1,Ta)
betahat_a <- solve(t(Xa) %*% Xa) %*% t(Xa) %*% Ya

# regression residual
e_a <- Ya - Xa %*% betahat_a

# plot
(scatter7_1 <- plot(Ta,e_a,main="residuals pattern",xlab = "T",ylab = "regression residual"))

# We could see the residual of the first half is different from the residual from the second half


#2


# define the function for the bias
BiasSimulator <- function(simulationSize, sampleSize, trueBeta) {
  OLSBiasGenerator <- function(sampleSize, trueBeta) {
    # data generation
    T <- 1:sampleSize
    e1 <- rnorm(n=sampleSize/2,mean=0,sd=2)
    e2 <- rnorm(n=sampleSize/2,mean=0,sd=8)
    e <- c(e1,e2)
    y <- trueBeta[1] + trueBeta[2] * T + e
    
    # coeffiecient esitmation
    X <- cbind(1,T)
    y <- matrix(y, ncol = 1)
    b.ols <- solve(t(X) %*% X) %*% t(X) %*% y
    b.ols %<>% as.vector()
    
    # Calculate bias
    biasBeta <- (trueBeta - b.ols) %>%
      matrix(ncol = 2) %>% 
      data.frame()
    # Set names
    names(biasBeta) <- c("interceptBias", "regressorBias")
    return(biasBeta)
  }
  
  # Run OLSBiasGenerator simulationSize times with given parameters
  simulation.dt <- lapply(
    X = 1:simulationSize,
    FUN = function(i) OLSBiasGenerator(sampleSize, trueBeta)) %>%
    # Bind together to output a data.frame
    bind_rows()
  
  return(simulation.dt)
}

# run the simulation for 10000 times
set.seed(12345)
sim10000 <- BiasSimulator(simulationSize=1e4, sampleSize=100,trueBeta=c(20,1.5))

# Check the results with a histogram for the slope
hist(sim10000[,2], breaks = 200,
                 main = "OLS  b unbiasedness of slope coefficient - sample size 10000",
                 xlab = "beta2 Bias")
abline(v = 0, col = "red", lwd = 6)

# We could see it's centered around 0.


#c


# run the simulation for 10000 times
set.seed(12345)
sim10000_1 <- BiasSimulator(simulationSize=1e4, sampleSize=100,trueBeta=c(20,1.5)) #N=100
sim10000_2 <- BiasSimulator(simulationSize=1e4, sampleSize=1e4,trueBeta=c(20,1.5)) #N=10000

# plot by overlapping
(hist7<-ggplot() + 
  geom_histogram(data = sim10000_1, aes(x = regressorBias, fill = "r")) +
  geom_histogram(data = sim10000_2, aes(x = regressorBias, fill = "g")) +
  scale_colour_manual(name ="sim.dt", values = c("r" = "red", "g"="green"), labels=c("r" = "N=100",  "g"="N=10000")) +
  scale_fill_manual(name ="Sample Size", values = c("r" = "red", "g"="green"), labels=c("r" = "N=100",  "g"="N=10000")))

# We can see the sample bias converges to zero as sample size increases.
