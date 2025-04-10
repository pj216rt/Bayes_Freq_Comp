library(tidyverse)

#interplay of Bayesian and Frequentist Analysis
#Code

#Example 2.1
#Data X_1, X_2,...,X_n ~ Poisson(\Theta) iid
#Estimate \Theta=X_bar
#using weighterd squared error.  More weight to \Theta when \Theta is small

#Constants in problem
C <- 0.2       #maximum expected loss
theta_0 <- 4   #guess
theta_U <- 9    #upper bound

##Compute sample size frequentist manner
n_point <- ceiling(sqrt(theta_0)/C)
n_worst <- ceiling(sqrt(theta_U)/C)

#Bayesian estimate
alpha <- 3
beta <- 1

prior_density <- function(theta) dgamma(theta, shape = alpha, rate = beta)
expected_sqrt_theta <- integrate(function(theta) sqrt(theta) * prior_density(theta), lower = 0, upper = Inf)$value

n_bayes <- ceiling(expected_sqrt_theta / C)

#Comparison
sample_size_df <- data.frame(
  Method = c("Point estimate", "Worst-case", "Bayesian (Gamma prior)"),
  SampleSize = c(n_point, n_worst, n_bayes)
)
print(sample_size_df)

#Exact output may change depending on seed


##Next portion discusses, an optimal design for analyzing \Theta depends on
##what the true \Theta is.  Need to account for uncertainty over the data that
##you'll observe and the true value of \Theta
##Bayesian expected power gives a much more holistic picture of how well a test
##will perform

#Poisson distribution assumption
freq.vs.bayes.power.calc <- function(theta_null, theta_alt, alpha, n,
                                     prior_a, prior_b){
  #crtical value calculation
  lambda_null <- n * theta_null
  crit_val <- qpois(1 - alpha, lambda_null)
  
  #power at thet_alt
  lambda_alt <- n * theta_alt
  classical_power <- 1 - ppois(crit_val, lambda_alt)
  
  #bayesian expected power
  #computes power at a given theta * prior prob at that given theta
  power_function <- function(theta) {
    lambda <- n * theta
    prob <- 1 - ppois(crit_val, lambda)  #power
    prob * dgamma(theta, shape = prior_a, rate = prior_b)  #weight
  }
  
  #lets us compute the expected power, integrates region of alternative hypothesis
  bayes_expected_power <- integrate(power_function, lower = 2, upper = Inf)$value
  
  power_df <- data.frame(
    Method = c("Classical", "Bayesian Expected Power"),
    Power = c(classical_power, bayes_expected_power)
  )
  
  return(power_df)
}


computation <- freq.vs.bayes.power.calc(theta_null = 2, theta_alt = 3, alpha = 0.05,
                                        n=10, prior_a = 5, prior_b = 2)


#Example 2.2
#basically recreating the coverage computations and plots from Brown, Cai, and DasGupta in 2001
#parameters
n <- 50
theta_grid <- seq(0.0, 1.0, length.out = 1000)
n_sim <- 10000

#need a function to compute CJ*
cj_star <- function(x, n, conf_level = 0.95) {
  alpha <- 1 - conf_level
  lower <- if (x == n) qbeta((alpha/2), (x+0.5), (n-x+0.5)) else qbeta((alpha/2), (x+0.5), (n-x+0.5))
  upper <- if (x == 0) qbeta((1-alpha/2), (x+0.5), (n-x+0.5)) else qbeta((1-alpha/2), (x+0.5), (n-x+0.5))
  
  #override lower and upper bounds if x==0 or x==n (edge cases)
  if (x == 0) lower <- 0
  if (x == n) upper <- 1
  
  #returns this interval as a vector
  c(lower, upper)
}

#actually do simulation
coverage <- sapply(theta_grid, function(theta) {
  x_vals <- rbinom(n_sim, size = n, prob = theta)
  
  #for each simulated set of s values, compute the CI using the cj star function
  intervals <- t(sapply(x_vals, cj_star, n = n))
  
  #check if the true theta falls inside the interval
  mean(theta >= intervals[, 1] & theta <= intervals[, 2])
})

#kernel smoothing process
#need to define a kernel smoothing function
a_theta <- function(theta, epsilon) {
  if (theta <= epsilon) return(1 - (2*epsilon))
  if (theta >= (1-epsilon)) return((1/epsilon) - 3 + (2*epsilon))
  return(((theta*(1 - theta)/(epsilon^2) - 1))*theta)
}


epsilon <- 0.05

#computing smoothed kernel coverage
smoothed_coverage <- sapply(theta_grid, function(t0) {
  a <- a_theta(t0, epsilon)
  b <- a_theta(1 - t0, epsilon)
  weights <- dbeta(theta_grid, a, b)
  sum(coverage*weights)/sum(weights)
})


#plot this
coverage_df <- data.frame(
  theta = theta_grid,
  Raw = coverage,
  Smoothed = smoothed_coverage
)

coverage_long <- coverage_df %>%
  pivot_longer(
    cols = c(Raw, Smoothed),
    names_to = "Type",
    values_to = "Coverage")

#plot this stuff
plot1 <- ggplot(coverage_long, aes(x = theta, y = Coverage, color = Type)) +
  geom_line(size = 1.2) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "gray40") +
  labs(
    title = "Raw vs. Smoothed Coverage of Modified Jeffreys Interval",
    subtitle = paste("Binomial(n = 50), epsilon =", epsilon),
    x = expression(theta),
    y = "Coverage Probability",
    color = "Coverage Type"
  ) +
  theme_minimal(base_size = 14)

plot1


#Example 3.1
#Medical diagnosis
comp.diagnosis.interval <- function(x, n, alpha=0.05, n_sim=10000){
  #input check
  if (length(x) != 3 || length(n) != 3) {
    stop("x and n must be vectors of length 3: (x0, x1, x2) and (n0, n1, n2)")
  }
  
  #need to sample from the p_i posterior distributions
  p_samples <- mapply(function(xi, ni) {
    rbeta(n_sim, shape1 = xi + 0.5, shape2 = ni - xi + 0.5)
  }, x, n)
  
  p0 <- p_samples[, 1]
  p1 <- p_samples[, 2]
  p2 <- p_samples[, 3]
  
  #compute theta
  theta_samples <- (p0*p1) / (p0*p1 + (1 - p0)*p2)
  
  lower <- quantile(theta_samples, probs = alpha / 2)
  upper <- quantile(theta_samples, probs = 1 - alpha / 2)
  
  return(list(
    lower = lower,
    upper = upper,
    credible_interval = c(lower, upper),
    theta_samples = theta_samples
  ))
}

x_vec <- c(10, 90, 10)
n_vec <- c(100, 100, 100)

result <- comp.diagnosis.interval(x = x_vec, n = n_vec)
print(result$credible_interval)
hist(result$theta_samples)


#exploring marginal distribbution for example 3.3
tau <- 2
sigma <- 1
n <- 1000000

mu <- rnorm(n, mean = 0, sd = sqrt(tau))
x <- rnorm(n, mean = mu, sd = sqrt(sigma))

#theoretical marginal: N(0, tau + sigma)
hist(x, breaks = 50, probability = TRUE, col = "skyblue",
     main = "Empirical Marginal Distribution of X",
     xlab = "X")

curve(dnorm(x, mean = 0, sd = sqrt(tau + sigma)),
      col = "red", lwd = 2, add = TRUE)
legend("topright", legend = c("Empirical", "Theoretical"),
       fill = c("skyblue", NA), border = NA, lty = c(NA, 1), col = c("skyblue", "red"))
