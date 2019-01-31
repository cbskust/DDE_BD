#BBD(): MCMC for A, B, alpha, and beta for birth-death process
#indata: input data set - trajctory of species 
#scale: scaling constant: for YFP model, we use 1/0.089
#nrepeat: number of iteration 
#tun.B, tun.al, tun.be: tunning constnats of MH algorithm for number of reaction r_ki, alpha, and beta respectively
# pri.A, pri.B, pri.al, pri.be: hyper prior parameters for A, B, alpha, and beta respectively. 
#                               All default settings are non-informative prior


#MCMC.summary(): plots and summary for MCMC results from BBD() function 
#result: input data set. must be same with object from BBD()
#nrepeat: number of iteration
#burn: number of burn-in steps
#thin: thining constant 
#A = T, B = F, delay = T: controlling display for MCMC simulation.  

BDD <- function(indata, scale = 1, B = 0, nrepeat, tun.B, tun.al, tun.be, 
                init = c(4, 0.5), pri.A = c(0.001, 0.001), pri.B = c(0.001, 0.001), 
                pri.al = c(0.001, 0.001), pri.be = c(0.001, 0.001)) {
  
  # Propensity for birth process with delay and death process : h_1(x,theta1) and h_2(x,theta2)
  H <- function(x, ki, th = theta) {
    G <- c(ki, x[1])
    return(th * G)
  }
  
  # Metropolis-Hastings step for number of reaction, r_ki, r_ki: number of k-th reaction and i-th time interval (k=1,2, i = 0,1,2,..., T-1)
  #      using Boys et al.'s block updatig method approximate version (dropping the Radon- Nikodym derivative dP/dQ)
  # returning updated number of reaction r_1i, r_2i 
  # We sample candiate r_1i^* from proposal distribution. This value, together with the population sizes at the ends of the interval,
  #      then determines r_2i^*.  
  # One choice of proposal is a (symmetric) discrete random walk in which 
  #      the current value is augmented by the difference (r_1i and r_1i^* ) 
  #      between two independent Poisson random variables with same mean (1+r_1i^2/b)
  # R: 2*1 vector, number of reaction at the previous iteration in MCMC, 
  # d: y_{i+1} - y_{i}, difference betwen observation at time i and time (i+1)
  # b: tunning constant for candidate 
  # H0: mean of Poisson distribution for r_ki H0 = (h(y_{i},theta)+h(y_{i+1},theta))/2; 1st equation Page 2 & 2nd Eq. on Page 3 of Supplement
  MH <- function(R, d, H0, b) {
    count = 0
    R1.m = R[1]                           # ri.m is r_1i at the previous step in MCMC iteration 
    z <- 1 + R1.m^2/b                     #(mean of proposal Poisson Dist. 1+r_1i^2/b, line 2 in P. 132 in Boys et al. paper)
    u <- rpois(1, z) - rpois(1, z)        # u: difference between two independent Poisson random variable  
    R1.star = R1.m + u                    # r_1i^* <- r_1i + u 
    R.star = cbind(R1.star, R1.star - d)  # return the candidate of number of reaction r^* = (r_1i^*, r_2i^*) 
      lambda.st <- 1 + R1.star^2/b        #(line 2 in P. 132 in Boys et al., 1+(r_1i^*)^2/b, using candidate r_1i^*)
      lambda <- 1 + R1.m^2/b              #(line 2 in P. 132 in Boys et al., 1+(r_1i)^2/b , using r_1i at current step in MCMC  )
      prop.st <- log(besselI(2 * lambda, abs(R1.star - R1.m), expon.scaled = T))   #The following two lines define the Skellam dist.
      prop <- log(besselI(2 * lambda.st, abs(R1.star - R1.m), expon.scaled = T))
      q.st <- sum(log(dpois(R.star, H0)+1e-300)) #(Poisson log likelihood in Eq. 5 of supplement for candidate r_i^*)
      q <- sum(log(dpois(R, H0)+1e-300))         #(Poisson log likelihood in Eq. 5 of supplement for r_i at current step)
      logMH <- q.st - q + prop - prop.st  #Log ccceptance ratio A (likelihood ratio scaled by proposal distribution) in the last Eq. on supplement Page 3. 
      if (!is.nan(logMH) && runif(1) < exp(logMH)) {
        R <- R.star; count = 1;
      }
    return(list(R=R, count=count))
  }
  
  # Sample the number of reaction.
  # using Boys et al.'s block updatimg method approximate version (dropping the Radon- Nikodym derivative dP/dQ)
  # This fucntion uses MH(), proposal(), and H()
  # R.in: (T-1)*2 matrix of rki, k=1,2, i=0,1,..,T-1, number of reactions at previous step in MCMC 
  #th: reaction constant, A and B
  #b; tunning constant for r_ki
  #dat: vecter for difference between y_{i} and y_{i+1}
  # block updating method updates r_1i, r_2i for each i-th time interval. This function repeats T-1 times and 
  #returns the matrix of the updated number of reactions r_ki 
  R.gen <- function(R.in, th, ki, b, dat, dif = diff) {
    tun <- b
    T1 <- nrow(R.in)
    res <- matrix(0, nrow = nrow(R.in), ncol = ncol(R.in))
    count = matrix(0, nrow = T1, ncol = 1)
    for (i in 1:T1) {
      Rin <- R.in[i, ]        # r_1k, r_2k
      D <- dif[i]             # y_{i+1} - y_{i}
      k <- ki[i]              # cumulative of kappa sum_(m=1)^i kappa(Delta,m) in Eq (5) & (6) in Sup. or gamma_k(m,Delta) in page 3.
      H.i <- (H(dat[i], k, th) + H(dat[i + 1], k, th))/2   # H0 = (h(y_{i},theta)+h(y_{i+1},theta))/2
      update = MH(R = Rin, d = D, H0 = H.i, b = tun)
      res[i,] = update$R 
      count[i,] = update$count
    }
    return(list(res=res, count=count))
  }
  
  
  # This function calculates kappa Eq. (5) and Eq. (6) in Supplement and gamma_k(m, Delta) in page 3.
  # Returning the cumulative sum of kappa: sum_{m=0}^{i}kappa(delta,m) in Eq. (6) in Supplement. 
  # P: shape parameter alpha and rate parameter of beta in gamma delay distribution. 
  KI <- function(P,ti=time) {
    a <- P[1]             #shape parameter alpha of gamma distribution in Eq. (5) and Eq. (6) in Supplement 
    b <- P[2]             #rate parameter beta of gamma distribution in Eq. (5) and Eq. (6) in Supplement .
    kappa <- rep(0, ti)
    for (m in 2:ti) {                   #discretized integral kappa(delta, m) of Eq. (5) and Eq. (6) in Supplement 
      kappa[m] = pgamma(m,a, rate = b) - pgamma((m-1), a, rate = b)
    }
    kappa[1] = pgamma(1, a, rate = b)
    k.j <- rep(1, ti)    # cumulative sum of kappa until i
    for (m in 1:ti) {
      k.all <- 0
      for (j in 1:m) k.all <- k.all + kappa[j]
      k.j[m] <- k.all
    }
    return(k.j)
  }
  
  
  
  
  # Metropolis-Hastings step for shape parameter alpha of gamma delay distribution
  # MH with random walk chain is used in the function 
  # contional posterior; alpha | A , beta, r_1i, i=0,1,...,(T-1)
  # p =(alpha, beta): parameter value at the current step in MCMC
  # Ri: T by 1 vertor (r_{11}, r_{12}, ..., r_{i,T-1} )
  # A ; A in Eq. (5) in Supplement
  # This function returns upadated alpha & current value of beta
  # refering Chib and Greenburg (1995) & son and Oh (2006)
  # using non informative gamma prior for alpha
  MH.P.a <- function(P, Ri, A, tun) {
    a.m <- P[1]
    b.m <- P[2]
    a.star <- a.m + rnorm(1, 0, tun[1])
    while (a.star < 0) { 
      a.star <- a.m + rnorm(1, 0, tun[1])
    }
    P.star <- c(a.star, b.m)
    KI.star <- KI(P.star) # calculating kappa using candidate of alpha & current beta
    KI.m <- KI(P)         # calculating kappa using current alpha & beta
    l.lik.st <- 0
    l.lik <- 0
    for (i in 1:v.num) {
      Rii <- Ri[, 1, i]
      l.lik.st <- l.lik.st + sum(Rii * log(KI.star + 1e-300)) - A * (sum(KI.star))
      l.lik <- l.lik + sum(Rii * log(KI.m + 1e-300)) - A * (sum(KI.m))
    }
    l.prior.st <- dgamma(a.star, pri.al[1], pri.al[2], log = TRUE)
    l.prior <- dgamma(a.m, pri.al[1], pri.al[2], log = TRUE)
    logMH <- (l.lik.st - l.lik + l.prior.st - l.prior 
              + pnorm(a.m, 0,  tun[1], log.p = T) - pnorm(a.star, 0, tun[1], log.p = T))
    if (!is.nan(logMH) && runif(1) < exp(logMH)) 
      P <- P.star
    return(P)
  }
  
  # Metropolis-Hastings step for rate parameter beta of gamma delay distribution
  # MH with random walk chain is used in the function 
  # contional posterior: beta|A, alpha, r_1i, i=0,1,...,(T-1)
  # p =(alpha, beta): parameter value at the current step in MCMC
  # Ri: T by 1 vertor (r_{11}, r_{12}, ..., r_{i,T-1} )
  # A ; A in Eq. (5) in Supplement
  # This function returns upadated beta & current value of alpha 
  # refering Chib and Greenburg (1995) & Son and Oh (2006)
  # using non informative gamma prior for beta
  MH.P.b <- function(P, Ri, A, tun) {
    a.m <- P[1]
    b.m <- P[2]
    b.star <- rgamma(1, b.m * tun, rate = tun)
    P.star <- c(a.m, b.star)
    KI.star <- KI(P.star)      # calculating kappa using current alpha & candidate of beta
    KI.m <- KI(P)              # calculating kappa using current alpha & beta
    l.lik.st <- 0
    l.lik <- 0
    for (i in 1:v.num) {
      Rii <- Ri[, 1, i]
      l.lik.st <- l.lik.st + sum(Rii * log(KI.star)) - A * (sum(KI.star))
      l.lik <- l.lik + sum(Rii * log(KI.m)) - A * (sum(KI.m))
    }
    l.prior.st <- dgamma(b.star, pri.be[1], pri.be[2], log = TRUE)
    l.prior <- dgamma(b.m, pri.be[1], pri.be[2], log = TRUE)
    l.prop.st <- dgamma(b.star, b.m * tun, rate = tun, log = T)
    l.prop <- dgamma(b.m, b.star * tun, rate = tun, log = T)
    logMH <- l.lik.st - l.lik + l.prior.st - l.prior - l.prop.st + l.prop
    if (!is.nan(logMH) && runif(1) < exp(logMH)) 
      P <- P.star
    return(P)
  }
  
  ################################################ iteration start!!!
  
  data <- as.matrix(round(indata * scale))
  time <- nrow(data) - 1
  v.num <- ncol(data)
  
  theta_out <- matrix(0, nrow = nrepeat, ncol = 4)
  count_R <- matrix(0, nrow = time, ncol = v.num)
  count_al <- 0
  count_be <- 0
  
  # setting  initial values
  g_1 <- time
  RS <- matrix(0, nrow = 2, ncol = v.num)
  g <- matrix(0, nrow = 2, ncol = v.num)
  RR <- array(0, dim = c(time, 2, v.num))
  diff <- matrix(0, nrow = time, ncol = v.num)
  for (j in 1:v.num) {
    g_2 <- sum(data[, j]) - (data[1, j] + data[time + 1, j])/2
    g[, j] <- c(g_1, g_2)
    for (i in 1:time) diff[i, j] <- data[i + 1, j] - data[i, j]
    for (i in 1:time) {
      RR[i, 1, j] <- max(diff[i, j], 0)
      RR[i, 2, j] <- max(-diff[i, j], 0)
    }
    RS[, j] <- apply(RR[, , j], 2, sum)
  }
  g <- apply(g, 1, sum)
  theta1 <- rgamma(1, shape = sum(RS[1, ]) + pri.A[1], rate = (g[1]) + pri.A[2])
  if (B == 0) {
    theta2 <- rgamma(1, shape = sum(RS[2, ]) + pri.B[1], rate = (g[2]) + pri.B[2])
  } else {
    theta2 <- B
  }
  
  theta <- c(theta1, theta2)
  p <- init
  K.i <- KI(p)
  
  for (rep in 1:nrepeat) {
    for (j in 1:v.num) {
      R.update <- R.gen(RR[ , ,j], th = theta, ki = K.i, b = tun.B, dat = data[ ,j], dif=diff[ ,j])
      RR[ , , j] <- R.update$res
      count_R[ ,j] <- count_R[ ,j] + R.update$count 
      RS[ ,j] <- apply(RR[ , ,j], 2, sum)  
    }
    g_11 <- sum(K.i)
    g[1] <- g_11 * v.num
    theta1 <- rgamma(1, shape = sum(RS[1, ]) + pri.A[1], rate = g[1] + pri.A[2])
    if (B == 0) {
      theta2 <- rgamma(1, shape = sum(RS[2, ]) + pri.B[1], rate = g[2] + pri.B[2])
    } else {
      theta2 <- B
    }
    theta <- c(theta1, theta2)
    p.update <- MH.P.a(p, RR, theta[1], tun = tun.al)
    if (p.update[1] != p[1]) count_al <- count_al + 1
    p <- p.update
    p.update <- MH.P.b(p, RR, theta[1], tun = tun.be)
    if (p.update[2] != p[2]) count_be <- count_be + 1
    p <- p.update
    K.i <- KI(p)
    theta_out[rep, 1:2] <- theta
    theta_out[rep, 3:4] <- p
  }
  print(summary(count_R/nrepeat))
  cat("Acceptance rato of alpha: ", count_al/nrepeat, "\n")
  cat("Acceptance rato of beta:  ", count_be/nrepeat, "\n")
  colnames(theta_out) <- c("A", "B", "alpha", "beta")
  return(theta_out)
}


MCMC.summary <- function(result, nrepeat, burn, thin, A = T, B = T, delay = T) {
  sel <- seq(0, nrepeat - burn, by = thin)
  sel <- sel[-1] + burn
  theta1 <- result[sel, 1]
  theta2 <- result[sel, 2]
  al <- result[sel, 3]
  be <- result[sel, 4]
  mean.delay <- (al)/(be)
  var.delay <- (al)/(be^2)
  
  if (A == T) {
    par(oma = c(0, 0, 4, 0), mar = c(4, 4, 1, 1))
    mat <- matrix(c(1, 1, 2, 3), 2, 2, byrow = T)
    layout(mat)
    plot(theta1, type = "l", main = "", xlab = "iteration", ylab = "A")
    acf(theta1, main = "")
    plot(density(theta1), main = "", xlab = "A")
    mtext(side = 3, line = 1, outer = T, text = "Summary of death rate (A) ", 
          cex = 1.5)
  }
  if (B == T) {
    par(oma = c(0, 0, 4, 0), mar = c(4, 4, 1, 1))
    mat <- matrix(c(1, 1, 2, 3), 2, 2, byrow = T)
    layout(mat)
    plot(theta2, type = "l", main = "", xlab = "iteration", ylab = "B")
    acf(theta2, main = "")
    plot(density(theta2), main = "", xlab = "B")
    mtext(side = 3, line = 1, outer = T, text = "Summary of death rate (B)", 
          cex = 1.5)
  }
  if (delay == T) {
    par(oma = c(0, 0, 4, 0), mar = c(4, 4, 1, 1))
    mat <- matrix(c(1, 1, 2, 3), 2, 2, byrow = T)
    layout(mat)
    plot(mean.delay, type = "l", main = "", xlab = "iteration", ylab = "Delay time")
    acf(mean.delay, main = "")
    plot(density(mean.delay), main = "", xlab = "Delay time")
    mtext(side = 3, line = 1, outer = T, text = "Summary of delay time", 
          cex = 1.5)
  }
  
  result.full <- cbind(theta1, theta2, mean.delay, var.delay, al, be)
  colnames(result.full) <- c("A", "B", "mean delay", "Var. delay", "alpha", 
                             "beta")
  print(summary(result.full))
  return(result.full)
}

################################################
# Simulation start!!!
################################################

load("YFP1.rda") # YFP experiment data set 1 (40 trajectories)
load("YFP2.rda") # YFP experiment data set 2 (29 trajectories)
load("YFP1_short.rda") # trancated at T=22 of YFP1
load("YFP2_short.rda") # trancated at T=22 of YFP2

#example of YFP1_short data set. using 40 frajectories 
short1 <-BDD(indata=YFP1_short,scale=1/0.089, B=0.015,nrepeat = 110000, tun.B=300, tun.al = 0.05, tun.be = 7000)
short11<-MCMC.summary(result = short1, nrepeat=110000, burn=10000, thin=100,T,F,T)

  #example of YFP1_short data set. using only trajectory #5  
short1 <-BDD(indata=YFP1_short[,5],scale=1/0.089, B=0.015,nrepeat = 110000, tun.B=300, tun.al = 1.3, tun.be = 50)
short11<-MCMC.summary(result = short1, nrepeat=110000, burn=10000, thin=100,T,F,T)




