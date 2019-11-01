func_simvals.ith <- function(rep, data){
  ordered <- sort(data)
  hat.sigma2 <- var(ordered)
  n.0 <- length(data)
  n.1 <- n.0 - 1
  k <- 1
  y <- ordered
  X <-c(rep(1, n.0)) 
  b <- (1/(n.0)) %*% t(X)%*%y
  avg.simvals.ith <- c()
  for (i.th in 1:n.0){
    print(i.th)
    X.tilde <- X
    X.tilde[i.th] <- 0
    E.y.tilde <- X.tilde %*% b
    M = t(X)%*%X + t(X.tilde)%*%X.tilde
    I = diag(n.0)
    H <- (I - X.tilde %*% (1/M) %*% t(X.tilde) ) / hat.sigma2
    cov.y.tilde <- ((n.1-k)/(n.1-k-2)) * inv(H)
    simvals <- rmvt(rep, E.y.tilde, cov.y.tilde, df=n.1-k)
    ordered.simvals <- apply(simvals, 1, sort) #ordered.simvals
    avg.ordered.simvals <- apply(ordered.simvals, 1, mean) #avg.ordered.simvals
    avg.simval.ith <- avg.ordered.simvals[i.th]
    avg.simvals.ith <- cbind(avg.simvals.ith, avg.simval.ith)
  }
  avg.simvals.ith
}

```
```{r, echo=FALSE, results='hide' }
set.seed(030202) #Homeowrk 03 | Problem 02 | Part 02
#plant.sim.10    <- func_simvals.ith(10000,   plants) #run before hand, takes about 15 sec
plant.sim.10 <- c(-56.471258,
                  -35.582806,
                  -21.728706,
                  -11.742245,
                  -3.014056,
                  5.209234,
                  12.456428,
                  19.608197,
                  26.880061,
                  34.464025,
                  41.925134,
                  51.107852,
                  61.070947,
                  74.436406,
                  95.816467)
#plant.sim.100   <- func_simvals.ith(100000,  plants) #run before hand, takes about 1-2 min
plant.sim.100 <- c(-56.847186,
                   -35.265482,
                   -22.043177,
                   -11.702150,
                   -2.987761,
                   4.941145,
                   12.420899,
                   19.609107,
                   26.844291,
                   34.316421,
                   42.167696,
                   50.755623,
                   61.064745,
                   74.224380,
                   95.455229)


func_simdensityplots <- function(data, i.th){
  ordered <- sort(data)
  
  hat.sigma2 <- var(ordered)
  n.0 <- length(data)
  n.1 <- n.0 - 1
  k <- 1
  y <- ordered
  X <-c(rep(1, n.0)) 
  b <- (1/(n.0)) %*% t(X)%*%y
  
  simvals.grp <- c()
  it.grp <- c()
  percent <- c()
  it.vals <- c(10000, 100000)
  
  for (it in it.vals){
    X.tilde <- X
    X.tilde[i.th] <- 0
    E.y.tilde <- X.tilde %*% b
    M = t(X)%*%X + t(X.tilde)%*%X.tilde
    I = diag(n.0)
    H <- (I - X.tilde %*% (1/M) %*% t(X.tilde) ) / hat.sigma2
    cov.y.tilde <- ((n.1-k)/(n.1-k-2)) * inv(H)
    simvals <- apply(rmvt(it, E.y.tilde, cov.y.tilde, df=n.1-k), 1, sort)[i.th,]
    p.it <- length(simvals[simvals >= ordered[i.th]]) / length(simvals)
    it.vec <- c(rep(it, length(simvals)))
    simvals.grp <- c(simvals.grp, simvals)
    it.grp <- c(it.grp, it.vec)
    simulations <- cbind(it.grp, simvals.grp)
    percent <- c(percent, p.it )
  }
  
  simulations.df <- as.data.frame(simulations)
  colnames(simulations.df)
  
  
  dat_text <- data.frame(
    label = c(
      paste("P(X > i.th | Y.-i) =", percent[1])
      ,paste("P(X > i.th | Y.-i) =", percent[2])
      #            ,paste("P(X > i.th | Y.-i) =", decimal(percent[3], dec))
    ),
    it.grp   = it.vals)
  
  plot <- ggplot(simulations.df, aes(x=simvals.grp)) + geom_histogram(aes(y=..density..), color="grey", alpha=0.4,) +
    facet_wrap(~it.grp, nrow=1) + 
    geom_vline(xintercept=ordered[i.th], color="red", linetype="dashed")+
    ggtitle(paste("Simulated Distribution of the", i.th, "ordered value, based on plant differences missing the", i.th, "value (", ordered[i.th], ")"))+
    xlab("") + #xlim(-600, 400)+ylim(0, 0.009)+
    geom_text( data    = dat_text,  
               mapping = aes(x = -Inf, y = -Inf, label = label)
               ,hjust   = -0,vjust   = -1, fontface="bold"
    )
  return(plot)
}