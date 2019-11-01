plants <- c(49, -67, 8, 16, 6, 23, 28, 41, 14, 29, 56, 24, 75, 60, -48) 
func_simvals <- function(rep, data){
  ordered <- sort(data)
  n <- length(data)
  df <- n - 1
  order.i <- c()
  for (i in 1:n){
    order.i <- cbind(order.i, ordered[-i])
  }
  xbar.i <- apply(order.i, 2, mean)
  sd.i <- apply(order.i, 2, sd)
  simval <-  matrix(rep(0, rep*n), nrow=rep, ncol=n)
  for (i in 1:n){ 
    bigmatrix<- matrix(rnorm(n*rep, mean=xbar.i[i], sd=sd.i[i]), ncol=rep, nrow=n) ### change to mvt 
    vec <- apply(bigmatrix, 2, sort)[i,]
    simval[,i] <- vec
  }
  est.simval <- apply(simval, 2, mean)
  est.simval
}


set.seed(030202) #Homeowrk 03 | Problem 02 | Part 02
plant.sim.10    <- func_simvals(10000, plants)
#plant.sim.100  <- func_simvals(100000, plants) #run before hand, takes about 2-3 min
plant.sim.100  <- c(-24.773063,
                    -16.368983,
                    -14.899603,
                    -5.963294,
                    1.274501,
                    8.127543, 
                    14.317594,
                    20.743906,
                    26.928987,
                    33.490793,
                    39.545905,
                    46.378990,
                    54.311824,
                    64.966700,
                    79.486059)
#plant.sim.1000 <- func_simvals(1000000, plants) #run before hand, takes about 13-14 min
plant.sim.1000 <- c(-24.761139,
                    -16.333066,
                    -14.903300, 
                    -6.034726,
                    1.264649,
                    8.135469,
                    14.297283,
                    20.722848,
                    26.900612,
                    33.455976,
                    39.492197,
                    46.333023,
                    54.318730,
                    64.966895,
                    79.474176)

data <- plants
func_simdensityplots <- function(data, i.th){
  ordered <- sort(data)
  n <- length(data)
  df <- n - 1
  order.i <-  ordered[-i.th]
  xbar.i <- mean(order.i)
  sd.i <- sd(order.i)
  simulations <- c()
  percent <- c()
  it.vals <- c(10000, 100000, 1000000)
  for (it in it.vals){
    it <- c(rep(it, it))
    simvals.it <- rnorm(n*it, mean=xbar.i, sd=sd.i)
    p.it <- length(simvals.it[simvals.it >= ordered[i.th]]) / length(simvals.it)
    simvals <- cbind(it, simvals.it)
    simulations <- rbind(simulations, simvals)
    percent <- rbind(percent, p.it)
  }
  simulations.df <- as.data.frame(simulations)
  
  dat_text <- data.frame(
    label = c(
      paste("P(X > i.th | Y.-i) =", decimal(percent[1], dec))
      ,paste("P(X > i.th | Y.-i) =", decimal(percent[2], dec))
      ,paste("P(X > i.th | Y.-i) =", decimal(percent[3], dec))
    ),
    it   = it.vals)
  
  plot <- ggplot(simulations.df, aes(x=simvals.it)) + geom_histogram(aes(y=..density..), color="grey", alpha=0.4,) +
    facet_wrap(~it, nrow=1) + 
    geom_vline(xintercept=ordered[i.th], color="red", linetype="dashed")+
    ggtitle(paste("Simulated Distribution based on Plant Differences, missing the", i.th, "value (", ordered[i.th], ")"))+
    xlab("") + #ylim(0, 0.015) + xlim(-150, 250)
    geom_text( data    = dat_text,  
               mapping = aes(x = -150, y = 0, label = label),
               hjust   = -0,vjust   = -1, fontface="bold")
  
  return(plot)
}