#function to test how mu_K varies between 2 and 100 species with 
#different values of mu_K
#mu is a contant, 2 
#N is a series, 2:100
#K is a constant, 40 or 80
#mu_K is a vector of length 7, c(0.1, 0.5, 1, 2, 3, 4, 5, 7.5, 10, 20, 50)

test_muK <- function(mu, N, K, mu_K){
 
  table <- matrix(nrow = 100, ncol = 13)  
  table[1,1] <- "N" 
  table[1,2:12] <- c("mu_K = 0.1", "mu_K = 0.5",  "mu_K = 1", 
                     "mu_K = 2", "mu_K = 3", "mu_K = 4", "mu_K = 5",
                    "mu_K = 7.5,", "mu_K = 10", "mu_K = 20", "mu_K = 50")
  table[1,13] <- "original"
  table[2:100,1] <- N
  
  for (i in N)
{
  ext_rate <- list()
  ext_rate[[i]] <- max(c(mu * (mu_K[1]/mu)^(i/K)),0,na.rm = T)
  table[i,2]<-ext_rate[[i]]
}
    for (i in N)
    {
      ext_rate <- list()
      ext_rate[[i]] <- max(c(mu * (mu_K[2]/mu)^(i/K)),0,na.rm = T)
      table[i,3]<-ext_rate[[i]]
    }
  for (i in N)
  {
    ext_rate <- list()
    ext_rate[[i]] <- max(c(mu * (mu_K[3]/mu)^(i/K)),0,na.rm = T)
    table[i,4]<-ext_rate[[i]]
  }
  for (i in N)
  {
    ext_rate <- list()
    ext_rate[[i]] <- max(c(mu * (mu_K[4]/mu)^(i/K)),0,na.rm = T)
    table[i,5]<-ext_rate[[i]]
  }
  for (i in N)
  {
    ext_rate <- list()
    ext_rate[[i]] <- max(c(mu * (mu_K[5]/mu)^(i/K)),0,na.rm = T)
    table[i,6]<-ext_rate[[i]]
  }
  for (i in N)
  {
    ext_rate <- list()
    ext_rate[[i]] <- max(c(mu * (mu_K[6]/mu)^(i/K)),0,na.rm = T)
    table[i,7]<-ext_rate[[i]]
  }
  for (i in N)
  {
    ext_rate <- list()
    ext_rate[[i]] <- max(c(mu * (mu_K[7]/mu)^(i/K)),0,na.rm = T)
    table[i,8]<-ext_rate[[i]]
  }
  for (i in N)
  {
    ext_rate <- list()
    ext_rate[[i]] <- max(c(mu * (mu_K[8]/mu)^(i/K)),0,na.rm = T)
    table[i,9]<-ext_rate[[i]]
  }
  for (i in N)
  {
    ext_rate <- list()
    ext_rate[[i]] <- max(c(mu * (mu_K[9]/mu)^(i/K)),0,na.rm = T)
    table[i,10]<-ext_rate[[i]]
  }
  for (i in N)
  {
    ext_rate <- list()
    ext_rate[[i]] <- max(c(mu * (mu_K[10]/mu)^(i/K)),0,na.rm = T)
    table[i,11]<-ext_rate[[i]]
  }
  for (i in N)
  {
    ext_rate <- list()
    ext_rate[[i]] <- max(c(mu * (mu_K[11]/mu)^(i/K)),0,na.rm = T)
    table[i,12]<-ext_rate[[i]]
  }
  for (i in N)
  {
    ext_rate <- list()
    ext_rate[[i]] <- max(c(mu * i),0,na.rm = T)
    table[i,13]<-ext_rate[[i]]
  }
  print(table)
  plot(table[,1],table[,13],col='darkred')
  points(table[,1],table[,6],col='red')
  points(table[,1],table[,7],col='darkorange')
  points(table[,1],table[,8],col='orange')
  points(table[,1],table[,9],col='yellow')
  points(table[,1],table[,10],col='green')
  points(table[,1],table[,11],col='blue')
  points(table[,1],table[,12],col='darkblue')
}

#mu is a contant, 2 
#NK is the fraction N/K is a series, NK <- seq(from = 0.05, to = 5, by = 0.1)
#mu_K is a vector of length 3, c(2.5, 3, 4)
test_ratioNK <- function(mu, N, K, mu_K){
  
  table <- matrix(nrow = 51, ncol = 4)  
  table[1,1] <- "NK" 
  table[1,2:4] <- c("mu_K = 2.5", "mu_K = 3",  "mu_K = 4")
  table[2:51,1] <- NK
  
  for (i in NK)
  {
    ext_rate <- list()
    ext_rate[[i]] <- max(c(mu * (mu_K[1]/mu)^i),0,na.rm = T)
    table[i,2]<-ext_rate[[i]]
  }
  
}
