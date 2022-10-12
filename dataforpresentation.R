require(DescTools)
set.seed(37)

# Code for for 3Min UROP Presentation 12/10/22
# Loki Dunn


# ------
# Implimentation of Lee's 1997 algorithm

getmat <- function(q1, q2){
  # a function which produces matrix for use by categoricalgen
  # based on nature of maginals q1, q2 parsed to it
  n1 <- length(q1)
  n2 <- length(q2)
  
  # Defining matrix A describing relations imposed by marginals
  
  A <- matrix(rep(0, (n1 + n2 - 1)*(n1*n2)), ncol = n1*n2)
  
  A[1,] <- rep(1, n1*n2)
  
  for(i in 1:(n1-1)){ # fill rows 2 to n1 
    A[i+1,(((i-1)*n2)+1):(i*n2)] <- rep(1,n2)
  }
  for(i in 1:(n2-1)){
    for(j in 0:(n1-1)){
      A[n1+i,i+j*n2] <- 1
    }
  }
  
  return(A)
}

maxvertex <- function(q1, q2, errs = FALSE){
  # A function accepting marginal vectors q1 and q2 which finds the table
  # with these marginals with a maximal Uncertainty Coefficient
  # This is done via an exhaustive vertex search
  best <- 0
  stored <- 0
  
  n1 <- length(q1)
  n2 <- length(q2)
  
  A <- getmat(q1,q2) # Obtaining A via setup func
  b <- c(1,q1[-n1],q2[-n2])
  
  # Iterating across all maximal square submatrices of A
  combos <- combn(n1*n2, (n1*n2) - (n1 + n2 - 1), simplify = FALSE) # as a list
  for(cols in combos){
    A_s <- A[,-cols] # excluding columns according to x
    if(det(A_s) != 0){
      tempx <- solve(A_s) %*% b # solving for solutions in remaining columns
      if(all(tempx >= 0)){ # check solution non-negative
        
        x <- rep(0, n1*n2)
        pos <- 1
        # filling x with 0s and entries from tempx to complete found vector
        
        for(m in 1:(n1*n2)){
          if(m %in% cols){
            x[m] <- 0
          } else {
            x[m] <- tempx[pos]
            pos <- pos+1
          }
        }
        matx <- matrix(x, byrow = TRUE, ncol = n2)
        if(!is.nan(uncert(matx))){
          U <- uncert(matx)
          if(U > best){ # comparing against best so far
            best <- U
            stored <- matx
          }
        }
      }
    }
  }
  if(errs == TRUE){
    # optional error messages
    print(paste("Max U achieved:", best))
  }
  return(stored)
}

uncert <- function(probs, nudge = 1e-10){
  # returns uncertainty coefficient calculated exactly from a
  # table of probabilities given in matrix form
  n1 <- nrow(probs) # number of values in image of X1
  n2 <- ncol(probs) # " " " " " " X2
  
  entX1 <- 0 # iinitialising marginal and conditional
  entX2 <- 0 # information entropies to be calculated
  entX1X2 <- 0
  entX2X1 <- 0
  
  for(i in 1:n1){ # first finding entX1 and entX2X1
    Px1 <- sum(probs[i,]) # prob(X1 = i)
    entX1 <- entX1 - (Px1*log(Px1 + nudge))
    for(j in 1:n2){
      Px1x2 <- probs[i,j]
      entX2X1 <- entX2X1 - (Px1x2*log((Px1x2) + nudge))
    }}
  for(i in 1:n2){ # repeating to find entX2 and entX1X2
    Px2 <- sum(probs[,i])
    entX2 <- entX2 - (Px2*log(Px2 + nudge))
  }
  Ux1x2 <- 2*(entX1 + entX2 - entX2X1)/(entX1 + entX2)
  #Ux2x1 <- (entX2 + entX1 - entX2X1)/entX2
  return(Ux1x2)
}

TupleCode <- function(j,n1,n2){
  # returns tuple corresponding to row wise ordering
  # of elements in table
  x1 <- 1
  x2 <- 1
  
  while(j > n2){
    j <- j - n2
    x1 <- x1+1
  }
  x2 <- j
  
  return(c(x1,x2))
}

InvSample <- function(probs,N){
  # accepts probs a n1 by n2 matrix and N a number of samples to generate
  res <- data.frame(matrix(ncol = 2, nrow = 0))
  p <- as.numeric(t(probs))
  n1 = nrow(probs)
  n2 = ncol(probs)
  z <- c(0)
  for(i in 1:(n1*n2)){
    z <- append(z, z[length(z)] + p[i])
  }
  for(i in 1:N){
    obs <- runif(1)
    for(j in 1:(n1*n2)){
      if(obs >= z[j] & obs < z[j+1]){
        res[nrow(res)+1,] <- TupleCode(j, n1, n2)
      }
    }
  }
  return(res)
}

NomGen <- function(Pi1, Pi2, dU, N){
  # Pi1/2 - X1/2 marginal distribution vectors
  # dU - desired value of U
  # N - number of observations to generate
  n1 <- length(Pi1)
  n2 <- length(Pi2)
  
  indtab <- Pi1 %*% t(Pi2) # table corresponding to independence
  maxtab <- maxvertex(Pi1,Pi2) # maximal table
  
  if(dU == 0){ # in case of indepence
    return(InvSample(indtab, N))
  }
  
  f <- function(t){ # function sending t to U(t*ind + (1-t)*max)
    probs <- t*indtab + (1-t)*maxtab
    U <- uncert(probs)
    return(U - dU)
  }
  
  lambda <- unlist(uniroot(f, c(0,1))["root"])
  tab <- (lambda*indtab) + (1-lambda)*maxtab # table with desired tau value
  return(InvSample(tab, N))
}

countup <- function(data, n1, n2){
  # returns a matrix of counts for each tuple for a 2D categorical vector
  counts <- rep(0, n1*n2)
  for(i in 1:nrow(data)){
    x1 <- data[i,1]
    x2 <- data[i,2]
    j <- n2*(x1-1) + x2
    counts[j] <- counts[j] + 1 # increment relevant count
  }
  return(matrix(counts, byrow = TRUE, ncol = n2))
}

# ------

# DATA FOR PRESENTATION START

# wales data
A <- matrix(c(199,59,142,114,120,43,28,98,100), nrow = 3)
A <- A / 903

# U-value for data
U0 <- TheilUncertaintyFromProbs(A)

# sample distributions of X = destination and Y = origin
p1 <- colSums(A)
p2 <- rowSums(A)

# simulating (X,Y) data
generated <- NomGen(p1,p2,U0, 903)

# as a labelled data frame
data <- data.frame(data = countup(generated, 3, 3))
colnames(data) <- c("A","B","C")
rownames(data) <- c("Wales", "UK exc. Wales", "International")

# Uncertainty coefficient for simulated data
UncertCoef(generated[,1], generated[,2], direction = "symmetric")

# repeating to observe distribution of U-values
res <- c()
for(i in 1:1000){
  generated <- NomGen(p1,p2,U0, 903)
  res <- append(res, UncertCoef(generated[,1], generated[,2], direction = "symmetric"))
}

hist(res, breaks = 12)
abline(v=U0,col="red",lwd=2)

# exporting as svg

# svg("IMAGE.svg", width = 327, height = 356)
# par(bg=NA)
# hist(res, breaks = 20, freq = FALSE, main = "Histogram of Observed U values", xlab = "U")
# abline(v=U0,col="red",lwd=2)
# mean(res)
# U0
# 100*(mean(res) - U0)/U0
# dev.off()



