# PREREQUISITES
require(DescTools)
require(lpSolve)
getmat <- function(q1, q2){
  # a function which produces matrix for use by catgen
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
          if(U >= 1){
            return(matx)
          }
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
NomGenTab <- function(Pi1, Pi2, dU){
  # Pi1/2 - X1/2 marginal distribution vectors
  # dU - desired value of U
  # N - number of observations to generate
  n1 <- length(Pi1)
  n2 <- length(Pi2)
  
  indtab <- Pi1 %*% t(Pi2) # table corresponding to independence
  maxtab <- maxvertex(Pi1,Pi2) # maximal table
  
  if(dU == 0){ # in case of indepence
    return(indtab)
  }
  
  f <- function(t){ # function sending t to U(t*ind + (1-t)*max)
    probs <- t*indtab + (1-t)*maxtab
    U <- uncert(probs)
    return(U - dU)
  }
  
  lambda <- unlist(uniroot(f, c(0,1))["root"])
  tab <- (lambda*indtab) + (1-lambda)*maxtab # table with desired tau value
  return(tab)
}
shift <- function(vec, N = 1){
  # shifts entries in a vector by N places to the right
  if(N == 0){
    return(vec)
  }
  n <- length(vec)
  v1 <- vec
  v2 <- rep(0,n)
  for(i in 1:N){
    v2 <- rep(0,n)
    v2[1] = v1[n]
    for(j in 1:(n-1)){
      v2[j+1] <- v1[j]
    }
    v1 <- v2
  }
  return(v2)
}

# NomGen3Var is a function which accepts 3 marginal distributions,
# 2 values for association and a number of observations

# and generates observations of the edsired vector (X,Y,Z)

# START

# a function to obtain matrix A used in linear programming
get3VarMat <- function(n1,n2,n3){
  # instantiating empty matrix of correct dimensions
  A <- matrix(ncol = n1*n2*n3, nrow = (n1 + n2 + n3 - 2) + (n1-1)*(n2-1) + (n1-1)*(n3-1))
  
  # filling entries for each section of RHS vector
  A[1,] <- rep(1,ncol(A))
  for(i in 1:(n1-1)){ # first marginal
    A[1+i,] <- shift(rep((c(1, rep(0, n1-1))),n2*n3), (i-1))
  }
  for(i in 1:(n2-1)){ # second marginal
    A[n1+i,] <- shift(rep(c(rep(1,n1), rep(rep(0,n1), n2-1)), n3), n1*(i-1))
  }
  for(i in 1:(n3-1)){ # third marginal
    A[n1 + n2 -1 + i,] <- shift(c(rep(1,n1*n2), rep(0, n1*n2*(n3-1))), n1*n2*(i-1))
  }
  for(i in 1:((n1-1)*(n2-1))){ # first table
    A[n1 + n2 + n3 - 2 + i, ] <- rep(0,n1*n2*n3) # first fill row with zeros
    val <- i
    xi <- 0 # then find which probability row corresponds to
    xj <- 1
    while(val > (n1-1)){
      val <- val-(n1-1)
      xj <- xj + 1
    }
    xi <- val
    for(l in 0:(n3-1)){
      A[n1 + n2 + n3 - 2 + i, xi + n1*(xj-1) + l*(n1*n2)] <- 1
    }
  }
  for(i in 1:((n1-1)*(n3-1))){ # second
    A[n1 + n2 + n3 - 2 + i + (n1-1)*(n2-1), ] <- rep(0,n1*n2*n3) # first fill row with zeros
    val <- i
    xi <- 0 # then find which probability row corresponds to
    xj <- 1
    while(val > (n1-1)){
      val <- val-(n1-1)
      xj <- xj + 1
    }
    xi <- val
    for(l in 0:(n2-1)){
      A[n1 + n2 + n3 - 2 + i + (n1-1)*(n2-1), xi + (n1*n2)*(xj-1) + l*(n1)] <- 1
    }
  }
  return(A)
}

# a function which samples N observations from a 3D table
Var3Sample <- function(probs, n1, n2, n3, N){
  # accepts properly ordered vector of probabilities and
  # returns observations row-wise as a dataframe
  
  # get list of samples as corresponding place in prob ordering
  codes <- sample(1:length(probs), size = N, prob = probs, replace = TRUE)
  
  res <- data.frame(matrix(ncol = 3, nrow = 0))
  
  for(i in 1:length(codes)){
    val <- codes[i]
    tup <- c(0,1,1) # initialising tuple
    while(val > n1*n2){
      val <- val - n1*n2
      tup[3] <- tup[3] + 1
    }
    while(val > n1){
      val <- val - n1
      tup[2] <- tup[2] +1
    }
    tup[1] <- val
    res[nrow(res)+1,] <- tup
  }
  
  return(res)
}

# engine function
NomGen3Var <- function(q1,q2,q3, U,N, givetable = FALSE){
  # accepts 3 marginal distributions as vectors
  # as well as U a vector (U(X1,X2),U(X1,X3))
  # outputs a 3d table of probabilities
  n1 <- length(q1)
  n2 <- length(q2)
  n3 <- length(q3)
  
  Tab1 <- NomGenTab(q1, q2, U[1])
  Tab2 <- NomGenTab(q1, q3, U[2])
  
  A <- get3VarMat(n1,n2,n3) # matrix describing relations on table entries
  
  b <- c(1, q1[-n1], q2[-n2], q3[-n3], as.numeric(Tab1[-n1,-n2]),as.numeric(Tab2[-n1,-n3]))
  
  result <- lp("max", rep(1,n1*n2*n3), A, rep("=", length(b)), b)
  if(result["status"] == 2){
    stop("Correlation values not coherent")
  }
  x <- as.numeric(unlist(result["solution"]))
  
  res <- list()
  for(i in 1:n3){
    level <- matrix(x[(1+(i-1)*(n1*n2)):(i*(n1*n2))], nrow = n1, ncol = n2)
    res[[i]] <- level
  }
  if(givetable == TRUE){
    print(res)
  }
  
  return(Var3Sample(x, n1, n2, n3, N))
  
}


