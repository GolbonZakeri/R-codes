#Craig Brinkerhoff
#MIE 620 Revised Simplex Method Programming Project
#11 December 2019

#main functions for part c----------------------------------------------------
fullrsm <- function(m,n,costVec,A,b){
  # Solves a linear program using Gauss-Jordon updates
  # Assumes standard computational form
  # Performs a Phase I procedure starting from an artificial basis
  # Input:
    # m,n = number of constraints and variables
    # costVec = nx1 cost vector
    # A = mxn constraint matrix
    # b = mx1 rhs vector
    # Output:
      # result = 1 if problem optimal, 0 if infeasible, -1 if unbounded
      # z = objective function value
      # xB = nx1 solution vector
      # pi = mx1 dual vector
      # Status = 1xm vector of what variables are in the optimal bfs and which are set to 0
  
  #set up phase I implicitly---------------------------------------------------------------
  intialProblem <- phaseIsetup(m, n, A)
  
  constraints <- intialProblem$Constraints
  cost <- intialProblem$Cost
  basic <- intialProblem$Basic
  status <- intialProblem$Status
  B <- intialProblem$B
  Binv <- intialProblem$Binv
  
  phase1 <- 1
  
  #Start simplex---------------------------------------------------------------------------
  iterator <- 2
  while(iterator == 2){
    duals <- findDuals(cost, basic, Binv)
    enterVar <-fullfindEV(n, cost, constraints, status, duals, phase1) #get entering variable
    s <- enterVar$s
    Cs = cost[s]
    piAs <- duals %*% constraints[,s]

    #Optimality test
    if((Cs - piAs) >= 0) {
      print('Solution is optimal')
      xB <- findxB(Binv, b)
      
      z <- cbind(xB, basic)
      z <- sum(cost[z[,2]] * z[,1])
      
      output <- list('Result'= 1, 'Z'= z, 'BFS'= xB, 'Duals'= duals, 'Status' = status)
      
      #check infeasiblity of phase I problem
      if (phase1 == 1){
        if (output$`Z` > 0){
          print(paste("Problem is infeasible with Z = ", output$`Z`, sep = ''))
          output <- list('Result'= 0, 'Z'= output$`Z`, 'BFS'= NA, 'Duals'= NA, 'Status' = NA)
          return(output) #return incorrect solution as problem is infeasible
        }
        
        #if feasible, intialize phase II
        else{
          print(paste('Problem is feasible with Z = ', output$`Z`, '. PhaseII starting w/ basis: ', sep=''))
          print(basic)
          
          #update status to non-artifical problem
          status <- rep(0,ncol(A))
          for (i in seq_along(status)){
            status[i] <- ifelse(i %in% basic, which(basic == i), 0)
          }
         
          constraints = A #revert back to non-artifical problem
          cost = costVec #revert back to non-artifical problem
          
          iterator = 2
          phase1 = 0
        }
        next #begin phase II
      }
      
      #if phase II, return optimal solution
      print('Optimization Complete')
      return(output)
    }
    
    #if not optimal or unbounded, run another iteration of algorithm
    else{
      BinvAs <- findBinvAs(Binv, A, s)
      xB <- findxB(Binv, b)
      leavingVar <- fullfindLV(m, xB, BinvAs, phase1, basic) #get leaving variable

      #extended leaving variable criterion
      if (phase1 == 0 && any(basic > n)) {
        if (BinvAs != 0){
          leavingVar$minratio = 0 #if artifical in basis, set minratio to 0 and set artifical variable to leaving variable
          leavingVar$r <- which(basic > m)
        }
      }
      r <- leavingVar$r
      
      #check for unboundedness
      if (all(leavingVar$minratio == 0)) {
        print('Problem is unbounded')
        result <- -1
        output <- list('Result'= result, 'Z'= NA, 'BFS'= xB, 'Duals'= duals, 'Status' = status)
        return(output)
      }
      
      else{
      cB <- cost[status != 0]
      update <- fullupdate(m,cost,s,r,BinvAs, phase1, status,basic,cB,Binv,xB)
      basic <- update$basic
      status <- update$status
      Binv <- update$Binv
      xB <- xB
        
      print('Basic Variables: ')
      print(basic)
      print('Iteration complete')
      }
    }
  }
}

#generic functions---------------------------------------------------------
inverse <- function(m, x){
  #Inverts a matrix using pivoting operations
  # Input:
    # m = number of constraints
    # x = mxn matrix to invert
  # Output:
    # output = mxn matrix solution to inversion
  
  I <- diag(nrow(x))
  matrix <- cbind(x,I)
  
  matrix[1,] <- matrix[1,]/matrix[1,1]
  
  i <- 2
  while (i < m+1) {
    j <- i
    while (j < m+1) {
      matrix[j, ] <- matrix[j, ] - matrix[i-1, ] * matrix[j, i-1]
      j <- j+1
    }
    while (matrix[i,i] == 0) {
      matrix <- rbind(matrix[-i,],matrix[i,])
    }
    matrix[i,] <- matrix[i,]/matrix[i,i]
    i <- i+1
  }
  for (i in m:2){
    for (j in i:2-1) {
      matrix[j, ] <- matrix[j, ] - matrix[i, ] * matrix[j, i]
    }
  }
  drop <- 1:(ncol(matrix)-m)
  output <- matrix[,-drop]
  return(output)
}

findBinvAs <- function(Binv, A, s){
  #multiplies current inverse B matrix against current As vector
    # Input:
    # A = mxn constraint matrix
  # Output:
    # BinvAs = nx1 vector solution to multiplication
  As <- A[,s]
  
  output <- Binv %*% As
  return(output)
}

findDuals <- function(c, basicvars, invB) {
  #Calculates the dual vector
  # Input:
    # c = 1xm cost vector
    # basicvars = indices of variables in the basis
    # invB = mxn matrix of the inverse of the current basis
  # Output:
    # output = duals vector
  
  costB <- c[basicvars]
  pi <- t(costB) %*% invB
  
  return(pi)
}

findxB <- function(Binv, b){
  #Finds current bfs solution
  # Input:
    # b = nx1 solution vector
    # invB = mxn matrix of the inverse of the current basis
  # Output:
    # output = current bfs vector
  
  output <- Binv%*%b
  return(output)
}

phaseIsetup <- function(m, n, A){
  #sets up phase I problem given constraint matrix A for problem
  # Input:
    # m, n = number of constraints, variables
    # A = constraint matrix
  # Output:
    # Constraints = mxn constraint matrix
    # Costs = nx1 cost vector
    # Basic = 1xm vector of indices in A of basic variables
    # Status = 1xn vector of indices of basic varaibles and 0 for non-basic variables
    # B = constraint matrix for variables in basic
    # Binv = inverse matrix for B
  
  constraints <- cbind(A, diag(nrow(A)))
  cost <- c(rep(0, n), rep(1, m))
    
  basic <- rep(0, (n+m)-(n)) #intialized as n+1, ... n+m
  for (i in (n+1):(n+m)){
    basic[i] = i
  }
  basic <- basic[-c(1:n)]
    
  status <- rep(0, ncol(constraints))
  for (i in seq_along(status)){
    status[i] <- ifelse(i %in% basic, which(basic == i), 0)
  }
    
  B <- constraints[,basic]
  Binv <- inverse(m, B)
    
  output <- list('Constraints' = constraints, 'Cost' = cost, 'Basic' = basic, 'Status' = status, 'B' = B, 'Binv' = Binv)
  return(output)
  
  #start up simplex------------------------------------------
  print('Starting up Phase I...')
}

#necessary functions for part c-----------------------------------------------
fullfindEV <- function(n,cost, A, varstatus, pi, phaseI){
  # Returns the index of the entering variable and it's reduced cost,
  # or returns 0 if no entering variable exists
  # Per part C requirements, DOES NOT compute the reduced costs as a vector
  
  # Input:
    # n = number of variables
    # cost = nx1 cost vector
    # A = mxn constraint matrix
    # varstatus = 1xn vector, varstatus(i) = position in basis of variable i,
      # or 0 if variable i is nonbasic
    # pi = duals vector
    # phase1 = boolean, phase1 = true if Phase 1, or false otherwise
  # Output:
    # s = index of the entering variable
    # minrc = reduced cost of the entering variable

  cN <- cost[varstatus == 0] #N cost vector
  N <- A[,varstatus == 0] #current N matrix
  temp <- (t(cN) - (pi%*%N))

  indexN <- which((t(cN) - (pi%*%N))==min((t(cN) - (pi%*%N))))
  indexN <- indexN[length(indexN)] #in case of tie of reduced cost, take highest coefficient variable as tie breaker
  #get s as minrc index in A
  ticker <- 0
  for(i in seq_along(varstatus)){
    if (varstatus[i] == 0) {
      ticker = ticker + 1
      if (ticker == indexN) {
        s <- i
      }
    }
  }
  
  #Identify if there are multiple optimal solutions.  Simplex will still find an optimal extreme point.
    #NOTE: This will always be triggered in phase 1.
  if (any((t(cN) - (pi%*%N))==0)){
    print('Infinite optimal solutions exist. Solver will find one.')
  }

  #return 0 if no entering variable exists
  if (all((t(cN) - (pi%*%N)) >= 0) == 1) {
    output <- list('s' = s, 'minrc' = 0)
    return(output)
  }
  else {
    output <- list('s' = s, 'minrc' = min((t(cN) - (pi%*%N))))
    return(output)
  }
}

fullfindLV <- function(m, xB, BinvAs, phase1, basic){
  #Returns the position in the basis of the leaving variable,
    # or returns 0 if no leaving variable exists
    # can implicitly simulate presence of artifical variables via the extended leaving variable criterion
  
  # Input:
    # m = number of constraints
    # xB = mx1 basic variable vector
    # BinvAs = mx1 vector of Binv*As
    # phase1 = boolean, phase1 = true if Phase 1, or false otherwise
    # basic = 1xm vector of indices of basic variables
    # Output:
      # r = position in the basis of the leaving variable
      # minratio = minimum ratio from ratio test
  
  finite <- is.finite(xB/BinvAs)
  minratio <- ifelse(finite == 1, xB/BinvAs, xB) #if equals infinity (i.e. number/0), set to number
  minratio <- ifelse(BinvAs <= 0, -9998, minratio) #if BinvAs is <= 0, ignore
  if (all(minratio == -9998)) { #if no leaving variable exists, set minratio to 0
    output <- list('r' = NA, 'minratio' = 0)
    return (output)
  }
  minimum <- min(minratio[minratio != -9998])
  r <- which(minratio == minimum)
  
  output <- list('r' = r, 'minratio' = minratio)
  return(output)
}

fullupdate <- function(m,cost,s,r,BinvAs,phase1,status,basic,cB,Binv,xB) {
  # Updates the basis representation using Gauss-Jordan Pivoting
  
  # Input:
    # m = number of constraints
    # cost = nx1 cost vector
    # s = index of entering variable
    # r = position in the basis of the leaving variable
    # BinvAs = mx1 Binv*As vector
    # phase1 = boolean, phase1 = true if Phase 1, or false otherwise
    # status = 1xn vector, status(i) = position in basis of variable i,
      # or 0 if variable i is nonbasic
    # basic = 1xm vector of indices of basic variables
    # cB = mx1 basic cost vector
    # Binv = mxm basis inverse matrix
    # xB = mx1 basic variable vector
    # Output:
      # status = 1xn updated status vector
      # basic = 1xm updated basic vector
      # cB = mx1 updated basic cost vector
      # Binv = mxm updated basis inverse matrix
      # xB = mx1 updated basic variable vector
  
  #form temporary matrix x for pivoting
  x = cbind(xB, Binv, BinvAs)
  pivot <- ncol(x) #pivot column
  
  #check if trying to pivot on a zero
  if (x[r,pivot] == 0) {
    stop(print("Can't pivot on a zero!  Something's wrong!"))
  }
  
  #Gauss-Jordan Pivot operation
  x[r,] = x[r,]/x[r,pivot]#perform pivot on r row
  for (i in 1:nrow(x)){
    if (i == r) {next}
    x[i,] = x[i,] - x[i,pivot] * x[r,]/x[r,pivot]
  }

  #extract new inverse basis B
  newBinv <- x[,-1]
  newBinv <- newBinv[,-ncol(newBinv)]
  
  #extract new Xb and basic variables
  newxB <- x[,1]
  newBasicvars <- basic
  newBasicvars[r] = s
  
  #update status vector
  newVarstatus <- status
  for (i in seq_along(newVarstatus)){
    newVarstatus[i] <- ifelse(i %in% newBasicvars, which(newBasicvars == i), 0)
  }
  
  #extract new basic cost vector
  newcB <- cost[newBasicvars]
  
  output <- list('status' = newVarstatus, 'basic' = newBasicvars, 'cB' = newcB, 'Binv' = newBinv, 'xB' = newxB)
  
  #error handling
  if (any(output$`xB` < 0)) {
    print('BFS can not be negative! Something is wrong.')
    return(NA)
  }

  else {
      return(output)
    }
}

#run simplex-----------------------------------------------------------
#Toy problems--------------------------------------------------------------------------
#Optimal
problem1 <- list('inputcost' = c(-25,-20,0,0), 
                 'constraints' =  rbind(c(20,12,1,0),c(5,5,0,1)), 
                 'b' = c(2000, 540), 
                 'm' = 2, 
                 'n' = 4)

problem2 <- list('inputcost' = c(1,2,3,0,0,0), 
                 'constraints' =  rbind(c(1,-1,1,1,0,0), c(1,1,0,0,-1,0), c(0,0,1,0,0,1)), 
                 'b' = c(4, 0, 6), 
                 'm' = 3, 
                 'n' = 6)

problem3 <- list('inputcost' = c(-3,-2,0,0), 
                 'constraints' =  rbind(c(2,3,1,0),c(2,1,0,1)), 
                 'b' = c(12,8), 
                 'm' = 2, 
                 'n' = 4)

#infinte/multiple optimal solutions
problem4 <- list('inputcost' = c(-1,-2,0,0), 
                 'constraints' =  rbind(c(2,4,1,0),c(4,3,0,1)), 
                 'b' = c(12,16), 
                 'm' = 2, 
                 'n' = 4)

#infeasible
problem5 <-  list('inputcost' = c(-1,1,-1), 
                  'constraints' =  rbind(c(2,-1,-2), c(-2,3,1), c(1,-1,-1)), 
                  'b' = c(4, 5, 1), 
                  'm' = 3, 
                  'n' = 3)

#unbounded
problem6 <- list('inputcost' = c(-2,-1,0,0), 
                 'constraints' =  rbind(c(1,-1,1,0), c(2,-1,0,1)), 
                 'b' = c(10,40), 
                 'm' = 2, 
                 'n' = 4)

problem7 <- list('inputcost' = c(-2,-1), 
                  'constraints' =  rbind(c(1,-1), c(2,-1)), 
                  'b' = c(10,40), 
                  'm' = 2, 
                  'n' = 2)

resultFULL <- fullrsm(problem4$m,problem4$n,problem4$inputcost, problem4$constraints, problem4$b)