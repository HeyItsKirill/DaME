#library(dplyr)
#library(Rcpp)
#library(prodlim)
#sourceCpp("helpFunction.cpp")

#' @title Multivariate Dabrowska estimator
#' 
#' @description Calculates the k-variate survival function. 
#' 
#' @param data An n x 2k data frame with the following items:
#' \describe{
#' \item{X1, X2, ..., Xk:}{Survival times for k dimensions.}
#' \item{Delta1, Delta2, ..., Deltak:}{Censoring indicators for k dimensions.}
#' }  
#' @param k The number of dimensions.
#' @param ... Additional arguments to be passed to the function. 
#' 
#' @return A data frame containing univariate, bivariate, and trivariate survival statistics.
#' 
#' @examples 
#' data <- dabrowska(genClaytonk(10, 1, c(0.1, 0.2, 0.3)), 3)
#' 
dabrowska = function(data, k,...){
  
  if(k == 2){
    colnames(data) = c("t1","t2","delta1","delta2")
  }else if(k == 3){
    colnames(data) = c("t1","t2","t3","delta1","delta2","delta3")
  }
  
  ################### unique failure times ###############
  
  
  bivariate = function(data2){
    
    colnames(data2) = c("t1","t2","delta1","delta2")
    
    
    t1.u = sort(unique(c(0,data2[data2$delta1 == 1,]$t1)))
    t2.u = sort(unique(c(0,data2[data2$delta2 == 1,]$t2)))
    
    
    df2 = data.frame(
      expand.grid(t1 = t1.u, t2 = t2.u)
    )
    
    df2 <- cbind(df2,bivariate_hazards(df2,data2))
    
    ################### univariate estimators ###############
    
    t1.km = rev(sort(unique(c(1,prodlim(Hist(t1,delta1)~1, data = data2[,c("t1","delta1")])$surv))))
    
    if(data2$delta1[data2$t1 == max(t1.u)] == 0){
      t1.km = t1.km[-length(t2.km)]
    }
    
    
    t2.km = rev(sort(unique(c(1,prodlim(Hist(t2,delta2)~1, data = data2[,c("t2","delta2")])$surv))))
    if(data2$delta2[data2$t2 == max(t2.u)] == 0){
      t2.km = t2.km[-length(t2.km)]
    }
    
    df2 <- 
      bind_cols(
        t1 = t1.u,
        t1.km = t1.km
      ) |>
      right_join(
        df2, by = "t1"
      )
    
    
    df2 <- 
      bind_cols(
        t2 = t2.u,
        t2.km = t2.km
      ) |>
      right_join(
        df2, by = "t2"
      )
    
    
    ################### bivariate estimator ###############
    
    df2$s.hat = df2$t1.km*df2$t2.km*df2$prod.odds
    
    ################### return ###############
    
    return(df2[,c("t1","t2","prod.odds","s.hat")])
  }
  
  
  if(k == 2){
    
    bivariate(data)
    
  }else if(k == 3){
    
    t1.u = sort(unique(c(0,data[data$delta1 == 1,]$t1)))
    t2.u = sort(unique(c(0,data[data$delta2 == 1,]$t2)))
    t3.u = sort(unique(c(0,data[data$delta3 == 1,]$t3)))
    
    df = data.frame(
      expand.grid(t1 = t1.u, t2 = t2.u, t3 = t3.u)
    )
    
    
    #univariate estimators
    
    t1.km = rev(sort(unique(c(1,prodlim(Hist(t1,delta1)~1, data = data[,c("t1","delta1")])$surv))))
    
    if(data$delta1[data$t1 == max(t1.u)] == 0){
      t1.km = t1.km[-length(t2.km)]
    }
    
    
    t2.km = rev(sort(unique(c(1,prodlim(Hist(t2,delta2)~1, data = data[,c("t2","delta2")])$surv))))
    if(data$delta2[data$t2 == max(t2.u)] == 0){
      t2.km = t2.km[-length(t2.km)]
    }
    
    t3.km = rev(sort(unique(c(1,prodlim(Hist(t3,delta3)~1, data = data[,c("t3","delta3")])$surv))))
    if(data$delta3[data$t3 == max(t3.u)] == 0){
      t3.km = t3.km[-length(t3.km)]
    }
    
    df <- 
      bind_cols(
        t1 = t1.u,
        t1.km = t1.km
      ) |>
      right_join(
        df, by = "t1"
      )
    
    
    df <- 
      bind_cols(
        t2 = t2.u,
        t2.km = t2.km
      ) |>
      right_join(
        df, by = "t2"
      )
    
    
    df <- 
      bind_cols(
        t3 = t3.u,
        t3.km = t3.km
      ) |>
      right_join(
        df, by = "t3"
      )
    
    
    #Hazards/odds measure and bivariate KM
    df <- cbind(df,trivariate_hazards(df,data))
    
    shat.110 = bivariate(data2 = data[,c("t1","t2","delta1","delta2")])  |> rename("prod.odds.12" = "prod.odds","s.hat.12" = "s.hat")
    shat.101 = bivariate(data2 = data[,c("t1","t3","delta1","delta3")])  |> rename("t3" = "t2","prod.odds.13" = "prod.odds","s.hat.13" = "s.hat")
    shat.011 = bivariate(data2 = data[,c("t2","t3","delta2","delta3")])  |> rename("t2" = "t1","t3" = "t2","prod.odds.23" = "prod.odds","s.hat.23" = "s.hat")
    
    
    df <- df |>
      left_join(shat.110, by = c("t1" = "t1","t2" = "t2")) |>
      left_join(shat.101, by = c("t1" = "t1","t3" = "t3")) |>
      left_join(shat.011, by = c("t2" = "t2","t3" = "t3")) 
    
    
    df[is.na(df)] = 0
    
    
    ################### trivariate estimator ###############
    
    df$s.hat = df$t1.km*df$t2.km*df$t3.km*df$prod.odds*df$prod.odds.12*df$prod.odds.13*df$prod.odds.23*df$prod.odds
    
    ################### return ###############
    
    
    return(df)
    
  }
  
}
