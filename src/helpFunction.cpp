#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
DataFrame bivariate_hazards(DataFrame df, DataFrame data) {
  int n = df.nrows();
  
  NumericVector lambda10(n);
  NumericVector lambda01(n);
  NumericVector lambda11(n);
  NumericVector odds(n);
  NumericVector prodOdds(n);
  
  NumericVector df_t1 = df[0];
  NumericVector df_t2 = df[1];
  
  NumericVector data_t1 = data[0];
  NumericVector data_t2 = data[1];
  IntegerVector data_delta1 = data[2];
  IntegerVector data_delta2 = data[3];
  
  for (int i = 0; i < n; ++i) {
    double u = df_t1[i];
    double v = df_t2[i];
    
    int y = 0;
    int dj1 = 0;
    int dj2 = 0;
    int d1 = 0;
    
    for (int j = 0; j < data.nrows(); ++j) {
      double t1 = data_t1[j];
      double t2 = data_t2[j];
      int delta1 = data_delta1[j];
      int delta2 = data_delta2[j];
      
      if (t1 >= u && t2 >= v) {
        y++;
        if (t1 == u && t2 >= v && delta1 == 1) {
          dj1++;
        }
        if (t2 == v && t1 >= u && delta2 == 1) {
          dj2++;
        }
        if (t1 == u && t2 == v && delta1 == 1 && delta2 == 1) {
          d1++;
        }
      }
    }
    
    lambda10[i] = (y > 0) ? static_cast<double>(dj1) / y : 0.0;
    lambda01[i] = (y > 0) ? static_cast<double>(dj2) / y : 0.0;
    lambda11[i] = (y > 0) ? static_cast<double>(d1) / y : 0.0;
    
    // Calculate odds
    odds[i] = (lambda10[i] * lambda01[i] - lambda11[i]) / ((1 - lambda10[i]) * (1 - lambda01[i]));
    odds[i] = 1 - odds[i];
  }
  
  for (int i = 0; i < n; ++i) {
    double u = df_t1[i];
    double v = df_t2[i];
    double productOdds = 1.0;
    
    for (int j = 0; j < n; ++j) {
      if (df_t1[j] <= u && df_t2[j] <= v) {
        productOdds *= odds[j];
      }
    }
    
    prodOdds[i] = productOdds;
  }
  
  DataFrame result = DataFrame::create(
    Named("lambda.10") = lambda10,
    Named("lambda.01") = lambda01,
    Named("lambda.11") = lambda11,
    Named("odds") = odds,
    Named("prod.odds") = prodOdds
  );
  
  return result;
}


// [[Rcpp::export]]
DataFrame trivariate_hazards(DataFrame df, DataFrame data) {
  int n = df.nrows();
  
  NumericVector lambda100(n);
  NumericVector lambda010(n);
  NumericVector lambda001(n);

  NumericVector lambda110(n);
  NumericVector lambda101(n);
  NumericVector lambda011(n);
  
  NumericVector lambda111(n);
  
  NumericVector pfija(n);
  NumericVector p1(n);
  NumericVector p2(n);
  NumericVector odds(n);
  NumericVector prodOdds(n);
  
  NumericVector df_t1 = df[4];
  NumericVector df_t2 = df[2];
  NumericVector df_t3 = df[0];
  
  NumericVector data_t1 = data["t1"];
  NumericVector data_t2 = data["t2"];
  NumericVector data_t3 = data["t3"];
  IntegerVector data_delta1 = data["delta1"];
  IntegerVector data_delta2 = data["delta2"];
  IntegerVector data_delta3 = data["delta3"];
  
  for (int i = 0; i < n; ++i) {
    double u = df_t1[i];
    double v = df_t2[i];
    double k = df_t3[i];
    
    int y = 0;
    int d1 = 0;
    int d2 = 0;
    int d3 = 0;
    int d1_double = 0;
    int d2_double = 0;
    int d3_double = 0;
    int d_triple = 0;
  
    for (int j = 0; j < data.nrows(); ++j) {
      double t1 = data_t1[j];
      double t2 = data_t2[j];
      double t3 = data_t3[j];
      int delta1 = data_delta1[j];
      int delta2 = data_delta2[j];
      int delta3 = data_delta3[j];

      if (t1 >= u && t2 >= v && t3 >= k) {
        y++;
        if (t1 == u && t2 >= v && t3 >= k && delta1 == 1) {
          d1++;
        }
        if (t1 >= u && t2 == v && t3 >= k && delta2 == 1) {
          d2++;
        }
        if (t1 >= u && t2 >= v && t3 == k && delta3 == 1) {
          d3++;
        }
        if (t1 == u && t2 == v && t3 >= k && delta1 == 1 && delta2 == 1) {
          d1_double++;
        }
        if (t1 == u && t2 >= v && t3 == k && delta1 == 1 && delta3 == 1) {
          d2_double++;
        }
        if (t1 >= u && t2 == v && t3 == k && delta2 == 1 && delta3 == 1) {
          d3_double++;
        }
        if (t1 == u && t2 == v && t3 == k && delta1 == 1 && delta2 == 1 && delta3 == 1) {
          d_triple++;
        }
      }
    }
    
    
    lambda100[i] = (y > 0) ? static_cast<double>(d1) / y : 0.0;
    lambda010[i] = (y > 0) ? static_cast<double>(d2) / y : 0.0;
    lambda001[i] = (y > 0) ? static_cast<double>(d3) / y : 0.0;
    
    lambda110[i] = (y > 0) ? static_cast<double>(d1_double) / y : 0.0;
    lambda101[i] = (y > 0) ? static_cast<double>(d2_double) / y : 0.0;
    lambda011[i] = (y > 0) ? static_cast<double>(d3_double) / y : 0.0;
    
    lambda111[i] = (y > 0) ? static_cast<double>(d_triple) / y : 0.0;
    
    
    // Calculate odds
    pfija[i] = 1 - lambda100[i] - lambda010[i] - lambda001[i] + lambda110[i] + lambda101[i] + lambda011[i] - lambda111[i];
    p1[i] = (1 - lambda100[i]) * (1 - lambda010[i]) * (1 - lambda001[i]);
    p2[i] = (1 - lambda100[i] - lambda010[i] + lambda110[i]) * (1 - lambda100[i] - lambda001[i] + lambda101[i]) * (1 - lambda010[i] - lambda001[i] + lambda011[i]);
    odds[i] = (p1[i] / p2[i]) * pfija[i];
    
  }
  
  for (int i = 0; i < n; ++i) {
    double u = df_t1[i];
    double v = df_t2[i];
    double k = df_t3[i];
    double productOdds = 1.0;
    
    for (int j = 0; j < n; ++j) {
      if (df_t1[j] <= u && df_t2[j] <= v && df_t3[j] <= k) {
        productOdds *= odds[j];
      }
    }
    
    prodOdds[i] = productOdds;
  }

  DataFrame result = DataFrame::create(
    Named("lambda.100") = lambda100,
    Named("lambda.010") = lambda010,
    Named("lambda.001") = lambda001,
    Named("lambda.110") = lambda110,
    Named("lambda.101") = lambda101,
    Named("lambda.011") = lambda011,
    Named("lambda.111") = lambda111,
    Named("prod.odds") = prodOdds
  );
 
  return result;
}


