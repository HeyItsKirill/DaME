#include <Rcpp.h>
#include <map>
#include <string>
#include <vector>
using namespace Rcpp;
//This is a slice function for conducting column slice in Rcpp
DataFrame sliceDataFrame(DataFrame df, int start, int end) {
  if (start < 0) {
    start = 0;
  }
  if (end > df.size()) {
    end = df.size();
  }
  if (start >= end) {
    stop("start must be less than end");
  }
  List slices(end - start);
  for (int i = start; i < end; ++i) {
    slices[i - start] = df[i];
  }
  if (!Rf_isNull(df.attr("names"))) {
    CharacterVector names = df.names();
    slices.attr("names") = names[Range(start, end - 1)];
  }
  
  return DataFrame(slices);
}


//[[Rcpp::export]]
//Transform the shape of DataFrame for selecting data by rows 
std::vector<NumericVector> getRowsAsNumericVectors(DataFrame df){
  int nrows = df.nrows();
  int ncols = df.size();
  std::vector<NumericVector> rows(nrows);
  for (int i = 0; i < nrows; ++i){
    rows[i] = NumericVector(ncols);
    for (int j = 0; j < ncols; ++j){
      switch(TYPEOF(df[j])){
      case REALSXP: case INTSXP: {
        rows[i][j] = as<NumericVector>(df[j])[i];
        break;
      }
      default: {
        stop("All columns must be numeric or integer");
      }
      }
    }
  }
  return rows;
}


//Compare two numeric vectors if for each element, vec1 >= vec2, return True 
bool compareNumericVectors(NumericVector vec1, NumericVector vec2) {
  if (vec1.size() != vec2.size()) {
    stop("Vectors must be of the same length");
  }
  
  for (int i = 0; i < vec1.size(); ++i) {
    if (vec1[i] < vec2[i]) {
      return false;
    }
  }
  
  return true;  
}


// Compare wheater two NumircVector are equal
bool areVectorsEqual(NumericVector vec1, NumericVector vec2) {
  if (vec1.size() != vec2.size()) {
    return false;
  }
  for (int i = 0; i < vec1.size(); ++i) {
    if (vec1[i] != vec2[i]) {
      return false;
    }
  }
  
  return true;
}


// Merge two dataframe by columns
DataFrame mergeDataFrames(DataFrame df1, DataFrame df2) {
  if (df1.nrows() != df2.nrows()) {
    stop("DataFrames do not have the same number of rows");
  }
  int nCols = df1.size() + df2.size();
  List merged(nCols);
  for (int i = 0; i < df1.size(); ++i) {
    merged[i] = df1[i];
  }
  for (int i = 0; i < df2.size(); ++i) {
    merged[i + df1.size()] = df2[i];
  }
  return DataFrame::create(merged);
}


//Generate all possible (1, 0, ...,1, 0) for k 1s in n dimensional vector
std::vector<std::string> generateCombinationsAsString(int n, int k) {
  std::vector<int> v(n, 0);
  std::vector<std::string> allCombinations;
  std::fill(v.begin(), v.begin() + k, 1);
  do {
    std::string combination;
    for (int i : v) {
      combination += std::to_string(i);
    }
    allCombinations.push_back(combination);
  } while (std::prev_permutation(v.begin(), v.end()));
  
  return allCombinations;
}

//Transform binary string in to a vector that record loaction of 1
std::vector<int> binaryStringToVector(const std::string& binaryStr){
  std::vector<int> result;
  for (int i = 0; i < binaryStr.size(); ++i) {
    if (binaryStr[i] == '1') {
      result.push_back(i);
    }
  }
  return result;
}

//Based on the location generated above, select columns of given DataFrame and construct a sub DataFrame
DataFrame selectColumnsCpp(DataFrame df, std::vector<int> columns) {
  List selected(columns.size());
  for (size_t i = 0; i < columns.size(); ++i) {
    int colIndex = columns[i]; 
    if (colIndex < 0 || colIndex >= df.size()) {
      stop("Column index out of range");
    }
    selected[i] = df[colIndex];
  }
  return DataFrame::create(selected);
}


//similarly get sub vector given the location vector and a numeric vector
NumericVector getSubVector(NumericVector vec, const std::vector<int> location) {
  NumericVector subVec(location.size());
  for (size_t i = 0; i < location.size(); ++i) {
    int index = location[i] - 1;
    if (index >= 0 && index < vec.size()) {
      subVec[i] = vec[index];
    } else {
      stop("Index out of bounds");
    }
  }
  
  return subVec;
}


// Dabrowska Estimator for k dimension
DataFrame DabrowskaEstimator(DataFrame df, DataFrame data, int k){
  int n = df.nrows();
  NumericVector ones(k, 1.0);
  NumericVector DabrowskaEst(n);
  DataFrame df_t = sliceDataFrame(df, 0, k);
  DataFrame data_t = sliceDataFrame(data, 0, k);
  DataFrame df_delta = sliceDataFrame(df, df.ncol() - k + 1, df.ncol());
  DataFrame data_delta = sliceDataFrame(data, data.ncol() - k + 1, data.ncol());
  std::vector<NumericVector> row_df_t = getRowsAsNumericVectors(df_t);
  std::vector<NumericVector> row_data_t = getRowsAsNumericVectors(data_t);
  std::vector<NumericVector> row_df_delta = getRowsAsNumericVectors(df_delta);
  std::vector<NumericVector> row_data_delta = getRowsAsNumericVectors(data_delta);
  // Dimension 1 as a starter for recursion
  if (k==1){
    NumericVector lambda(n);
    for (int i = 0; i < n; ++i){
       NumericVector u = row_df_t[i];
       int y = 0;
       int d = 0;
       for(int j = 0; j < data.nrows(); ++j){
         NumericVector t = row_data_t[i];
         NumericVector delta = row_data_delta[i];
         if (compareNumericVectors(t, u)){
           y++;
           if(areVectorsEqual(delta, ones)){
             d++;
           }
         }
       }
       lambda[i] = (y > 0) ? static_cast<double>(d) / y : 0.0;
    }
    for(int i = 0; i < n; ++i){
      double survival = 1;
      NumericVector u = row_df_t[i];
      for(int j = 0; j < n; ++j){
        NumericVector v = row_df_t[j];
        if (compareNumericVectors(u, v)){
          survival *= (1 - lambda[j]);
        }
      }
      DabrowskaEst[i] = survival;
    }
    DataFrame result = DataFrame::create(
      Named("Dabrowska") = DabrowskaEst
    );
    return result;
  }
  
  else{
    //Using dfMap to store lower cases of Dabrowska Estimator and use binary string representation as keys of it
    std::map<std::string, Rcpp::DataFrame> dfMap;
    for (int i = 0; i < k; ++i){
      std::vector<std::string> combinations = generateCombinationsAsString(k, i);
      for (size_t j = 0; j < combinations.size(); ++j){
        std::vector<int> location = binaryStringToVector(combinations[j]);
        DataFrame sub_df_t = selectColumnsCpp(df_t, location);
        DataFrame sub_df_delta = selectColumnsCpp(df_delta, location);
        DataFrame sub_data_t = selectColumnsCpp(data_t, location);
        DataFrame sub_data_delta = selectColumnsCpp(data_delta, location);
        DataFrame sub_df = mergeDataFrames(sub_df_t, sub_df_delta);
        DataFrame sub_data = mergeDataFrames(sub_data_t, sub_data_delta);
        DataFrame sub_Dabrowska = DabrowskaEstimator(sub_df, sub_data, i);
        dfMap[combinations[j]] = mergeDataFrames(sub_df, sub_Dabrowska);
      }
    }
    //After store the lower cases for each point, extract lower cases terms to compute Dabrowska Estimators
    for (int i = 0; i < n; ++i){
      double prod = 1;
      //The computation part is unfinished. Ideally after this part, a DataFrame which documents k-dimensional estimator
      //and its time variables will be output
      NumericVector u = row_df_t[i];
        for (int x = 0; x < k; ++x){
          std::vector<std::string> combinations = generateCombinationsAsString(k, i);
          for (size_t y = 0; y < combinations.size(); ++k){
            std::vector<int> location = binaryStringToVector(combinations[y]);
            NumericVector sub_u = getSubVector(u, location);
            DataFrame comparing = dfMap[combinations[y]];
            DataFrame comparing_t = sliceDataFrame(df, 0, x);
            std::vector<NumericVector> row_comparing_t = getRowsAsNumericVectors(comparing_t);
            for (int j = 0; j < n; ++j){
              NumericVector v = row_comparing_t[j];
              if (compareNumericVectors(sub_u, v)){
                //....more computations needs to be done
              }
            }
          }
        }
    }
    return DataFrame::create();
  }
}