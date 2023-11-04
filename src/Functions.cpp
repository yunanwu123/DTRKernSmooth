// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(cpp11)]]
#include <RcppEigen.h>
#include <vector>
#include <iostream>
#include <Eigen/Dense>

using namespace Rcpp;
using namespace std;
using namespace Eigen;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//


double eigenSD(const Eigen::VectorXd & inVec) {
  Eigen::VectorXd inVec_cent = inVec.array() - inVec.mean();
  int n = inVec.size();

  return sqrt(inVec_cent.array().square().sum()/(n-1));
}


Eigen::VectorXd eigenvecPercentiles(const Eigen::VectorXd& vec, const Eigen::VectorXd& percentiles) {
   int n = vec.size();
   int numPercentiles = percentiles.size();

   // Create a vector to store the percentile values
   Eigen::VectorXd results(numPercentiles);

   // Sort the vector
   Eigen::VectorXd sortedVec = vec;
   std::sort(sortedVec.data(), sortedVec.data() + n);

   // Calculate the requested percentiles
   for (int i = 0; i < numPercentiles; i++) {
     double id = (n - 1) * percentiles[i];
     int low = static_cast<int>(id);
     int high = low + 1;
     double weight = id - low;

     results[i] = (1.0 - weight) * sortedVec[low] + weight * sortedVec[high];
   }

   return results;
 }


double hbase(const Eigen::VectorXd & inVec, double err_tol){
  Eigen::VectorXd percentiles(2);
  percentiles << 0.25,0.75;
  Eigen::VectorXd quant = eigenvecPercentiles(inVec,percentiles);
  double iqr = quant[1] - quant[0];
  double stdA = eigenSD(inVec);
  double res = max(min(stdA,iqr/1.34),err_tol);
  return res;
}


Eigen::MatrixXd eigenPercentile(const Eigen::MatrixXd& inMat, const Eigen::VectorXd& perc) {
  int nvars = inMat.cols();

  Eigen::MatrixXd percentileMatrix(nvars,perc.size());

  for (int col = 0; col < nvars; col++) {
    Eigen::VectorXd columnVector = inMat.col(col);
    Eigen::VectorXd columnPercentiles = eigenvecPercentiles(columnVector, perc);
    percentileMatrix.row(col) = columnPercentiles;
  }

  return percentileMatrix;
}


Eigen::VectorXd kernel(const Eigen::VectorXd v, const int out) {
  Eigen::VectorXd a;
  if(out==1){
    a =  exp(-v.array().square()/2)/sqrt(2.0*M_PI);
  }  else if(out==2){
    a = 0.5 * v.unaryExpr([](double x) {return erfc(-x* M_SQRT1_2);});
  } else if(out==3){
    Eigen::VectorXd v5 = v*0.2;
    a = (v5.array().abs()<=1.0).select(105.0/64.0*(0.2-v5.array().square()+1.4*v5.array().pow(4)-0.6*v5.array().pow(6)),0.0);
  } else if(out==4){
    Eigen::VectorXd v5 = v*0.2;
    a = (v5.array()< -1.0).select(0.0,0.5+105.0/64.0*(v5.array()-5.0/3.0*v5.array().pow(3)+1.4*v5.array().pow(5)-3.0/7.0*v5.array().pow(7)));
    a = (v5.array()>1.0).select(1.0,a);
  }else if(out==5){
    Eigen::VectorXd v5 = v*0.2;
    a = (v5.array().abs()<=1.0).select(225.0/128.0*(0.2-14.0/15.0*v5.array().square()+21.0/25.0*v5.array().pow(4)),0.0);
  } else{
    Eigen::VectorXd v5 = v*0.2;
    a = (v5.array()< -1.0).select(0.0,0.5+225.0/128.0*(v5.array()-14.0/9.0*v5.array().pow(3)+21.0/25.0*v5.array().pow(5)));
    a = (v5.array()>1.0).select(1.0,a);
  }
  return a;
}

List obq(Eigen::MatrixXd X, Eigen::VectorXd y1, Eigen::VectorXd bet, const int out,
         const double phi, const double err_tol){
  int nData = X.rows();
  Eigen::VectorXd z = X*bet;
  double exp_ind;
  if(out==1){
    exp_ind = -0.2;
  }else if(out==3){
    exp_ind = -1.0/9.0;
  }else{
    exp_ind = -1.0/13.0;
  } 

  double h = 0.9*pow(nData,exp_ind)*hbase(z,err_tol);
  Eigen::VectorXd g1 = kernel(z/h,out);
  Eigen::VectorXd g2 = kernel(z/h,out+1);

  Eigen::VectorXd weight1 = y1.array()*g2.array() ;
  Eigen::VectorXd weight2 = y1.array()*g1.array();
  Eigen::VectorXd eta = bet + 0.5*X.transpose()*weight2/nData/phi/h;

  List output;
  output["eta"] = eta;
  output["value"] = weight1.mean();
  return output;
}

List LAMM(Eigen::MatrixXd X, Eigen::VectorXd y1, Eigen::VectorXd bet,
          const char* kn, const double phi0, double phi, const double gamma,
          const double err_tol){
  phi = max(phi0, phi/gamma);
  int out; 
  if(strcmp(kn, "normal")){
    out = 1;
  }else if(strcmp(kn, "poly1")){
    out = 3;
  }else{
    out = 5;
  } 
   
  List output0 = obq(X,y1, bet,out,phi,err_tol);
  Eigen::VectorXd bet0 = output0["eta"];
  double value0 = output0["value"];
  double value = value0 - 1;
  double delta = 0;
  while(value0 - value > phi*delta){
    delta = (bet0-bet).squaredNorm();
    value = value0;
    bet = bet0;
    phi = phi*gamma;
    output0 = obq(X,y1,bet,out,phi,err_tol);
    bet0 = output0["eta"];
    value0 = output0["value"];
  }

  List output;
  output["bet"] = bet;
  output["phi"] = phi;
  return output;
}


// [[Rcpp::export]]
double obj_value_C(Eigen::MatrixXd X, Eigen::VectorXd y, Eigen::VectorXd a,
                   const Eigen::VectorXd & m1, const Eigen::VectorXd & m0,
                   const Eigen::VectorXd & eta,const Eigen::VectorXd & prob) {
  Eigen::VectorXd g = Eigen::VectorXd::Zero(y.size());
  g = ((X*eta).array()>0).select(1,g);

  Eigen::VectorXd weighted_y1 = (2*a.array()-1)*y.array() - m1.array()*(a.array()-prob.array())+m0.array()*(1-a.array()-prob.array());
  Eigen::VectorXd weighted_y2 = y.array()*(1-a.array())-(1-a.array()-prob.array())*m0.array();
  Eigen::VectorXd c =  weighted_y1.array()*g.array() + weighted_y2.array();
  double value = (c.array()/prob.array()).mean();
  return value;
}

// [[Rcpp::export]]
List Smooth_C(Eigen::MatrixXd X, Eigen::VectorXd y, Eigen::VectorXd a, bool intercept,
              Eigen::VectorXd initial, const Eigen::VectorXd & prob, const char* kn,
              const Eigen::VectorXd & m1, const Eigen::VectorXd & m0, double phi,
              const double gamma, const double err_tol, const int iter_tol){

  Eigen::VectorXd weighted_y = (2*a.array()-1)*y.array() - m1.array()*(a.array()-prob.array())+m0.array()*(1-a.array()-prob.array());
  weighted_y = weighted_y.cwiseQuotient(prob);
  int p = X.cols();
  int n = X.rows();
  Eigen::MatrixXd newX(n,initial.size());
  if(intercept){
    newX << Eigen::MatrixXd::Ones(n, 1), X;
  }else{
    newX = X;
  }

  List L = LAMM(newX,weighted_y,initial,kn,phi,phi,gamma,err_tol);

  Eigen::VectorXd bet0 = L["bet"];
  double phi1 = L["phi"];
  int i = 1;
  while(((initial-bet0).norm()>err_tol)&(i<=iter_tol)){
    initial = bet0;
    L = LAMM(newX,weighted_y,initial,kn,phi,phi1,gamma,err_tol);
    bet0 = L["bet"];
    phi1 = L["phi"];
    i += 1;
  }

  Eigen::VectorXd beta(p+1);
  if(intercept){
    beta << bet0;
  }else{
    beta << 0, bet0;
  }
  beta = beta/abs(beta(1));
  double value = obj_value_C(newX,y,a,m0,m1,bet0,prob);

  Eigen::VectorXd opt_trt = Eigen::VectorXd::Zero(y.size());
  opt_trt = ((newX*bet0).array()>0).select(1,opt_trt);

  bool converge = (initial-bet0).norm()<=err_tol;  

  List output;
  output["X"] = X;
  output["y"] = y;
  output["a"] = a;
  output["intercept"] = intercept;
  output["prob"] = prob;
  output["m0"] = m0;
  output["m1"] = m1;
  output["kernel"] = kn;
  output["beta_smooth"] = beta;
  output["opt_treatment"] = opt_trt;
  output["value_smooth"] = value;
  output["converge"] = converge;
  output["iter_num"] = i;
  return output;
}

// [[Rcpp::export]]
List Boots_C(Eigen::MatrixXd X, Eigen::VectorXd y, Eigen::VectorXd a, bool intercept,
             Eigen::VectorXd initial, const Eigen::VectorXd & prob, const char* kn,
             const Eigen::VectorXd & m1, const Eigen::VectorXd & m0,
             Eigen::MatrixXd weights, const double alpha, double phi,
             const double gamma, const double err_tol, const int iter_tol){
  List smooth_est = Smooth_C(X,y,a,intercept,initial,prob,kn,m1,m0,phi,gamma,err_tol,iter_tol);
  Eigen::VectorXd beta_est = smooth_est["beta_smooth"];
  double beta_value = smooth_est["value_smooth"];

  const int B = weights.cols();
  Eigen::VectorXd value_boots(B);
  Eigen::MatrixXd Beta_boots(B,X.cols()+1);
  Eigen::VectorXd weights_col_i, weighted_y, weighted_m0, weighted_m1, beta_boots;
  List smooth_boots, output;
  for(int i=0; i<B; i++){
    weights_col_i = weights.col(i);
    weighted_y = y.array()*weights_col_i.array();
    weighted_m0 = m0.array()*weights_col_i.array();
    weighted_m1 = m1.array()*weights_col_i.array();
    smooth_boots = Smooth_C(X,weighted_y,a,intercept,initial,prob,kn,weighted_m1,weighted_m0,phi,gamma,err_tol,iter_tol);
    beta_boots = smooth_boots["beta_smooth"] ;
    Beta_boots.row(i) = beta_boots ;
    value_boots(i) =  obj_value_C(X,weighted_y,a,weighted_m1,weighted_m0,beta_est,prob);
  }
  Eigen::VectorXd perc(2);
  perc << alpha*0.5,1.0-alpha*0.5;
  Eigen::MatrixXd Beta_bound = eigenPercentile(Beta_boots,perc);
  Eigen::MatrixXd value_bound = eigenPercentile(value_boots,perc);
  Eigen::MatrixXd Beta_CI(Beta_bound.rows(),2);
  Eigen::VectorXd value_CI(2);

  Beta_CI << 2*beta_est.array()-Beta_bound.col(1).array(),2*beta_est.array()-Beta_bound.col(0).array();
  value_CI << 2*beta_value-value_bound(1), 2*beta_value-value_bound(0);

  output["alpha"] = alpha;
  output["B"] = B;
  output["smooth_est"] = smooth_est;
  output["Beta_CI"] = Beta_CI;
  output["value_CI"] = value_CI;
  return output;
}
