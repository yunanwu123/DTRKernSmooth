// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]
// [[Rcpp::plugins(cpp11)]]
#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <vector>
#include <iostream>
#include <Eigen/Dense>

using namespace Rcpp;
using namespace std;
using namespace arma;
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
  arma::vec A(inVec.data(), inVec.rows());
  return arma::stddev(A);
}

double eigenIQR(const Eigen::VectorXd & inVec){
  arma::vec A(inVec.data(), inVec.rows());
  arma::vec P = { 0.25, 0.75 };
  arma::vec Q = quantile(A, P);
  double res = Q(1)-Q(0);
  return res;
}

Eigen::MatrixXd eigenPercentile(const Eigen::MatrixXd & inMat,
                                const Eigen::VectorXd & perc){
  const int nData = inMat.rows();
  const int nvars = inMat.cols();
  const int nperc = perc.size();
  
  arma::mat A =  arma::mat(inMat.data(),nData,nvars);
  arma::vec perc_arma =  arma::vec(perc.data(), nperc);
  arma::mat Q_arma =  quantile(A,perc_arma,0);
  Eigen::MatrixXd Q =  Eigen::Map<Eigen::MatrixXd>(Q_arma.memptr(),nperc,nvars);
  
  return Q.transpose();
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
  } else{
    Eigen::VectorXd v5 = v*0.2;
    a = (v5.array()< -1.0).select(0.0,0.5+105.0/64.0*(v5.array()-5.0/3.0*v5.array().pow(3)+1.4*v5.array().pow(5)-3.0/7.0*v5.array().pow(7)));
    a = (v5.array()>1.0).select(1.0,a);
  }
  return a;
}
 
Eigen::VectorXd obq(Eigen::MatrixXd X, Eigen::VectorXd y1, Eigen::VectorXd y2,
                    Eigen::VectorXd prob, Eigen::VectorXd bet, const int out, 
                    const double phi, const double err_tol){
  int nData = X.rows();
  Eigen::VectorXd z = X*bet;
  double exp_ind = out == 1?-0.2:-1.0/9.0;
  
  double h = 0.9*pow(nData,exp_ind)*max(min(eigenSD(z),eigenIQR(z)/1.34),err_tol);
  Eigen::VectorXd g1 = kernel(z/h,out);
  Eigen::VectorXd g2 = kernel(z/h,out+1); 
  
  Eigen::VectorXd weight1 = y1.array()*g2.array() + y2.array(); 
  Eigen::VectorXd weight2 = y1.array()*g1.array()/h;
  Eigen::VectorXd weight3 = (2*prob.array()-1)*g2.array() + 1 - prob.array();
  Eigen::VectorXd weight4 = (2*prob.array()-1)*g1.array()/h;
  Eigen::VectorXd weight = (weight2.array()*weight3.array()-weight1.array()*weight4.array())/weight3.array().square();
  Eigen::VectorXd eta = bet + 0.5*X.transpose()*weight/nData/phi; 
  return eta;
}

bool FPhi(Eigen::MatrixXd X, Eigen::VectorXd y1, Eigen::VectorXd y2, 
          Eigen::VectorXd prob, Eigen::VectorXd bet0, Eigen::VectorXd bet1, 
          const char* kn, const double phi, const double err_tol){
  int nData = y1.size();
  Eigen::VectorXd z0 = X*bet0;
  Eigen::VectorXd z1 = X*bet1;
  double exp_ind;
  int out;
  
  if(strcmp(kn, "normal")){
    exp_ind = -0.2;
    out = 1;
  } else{
    exp_ind = -1.0/9.0;
    out = 3;
  }
  
  double h0 = 0.9*pow(nData,exp_ind)*max(min(eigenSD(z0),eigenIQR(z0)/1.34),err_tol);
  double h1 = 0.9*pow(nData,exp_ind)*max(min(eigenSD(z1),eigenIQR(z1)/1.34),err_tol);
  Eigen::VectorXd g0 = kernel(z0/h0,out+1);
  Eigen::VectorXd g1 = kernel(z1/h1,out+1); 
  
  Eigen::VectorXd value1 = (y1.array()*g0.array() + y2.array())/((2*prob.array()-1)*g0.array() + 1 - prob.array()); 
  Eigen::VectorXd value2 = (y1.array()*g1.array() + y2.array())/((2*prob.array()-1)*g1.array() + 1 - prob.array());
  
  double f0 = (value1.array()-value2.array()).mean();
  double f1 = phi*(bet0-bet1).squaredNorm();
  
  return f0>f1;
}


List LAMM(Eigen::MatrixXd X, Eigen::VectorXd y1, Eigen::VectorXd y2, 
          Eigen::VectorXd prob, Eigen::VectorXd bet, const char* kn, const double phi0,
          double phi, const double gamma, const double err_tol){
  phi = max(phi0, phi/gamma);
  int out = strcmp(kn, "normal")?1:3;
  Eigen::VectorXd bet0 = obq(X,y1,y2,prob,bet,out,phi,err_tol);
  while(FPhi(X,y1,y2,prob,bet0,bet,kn,phi,err_tol)){
    bet = bet0;
    phi = phi*gamma;
    bet0 = obq(X,y1,y2,prob,bet,out,phi,err_tol);
  }
  
  List output;
  output["bet"] = bet0;
  output["phi"] = phi;
  return output;
}


// [[Rcpp::export]]
double obj_value_C(Eigen::MatrixXd X, Eigen::VectorXd y, Eigen::VectorXd a,
                   const Eigen::VectorXd & m1, const Eigen::VectorXd & m0,
                   const Eigen::VectorXd & eta,const Eigen::VectorXd & prob) {
  Eigen::VectorXd g = Eigen::VectorXd::Zero(y.size());
  g = ((X*eta).array()>0).select(1,g); 
  Eigen::VectorXd weighted_y1 = y.array()*(2*a.array()-1)-(a-prob).array()*(m1+m0).array();
  Eigen::VectorXd weighted_y2 = y.array()*(1-a.array())+(a-prob).array()*m0.array();
  Eigen::VectorXd c =  weighted_y1.array()*g.array() + weighted_y2.array();
  Eigen::VectorXd weights = a.array()*prob.array()+(1-a.array())*(1-prob.array());
  double value = (c.array()/weights.array()).mean();
  return value;
}

// [[Rcpp::export]]
List Smooth_C(Eigen::MatrixXd X, Eigen::VectorXd y, Eigen::VectorXd a,
              Eigen::VectorXd initial, const Eigen::VectorXd & prob, const char* kn,
              const Eigen::VectorXd & m1, const Eigen::VectorXd & m0, double phi, 
              const double gamma, const double err_tol, const int iter_tol){
  
  Eigen::VectorXd weighted_y1 = y.array()*(2*a.array()-1)-(a-prob).array()*(m1+m0).array();
  Eigen::VectorXd weighted_y2 = y.array()*(1-a.array())+(a-prob).array()*m0.array();//(y.array()-y.mean())/weights.array();
  List L = LAMM(X,weighted_y1,weighted_y2,prob,initial,kn,phi,phi,gamma,err_tol);
  
  Eigen::VectorXd bet0 = L["bet"];
  double phi1 = L["phi"];
  int i = 1;
  while(((initial-bet0).squaredNorm()>err_tol)&(i<=iter_tol)){
    initial = bet0;
    L = LAMM(X,weighted_y1,weighted_y2,prob,initial,kn,phi,phi1,gamma,err_tol);
    bet0 = L["bet"];
    phi1 = L["phi"];
    i += 1;
  }
  bet0 = bet0/abs(bet0(1));
  double value = obj_value_C(X,y,a,m0,m1,bet0,prob);
  
  Eigen::VectorXd opt_trt = Eigen::VectorXd::Zero(y.size());
  opt_trt = ((X*bet0).array()>0).select(1,opt_trt);
  
  bool converge = (initial/abs(initial(1))-bet0).squaredNorm()<=err_tol; 
  
  List output;
  output["X"] = X;
  output["y"] = y;
  output["a"] = a;
  output["prob"] = prob;
  output["m0"] = m0;
  output["m1"] = m1;
  output["kernel"] = kn;
  output["beta_smooth"] = bet0;
  output["opt_treatment"] = opt_trt;
  output["value_smooth"] = value;
  output["converge"] = converge;
  output["iter_num"] = i; 
  return output;
}

// [[Rcpp::export]]
List Boots_C(Eigen::MatrixXd X, Eigen::VectorXd y, Eigen::VectorXd a,
             Eigen::VectorXd initial, const Eigen::VectorXd & prob, const char* kn,
             const Eigen::VectorXd & m1, const Eigen::VectorXd & m0,
             Eigen::MatrixXd weights, const double alpha, double phi,
             const double gamma, const double err_tol, const int iter_tol){
  List smooth_est = Smooth_C(X,y,a,initial,prob,kn,m1,m0,phi,gamma,err_tol,iter_tol);
  Eigen::VectorXd beta_est = smooth_est["beta_smooth"];
  double beta_value = smooth_est["value_smooth"];
  
  const int B = weights.cols();
  Eigen::VectorXd value_boots(B);
  Eigen::MatrixXd Beta_boots(B,X.cols());
  Eigen::VectorXd weights_col_i, weighted_y, weighted_m0, weighted_m1, beta_boots;
  List smooth_boots, output;
  for(int i=0; i<B; i++){
    weights_col_i = weights.col(i);
    weighted_y = y.array()*weights_col_i.array();
    weighted_m0 = m0.array()*weights_col_i.array();
    weighted_m1 = m1.array()*weights_col_i.array();
    smooth_boots = Smooth_C(X,weighted_y,a,initial,prob,kn,weighted_m1,weighted_m0,phi,gamma,err_tol,iter_tol);
    beta_boots = smooth_boots["beta_smooth"] ;
    Beta_boots.row(i) = beta_boots ;
    value_boots(i) =  obj_value_C(X,weighted_y,a,weighted_m1,weighted_m0,beta_est,prob);
  }
  Eigen::VectorXd perc(2);
  perc << alpha*0.5,1.0-alpha*0.5;
  Eigen::MatrixXd Beta_bound = eigenPercentile(Beta_boots,perc);
  Eigen::MatrixXd value_bound = eigenPercentile(value_boots,perc);
  Beta_bound << 2*beta_est.array()-Beta_bound.col(1).array(),2*beta_est.array()-Beta_bound.col(0).array();
  value_bound << 2*beta_value-value_bound(1), 2*beta_value-value_bound(0);
  
  output["alpha"] = alpha;
  output["B"] = B;
  output["smooth_est"] = smooth_est;
  output["Beta_CI"] = Beta_bound;
  output["value_CI"] = value_bound;
  return output;
}