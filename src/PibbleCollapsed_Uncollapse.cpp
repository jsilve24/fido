#include <fido.h>
#include <Rcpp/Benchmark/Timer.h>
#include <boost/random/mersenne_twister.hpp>

#ifdef FIDO_USE_PARALLEL
#include <omp.h>
#endif 

using namespace Rcpp;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::ArrayXXd;
using Eigen::Map;
using Eigen::Lower;

//Eta should be array with dim [D-1, N, iter]


//' Uncollapse output from optimPibbleCollapsed to full pibble Model
//' 
//' See details for model. Should likely be called following 
//' \code{\link{optimPibbleCollapsed}}. Notation: \code{N} is number of samples,
//' \code{D} is number of multinomial categories, \code{Q} is number
//' of covariates, \code{iter} is the number of samples of \code{eta} (e.g., 
//' the parameter \code{n_samples} in the function \code{optimPibbleCollapsed})
//' 
//' @param eta array of dimension (D-1) x N x iter (e.g., \code{Pars} output of 
//'   function optimPibbleCollapsed)
//' @param X matrix of covariates of dimension Q x N
//' @param Theta matrix of prior mean of dimension (D-1) x Q
//' @param Gamma covariance matrix of dimension Q x Q
//' @param Xi covariance matrix of dimension (D-1) x (D-1)
//' @param upsilon scalar (must be > D) degrees of freedom for InvWishart prior
//' @param ret_mean if true then uses posterior mean of Lambda and Sigma 
//'   corresponding to each sample of eta rather than sampling from 
//'   posterior of Lambda and Sigma (useful if Laplace approximation
//'   is not used (or fails) in optimPibbleCollapsed)
//' @param seed seed to use for random number generation 
//' @param ncores (default:-1) number of cores to use, if ncores==-1 then 
//' uses default from OpenMP typically to use all available cores. 
//'  
//' @details Notation: Let \eqn{Z_j}{Z_j} denote the J-th row of a matrix Z.
//' While the collapsed model is given by:
//'    \deqn{Y_j \sim \text{Multinomial}(\pi_j)}{Y_j \sim Multinomial(Pi_j)}
//'    \deqn{\pi_j = \Phi^{-1}(\eta_j)}{Pi_j = Phi^{-1}(Eta_j)}
//'    \deqn{\eta \sim T_{D-1, N}(\upsilon, \Theta X, K, A)}{Eta \sim T_{D-1, N}(upsilon, Theta*X, K, A)}
//' Where \eqn{A = I_N + X \Gamma X'}{A = I_N + X * Gamma * X'}, \eqn{K=\Xi}{K = Xi} is a (D-1)x(D-1) covariance 
//' matrix, \eqn{\Gamma}{Gamma} is a Q x Q covariance matrix, and \eqn{\Phi^{-1}}{Phi^(-1)} is ALRInv_D 
//' transform. 
//' 
//' The uncollapsed model (Full pibble model) is given by:
//'    \deqn{Y_j \sim \text{Multinomial}(\pi_j)}{Y_j \sim Multinomial(Pi_j)}
//'    \deqn{\pi_j = \Phi^{-1}(\eta_j)}{Pi_j = Phi^(-1)(Eta_j)}
//'    \deqn{\eta \sim MN_{D-1 \times N}(\Lambda X, \Sigma, I_N)}{Eta \sim MN_(D-1 x N)(Lambda*X, Sigma, I_N)}
//'    \deqn{\Lambda \sim MN_{D-1 x Q}(\Theta, \Sigma, \Gamma)}{Lambda \sim MN_(D-1 x Q)(Theta, Sigma, Gamma)}
//'    \deqn{\Sigma \sim InvWish(\upsilon, \Xi)}{Sigma \sim InvWish(upsilon, Xi)}
//' This function provides a means of sampling from the posterior distribution of 
//' \code{Lambda} and \code{Sigma} given posterior samples of \code{Eta} from 
//' the collapsed model. 
//' @return List with components 
//' 1. Lambda Array of dimension (D-1) x Q x iter (posterior samples)
//' 2. Sigma Array of dimension (D-1) x (D-1) x iter (posterior samples)
//' 3. The number of cores used
//' 4. Timer
//' @export
//' @md
//' @seealso \code{\link{optimPibbleCollapsed}}
//' @references JD Silverman K Roche, ZC Holmes, LA David, S Mukherjee. 
//'   Bayesian Multinomial Logistic Normal Models through Marginally Latent Matrix-T Processes. 
//'   2019, arXiv e-prints, arXiv:1903.11695
//' @examples
//' sim <- pibble_sim()
//' 
//' # Fit model for eta
//' fit <- optimPibbleCollapsed(sim$Y, sim$upsilon, sim$Theta%*%sim$X, sim$KInv, 
//'                              sim$AInv, random_pibble_init(sim$Y))  
//' 
//' # Finally obtain samples from Lambda and Sigma
//' fit2 <- uncollapsePibble(fit$Samples, sim$X, sim$Theta, 
//'                                    sim$Gamma, sim$Xi, sim$upsilon, 
//'                                    seed=2849)
// [[Rcpp::export]]
List uncollapsePibble(const Eigen::Map<Eigen::VectorXd> eta, // note this is essentially eta
                    const Eigen::Map<Eigen::MatrixXd> X, 
                    const Eigen::Map<Eigen::MatrixXd> Theta,
                    const Eigen::Map<Eigen::MatrixXd> Gamma, 
                    const Eigen::Map<Eigen::MatrixXd> Xi, 
                    const double upsilon, 
                    long seed, 
                    bool ret_mean = false, 
                    int ncores=-1){
  #ifdef FIDO_USE_PARALLEL
    Eigen::initParallel();
    if (ncores > 0) Eigen::setNbThreads(ncores);
    if (ncores > 0) {
      omp_set_num_threads(ncores);
    } else {
      omp_set_num_threads(omp_get_max_threads());
    }
  #endif 
  Timer timer;
  timer.step("Overall_start");
  List out(4);
  out.names() = CharacterVector::create("Lambda", "Sigma", "Timer", "NoCores");
  int Q = Gamma.rows();
  int D = Xi.rows()+1;
  int N = X.cols();
  int iter = eta.size()/(N*(D-1)); // assumes result is an integer !!!
  double upsilonN = upsilon + N;
  const MatrixXd GammaInv(Gamma.lu().inverse());
  const MatrixXd GammaInvN(GammaInv + X*X.transpose());
  const MatrixXd GammaN(GammaInvN.lu().inverse());
  const MatrixXd LGammaN(GammaN.llt().matrixL());
  //const Map<const MatrixXd> Eta(NULL);
  const MatrixXd ThetaGammaInvGammaN(Theta*GammaInv*GammaN);
  const MatrixXd XTGammaN(X.transpose()*GammaN);

  // Storage for output
  MatrixXd LambdaDraw0((D-1)*Q, iter);
  MatrixXd SigmaDraw0((D-1)*(D-1), iter);
  
  //iterate over all draws of eta - embarrassingly parallel with parallel rng
  #ifdef FIDO_USE_PARALLEL
    Eigen::setNbThreads(1);
    //Rcout << "thread: "<< omp_get_max_threads() << std::endl;
  #endif 
  #pragma omp parallel shared(D, N, Q, LambdaDraw0, SigmaDraw0)
  {
    boost::random::mt19937 rng(seed);
  #ifdef FIDO_USE_PARALLEL
    rng.discard(omp_get_thread_num()*iter);
  #endif 
  // storage for computation
  MatrixXd LambdaN(D-1, Q);
  MatrixXd XiN(D-1, D-1);
  MatrixXd LSigmaDraw(D-1, D-1);
  MatrixXd ELambda(D-1, Q);
  MatrixXd EEta(D-1, N);
  #pragma omp for 
  for (int i=0; i < iter; i++){
    //R_CheckUserInterrupt();
    const Map<const MatrixXd> Eta(&eta(i*N*(D-1)),D-1, N);
    LambdaN.noalias() = Eta*XTGammaN+ThetaGammaInvGammaN;
    ELambda = LambdaN-Theta;
    EEta.noalias() = Eta-LambdaN*X;
    XiN.noalias() = Xi+ EEta*EEta.transpose() + ELambda*GammaInv*ELambda.transpose();
    
    if (ret_mean){
      Map<VectorXd> LambdaNVec(LambdaN.data(), LambdaN.size());
      Map<VectorXd> XiNVec(XiN.data(), XiN.size());
      LambdaDraw0.col(i) = LambdaNVec;
      SigmaDraw0.col(i) = (upsilonN-D)*XiNVec; // mean of inverse wishart
    } else {
      // Draw Random Component
      rInvWishRevCholesky_thread_inplace(LSigmaDraw, upsilonN, XiN, rng);
      // Note: Below is valid even though LSigmaDraw is reverse cholesky factor
      Eigen::Ref<VectorXd> LambdaDraw_tmp = LambdaDraw0.col(i);
      Eigen::Map<MatrixXd> LambdaDraw(LambdaDraw_tmp.data(), D-1, Q);
      rMatNormalCholesky_thread_inplace(LambdaDraw, LambdaN, LSigmaDraw, 
                                        LGammaN.matrix(), rng);

      Eigen::Ref<VectorXd> SigmaDraw_tmp = SigmaDraw0.col(i);
      Eigen::Map<MatrixXd> SigmaDraw_tosquare(SigmaDraw_tmp.data(), D-1, D-1);
      SigmaDraw_tosquare.noalias() = LSigmaDraw*LSigmaDraw.transpose();
      
    }
  }
  }
  #ifdef FIDO_USE_PARALLEL
  if (ncores > 0){
    Eigen::setNbThreads(ncores);
  } else {
    Eigen::setNbThreads(omp_get_max_threads());  
  }
  #endif 
  int n_coresUsed = Eigen::nbThreads();

  IntegerVector dLambda = IntegerVector::create(D-1, Q, iter);
  IntegerVector dSigma = IntegerVector::create(D-1, D-1, iter);
  NumericVector nvLambda = wrap(LambdaDraw0);
  NumericVector nvSigma = wrap(SigmaDraw0);
  nvLambda.attr("dim") = dLambda;
  nvSigma.attr("dim") = dSigma;
  out[0] = nvLambda;
  out[1] = nvSigma;
  out[3] = n_coresUsed;
  timer.step("Overall_stop");
  NumericVector t(timer);
  out[2] = timer;
  return out;
}

// A few functions for testing MatDist Functions
// [[Rcpp::export]]
Eigen::MatrixXd rMatNormalCholesky_test(Eigen::MatrixXd M, 
                                        Eigen::MatrixXd LU, 
                                        Eigen::MatrixXd LV, 
                                        int discard){
  boost::random::mt19937 rng;
  rng.discard(discard);
  MatrixXd res(M.rows(), M.cols());
  rMatNormalCholesky_thread_inplace(res, M, LU, LV, rng);
  return res;
}

// [[Rcpp::export]]
Eigen::MatrixXd rInvWishRevCholesky_test(int v, Eigen::MatrixXd Psi){
  return rInvWishRevCholesky(v, Psi);
}

// [[Rcpp::export]]
Eigen::MatrixXd rInvWishRevCholesky_thread_test(int v, Eigen::MatrixXd Psi, 
                                                        int discard){
  boost::random::mt19937 rng;
  rng.discard(discard);
  return rInvWishRevCholesky_thread(v, Psi, rng);
}

// [[Rcpp::export]]
Eigen::MatrixXd rInvWishRevCholesky_thread_inplace_test(int v, Eigen::MatrixXd Psi, 
                                                        int discard){
  boost::random::mt19937 rng;
  rng.discard(discard);
  MatrixXd res(Psi.rows(), Psi.cols());
  rInvWishRevCholesky_thread_inplace(res, v, Psi, rng);
  return res;
}

// [[Rcpp::export]]
Eigen::MatrixXd rMatUnitNormal_test1(int n, int m){
  MatrixXd X(n,m);
  boost::random::mt19937 rng;
  fillUnitNormal_thread(X, rng);
  return X;
}

// [[Rcpp::export]]
Eigen::MatrixXd rMatUnitNormal_test2(int n){
  VectorXd X(n);
  fillUnitNormal(X);
  return X;
}
