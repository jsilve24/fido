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


//' Uncollapse output from optimPibbleCollapsed to full pibble Model when Sigma is known
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
//' @param GammaComb summed covariance matrix across additive components of dimension Q x Q.
//' @param Xi covariance matrix of dimension (D-1) x (D-1)
//' @param sigma known covariance matrix of dimension (D-1) x (D-1) x iter
//' @param upsilon scalar (must be > D) degrees of freedom for InvWishart prior
//' @param ret_mean if true then uses posterior mean of Lambda and Sigma 
//'   corresponding to each sample of eta rather than sampling from 
//'   posterior of Lambda and Sigma (useful if Laplace approximation
//'   is not used (or fails) in optimPibbleCollapsed)
//' @param seed seed to use for random number generation
//' @param linear Boolean. Is this for a linear parameter? 
//' @param ncores (default:-1) number of cores to use, if ncores==-1 then 
//' uses default from OpenMP typically to use all available cores. 
//'  
//' @details Notation: Let \eqn{Z_j} denote the J-th row of a matrix Z.
//' While the collapsed model is given by:
//'    \deqn{Y_j \sim Multinomial(\pi_j)}{Y_j \sim Multinomial(Pi_j)}
//'    \deqn{\pi_j = \Phi^{-1}(\eta_j)}{Pi_j = Phi^(-1)(Eta_j)}
//'    \deqn{\eta \sim T_{D-1, N}(\upsilon, \Theta X, K, A)}{Eta \sim T_(D-1, N)(upsilon, Theta*X, K, A)}
//' Where \eqn{A = I_N + X \Gamma X'}{A = I_N + X * Gamma * X'}, \eqn{K = \Xi}{K = Xi} is a (D-1)x(D-1) covariance 
//' matrix, \eqn{\Gamma}{Gamma} is a Q x Q covariance matrix, and \eqn{\Phi^{-1}}{Phi^(-1)} is ALRInv_D 
//' transform. 
//' 
//' The uncollapsed model (Full pibble model) is given by:
//'    \deqn{Y_j \sim Multinomial(\pi_j)}{Y_j \sim Multinomial(Pi_j)}
//'    \deqn{\pi_j = \Phi^{-1}(\eta_j)}{Pi_j = Phi^(-1)(Eta_j)}
//'    \deqn{\eta \sim MN_{D-1 \times N}(\Lambda X, \Sigma, I_N)}{Eta \sim MN_{D-1 \times N}(Lambda*X, Sigma, I_N)}
//'    \deqn{\Lambda \sim MN_{D-1 \times Q}(\Theta, \Sigma, \Gamma)}{Lambda \sim MN_{D-1 x Q}(Theta, Sigma, Gamma)}
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
// [[Rcpp::export]]
List uncollapsePibble_sigmaKnown(const Eigen::Map<Eigen::VectorXd> eta, // note this is essentially eta
                    const Eigen::Map<Eigen::MatrixXd> X, 
                    const Eigen::Map<Eigen::MatrixXd> Theta,
                    const Eigen::Map<Eigen::MatrixXd> Gamma, 
                    const Eigen::Map<Eigen::MatrixXd> GammaComb,
                    const Eigen::Map<Eigen::MatrixXd> Xi, 
                    const Eigen::Map<Eigen::VectorXd> sigma,
                    const double upsilon, 
                    long seed, 
                    bool ret_mean = false, 
                    bool linear = false,
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
  List out(3);
  out.names() = CharacterVector::create("Lambda", "Timer", "NoCores");
  int Q = Gamma.rows();
  int D = Xi.rows()+1;
  int N = X.cols();
  int iter = eta.size()/(N*(D-1)); // assumes result is an integer !!!
  const MatrixXd Xt(X.transpose());
  const MatrixXd GammaInv(Gamma.lu().inverse());
  const MatrixXd GammaCombInv(GammaComb.lu().inverse());
  const MatrixXd XGammaCombInvXt(X*GammaCombInv*Xt);
  const MatrixXd GammaInvSum(GammaInv+XGammaCombInvXt);
  const MatrixXd GammaInvSumInv(GammaInvSum.lu().inverse());
  const MatrixXd LGammaInv(GammaInvSumInv.llt().matrixL());
  

  // Storage for output
  MatrixXd LambdaDraw0((D-1)*Q, iter);

  //iterate over all draws of eta - embarrassingly parallel with parallel rng
  #ifdef FIDO_USE_PARALLEL
    Eigen::setNbThreads(1);
    //Rcout << "thread: "<< omp_get_max_threads() << std::endl;
  #endif 
  #pragma omp parallel shared(D, N, Q, LambdaDraw0)
  {
    boost::random::mt19937 rng(seed);
  #ifdef FIDO_USE_PARALLEL
    rng.discard(omp_get_thread_num()*iter);
  #endif 
  // storage for computation
  MatrixXd LambdaN(D-1, Q);
  #pragma omp for 
  for (int i=0; i < iter; i++){
    //R_CheckUserInterrupt();
    const Map<const MatrixXd> Eta(&eta(i*N*(D-1)),D-1, N);
    const Map<const MatrixXd> Sigma(&sigma(i*(D-1)*(D-1)), D-1, D-1);
    const MatrixXd LSigma(Sigma.llt().matrixL());
    
    
    LambdaN.noalias() = (Eta*GammaCombInv*Xt*GammaInvSumInv + Theta*GammaInv*GammaInvSumInv);
    
    if (ret_mean){
      Map<VectorXd> LambdaNVec(LambdaN.data(), LambdaN.size());
      LambdaDraw0.col(i) = LambdaNVec;
    } else {
      // Note: Below is valid even though LSigmaDraw is reverse cholesky factor
      Eigen::Ref<VectorXd> LambdaDraw_tmp = LambdaDraw0.col(i);
      Eigen::Map<MatrixXd> LambdaDraw(LambdaDraw_tmp.data(), D-1, Q);
      rMatNormalCholesky_thread_inplace(LambdaDraw, LambdaN, LSigma.matrix(), 
                                          LGammaInv.matrix(), rng);
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
  NumericVector nvLambda = wrap(LambdaDraw0);
  nvLambda.attr("dim") = dLambda;
  out[0] = nvLambda;
  out[2] = n_coresUsed;
  timer.step("Overall_stop");
  NumericVector t(timer);
  out[1] = timer;
  return out;
}
