# EA_manuscript_simulations.R -- standalone reproduction of the manuscript
# Model: X ~ N_d(mu,Sigma), A ~ W_d(n,Sigma), independently; H0: mu=0.
# n is the Wishart degrees of freedom. All signal vectors are means of X.
# QCQ/ZWS additionally use the parameter-free Gaussian residual orientation.
# The n+1 Gaussian rows below construct that joint experiment; there is no
# second study-size parameter. Scores retain every validated normalization.
#
# Source this file to define functions without starting any simulation.
# Examples (full budgets are expensive):
#   source("EA_manuscript_simulations.R")
#   manuscript_settings("original")
#   run_manuscript("original", budget="quick", group_ids="cor_0.1")
#   run_manuscript("confirmation", budget="full", workers=2)
#   source_evaluation("outputs", "manuscript_results/source")
# New input values: open EA_explore_settings.R and edit the six input vectors.
# Individual existing panels: EA_individual_figures.R (12 runnable blocks).
#   plot_manuscript_figure(1, "EA_manuscript_results")  # RStudio Plots pane
#   plot_manuscript_figure(1, "EA_manuscript_results", "main_known_rho.png")
# No package installation. Optional Rcpp/C++ accelerates the SAME formulas;
# base R is the portable fallback. Source-evaluation reuses validated results.
# The estimated panel uses oracle null banks, not general feasible calibration.

manuscript_methods <- c("EA","Z-Max","Z-TopK","IHC","Feng-COM",
                         "QCQ","ZWS-{2,6}","Min-p","Quadratic")
.ea_runtime <- new.env(parent=baseenv())
.ea_runtime$raw_cpp <- .ea_runtime$precision_cpp <- FALSE
.ea_runtime$backend <- "base"
.ea_runtime$source_file <- local({
  paths <- unlist(lapply(sys.frames(),function(f)f$ofile),use.names=FALSE)
  if(!length(paths))paths <- sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE))
  if(length(paths))normalizePath(tail(paths,1L),mustWork=FALSE) else NULL
})

# 1. Evidence scores. Internal transforms precede complete-score calibration.
log_two_tail <- function(x,df=Inf) if(is.infinite(df))
  log(2)+pnorm(-abs(x),log.p=TRUE) else log(2)+pt(-abs(x),df,log.p=TRUE)
hc_score <- function(p,lo,hi) {
  p<-sort(p);j<-seq_along(p);eligible<-p>=lo & p<=hi
  if(!any(eligible))return(-Inf)
  max(sqrt(length(p))*(j[eligible]/length(p)-p[eligible])/
        sqrt(p[eligible]*(1-p[eligible])))
}
ea_logscore <- function(lp,lq)
  -log(length(lp))-min(pmin(log(2)+pmin(lp,lq),pmax(lp,lq)))
make_context <- function(Sigma,n) {
  stopifnot(is.matrix(Sigma),n>=nrow(Sigma),n==as.integer(n))
  if(any(!is.finite(Sigma))||max(abs(Sigma-t(Sigma)))>1e-12)stop("Sigma must be a finite symmetric covariance matrix")
  d<-nrow(Sigma);L<-t(chol(Sigma));Omega<-chol2inv(chol(Sigma))
  sizes<-attr(Sigma,"block_sizes");rhos<-attr(Sigma,"block_rhos")
  if(!is.null(sizes)) {
    reference<-diag(d);start<-1L
    if(sum(sizes)!=d||length(rhos)!=length(sizes))sizes<-NULL else {
      for(j in seq_along(sizes)){ii<-start:(start+sizes[j]-1L)
        reference[ii,ii]<-rhos[j];diag(reference)[ii]<-1;start<-start+sizes[j]}
      if(max(abs(Sigma-reference))>1e-14)sizes<-NULL
    }
  }
  U0<-solve(L);U<-U0;b<-max(1,floor(log(d)))
  U[row(U)-col(U)>=b]<-0;U<-sweep(U,2,sqrt(colSums(U^2)),"/")
  list(Sigma=Sigma,n=n,d=d,L=L,Omega=Omega,IHC=crossprod(U,U0),b=b,
       block_sizes=sizes,block_rhos=if(is.null(sizes))NULL else rhos)
}
known_scores <- function(X,ctx) {
  d<-ctx$d;n<-ctx$n;m<-X/sqrt(diag(ctx$Sigma))
  z<-drop(ctx$Omega%*%X)/sqrt(diag(ctx$Omega));w<-sort(z*z,decreasing=TRUE)
  k<-seq_len(min(40,d,n));lp<-log_two_tail(m);lq<-log_two_tail(z)
  c(EA=ea_logscore(lp,lq),`Z-Max`=max(w),
    `Z-TopK`=max(sqrt((1-k/(n+1))/(2*k))*(cumsum(w)[k]-k)),
    IHC=hc_score(2*pnorm(-abs(ctx$IHC%*%X)),1/d,.5)/sqrt(2*ctx$b-1),
    `Min-p`=-min(lp),Quadratic=sum(X*(ctx$Omega%*%X)))
}
unknown_prepare <- function(A,n) {
  list(A=A,Ai=chol2inv(chol(A)),Gamma=nodewise_lasso(A,n)$Gamma,n=n,d=nrow(A))
}
unknown_scores <- function(X,prep) {
  n<-prep$n;d<-prep$d;u<-drop(prep$Ai%*%X);q<-sum(X*u)
  den<-diag(prep$Ai)*(1+q)-u*u
  if(any(den<=0))stop("Undefined EA adjusted Student denominator")
  lp<-log_two_tail(sqrt(n)*X/sqrt(diag(prep$A)),n)
  lq<-log_two_tail(sqrt(n-d+1)*u/sqrt(den),n-d+1)
  z<-drop(prep$Gamma%*%X)/sqrt(diag(prep$Gamma));w<-sort(z*z,decreasing=TRUE)
  k<-seq_len(min(40,d,n))
  c(EA=ea_logscore(lp,lq),`Z-Max`=max(w),
    `Z-TopK`=max(sqrt((1-k/(n+1))/(2*k))*(cumsum(w)[k]-k)),
    `Min-p`=-min(lp),Quadratic=n*q)
}

# 2. Raw moment preparation: same centered ancillary rows for every method.


.raw_e2 <- function(p1, p2) (p1 * p1 - p2) / 2

# Degree-six Newton identity. Inputs can be coordinate vectors or d-by-d
# matrices of power sums. This exactly replaces enumeration of choose(n+1,6)
# subsets; its six matrix power sums are crossprod(C^r), r=1,...,6.
.raw_e6 <- function(p1, p2, p3, p4, p5, p6) {
  p1sq <- p1 * p1
  p2sq <- p2 * p2
  (p1sq^3 - 15 * p1sq^2 * p2 + 45 * p1sq * p2sq - 15 * p2sq * p2 +
     40 * p1sq * p1 * p3 - 120 * p1 * p2 * p3 + 40 * p3 * p3 -
     90 * p1sq * p4 + 90 * p2 * p4 + 144 * p1 * p5 - 120 * p6) / 720
}

.raw_assert_positive <- function(x, label) {
  if (length(x) != 1L || !is.finite(x) || x <= 0) {
    stop(sprintf("Undefined %s (%s); do not silently truncate the variance.",
                 label, format(x)), call. = FALSE)
  }
  x
}

raw_prepare <- function(Y0) {
  Y0<-as.matrix(Y0);storage.mode(Y0)<-"double"
  rows<-nrow(Y0);n<-rows-1L;d<-ncol(Y0)
  if(rows<6L||d<3L||any(!is.finite(Y0)))stop("Need finite rows, n>=5 and d>=3")
  mean0<-colMeans(Y0);C<-sweep(Y0,2,mean0,"-")
  scl<-sqrt(max(colSums(C*C)/n));.raw_assert_positive(scl,"common scale")
  C<-C/scl;mean0<-mean0/scl
  if(.ea_runtime$raw_cpp) {
    fused<-.ea_runtime$raw_pair_moments_cpp(C)
    A<-fused$A;moment<-fused$moment;V2sum<-fused$V2sum;V6sum<-fused$V6sum
  } else {
    moment<-matrix(NA_real_,7L,d);moment[1L,]<-rows
    power<-C;pairs<-vector("list",6L)
    for(r in 1:6){moment[r+1L,]<-colSums(power);pairs[[r]]<-crossprod(power)
      if(r<6L)power<-power*C}
    A<-pairs[[1L]];V2sum<-sum(.raw_e2(pairs[[1L]],pairs[[2L]]))
    V6sum<-sum(.raw_e6(pairs[[1L]],pairs[[2L]],pairs[[3L]],pairs[[4L]],pairs[[5L]],pairs[[6L]]))
  }
  diagA<-diag(A);if(any(!is.finite(diagA))||any(diagA<=0))stop("Undefined scatter diagonal")
  .raw_assert_positive(V2sum,"ZWS-2 centered variance")
  .raw_assert_positive(V6sum,"ZWS-6 centered variance")
  trR2<-sum((A/sqrt(outer(diagA,diagA)))^2)
  FengVar<-2*(trR2-d*d/n);.raw_assert_positive(FengVar,"Feng-SUM variance")
  ai<-rowSums(C*C);qcq_v<-colSums(C*ai)
  J0<-(n*n*sum(A*A)-rows*rows*sum(ai*ai)+sum(ai)^2)/(n-1)^2
  list(n=n,rows=rows,d=d,scale=scl,mean0=mean0,moment=moment,A=A,
       diagA=diagA,zws_norm2=rows/sqrt(2*choose(rows,2)*V2sum),
       zws_norm6=rows^3/sqrt(720*choose(rows,6)*V6sum),
       FengSD=sqrt(FengVar),qcq_v=qcq_v,qcq_J0=J0)
}


.raw_shifted_moments <- function(moment, m) {
  out <- matrix(0, 6L, length(m))
  for (r in 1:6) {
    for (k in 0:r) {
      out[r, ] <- out[r, ] + choose(r, k) * m^(r - k) * moment[k + 1L, ]
    }
  }
  out
}

raw_from_prepared <- function(prep, mu = NULL) {
  rows <- prep$rows
  d <- prep$d
  if (is.null(mu)) mu <- numeric(d)
  if (length(mu) != d || any(!is.finite(mu))) {
    stop("mu must be a finite vector of length ncol(Y0).", call. = FALSE)
  }
  m <- prep$mean0 + mu / (sqrt(rows) * prep$scale)
  moments <- .raw_shifted_moments(prep$moment, m)
  U2sum <- sum(.raw_e2(moments[1L, ], moments[2L, ]))
  U6sum <- sum(.raw_e6(moments[1L, ], moments[2L, ], moments[3L, ],
                      moments[4L, ], moments[5L, ], moments[6L, ]))
  ZWS2 <- prep$zws_norm2 * U2sum
  ZWS6 <- prep$zws_norm6 * U6sum
  pz <- min(pnorm(ZWS2, lower.tail = FALSE),
            pnorm(ZWS6, lower.tail = FALSE))
  ZWScomb <- -pz * (2 - pz)  # exactly -[1-(1-min(p2,p6))^2]

  # S=A/rows, so rows*m_j^2/S_jj = rows^2*m_j^2/A_jj.
  standardized2 <- rows^2 * m * m / prep$diagA
  FengSum <- (sum(standardized2) - (rows - 1) * d / (rows - 3)) / prep$FengSD
  FengMax <- max(standardized2)
  gumbel_argument <- FengMax - 2 * log(d) + log(log(d))
  pFengMax <- -expm1(-exp(-gumbel_argument / 2) / sqrt(pi))
  pFengSum <- pnorm(FengSum, lower.tail = FALSE)
  pf <- min(pFengSum, pFengMax)
  FengComb <- -pf * (2 - pf)

  # J(m)=J(0)-2N*m'v/(rows-2)^2+rows*m'A*m/(rows-2), derived by
  # summing the leave-two-out expression over ordered i != j.
  mAm <- sum(m * as.vector(prep$A %*% m))
  J <- prep$qcq_J0 - 2 * rows * sum(m * prep$qcq_v) / (rows - 2)^2 +
    rows * mAm / (rows - 2)
  .raw_assert_positive(J, "QCQ leave-two-out variance")
  qcq_numerator <- rows * (rows - 1) * sum(m * m) - sum(prep$diagA)
  T2 <- qcq_numerator / sqrt(2 * J)
  t2 <- standardized2 * (rows - 1) / rows
  p <- sort(2 * pnorm(-sqrt(t2)))
  eligible <- p >= log(d) / d & p <= d^(-1/4)
  H <- if (any(eligible)) {
    jj <- seq_len(d)[eligible]
    pp <- p[eligible]
    max(sqrt(d) * (jj / d - pp) / sqrt(pp * (1 - pp)))
  } else -Inf
  Einf <- as.numeric(max(t2) > 2 * log(d) + 3 * log(log(d)))
  threshold_sq <- 6 * log(log(d)) + log(log(log(d)))
  if(!is.finite(threshold_sq)||threshold_sq<0)stop('QCQ enhancement threshold is undefined for this dimension')
  EH <- as.numeric(H > sqrt(threshold_sq))
  QCQ <- T2 + sqrt(d) / 2 * (Einf + EH)
  c("Feng-SUM" = FengSum, "Feng-MAX" = FengMax,
    "Feng-COM" = FengComb, "QCQ" = QCQ,
    "ZWS-2" = ZWS2, "ZWS-6" = ZWS6, "ZWS-{2,6}" = ZWScomb)
}
raw_scores <- function(Y) raw_from_prepared(raw_prepare(Y))



# 3. Optional fused moment kernel (portable scalar path; ARM SIMD when available).


compile_raw_fast <- function(verbose = FALSE) {
  if (!requireNamespace("Rcpp", quietly = TRUE)) return(FALSE)
  if (exists("raw_pair_moments_cpp", envir=.ea_runtime, inherits=FALSE)) return(TRUE)
  Rcpp::cppFunction(code = '
Rcpp::List raw_pair_moments_cpp(Rcpp::NumericMatrix C) {
  const int rows=C.nrow(), d=C.ncol();
  Rcpp::NumericMatrix A(d,d), moments(7,d);
  const double * cp=C.begin();
  for(int j=0; j<d; ++j) {
    const double *cj=cp+j*rows;
    double s1=0,s2=0,s3=0,s4=0,s5=0,s6=0;
    for(int i=0; i<rows; ++i) {
      const double x=cj[i], x2=x*x, x3=x2*x;
      s1+=x; s2+=x2; s3+=x3; s4+=x2*x2;
      s5+=x3*x2; s6+=x3*x3;
    }
    moments(0,j)=rows; moments(1,j)=s1; moments(2,j)=s2;
    moments(3,j)=s3; moments(4,j)=s4; moments(5,j)=s5;
    moments(6,j)=s6;
  }
  double v2=0,v6=0;
  for(int j=0; j<d; ++j) {
    const double *cj=cp+j*rows;
    for(int k=0; k<=j; ++k) {
      const double *ck=cp+k*rows;
      double s1=0,s2=0,s3=0,s4=0,s5=0,s6=0;
      int i=0;
#ifdef __aarch64__
      float64x2_t t1=vdupq_n_f64(0),t2=vdupq_n_f64(0),t3=vdupq_n_f64(0),
        t4=vdupq_n_f64(0),t5=vdupq_n_f64(0),t6=vdupq_n_f64(0);
      for(; i+1<rows; i+=2) {
        const float64x2_t x=vmulq_f64(vld1q_f64(cj+i),vld1q_f64(ck+i));
        const float64x2_t x2=vmulq_f64(x,x), x3=vmulq_f64(x2,x);
        t1=vaddq_f64(t1,x); t2=vaddq_f64(t2,x2); t3=vaddq_f64(t3,x3);
        t4=vaddq_f64(t4,vmulq_f64(x2,x2));
        t5=vaddq_f64(t5,vmulq_f64(x3,x2));
        t6=vaddq_f64(t6,vmulq_f64(x3,x3));
      }
      s1=vaddvq_f64(t1);s2=vaddvq_f64(t2);s3=vaddvq_f64(t3);
      s4=vaddvq_f64(t4);s5=vaddvq_f64(t5);s6=vaddvq_f64(t6);
#endif
      for(; i<rows; ++i) {
        const double x=cj[i]*ck[i], x2=x*x, x3=x2*x;
        s1+=x; s2+=x2; s3+=x3; s4+=x2*x2;
        s5+=x3*x2; s6+=x3*x3;
      }
      A(j,k)=A(k,j)=s1;
      const double p1sq=s1*s1, p2sq=s2*s2;
      const double e2=(p1sq-s2)/2;
      const double e6=(p1sq*p1sq*p1sq-15*p1sq*p1sq*s2+
        45*p1sq*p2sq-15*p2sq*s2+40*p1sq*s1*s3-
        120*s1*s2*s3+40*s3*s3-90*p1sq*s4+
        90*s2*s4+144*s1*s5-120*s6)/720;
      const double mult=(j==k)?1:2;
      v2+=mult*e2; v6+=mult*e6;
    }
  }
  return Rcpp::List::create(Rcpp::Named("A")=A,
    Rcpp::Named("moment")=moments, Rcpp::Named("V2sum")=v2,
    Rcpp::Named("V6sum")=v6);
}', includes = '#ifdef __aarch64__\n#include <arm_neon.h>\n#endif',
    env = .ea_runtime, verbose = verbose)
  TRUE
}


# 4. Ordinary nodewise Lasso; full KKT checks on every fitted row.


precision_engine_compile <- function() {
  if (exists("nodewise_lasso_cpp", .ea_runtime, inherits = FALSE)) return(invisible(TRUE))
  if (!requireNamespace("Rcpp", quietly = TRUE)) return(invisible(FALSE))
  Rcpp::cppFunction(env = .ea_runtime, plugins = "cpp11", code = '
  Rcpp::List nodewise_lasso_cpp(Rcpp::NumericMatrix S, double lambda,
       double tol = 1e-10, int maxit = 100000,
       Rcpp::Nullable<Rcpp::NumericMatrix> warm = R_NilValue) {
    int d = S.nrow();
    if (S.ncol() != d) Rcpp::stop("S must be square");
    Rcpp::NumericMatrix Beta(d,d), Gamma(d,d);
    Rcpp::NumericVector tau2(d), kkt(d);
    Rcpp::IntegerVector iterations(d);
    if (warm.isNotNull()) {
      Rcpp::NumericMatrix W(warm);
      if(W.nrow()!=d || W.ncol()!=d) Rcpp::stop("warm must be d by d");
      std::copy(W.begin(), W.end(), Beta.begin());
    }
    std::vector<double> b(d), r(d), diag(d);
    for(int k=0;k<d;++k) {
      diag[k]=S(k,k);
      if(!R_finite(diag[k]) || diag[k]<=0) Rcpp::stop("S has nonpositive diagonal");
    }
    for(int j=0;j<d;++j) {
      std::vector<int> active;
      std::vector<unsigned char> isactive(d,0);
      for(int k=0;k<d;++k) b[k]=(k==j?0.0:Beta(j,k));
      for(int k=0;k<d;++k) r[k]=S(k,j);
      for(int l=0;l<d;++l) if(b[l]!=0)
        for(int k=0;k<d;++k) r[k]-=S(k,l)*b[l];
      for(int k=0;k<d;++k) {
        if(k!=j && (b[k]!=0 || std::abs(r[k])>lambda+tol)) {
          active.push_back(k); isactive[k]=1;
        }
      }
      int cycles=0; double violation=R_PosInf;
      while(cycles<maxit) {
        ++cycles; double delta=0;
        for(size_t ai=0;ai<active.size();++ai) {
          int k=active[ai];
          double v=r[k]+diag[k]*b[k];
          double fresh=(v>lambda?(v-lambda)/diag[k]:(v< -lambda?(v+lambda)/diag[k]:0.0));
          double change=fresh-b[k];
          delta=std::max(delta,std::abs(change));
          if(change!=0) {
            b[k]=fresh;
            for(int l=0;l<d;++l) r[l]-=S(l,k)*change;
          }
        }
        // Once a support and its signs stabilize, solve its quadratic exactly.
        // Accept only a sign-preserving candidate, then verify ALL KKT equations.
        // This accelerates correlated blocks without changing the Lasso objective.
        if(cycles%10==0 && delta>tol*.1) {
          std::vector<int> support;
          for(size_t ai=0;ai<active.size();++ai)
            if(b[active[ai]]!=0) support.push_back(active[ai]);
          int m=support.size();
          if(m>1) {
            std::vector<double> H(m*m), rhs(m), cand(m);
            for(int p=0;p<m;++p) {
              rhs[p]=S(support[p],j)-std::copysign(lambda,b[support[p]]);
              for(int q=0;q<m;++q) H[p*m+q]=S(support[p],support[q]);
            }
            bool ok=true;
            // Cholesky factorization of the active principal submatrix.
            for(int p=0;p<m && ok;++p) {
              for(int q=0;q<=p;++q) {
                double val=H[p*m+q];
                for(int h=0;h<q;++h) val-=H[p*m+h]*H[q*m+h];
                if(p==q) {
                  if(val<=1e-14*diag[support[p]]) { ok=false; break; }
                  H[p*m+q]=std::sqrt(val);
                } else H[p*m+q]=val/H[q*m+q];
              }
            }
            if(ok) {
              for(int p=0;p<m;++p) {
                double val=rhs[p];
                for(int q=0;q<p;++q) val-=H[p*m+q]*cand[q];
                cand[p]=val/H[p*m+p];
              }
              for(int p=m-1;p>=0;--p) {
                double val=cand[p];
                for(int q=p+1;q<m;++q) val-=H[q*m+p]*cand[q];
                cand[p]=val/H[p*m+p];
                if(cand[p]*b[support[p]]<=0 || !R_finite(cand[p])) ok=false;
              }
              if(ok) {
                for(int p=0;p<m;++p) b[support[p]]=cand[p];
                delta=0;
              }
            }
          }
        }
        // Examine every excluded coordinate only after restricted convergence.
        if(delta>tol*.1 && !active.empty()) continue;
        violation=0; bool added=false;
        // Recompute residual to remove floating-point accumulation.
        for(int k=0;k<d;++k) {
          double v=S(k,j);
          for(size_t ai=0;ai<active.size();++ai) {
            int l=active[ai]; if(b[l]!=0) v-=S(k,l)*b[l];
          }
          r[k]=v;
        }
        for(int k=0;k<d;++k) if(k!=j) {
          double v=(b[k]==0?std::max(std::abs(r[k])-lambda,0.0):
             std::abs(r[k]-std::copysign(lambda,b[k])));
          violation=std::max(violation,v);
          if(!isactive[k] && v>tol) {
            active.push_back(k); isactive[k]=1; added=true;
          }
        }
        if(violation<=tol && !added) break;
      }
      if(violation>tol || !R_finite(violation))
        Rcpp::stop("Nodewise Lasso failed full KKT convergence at node %d",j+1);
      double linear=0, quadratic=0, l1=0;
      for(int k=0;k<d;++k) if(b[k]!=0) {
        linear+=b[k]*S(k,j); l1+=std::abs(b[k]);
        for(int l=0;l<d;++l) if(b[l]!=0) quadratic+=b[k]*S(k,l)*b[l];
      }
      double t=diag[j]-2*linear+quadratic+lambda*l1;
      if(!R_finite(t) || t<=0) Rcpp::stop("Nonpositive nodewise residual variance");
      tau2[j]=t; kkt[j]=violation; iterations[j]=cycles;
      for(int k=0;k<d;++k) {
        Beta(j,k)=b[k]; Gamma(j,k)=(k==j?1.0:-b[k])/t;
      }
    }
    return Rcpp::List::create(Rcpp::Named("Gamma")=Gamma,
      Rcpp::Named("beta")=Beta,Rcpp::Named("tau2")=tau2,
      Rcpp::Named("kkt")=kkt,Rcpp::Named("iterations")=iterations);
  }')
  invisible(TRUE)
}

nodewise_lasso_reference <- function(A, n, lambda = sqrt(2 * log(nrow(A)) / (n+1)),
                                     tol = 1e-10, maxit = 100000L, warm = NULL) {
  S <- A / (n+1); d <- nrow(S)
  beta <- if (is.null(warm)) matrix(0, d, d) else warm
  Gamma <- matrix(0, d, d); tau2 <- numeric(d); kkt <- numeric(d)
  iterations <- integer(d)
  for (j in seq_len(d)) {
    ii <- setdiff(seq_len(d), j); b <- beta[j, ii]
    if (!length(ii)) { tau2[j] <- S[j,j]; Gamma[j,j] <- 1/tau2[j]; next }
    G <- S[ii, ii, drop=FALSE]; s <- S[ii,j]; gd <- diag(G)
    for (it in seq_len(maxit)) {
      for (k in seq_along(ii)) {
        val <- s[k] - sum(G[k,] * b) + gd[k] * b[k]
        b[k] <- sign(val) * max(abs(val)-lambda, 0) / gd[k]
      }
      g <- drop(G %*% b) - s
      violation <- ifelse(b==0, pmax(abs(g)-lambda,0), abs(g+lambda*sign(b)))
      if (max(violation) <= tol) break
    }
    if(max(violation)>tol) stop("Nodewise Lasso failed full KKT convergence")
    beta[j,ii] <- b
    tau2[j] <- S[j,j] - 2*sum(b*s) + sum(b * (G %*% b)) + lambda*sum(abs(b))
    if (!is.finite(tau2[j]) || tau2[j]<=0) stop("Nonpositive nodewise residual variance")
    Gamma[j,j] <- 1/tau2[j]; Gamma[j,ii] <- -b/tau2[j]
    kkt[j] <- max(violation); iterations[j] <- it
  }
  list(Gamma=Gamma, beta=beta, tau2=tau2, kkt=kkt, iterations=iterations)
}

nodewise_lasso <- function(A, n, lambda = sqrt(2 * log(nrow(A)) / (n+1)),
                           tol = 1e-10, maxit = 100000L, warm = NULL,
                           use_cpp = TRUE) {
  if (!is.matrix(A) || nrow(A)!=ncol(A) || any(!is.finite(A))) stop("Invalid scatter matrix")
  if (n<nrow(A) || lambda<0 || tol<=0) stop("Invalid n, lambda, or tolerance")
  if (use_cpp && .ea_runtime$precision_cpp)
    .ea_runtime$nodewise_lasso_cpp(A/(n+1), lambda, tol, as.integer(maxit), warm)
  else nodewise_lasso_reference(A,n,lambda,tol,maxit,warm)
}

# Optional compilation never changes the RNG stream and never installs packages.
enable_acceleration <- function(backend=c("auto","base","cpp"),verbose=FALSE) {
  backend<-match.arg(backend)
  .ea_runtime$raw_cpp<-.ea_runtime$precision_cpp<-FALSE
  if(backend!="base") {
    attempt<-function(fun)tryCatch(isTRUE(fun()),error=function(e){
      message("Optional C++ unavailable; using base R: ",conditionMessage(e));FALSE})
    .ea_runtime$raw_cpp<-attempt(function()compile_raw_fast(verbose))
    .ea_runtime$precision_cpp<-attempt(precision_engine_compile)
    if(backend=="cpp"&&!all(.ea_runtime$raw_cpp,.ea_runtime$precision_cpp))
      stop("Explicit cpp backend requires Rcpp and a working C++ compiler")
  }
  .ea_runtime$backend<-paste0("raw_",if(.ea_runtime$raw_cpp)"cpp" else "base",
                              "+nodewise_",if(.ea_runtime$precision_cpp)"cpp" else "base")
  invisible(.ea_runtime$backend)
}

# 5. Frozen designs. Indices are permanent: selecting groups never reseeds them.
make_block_cov <- function(sizes,rhos) {
  stopifnot(length(sizes)==length(rhos));Sigma<-diag(sum(sizes));start<-1L
  for(j in seq_along(sizes)) {
    ii<-start:(start+sizes[j]-1L);Sigma[ii,ii]<-rhos[j];diag(Sigma)[ii]<-1
    start<-start+sizes[j]
  }
  chol(Sigma);attr(Sigma,"block_sizes")<-sizes;attr(Sigma,"block_rhos")<-rhos
  Sigma
}
manuscript_settings <- function(which=c("original","confirmation","curves","all")) {
  which<-match.arg(which);original<-primary<-extended<-list()
  for(rho in c(.1,.3,.5,.7,.9)) {
    aa<-list(list(id=paste0("rho=",rho),mu=c(2,rep(0,49)),design="rho",value=rho))
    if(rho==.9)for(a in c(1,1.5,2.5,3))aa[[length(aa)+1L]]<-
      list(id=paste0("a=",a),mu=c(a,rep(0,49)),design="amplitude",value=a)
    original[[paste0("cor_",rho)]]<-list(Sigma=make_block_cov(c(5,45),c(rho,.1)),n=249L,alts=aa)
  }
  for(d in c(50,100,150,200,250))original[[paste0("dim_",d)]]<-list(
    Sigma=make_block_cov(rep(10,d/10),c(.7,rep(.1,d/10-1))),n=3*d-1,
    alts=list(list(id=paste0("d=",d),mu=c(3,rep(.5,9),rep(0,d-10)),design="d",value=d)))
  aa<-list()
  for(s in c(1,2,5,10,20,50))for(layout in c("C","S"))for(signs in c("+","mixed")) {
    if((s==1||s==50)&&(layout=="S"||signs=="mixed"&&s==1))next
    ii<-if(layout=="C")seq_len(s) else as.vector(t(matrix(seq_len(50),nrow=5)))[seq_len(s)]
    mu<-numeric(50);mu[ii]<-1;if(signs=="mixed")mu[sort(ii)]<-rep(c(1,-1),length.out=s)
    aa[[length(aa)+1L]]<-list(id=paste(s,layout,signs,sep="_"),mu=mu,design="geometry",value=s)
  }
  original$geometry<-list(Sigma=make_block_cov(rep(5,10),rep(.9,10)),n=249L,alts=aa)
  for(d in c(20,50))for(rho in c(.9,.95)) {
    ts<-(1-sqrt(1-rho^2))/rho;mu<-c(3.5,3.5*ts,rep(0,d-2))
    aa<-list(list(id="matched",mu=mu,design="matched",value=ts))
    if(d==50&&rho==.95)for(eps in c(-.05,.05)) {
      mm<-mu;mm[2]<-3.5*(ts+eps)
      aa[[length(aa)+1L]]<-list(id=paste0("matched_t_offset_",eps),mu=mm,design="matched",value=ts+eps)
    }
    primary[[paste0("matched_d",d,"_r",rho,"_N",5*d)]]<-
      list(Sigma=make_block_cov(c(2,d-2),c(rho,0)),n=5*d-1,alts=aa)
  }
  ts<-(1-sqrt(1-.95^2))/.95;mu<-c(3.5,3.5*ts,rep(0,48))
  for(n in c(99,499))primary[[paste0("matched_d50_r0.95_N",n+1)]]<-list(
    Sigma=make_block_cov(c(2,48),c(.95,0)),n=n,
    alts=list(list(id="matched",mu=mu,design="matched",value=ts)))
  primary$coherent_d50_r0.95<-list(Sigma=make_block_cov(c(2,48),c(.95,0)),n=249L,
    alts=list(list(id="coherent",mu=c(3.5,3.5,rep(0,48)),design="coherent",value=3.5)))
  primary$identity_d50<-list(Sigma=diag(50),n=249L,
    alts=list(list(id="isolated",mu=c(3.5,rep(0,49)),design="control",value=3.5)))
  primary$isolated_d50_r0.95<-list(Sigma=make_block_cov(c(2,48),c(.95,0)),n=249L,
    alts=list(list(id="isolated",mu=c(1.25,rep(0,49)),design="control",value=1.25)))
  for(d in c(50,100))extended[[paste0("complementary_d",d)]]<-list(
    Sigma=make_block_cov(c(2,2,d-4),c(.95,.9,0)),n=5*d-1,
    alts=lapply(c(.6,.8),function(a)list(id=paste0("coherent_contrast_",a),
      mu=c(3.5,3.5,a,-a,rep(0,d-4)),design="complementary",value=a)))
  extended$matched_d100_r0.95<-list(Sigma=make_block_cov(c(2,98),c(.95,0)),n=499L,
    alts=list(list(id="matched",mu=c(3.75,3.75*ts,rep(0,98)),design="matched",value=ts)))
  register<-function(specs,stage,master,directory) {
    for(j in seq_along(specs))specs[[j]]<-c(specs[[j]],list(stage=stage,
      master_seed=master,index=j,source_dir=directory,source_id=names(specs)[j]))
    names(specs)<-vapply(specs,function(g){id<-g$source_id
      if(grepl("_N[0-9]+$",id))id<-sub("_N[0-9]+$",paste0("_n",g$n),id)
      id},"");specs
  }
  oo<-register(original,"original",20261011,"reproduction_full")
  cc<-c(register(primary,"confirmation",2026101031,"confirmation"),
        register(extended,"confirmation",2026101041,"confirmation_extended"))
  curves<-list()
  curve<-function(d,rho,amplitudes)list(
    Sigma=make_block_cov(c(2,d-2),c(rho,0)),n=5*d-1,
    alts=lapply(amplitudes,function(a)list(id=paste0("matched_a",a),
      mu=c(a,a*(1-sqrt(1-rho^2))/rho,rep(0,d-2)),
      design="amplitude",value=a)),rho=rho)
  for(rho in c(.5,.7,.8,.9,.95))curves[[paste0("main_rho",100*rho)]]<-
    curve(50,rho,if(rho==.95)c(2.5,3,3.5,4,4.5) else 3.5)
  for(d in c(20,100,150,200))curves[[paste0("main_d",d)]]<-curve(d,.95,3.5)
  mm<-register(curves,"curves",2026101407,NA_character_)
  switch(which,original=oo,confirmation=cc,curves=mm,all=c(oo,cc,mm))
}
settings_table <- function(settings=manuscript_settings("all")) {
  rows<-list();j<-0L
  for(id in names(settings)) {
    g<-settings[[id]];aa<-g$alts
    if(g$stage=="curves") {
      if(id=="main_rho95") {
        base<-aa[[which(vapply(aa,function(a)a$mu[1]==3.5,FALSE))]]
        ar<-base;ar$design<-"rho";ar$value<-.95
        ad<-base;ad$design<-"d";ad$value<-50
        aa<-c(aa,list(ar,ad))
      } else {
        aa[[1]]$design<-if(startsWith(id,"main_rho"))"rho" else "d"
        aa[[1]]$value<-if(startsWith(id,"main_rho"))g$rho else nrow(g$Sigma)
      }
    }
    if(id=="cor_0.9")aa<-c(aa,list(list(id="rho=0.9",mu=aa[[1]]$mu,design="amplitude",value=2)))
    for(al in aa){j<-j+1L;rows[[j]]<-data.frame(group=id,source_group=g$source_id,
      stage=g$stage,figure_set=if(g$stage=="curves")"main" else g$stage,
      scenario=al$id,sweep=al$design,x=al$value,n=g$n,d=nrow(g$Sigma),
      mean=paste(al$mu,collapse=","),seed_cal=g$master_seed+10000*g$index,
      seed_eval=g$master_seed+10000*g$index+1,stringsAsFactors=FALSE)}
  };do.call(rbind,rows)
}

# 6. Gaussian computation preserves the legacy symmetric block root and draw order.
# Centering these n+1 rows yields independent X,A plus the ancillary orientation.
gaussian_draw <- function(ctx) {
  rows<-ctx$n+1L;E<-matrix(rnorm(rows*ctx$d),rows,ctx$d)
  sizes<-ctx$block_sizes;rhos<-ctx$block_rhos
  if(is.null(sizes))return(E%*%t(ctx$L))
  start<-1L
  for(j in seq_along(sizes)) {
    b<-sizes[j];ii<-start:(start+b-1L);a<-sqrt(1-rhos[j]);c<-sqrt(1+(b-1)*rhos[j])
    block<-E[,ii,drop=FALSE];E[,ii]<-a*block+(c-a)*rowMeans(block);start<-start+b
  };E
}
prepare_draw <- function(Y0,ctx) {
  rp<-raw_prepare(Y0)
  list(raw=rp,X0=sqrt(ctx$n+1)*rp$mean0*rp$scale,
       unknown=unknown_prepare(rp$A*rp$scale^2,ctx$n),ctx=ctx)
}
score_draw <- function(prepared,alternatives) {
  do.call(rbind,lapply(alternatives,function(al){
    X<-prepared$X0+al$mu;rr<-raw_from_prepared(prepared$raw,al$mu)
    kn<-c(known_scores(X,prepared$ctx),rr)[manuscript_methods]
    un<-c(unknown_scores(X,prepared$unknown),rr)[setdiff(manuscript_methods,"IHC")]
    c(setNames(kn,paste0("known:",names(kn))),setNames(un,paste0("unknown:",names(un))))
  }))
}
rank_cutoff <- function(bank,alpha=.05) {
  B<-nrow(bank);k<-B+1-floor(alpha*(B+1))
  if(k>B)return(setNames(rep(Inf,ncol(bank)),colnames(bank)))
  setNames(vapply(seq_len(ncol(bank)),function(j)sort(bank[,j],partial=k)[k],0),colnames(bank))
}
wilson <- function(count,R,alpha=.05) {
  z<-qnorm(1-alpha/2);p<-count/R;a<-1+z*z/R
  c(low=(p+z*z/(2*R)-z*sqrt(p*(1-p)/R+z*z/(4*R*R)))/a,
    high=(p+z*z/(2*R)+z*sqrt(p*(1-p)/R+z*z/(4*R*R)))/a)
}
# Paired SE uses discordant decisions, not two independent binomial SEs.
paired_intervals <- function(delta,squared_difference,R,family_size=1L,alpha=.05) {
  se<-sqrt(pmax(0,(squared_difference-R*delta^2)/(R-1))/R)
  z<-qnorm(1-alpha/2);zs<-qnorm(1-alpha/(2*family_size))
  c(delta=delta,se=se,low=delta-z*se,high=delta+z*se,
    low_simultaneous=delta-zs*se,high_simultaneous=delta+zs*se)
}
write_result_tables <- function(result,outdir) {
  dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
  for(nm in c("power","paired","analytic_EA","settings"))
    write.csv(result[[nm]],file.path(outdir,paste0(nm,".csv")),row.names=FALSE)
  invisible(result)
}

# Parallel replicates: at most four workers plus the coordinating R process.
# Workers use local sockets, not internet services, and create no result/log files.
.ea_check_workers <- function(workers) {
  if(!is.numeric(workers)||is.complex(workers)||length(workers)!=1L||
     !is.finite(workers)||workers!=floor(workers)||workers<1L||workers>4L)
    stop('workers must be 1--4; four workers plus R stay within the five-core budget')
  as.integer(workers)
}
.ea_score_code <- function() {
  names <- c('log_two_tail','hc_score','ea_logscore','make_context','known_scores',
    'unknown_prepare','unknown_scores','raw_prepare','.raw_e2','.raw_e6',
    '.raw_assert_positive','.raw_shifted_moments','raw_from_prepared',
    'nodewise_lasso_reference','nodewise_lasso','gaussian_draw','prepare_draw','score_draw',
    'compile_raw_fast','precision_engine_compile')
  vapply(names,function(nm)paste(deparse(get(nm,envir=environment(.ea_score_code)),
    width.cutoff=500),collapse='\n'),'')
}
.ea_worker_init <- function(path,backend) {
  sys.source(path,envir=.GlobalEnv)
  .GlobalEnv$enable_acceleration(backend)
  list(backend=.GlobalEnv$.ea_runtime$backend,code=.GlobalEnv$.ea_score_code())
}
environment(.ea_worker_init) <- baseenv()
.ea_worker_context <- function(ctx,alternatives) {
  .GlobalEnv$.ea_runtime$parallel_ctx <- ctx
  .GlobalEnv$.ea_runtime$parallel_alts <- alternatives
  NULL
}
.ea_worker_score <- function(Y) {
  rt <- .GlobalEnv$.ea_runtime
  .GlobalEnv$score_draw(.GlobalEnv$prepare_draw(Y,rt$parallel_ctx),rt$parallel_alts)
}
.ea_make_cluster <- function(workers,backend) {
  workers <- .ea_check_workers(workers)
  path <- .ea_runtime$source_file
  if(is.null(path)||length(path)!=1L||!file.exists(path))
    stop('Cannot locate the simulation engine file. Rerun SETUP before parallel simulations.')
  # Limit native numerical libraries in each child to one thread.
  variables <- c('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS',
                 'VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS')
  previous <- Sys.getenv(variables,unset=NA_character_)
  on.exit({
    missing <- is.na(previous)
    if(any(missing))Sys.unsetenv(variables[missing])
    if(any(!missing))do.call(Sys.setenv,as.list(previous[!missing]))
  },add=TRUE)
  do.call(Sys.setenv,as.list(setNames(rep('1',length(variables)),variables)))
  cl <- tryCatch(parallel::makePSOCKcluster(workers,outfile='',rscript_args='--vanilla'),
    error=function(e)stop('Could not start local parallel workers: ',conditionMessage(e),
      '. Use workers <- 1 for serial execution.',call.=FALSE))
  ready <- FALSE
  on.exit(if(!ready)parallel::stopCluster(cl),add=TRUE)
  info <- parallel::clusterCall(cl,.ea_worker_init,normalizePath(path),backend)
  code <- .ea_score_code()
  if(any(!vapply(info,function(x)identical(x$code,code),FALSE)))
    stop('The saved engine differs from the loaded scoring functions. Save the file and rerun SETUP.')
  if(any(vapply(info,function(x)x$backend!=.ea_runtime$backend,FALSE))) {
    if(backend=='cpp')stop('Not all workers could initialize the requested C++ backend')
    enable_acceleration('base')
    parallel::clusterCall(cl,.ea_worker_init,normalizePath(path),'base')
  }
  ready <- TRUE
  cl
}
.ea_draw_batch <- function(count,ctx,alternatives,cluster=NULL) {
  if(is.null(cluster))return(lapply(seq_len(count),function(i)
    score_draw(prepare_draw(gaussian_draw(ctx),ctx),alternatives)))
  # All Gaussian draws are made in their original order in the coordinator.
  # Scoring uses no randomness, so worker count cannot change the seeded draws.
  draws <- lapply(seq_len(count),function(i)gaussian_draw(ctx))
  parallel::parLapplyLB(cluster,draws,'.ea_worker_score',chunk.size=1L)
}
.ea_batch_size <- function(ctx,cluster) {
  if(is.null(cluster))return(1L)
  # Bound queued Gaussian matrices to about 64 MiB in the coordinator.
  max(1L,min(2L*length(cluster),floor(64*1024^2/(8*(ctx$n+1)*ctx$d))))
}
.ea_progress <- function(id,phase,done,total,previous,progress) {
  if(progress && (done==total || floor(10*done/total)>floor(10*previous/total)))
    cat(id,phase,done,'/',total,'\n')
}

# Hash score/summary code as well as design, seeds, budgets and backend. A cache
# cannot silently survive a changed statistic, normalization or summary formula.
.manuscript_signature <- function() {
  functions<-c("log_two_tail","hc_score","ea_logscore","make_context","known_scores",
    "unknown_prepare","unknown_scores","raw_prepare",".raw_e2",".raw_e6",
    ".raw_shifted_moments","raw_from_prepared","compile_raw_fast",
    "precision_engine_compile","nodewise_lasso_reference","nodewise_lasso",
    "gaussian_draw","prepare_draw","score_draw","wilson","paired_intervals","run_setting",
    ".ea_draw_batch",".ea_batch_size",".ea_worker_context",".ea_worker_score")
  code<-vapply(functions,function(nm)paste(deparse(get(nm,
    envir=environment(.manuscript_signature)),width.cutoff=500),collapse="\n"),"")
  file<-tempfile();on.exit(unlink(file));writeLines(code,file)
  unname(tools::md5sum(file))
}
# 7. One simulation group. outdir=NULL keeps all results in memory.
# Supplying a directory explicitly enables completed-group caches.
# Independent calibration/evaluation seeds and fixed group indices are recorded.
run_setting <- function(setting,id,B=19999L,R=10000L,outdir=NULL,
                        backend=c("auto","base","cpp"),alpha=.05,
                        save_scores=FALSE,resume=TRUE,progress=TRUE,
                        workers=1L,cluster=NULL) {
  backend<-match.arg(backend);stopifnot(B>=1,R>=2,alpha>0,alpha<1,
    B==as.integer(B),R==as.integer(R))
  persist<-!is.null(outdir)
  if(persist && (!is.character(outdir)||length(outdir)!=1L||is.na(outdir)||!nzchar(outdir)))
    stop("outdir must be NULL (memory only) or one nonempty directory name")
  B<-as.integer(B);R<-as.integer(R)
  workers<-.ea_check_workers(if(is.null(cluster))workers else length(cluster))
  if(!grepl("^[A-Za-z0-9_.-]+$",id)||id%in%c(".",".."))stop("Use a simple group identifier")
  if(anyDuplicated(vapply(setting$alts,`[[`,"","id")))stop("Alternative ids must be unique")
  enable_acceleration(backend);ctx<-make_context(setting$Sigma,setting$n)
  if(is.null(cluster)&&workers>1L) {
    cluster<-.ea_make_cluster(workers,backend)
    on.exit(parallel::stopCluster(cluster),add=TRUE)
  }
  if(any(vapply(setting$alts,function(a)length(a$mu)!=ctx$d,FALSE)))stop("Mean dimension mismatch")
  seed_cal<-setting$master_seed+10000*setting$index;seed_eval<-seed_cal+1
  config<-list(version="canonical-nine-v1-qcq-corrected",
    implementation=if(persist).manuscript_signature() else NA_character_,setting=setting,
    B=B,R=R,alpha=alpha,seed_cal=seed_cal,seed_eval=seed_eval,
    backend=.ea_runtime$backend,save_scores=save_scores)
  if(persist) {
    directory<-file.path(outdir,"groups");dir.create(directory,recursive=TRUE,showWarnings=FALSE)
    file<-file.path(directory,paste0(id,"_result.rds"))
    if(resume&&file.exists(file)) {
      prior<-readRDS(file)
      if(!identical(prior$config,config))stop("Cache configuration differs for ",id,"; use a different output directory")
      return(prior)
    }
  }
  RNGkind("Mersenne-Twister","Inversion","Rejection")
  null<-list(id="null",mu=numeric(ctx$d));aa<-c(list(null),setting$alts)
  batch_size<-.ea_batch_size(ctx,cluster)
  if(progress)cat('Start',id,'d',ctx$d,'n',ctx$n,'B',B,'R',R,
    'simulation workers',if(is.null(cluster))1L else length(cluster),'\n')
  if(!is.null(cluster))parallel::clusterCall(cluster,'.ea_worker_context',ctx,list(null))
  set.seed(seed_cal);initial_cal<-.Random.seed;bank<-NULL
  for(first in seq.int(1L,B,by=batch_size)) {
    last<-min(B,first+batch_size-1L)
    batch<-.ea_draw_batch(last-first+1L,ctx,list(null),cluster)
    for(k in seq_along(batch)) {
      ss<-batch[[k]]
      if(is.null(bank))bank<-matrix(NA_real_,B,ncol(ss),dimnames=list(NULL,colnames(ss)))
      bank[first+k-1L,]<-ss[1,]
    }
    .ea_progress(id,'calibration',last,B,first-1L,progress)
  }
  final_cal<-.Random.seed;cut<-rank_cutoff(bank,alpha)
  ids<-vapply(aa,`[[`,"","id");labels<-colnames(bank);RNG_cal<-.Random.seed
  counts<-matrix(0,length(aa),length(labels),dimnames=list(ids,labels))
  discordant<-matrix(0,length(aa),length(labels),dimnames=list(ids,labels))
  analytic<-matrix(0,length(aa),2,dimnames=list(ids,c("known","unknown")))
  scores<-if(save_scores)array(NA_real_,c(R,length(aa),length(labels)),dimnames=list(NULL,ids,labels)) else NULL
  if(!is.null(cluster))parallel::clusterCall(cluster,'.ea_worker_context',ctx,aa)
  set.seed(seed_eval);initial_eval<-.Random.seed
  for(first in seq.int(1L,R,by=batch_size)) {
    last<-min(R,first+batch_size-1L)
    batch<-.ea_draw_batch(last-first+1L,ctx,aa,cluster)
    for(k in seq_along(batch)) {
      i<-first+k-1L;ss<-batch[[k]]
      reject<-sweep(ss,2,cut,">");counts<-counts+reject
      if(save_scores)scores[i,,]<-ss
      for(panel in c("known","unknown")) {
        ii<-which(startsWith(labels,paste0(panel,":")));ee<-match(paste0(panel,":EA"),labels)
        discordant[,ii]<-discordant[,ii]+(reject[,ii,drop=FALSE]-reject[,ee])^2
        analytic[,panel]<-analytic[,panel]+(ss[,ee]>=-log(alpha))
      }
    }
    .ea_progress(id,'evaluation',last,R,first-1L,progress)
  }
  power<-paired<-diagnostic<-list();j<-k<-h<-0L
  for(a in seq_along(aa))for(panel in c("known","unknown")) {
    ii<-which(startsWith(labels,paste0(panel,":")));ee<-match(paste0(panel,":EA"),labels)
    for(l in ii) {
      j<-j+1L;ci<-wilson(counts[a,l],R)
      power[[j]]<-data.frame(stage=setting$stage,group=id,scenario=ids[a],panel=panel,
        method=sub("^[^:]+:","",labels[l]),R=R,rejections=counts[a,l],power=counts[a,l]/R,
        mcse=sqrt((counts[a,l]/R)*(1-counts[a,l]/R)/R),low=ci[1],high=ci[2],
        B=B,alpha=alpha,n=setting$n,d=ctx$d,
        calibration=if(panel=="known")"known complete-score rank" else "oracle complete-score rank",row.names=NULL)
      if(l!=ee) {
        k<-k+1L;ci<-paired_intervals((counts[a,ee]-counts[a,l])/R,discordant[a,l],R)
        paired[[k]]<-data.frame(stage=setting$stage,group=id,scenario=ids[a],panel=panel,
          comparison=paste("EA minus",sub("^[^:]+:","",labels[l])),R=R,t(ci),family_size=1,row.names=NULL)
      }
    }
    h<-h+1L;ci<-wilson(analytic[a,panel],R)
    diagnostic[[h]]<-data.frame(group=id,scenario=ids[a],panel=panel,method="EA-analytic",
      R=R,power=analytic[a,panel]/R,low=ci[1],high=ci[2],row.names=NULL)
  }
  result<-list(config=config,power=do.call(rbind,power),paired=do.call(rbind,paired),
    analytic_EA=do.call(rbind,diagnostic),settings=settings_table(setNames(list(setting),id)),
    cutoff=cut,counts=counts,discordant=discordant,
    RNG=list(kind=RNGkind(),initial_cal=initial_cal,final_cal=final_cal,
             initial_eval=initial_eval,final_eval=.Random.seed),
    backend=.ea_runtime$backend,session=sessionInfo(),
    parallel=list(workers=if(is.null(cluster))0L else length(cluster),
      coordinator=1L,scheme='Original ordered Gaussian draws; parallel deterministic scoring'))
  if(save_scores){result$bank<-bank;result$scores<-scores}
  if(persist) {
    temporary<-paste0(file,".tmp-",Sys.getpid());saveRDS(result,temporary)
    if(!file.rename(temporary,file))stop("Could not save completed group cache")
  }
  if(progress)cat("Complete",id,"B",B,"R",R,"\n")
  result
}
# Simultaneous intervals are conditional on the realized calibration banks.
# Fresh/imported original comparisons are pointwise (family1); confirmation
# retains conservative legacy432 and curves195, including in all mode.
adjust_paired <- function(paired,family_size=NULL,alpha=.05) {
  if(is.null(family_size))family_size<-sum(paired$scenario!="null")
  family_size<-max(1L,family_size);z<-qnorm(1-alpha/(2*family_size))
  paired$low_simultaneous<-paired$delta-z*paired$se
  paired$high_simultaneous<-paired$delta+z*paired$se
  paired$family_size<-family_size;paired
}
run_manuscript <- function(which=c("original","confirmation","curves","all"),
                           budget=c("quick","full"),group_ids=NULL,
                           workers=2L,outdir="manuscript_results",backend="auto",
                           B=NULL,R=NULL,family_size=NULL,save_scores=FALSE,resume=TRUE) {
  which<-match.arg(which);budget<-match.arg(budget)
  if(is.null(B))B<-if(budget=="full")19999L else 199L
  if(is.null(R))R<-if(budget=="full")10000L else 200L
  settings<-manuscript_settings(which)
  if(!is.null(group_ids)){if(any(!group_ids%in%names(settings)))stop("Unknown group id")
    settings<-settings[group_ids]}
  if(.Platform$OS.type=="windows")workers<-1L
  enable_acceleration(backend)
  ans<-parallel::mclapply(seq_along(settings),function(j)
    run_setting(settings[[j]],names(settings)[j],B,R,outdir,backend,
                save_scores=save_scores,resume=resume),mc.cores=workers,mc.set.seed=FALSE)
  if(any(vapply(ans,inherits,FALSE,"try-error")))stop("A group failed; inspect its error")
  paired<-do.call(rbind,lapply(ans,`[[`,"paired"))
  if(is.null(family_size))paired<-do.call(rbind,lapply(unique(paired$stage),function(stage) {
    adjust_paired(paired[paired$stage==stage,],switch(stage,original=1L,confirmation=432L,curves=195L))
  })) else paired<-adjust_paired(paired,family_size)
  result<-list(power=do.call(rbind,lapply(ans,`[[`,"power")),
    paired=paired,
    analytic_EA=do.call(rbind,lapply(ans,`[[`,"analytic_EA")),settings=settings_table(settings))
  write_result_tables(result,outdir)
}

# 8. Import already-completed full results; this never simulates observations.
# Historical source_id keys locate the old files; public group names use n.
source_evaluation <- function(legacy_root="outputs",outdir="manuscript_results/source",
                              include_curves=TRUE) {
  settings<-manuscript_settings("all")
  if(!include_curves)settings<-settings[vapply(settings,function(g)g$stage!="curves",FALSE)]
  pp<-qq<-aa<-manifest<-list();used<-list();j<-0L
  fill<-function(x,fields){for(k in fields)if(!k%in%names(x))x[[k]]<-NA;x[,fields,drop=FALSE]}
  pfields<-c("stage","group","scenario","panel","method","R","rejections","power","mcse",
             "low","high","B","alpha","n","d","calibration")
  qfields<-c("stage","group","scenario","panel","comparison","R","delta","se","low","high",
             "low_simultaneous","high_simultaneous","family_size","bound_origin")
  for(id in names(settings)) {
    g<-settings[[id]];directory<-if(g$stage=="curves")"manuscript_curves_full" else g$source_dir
    folder<-file.path(legacy_root,directory)
    aggregate<-if(g$stage=="confirmation")file.path(legacy_root,"confirmation_all") else folder
    powerfile<-if(file.exists(file.path(aggregate,"power.csv")))file.path(aggregate,"power.csv") else
      file.path(folder,paste0(g$source_id,"_power.csv"))
    pairfile<-if(file.exists(file.path(aggregate,"paired.csv")))file.path(aggregate,"paired.csv") else
      file.path(folder,paste0(g$source_id,"_paired.csv"))
    if(!file.exists(powerfile)||!file.exists(pairfile)) {
      if(g$stage=="curves")next
      stop("Missing completed source results for ",id,": ",powerfile)
    }
    if(normalizePath(outdir,mustWork=FALSE)%in%normalizePath(c(folder,aggregate),mustWork=FALSE))
      stop("Source-eval output must differ from the source directory")
    p<-read.csv(powerfile,stringsAsFactors=FALSE);p<-p[p$group==g$source_id & p$method%in%manuscript_methods,]
    q<-read.csv(pairfile,stringsAsFactors=FALSE)
    q<-q[q$group==g$source_id & sub("^EA minus ","",q$comparison)%in%setdiff(manuscript_methods,"EA"),]
    if(!nrow(p)||!nrow(q))stop("No requested source rows for ",id)
    if(any(p$R!=10000L)||nrow(p)!=(length(g$alts)+1L)*17L)
      stop("Source-eval requires the completed full R=10000 nine-method rows for ",id)
    if(anyDuplicated(p[,c("scenario","panel","method")]))stop("Duplicate source power keys for ",id)
    bankfile<-file.path(folder,paste0(g$source_id,"_bank.rds"));bank_verified<-FALSE
    if(file.exists(bankfile)) {
      archived<-readRDS(bankfile)
      if(archived$B!=19999L||archived$seed!=g$master_seed+10000*g$index)
        stop("Source bank budget/seed differs from the frozen design for ",id)
      bank_verified<-TRUE;rm(archived)
    }
    if(any(p$panel=="unknown" & p$method=="IHC"))stop("Unexpected estimated IHC")
    p$group<-q$group<-id;p$stage<-q$stage<-g$stage
    p$B<-19999L;p$alpha<-.05;p$n<-g$n;p$d<-nrow(g$Sigma)
    p$mcse<-sqrt(p$power*(1-p$power)/p$R)
    p$calibration<-ifelse(p$panel=="known","known complete-score rank","oracle complete-score rank")
    q$R<-vapply(seq_len(nrow(q)),function(i)unique(p$R[p$scenario==q$scenario[i]&p$panel==q$panel[i]])[1],0)
    if(g$stage=="curves") {
      q<-adjust_paired(q,195L);q$bound_origin<-"requested 195-comparison curve family"
    } else if(g$stage=="confirmation") {
      # Use the complete existing 432 family, including omitted sensitivities.
      q<-adjust_paired(q,432L);q$bound_origin<-"conservative legacy 432-comparison family"
    } else {q<-adjust_paired(q,1L);q$bound_origin<-"legacy pointwise paired interval"}
    analyticfile<-file.path(aggregate,"analytic_EA.csv")
    if(file.exists(analyticfile)) {
      a<-read.csv(analyticfile,stringsAsFactors=FALSE);a<-a[a$group==g$source_id,]
      if(nrow(a)){a$group<-id;aa[[length(aa)+1L]]<-a}
    }
    j<-j+1L;pp[[j]]<-fill(p,pfields);qq[[j]]<-fill(q,qfields);used[[id]]<-g
    manifest[[j]]<-data.frame(group=id,source_group=g$source_id,stage=g$stage,n=g$n,
      B=19999L,R=paste(unique(p$R),collapse=","),seed_cal=g$master_seed+10000*g$index,
      seed_eval=g$master_seed+10000*g$index+1,
      power_file=powerfile,power_md5=unname(tools::md5sum(powerfile)),
      paired_file=pairfile,paired_md5=unname(tools::md5sum(pairfile)),
      bank_budget_and_seed_verified=bank_verified,
      evaluation="imported completed results; no fresh Monte Carlo",stringsAsFactors=FALSE)
  }
  if(!length(pp))stop("No completed results found")
  result<-list(power=do.call(rbind,pp),paired=do.call(rbind,qq),
    analytic_EA=if(length(aa))do.call(rbind,aa) else data.frame(),settings=settings_table(used))
  write_result_tables(result,outdir)
  write.csv(do.call(rbind,manifest),file.path(outdir,"source_manifest.csv"),row.names=FALSE)
  generate_tables(result,outdir);invisible(result)
}
generate_tables <- function(result,outdir="manuscript_results") {
  if(is.character(result))result<-list(power=read.csv(file.path(result,"power.csv")))
  p<-result$power;keys<-unique(p[,c("stage","group","scenario","panel","R")])
  for(method in manuscript_methods)keys[[method]]<-vapply(seq_len(nrow(keys)),function(i){
    z<-p[p$group==keys$group[i] & p$scenario==keys$scenario[i] &
           p$panel==keys$panel[i] & p$method==method,"power"]
    if(length(z))z[1] else NA_real_},0)
  dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
  write.csv(keys,file.path(outdir,"power_wide.csv"),row.names=FALSE)
  invisible(keys)
}
# Edit a design safely: copy one spec, edit n/Sigma/each mu together, and use
# a NEW output directory. E.g. g<-manuscript_settings('curves')$main_rho95;
# g$n<-499; g$alts[[3]]$mu[1:2]<-c(3.5,2.53); run_setting(g,'custom',...).
# The covariance dimension must equal every mean length and satisfy n>=d.
# Reusing an edited design in an existing cache directory is rejected.

# 9. Command line. --mode source-eval uses completed files; other modes draw.
main <- function(args=commandArgs(trailingOnly=TRUE)) {
  opts<-list(mode="help",budget="quick",out="manuscript_results",workers="2",
             backend="auto",source="outputs",groups=NULL,`save-figures`="false")
  if(length(args)%%2L)stop("Options are --name value pairs")
  if(length(args))for(i in seq(1L,length(args),by=2L)) {
    key<-sub("^--","",args[i]);if(!key%in%names(opts))stop("Unknown option: ",args[i])
    opts[[key]]<-args[i+1L]
  }
  if(opts$mode=="help") {
    cat("Modes: original, confirmation, curves, all, source-eval, figures\n",
        "Options: --budget quick|full --workers 2 --out DIR --backend auto|base|cpp\n",
        "         --groups id1,id2 --source LEGACY_OUTPUT_ROOT\n",
        "         --save-figures false|true (default false; optional PNG export)\n",
        "Full B=19999,R=10000. Quick B=199,R=200 is a smoke run only.\n",
        "Source-eval imports existing full CSVs and completed curves without simulation.\n")
    return(invisible(NULL))
  }
  if(!opts[["save-figures"]]%in%c("false","true"))
    stop("--save-figures must be false or true")
  if(opts$mode=="source-eval")result<-source_evaluation(opts$source,opts$out) else
  if(opts$mode=="figures") {
    result<-list(power=read.csv(file.path(opts$source,"power.csv")),
                 settings=read.csv(file.path(opts$source,"settings.csv")))
  } else result<-run_manuscript(opts$mode,opts$budget,
    if(is.null(opts$groups))NULL else strsplit(opts$groups,",",fixed=TRUE)[[1]],
    as.integer(opts$workers),opts$out,opts$backend)
  generate_tables(result,opts$out)
  if(opts[["save-figures"]]=="true" &&
     exists("generate_manuscript_figures",mode="function")) {
    budgets<-sort(unique(result$power$R))
    if(opts$mode%in%c("source-eval","figures"))plot_B<-19999L else
      plot_B<-if(opts$budget=="full")19999L else 199L
    generate_manuscript_figures(result$power,result$settings,file.path(opts$out,"figures"),
      B=plot_B,R=if(length(budgets)==1L)budgets else 10000L)
  }
  invisible(result)
}
# The optional figure module is inserted before this guard in the final file.


# 10. Publication figure helpers (no score calculation or simulation).
# Manuscript plotting only: no simulation, score calculation, fitting or cache mutation.
# All panels use the same nine-method style key. Figure captions contain the full
# textual legend; the PNGs contain axes/curves/Wilson whiskers only.

manuscript_method_styles <- function(panel = NULL, line_width = NULL, point_size = NULL) {
  z <- data.frame(
    method = c('EA','Z-Max','Z-TopK','IHC','Feng-COM','QCQ','ZWS-{2,6}',
               'Min-p','Quadratic'),
    color = c('#111111','#0072B2','#D55E00','#6F4C9B','#009E73',
              '#E69F00','#CC79A7','#56B4E9','#666666'),
    color_name = c('black','blue','vermilion','purple','green','orange',
                   'reddish purple','sky blue','gray'),
    pch = c(18L,17L,25L,16L,15L,8L,7L,1L,4L),
    point_type = c('filled diamond','filled upward triangle','filled downward triangle',
                   'filled circle','filled square','asterisk','crossed square',
                   'open circle','cross'),
    lty = c(1L,1L,1L,2L,2L,2L,3L,3L,4L),
    stringsAsFactors = FALSE)
  z$pgf_mark <- c('diamond*','triangle*','triangle*','*','square*',
                  'manuscript-asterisk8','manuscript-squarex','o','x')
  z$pgf_rotation <- c(0L,0L,180L,0L,0L,0L,0L,0L,0L)
  z$pgf_line_style <- c('solid','solid','solid','dashed','dashed','dashed',
                        'dotted','dotted','dashdotted')
  z$line_width <- ifelse(z$method=='EA',2.8,1.9)
  z$point_size <- ifelse(z$method=='EA',.76,.69)
  # An explicit width overrides every method curve at plotting time.
  if(!is.null(line_width)) {
    if(!is.numeric(line_width) || is.complex(line_width) || length(line_width)!=1L ||
       !is.finite(line_width) || line_width<=0)
      stop('line_width must be NULL or one positive finite number')
    z$line_width <- as.numeric(line_width)
  }
  # point_size is the marker cex scale; it does not change text sizes.
  if(!is.null(point_size)) {
    if(!is.numeric(point_size) || is.complex(point_size) || length(point_size)!=1L ||
       !is.finite(point_size) || point_size<=0)
      stop('point_size must be NULL or one positive finite number')
    z$point_size <- as.numeric(point_size)
  }
  z$caption_label <- paste0(z$method,' (',z$point_type,', ',z$color_name,')')
  if(!is.null(panel)) {
    panel <- match.arg(panel,c('known','unknown'))
    if(panel=='unknown') z <- z[z$method!='IHC',,drop=FALSE]
  }
  rownames(z) <- NULL
  z
}

# Existing original-design archive mapping, with canonical Wishart degrees n.
# The a=2 display point intentionally shares the rho=.9 source group/scenario.
# Users of another registry supply its matching settings data frame directly.
manuscript_original_sweep_settings <- function() {
  rhos <- c(.1,.3,.5,.7,.9)
  dims <- c(50L,100L,150L,200L,250L)
  amps <- c(1,1.5,2,2.5,3)
  rbind(
    data.frame(group=paste0('cor_',rhos),scenario=paste0('rho=',rhos),
      sweep='rho',x=rhos,n=249L,d=50L,stringsAsFactors=FALSE),
    data.frame(group=paste0('dim_',dims),scenario=paste0('d=',dims),
      sweep='d',x=dims,n=3L*dims-1L,d=dims,stringsAsFactors=FALSE),
    data.frame(group='cor_0.9',scenario=ifelse(amps==2,'rho=0.9',paste0('a=',amps)),
      sweep='amplitude',x=amps,n=249L,d=50L,stringsAsFactors=FALSE))
}

plot_manuscript_curves <- function(power, settings, panel, outpath = NULL,
    B = 19999L, alpha = .05, R = 10000L, show_intervals = TRUE,
    width = 3600L, height = 1350L, res = 450L, pointsize = 10,
    sweeps = c('rho','d','amplitude'), line_width = NULL, point_size = NULL) {
  panel <- match.arg(panel,c('known','unknown'))
  if(is.character(power) && length(power)==1L)
    power <- read.csv(power,stringsAsFactors=FALSE,check.names=FALSE)
  if(is.character(settings) && length(settings)==1L)
    settings <- read.csv(settings,stringsAsFactors=FALSE,check.names=FALSE)
  # Accept canonical manuscript registry names and existing archive names.
  if(!'scenario' %in% names(power) && 'alt' %in% names(power)) power$scenario <- power$alt
  if(!'scenario' %in% names(settings) && 'alt' %in% names(settings)) settings$scenario <- settings$alt
  if(!'sweep' %in% names(settings) && 'curve' %in% names(settings)) settings$sweep <- settings$curve
  if(!'R' %in% names(power)) power$R <- R
  req_power <- c('group','scenario','panel','method','R','power','low','high')
  req_settings <- c('group','scenario','sweep','x')
  if(!is.data.frame(power) || !all(req_power %in% names(power)))
    stop('power must be a long data frame/CSV with ',paste(req_power,collapse=', '))
  if(!is.data.frame(settings) || !all(req_settings %in% names(settings)))
    stop('settings must be a data frame/CSV with ',paste(req_settings,collapse=', '))
  save_png <- !is.null(outpath)
  if(save_png && (length(outpath)!=1L || is.na(outpath) ||
     !grepl('\\.png$',outpath,ignore.case=TRUE))) stop('outpath must end in .png')
  if(!save_png && !interactive() && grDevices::dev.cur()==1L)
    stop('For Rscript, supply outpath="figure.png"; outpath=NULL displays on an active graphics device.')
  if(anyDuplicated(power[,c('group','scenario','panel','method')]))
    stop('Duplicated group/scenario/panel/method score-summary keys')
  if(any(power$panel=='unknown' & power$method=='IHC' & is.finite(power$power)))
    stop('IHC is unavailable for estimated covariance; no invented variant may be plotted')
  aliases <- c(correlation='rho',dimension='d',a='amplitude')
  hit <- settings$sweep %in% names(aliases)
  settings$sweep[hit] <- unname(aliases[settings$sweep[hit]])
  sweeps <- as.character(sweeps)
  hit <- sweeps %in% names(aliases);sweeps[hit] <- unname(aliases[sweeps[hit]])
  if(!length(sweeps) || anyNA(sweeps) || anyDuplicated(sweeps) ||
     any(!sweeps %in% c('rho','d','amplitude')))
    stop('sweeps must select distinct rho, d or amplitude panels')
  settings <- settings[settings$sweep %in% sweeps,,drop=FALSE]
  if(!setequal(unique(settings$sweep),sweeps)) stop('settings are missing a requested sweep')
  if(anyDuplicated(settings[,c('sweep','x')])) stop('Each sweep/x display key must be unique')
  if(any(!is.finite(settings$x))) stop('Require finite sampled coordinates')
  if(xor('n' %in% names(settings),'d' %in% names(settings)))
    stop('When supplying degree metadata, include both n and d')
  if(all(c('n','d') %in% names(settings)) &&
     (any(!is.finite(settings$n)) || any(!is.finite(settings$d)) ||
      any(settings$n < settings$d) || any(settings$n!=as.integer(settings$n)) ||
      any(settings$d!=as.integer(settings$d))))
    stop('Require canonical integer degree metadata n>=d')
  settings$alt <- settings$scenario;settings$curve <- settings$sweep
  if(any(vapply(sweeps,function(s)sum(settings$sweep==s)<2L,FALSE)))
    stop('Each sweep requires at least two sampled points')
  if(!is.numeric(B) || length(B)!=1L || B<1 || B!=as.integer(B) ||
     !is.numeric(alpha) || length(alpha)!=1L || alpha<=0 || alpha>=1)
    stop('Invalid calibration B or alpha')
  settings <- settings[order(match(settings$sweep,sweeps),settings$x),,drop=FALSE]
  rownames(settings) <- NULL
  meta <- settings[,setdiff(names(settings),c('panel','method','R','rejections',
    'power','low','high','available','B','alpha')),drop=FALSE]
  meta$setting_order <- seq_len(nrow(meta))
  style <- manuscript_method_styles(line_width=line_width,point_size=point_size)
  plotted_styles <- manuscript_method_styles(panel,line_width=line_width,point_size=point_size)
  rows <- list();k <- 0L
  for(i in seq_len(nrow(meta))) {
    eligible <- power[power$group==meta$group[i] & power$scenario==meta$scenario[i] &
                     power$panel==panel & power$method %in% plotted_styles$method,,drop=FALSE]
    if(nrow(eligible)!=nrow(plotted_styles) ||
       !setequal(eligible$method,plotted_styles$method))
      stop('Incomplete requested method set at ',meta$group[i],' / ',meta$scenario[i])
    if(length(unique(eligible$R))!=1L || any(!is.finite(eligible$R)) || any(eligible$R<1))
      stop('Evaluation budget must be common within each setting')
    if('B' %in% names(eligible) && any(eligible$B!=B)) stop('B differs from score-summary metadata')
    for(method in style$method) {
      available <- !(panel=='unknown' && method=='IHC')
      if(available) {
        z <- eligible[eligible$method==method,,drop=FALSE]
        # Wilson bounds at 0/1 can miss the boundary by floating-point roundoff.
        # Keep the exported values; tolerate only this negligible plotting error.
        if(!all(is.finite(c(z$power,z$low,z$high))) || z$low < -1e-12 || z$high > 1+1e-12 ||
           z$low>z$power+1e-12 || z$power>z$high+1e-12)
          stop('Invalid power/Wilson summary at ',meta$group[i],' / ',method)
        rejections <- if('rejections' %in% names(z)) z$rejections else NA_integer_
        if(is.finite(rejections) && abs(z$power-rejections/z$R)>1e-12)
          stop('Rejection count and power disagree')
        stats <- data.frame(panel=panel,method=method,R=z$R,rejections=rejections,
          power=z$power,low=z$low,high=z$high,available=TRUE,B=B,alpha=alpha,
          stringsAsFactors=FALSE)
      } else {
        stats <- data.frame(panel=panel,method=method,R=eligible$R[1],
          rejections=NA_integer_,power=NA_real_,low=NA_real_,high=NA_real_,
          available=FALSE,B=B,alpha=alpha,stringsAsFactors=FALSE)
      }
      k <- k+1L;rows[[k]] <- cbind(meta[i,,drop=FALSE],stats)
    }
  }
  plot_data <- do.call(rbind,rows)
  data_path <- caption_path <- style_path <- NULL
  num <- function(x) paste(format(x,trim=TRUE,scientific=FALSE,digits=7),collapse=', ')
  Rvals <- sort(unique(plot_data$R))
  caption <- paste0(if(panel=='known') 'Known-covariance' else 'Estimated-covariance',
    ' power under complete-statistic common null calibration. ',
    if(identical(sweeps,c('rho','d','amplitude')))
      'From left to right, the panels vary correlation (rho), dimension (d), and signal amplitude (a). '
    else paste0('The ',if(length(sweeps)==1L)'panel varies 'else'panels vary ',
      paste(c(rho='correlation (rho)',d='dimension (d)',amplitude='signal amplitude (a)')[sweeps],collapse=', '),'. '),
    'The supplied settings table gives every covariance, mean vector and Wishart degree n. ',
    'B = ',format(B,big.mark=',',trim=TRUE),' independent null draws per bank; R = ',
    paste(format(Rvals,big.mark=',',trim=TRUE),collapse='/'),' independent evaluation draws per setting; alpha = ',num(alpha),'. ',
    if(panel=='unknown') 'The null bank uses the true covariance while statistics retain their estimators (oracle simulation calibration). ' else '',
    if(show_intervals) 'Whiskers are pointwise 95% Wilson intervals conditional on the realized bank and exclude between-bank uncertainty. ' else '',
    'Lines connect only sampled settings. ',paste(plotted_styles$caption_label,collapse='; '),'. ',
    if(panel=='unknown') 'IHC is defined only for known covariance and is not plotted here.' else '')
  if(save_png) {
    stem <- sub('\\.png$','',outpath,ignore.case=TRUE)
    dir.create(dirname(outpath),recursive=TRUE,showWarnings=FALSE)
    data_path <- paste0(stem,'_data.csv');caption_path <- paste0(stem,'_caption.txt')
    style_path <- paste0(stem,'_styles.csv')
    write.csv(plot_data,data_path,row.names=FALSE,na='NA')
    write.csv(style,style_path,row.names=FALSE)
    writeLines(caption,caption_path,useBytes=TRUE)
    grDevices::png(outpath,width=width,height=height,res=res,pointsize=pointsize,bg='white')
    on.exit(grDevices::dev.off(),add=TRUE)
  } else {
    if(grDevices::dev.cur()==1L)grDevices::dev.new()
    old_par <- graphics::par(no.readonly=TRUE)
    on.exit(graphics::par(old_par),add=TRUE)
  }
  {
    graphics::par(mfrow=c(1,length(sweeps)),oma=c(0,0,0,0),family='sans',las=1)
    for(s in sweeps) {
      first_panel <- s==sweeps[1L]
      z <- plot_data[plot_data$sweep==s & plot_data$available,,drop=FALSE]
      xx <- sort(unique(z$x));span <- diff(range(xx))
      vertical_ticks <- min(diff(xx))/span < .16
      graphics::par(mar=c(if(vertical_ticks)4.6 else 3.8,if(first_panel)4.0 else 3.0,.4,.5))
      graphics::plot(range(xx)+c(-1,1)*span*.035,c(-.015,1.015),type='n',axes=FALSE,
        xaxs='i',yaxs='i',xlab='',ylab='')
      graphics::abline(h=seq(0,1,.2),col='#e7e9ec',lty=1,lwd=.6)
      graphics::axis(1,at=xx,labels=format(xx,trim=TRUE,scientific=FALSE),
        las=if(vertical_ticks)2 else 1,cex.axis=.86,tcl=-.25)
      graphics::axis(2,at=seq(0,1,.2),labels=if(first_panel)format(seq(0,1,.2),trim=TRUE) else FALSE,
        cex.axis=.86,tcl=-.25)
      graphics::mtext(switch(s,rho='Correlations',d='Number of tests',amplitude='Signal strengths'),
        side=1,line=if(vertical_ticks)3.4 else 2.55,cex=1.02)
      if(first_panel) graphics::mtext('Simulated power',side=2,line=2.65,las=0,cex=1.0)
      graphics::box(bty='l',col='#343a40',lwd=.7)
      for(j in seq_len(nrow(plotted_styles))) {
        sty <- plotted_styles[j,];r <- z[z$method==sty$method,,drop=FALSE]
        r <- r[order(r$x),,drop=FALSE]
        if(show_intervals) {
          error_color <- grDevices::adjustcolor(sty$color,alpha.f=.48)
          graphics::segments(r$x,r$low,r$x,r$high,col=error_color,lwd=.65)
          graphics::segments(r$x-span*.007,r$low,r$x+span*.007,r$low,col=error_color,lwd=.65)
          graphics::segments(r$x-span*.007,r$high,r$x+span*.007,r$high,col=error_color,lwd=.65)
        }
        graphics::lines(r$x,r$power,col=sty$color,lty=sty$lty,lwd=sty$line_width)
        graphics::points(r$x,r$power,pch=sty$pch,col=sty$color,bg=sty$color,
          cex=sty$point_size,lwd=.85)
      }
    }
  }
  invisible(list(png=outpath,data_csv=data_path,caption_file=caption_path,
    style_csv=style_path,caption=caption,styles=style,plot_data=plot_data))
}

# Place these definitions literally in a standalone TikZ/PGFplots preamble.
# The two compound marks match R pch8 and pch7 without an external image file.
manuscript_pgf_marker_definitions <- function() {
  paste(c(
    '\\pgfdeclareplotmark{manuscript-asterisk8}{%',
    '\\pgfpathmoveto{\\pgfpoint{-\\pgfplotmarksize}{0pt}}%',
    '\\pgfpathlineto{\\pgfpoint{\\pgfplotmarksize}{0pt}}%',
    '\\pgfpathmoveto{\\pgfpoint{0pt}{-\\pgfplotmarksize}}%',
    '\\pgfpathlineto{\\pgfpoint{0pt}{\\pgfplotmarksize}}%',
    '\\pgfpathmoveto{\\pgfpoint{-.70710678\\pgfplotmarksize}{-.70710678\\pgfplotmarksize}}%',
    '\\pgfpathlineto{\\pgfpoint{.70710678\\pgfplotmarksize}{.70710678\\pgfplotmarksize}}%',
    '\\pgfpathmoveto{\\pgfpoint{-.70710678\\pgfplotmarksize}{.70710678\\pgfplotmarksize}}%',
    '\\pgfpathlineto{\\pgfpoint{.70710678\\pgfplotmarksize}{-.70710678\\pgfplotmarksize}}%',
    '\\pgfusepathqstroke}',
    '\\pgfdeclareplotmark{manuscript-squarex}{%',
    '\\pgfpathrectangle{\\pgfpoint{-\\pgfplotmarksize}{-\\pgfplotmarksize}}{\\pgfpoint{2\\pgfplotmarksize}{2\\pgfplotmarksize}}%',
    '\\pgfpathmoveto{\\pgfpoint{-\\pgfplotmarksize}{-\\pgfplotmarksize}}%',
    '\\pgfpathlineto{\\pgfpoint{\\pgfplotmarksize}{\\pgfplotmarksize}}%',
    '\\pgfpathmoveto{\\pgfpoint{-\\pgfplotmarksize}{\\pgfplotmarksize}}%',
    '\\pgfpathlineto{\\pgfpoint{\\pgfplotmarksize}{-\\pgfplotmarksize}}%',
    '\\pgfusepathqstroke}'),collapse='\n')
}

# settings may be list(main=...,original=...), a combined data frame carrying
# figure_set, or one settings data frame. Complete groups are plotted separately.
# Missing method/setting cells are reported and skipped; no partial curve is drawn.
generate_manuscript_figures <- function(power, settings, outdir,
    B=19999L, alpha=.05, R=10000L, show_intervals=TRUE,
    skip_incomplete=TRUE, line_width=NULL, point_size=NULL, ...) {
  if(is.character(settings) && length(settings)==1L)
    settings <- read.csv(settings,stringsAsFactors=FALSE,check.names=FALSE)
  if(is.data.frame(settings)) {
    sets <- if('figure_set' %in% names(settings)) split(settings,settings$figure_set)
      else list(curves=settings)
  } else if(is.list(settings) && length(settings)>0L) {
    sets <- settings
    if(is.null(names(sets)) || any(!nzchar(names(sets))))
      stop('Name settings lists, e.g. list(main=...,original=...)')
  } else stop('Supply a settings data frame/CSV or named list of such data frames')
  if(any(!grepl('^[A-Za-z0-9_-]+$',names(sets)))) stop('Unsafe figure-set name')
  dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
  results <- list();skipped <- list()
  for(name in names(sets)) {
    ss <- sets[[name]]
    if(is.character(ss) && length(ss)==1L) ss <- read.csv(ss,stringsAsFactors=FALSE,check.names=FALSE)
    curve_col <- if('sweep' %in% names(ss)) 'sweep' else if('curve' %in% names(ss)) 'curve' else NULL
    if(is.null(curve_col)) stop('Each figure set requires curve or sweep')
    values <- as.character(ss[[curve_col]])
    aliases <- c(correlation='rho',dimension='d',a='amplitude')
    hit <- values %in% names(aliases);values[hit] <- unname(aliases[values[hit]])
    ss <- ss[values %in% c('rho','d','amplitude'),,drop=FALSE]
    values <- values[values %in% c('rho','d','amplitude')]
    if(!nrow(ss)) { skipped[[name]] <- 'No three-sweep curve rows; geometry/confirmation tables are not curve figures';next }
    count <- table(factor(values,levels=c('rho','d','amplitude')))
    if(!all(count==5L)) {
      reason <- paste('Require all five points in each rho/d/amplitude sweep; found',paste(count,collapse='/'))
      if(!skip_incomplete) stop(reason)
      skipped[[name]] <- reason;next
    }
    for(panel in c('known','unknown')) {
    path <- file.path(outdir,paste0(name,'_',if(panel=='known')'known' else 'estimated','.png'))
    result <- tryCatch(plot_manuscript_curves(power,ss,panel,path,
      B=B,alpha=alpha,R=R,show_intervals=show_intervals,line_width=line_width,
      point_size=point_size,...),error=function(e)e)
    label <- paste(name,panel,sep='/')
    if(inherits(result,'error')) {
      if(!skip_incomplete || !startsWith(conditionMessage(result),'Incomplete requested method set'))
        stop(result)
      skipped[[label]] <- conditionMessage(result)
    } else results[[label]] <- result
    }
  }
  invisible(list(figures=results,skipped=skipped,
    styles=manuscript_method_styles(line_width=line_width,point_size=point_size),
                  pgf_marker_definitions=manuscript_pgf_marker_definitions()))
}

# 11. Individual panels. These functions use completed results; they do not
# simulate observations. See EA_individual_figures.R for 12 select-and-run blocks.
manuscript_figure_index <- function() {
  study <- rep(c('main','original'),each=6L)
  panel <- rep(rep(c('known','unknown'),each=3L),2L)
  sweep <- rep(c('rho','d','amplitude'),4L)
  data.frame(figure=seq_len(12L),study=study,
    covariance=ifelse(panel=='known','known','estimated'),panel=panel,sweep=sweep,
    filename=sprintf('%02d_%s_%s_%s.png',seq_len(12L),study,
      ifelse(panel=='known','known','estimated'),sweep),stringsAsFactors=FALSE)
}

# Supply a run_manuscript()/source_evaluation() result, or a directory containing
# power.csv and settings.csv. Loading files never reruns a simulation.
read_manuscript_plot_results <- function(results) {
  if(is.character(results) && length(results)==1L && !is.na(results)) {
    files <- file.path(results,c('power.csv','settings.csv'))
    if(!all(file.exists(files)))
      stop('Missing power.csv/settings.csv in ',results,
        '. Supply a completed results directory or a simulation result object.')
    results <- list(power=read.csv(files[1L],stringsAsFactors=FALSE),
      settings=read.csv(files[2L],stringsAsFactors=FALSE))
  }
  if(!is.list(results) || !is.data.frame(results$power) ||
     !is.data.frame(results$settings))
    stop('results must contain power/settings data frames, or name their directory')
  results
}

# General single-panel entry point. Omitting outpath displays the panel on the
# active graphics device (normally RStudio's Plots pane). Give a .png path to save.
plot_manuscript_panel <- function(power,settings,panel,sweep,outpath=NULL,
    width=1800L,height=1800L,res=450L,pointsize=10,...) {
  plot_manuscript_curves(power,settings,panel,outpath=outpath,
    width=width,height=height,res=res,pointsize=pointsize,sweeps=sweep,...)
}

# IDs 1:3 = main known rho/d/amplitude; 4:6 = main estimated;
# 7:9 = supplementary known; 10:12 = supplementary estimated.
plot_manuscript_figure <- function(figure,results,outpath=NULL,res=450L,...) {
  index <- manuscript_figure_index()
  if(!is.numeric(figure) || length(figure)!=1L || !is.finite(figure) ||
     figure!=as.integer(figure) || !figure %in% index$figure)
    stop('figure must be one integer from 1 to 12; see manuscript_figure_index()')
  entry <- index[index$figure==figure,,drop=FALSE]
  results <- read_manuscript_plot_results(results)
  settings <- results$settings
  if(!'sweep' %in% names(settings) && 'curve' %in% names(settings))
    settings$sweep <- settings$curve
  if(!'figure_set' %in% names(settings) && 'stage' %in% names(settings))
    settings$figure_set <- ifelse(settings$stage=='curves','main',settings$stage)
  if(!all(c('figure_set','sweep','group') %in% names(settings)))
    stop('settings must include figure_set (main/original), sweep and group')
  aliases <- c(correlation='rho',dimension='d',a='amplitude')
  hit <- settings$sweep %in% names(aliases)
  settings$sweep[hit] <- unname(aliases[settings$sweep[hit]])
  settings <- settings[settings$figure_set==entry$study &
    settings$sweep==entry$sweep,,drop=FALSE]
  if(nrow(settings)!=5L)
    stop('Figure ',figure,' requires its five completed ',entry$study,' ',
      entry$sweep,' settings; this result contains ',nrow(settings),'.')
  power <- results$power
  power <- power[power$group %in% settings$group & power$panel==entry$panel,,drop=FALSE]
  arguments <- list(power=power,settings=settings,panel=entry$panel,
    sweep=entry$sweep,outpath=outpath,res=res)
  dots <- list(...)
  # Infer full/quick budgets from the data, unless explicitly overridden.
  for(key in c('B','R','alpha'))if(!key %in% names(dots) && key %in% names(power)) {
    value <- unique(power[[key]])
    if(length(value)!=1L || !is.finite(value))stop('Inconsistent ',key,' for this figure')
    arguments[[key]] <- value
  }
  result <- do.call(plot_manuscript_panel,c(arguments,dots))
  result$figure <- entry
  invisible(result)
}

# Optional: save every panel separately, or just a chosen subset, e.g. figures=4:6.
generate_individual_figures <- function(results,outdir,figures=1:12,res=450L,...) {
  results <- read_manuscript_plot_results(results)
  index <- manuscript_figure_index()
  if(!is.numeric(figures) || !length(figures) || anyNA(figures) ||
     any(!is.finite(figures)) || any(figures!=as.integer(figures)) ||
     any(!figures %in% index$figure) || anyDuplicated(figures))
    stop('figures must select distinct integers from 1 to 12')
  answer <- lapply(figures,function(i)plot_manuscript_figure(i,results,
    outpath=file.path(outdir,index$filename[index$figure==i]),res=res,...))
  names(answer) <- as.character(figures)
  invisible(answer)
}

# 12. Editable five-point sweeps. See EA_explore_settings.R for select-and-run
# sections. The statistical engine above is shared with the fixed manuscript runs.
# n ALWAYS denotes Wishart degrees of freedom. The estimated panel retains
# oracle common-bank calibration, exactly as in the manuscript.

.custom_sweep_design <- function(sweep,values,study,d,n,rho,amplitude,dimension_n,seed) {
  scalar <- function(x,label,minimum=-Inf,integer=FALSE) {
    if(!is.numeric(x)||is.complex(x)||length(x)!=1L||!is.finite(x)||x<minimum||
       (integer&&(x!=floor(x)||x>.Machine$integer.max)))
      stop(label," must be one finite ",if(integer)"integer " else "number ",
           "at least ",minimum,call.=FALSE)
    if(integer)as.integer(x) else as.numeric(x)
  }
  if(!is.numeric(values)||is.complex(values)||!is.null(dim(values))||length(values)!=5L||
     any(!is.finite(values))||anyDuplicated(values))
    stop("values must contain exactly five finite, distinct numeric values",call.=FALSE)
  values<-as.numeric(values)
  minimum_d<-if(study=="main")4L else if(sweep=="d")10L else 5L
  if(sweep=="d") {
    if(any(values<minimum_d|values!=floor(values)|values>.Machine$integer.max))
      stop("Dimension values must be integers at least ",minimum_d,call.=FALSE)
    dimensions<-as.integer(values)
    degrees<-if(is.null(dimension_n))
      (if(study=="main")5 else 3)*dimensions-1 else dimension_n
    if(!is.numeric(degrees)||is.complex(degrees)||!length(degrees)%in%c(1L,5L))
      stop("dimension_n must be NULL, one integer, or five integers in the dimension order",call.=FALSE)
    if(length(degrees)==1L)degrees<-rep(degrees,5L)
  } else {
    d<-scalar(d,"d",minimum_d,TRUE);n<-scalar(n,"n",max(d,5L),TRUE)
    dimensions<-rep(d,5L);degrees<-rep(n,5L)
  }
  if(any(!is.finite(degrees)|degrees!=floor(degrees)|degrees< pmax(dimensions,5L)|
         degrees>.Machine$integer.max))
    stop("Each Wishart degree must be an integer n >= max(d,5)",call.=FALSE)
  degrees<-as.integer(degrees)
  if(sweep=="amplitude"&&any(values<0))stop("Signal strengths must be nonnegative",call.=FALSE)
  if(sweep=="rho")rhos<-values else {
    default_rho<-if(study=="main") .95 else if(sweep=="d") .7 else .9
    fixed_rho<-scalar(if(is.null(rho))default_rho else rho,"rho")
    rhos<-rep(fixed_rho,5L)
  }
  if(any(rhos<= -1|rhos>=1))stop("Correlations must satisfy -1 < rho < 1",call.=FALSE)
  if(study=="original"&&sweep!="d"&&any(rhos<= -1/4))
    stop("The original five-coordinate block requires -1/4 < rho < 1",call.=FALSE)
  if(study=="original"&&sweep=="d"&&any(rhos<= -1/9))
    stop("The original ten-coordinate block requires -1/9 < rho < 1",call.=FALSE)
  default_a<-if(study=="main")3.5 else if(sweep=="d")3 else 2
  amplitudes<-if(sweep=="amplitude")values else
    rep(scalar(if(is.null(amplitude))default_a else amplitude,"amplitude",0),5L)
  groups<-list();points<-list()
  group_count<-if(sweep=="amplitude")1L else 5L
  seed<-scalar(seed,"seed",0,TRUE)
  if(seed>.Machine$integer.max-10000*group_count-1)
    stop("seed is too large for the independent calibration/evaluation seed offsets",call.=FALSE)
  for(i in seq_len(5L)) {
    dd<-dimensions[i];rr<-rhos[i];a<-amplitudes[i]
    if(study=="main") {
      Sigma<-make_block_cov(c(2L,dd-2L),c(rr,0))
      # Stable at zero and for signed correlations: no subtraction cancellation.
      tstar<-rr/(1+sqrt(1-rr*rr));mu<-c(a,a*tstar,rep(0,dd-2L))
    } else if(sweep=="d") {
      # The first 10 coordinates form the signal block. Background blocks have
      # size 10, with a shorter LAST block when d is not a multiple of 10.
      remainder<-dd-10L;background<-rep(10L,remainder%/%10L)
      if(remainder%%10L)background<-c(background,remainder%%10L)
      Sigma<-make_block_cov(c(10L,background),c(rr,rep(.1,length(background))))
      mu<-c(a,rep(.5,9L),rep(0,dd-10L))
    } else {
      sizes<-c(5L,if(dd>5L)dd-5L else integer())
      Sigma<-make_block_cov(sizes,c(rr,if(dd>5L).1 else numeric()))
      mu<-c(a,rep(0,dd-1L))
    }
    # make_block_cov checks positive definiteness; these are complete mean vectors.
    index<-if(sweep=="amplitude")1L else i
    id<-paste("custom",study,sweep,index,sep="_")
    alt<-list(id=paste0(sweep,"_",i),mu=mu,design=sweep,value=values[i])
    if(is.null(groups[[id]]))groups[[id]]<-list(Sigma=Sigma,n=degrees[i],alts=list(),
      stage="custom",master_seed=seed,index=index,source_dir=NA_character_,source_id=id)
    groups[[id]]$alts<-c(groups[[id]]$alts,list(alt))
    points[[i]]<-data.frame(group=id,scenario=alt$id,value=values[i],d=dd,n=degrees[i],
      rho=rr,amplitude=a,seed_master=seed,seed_cal=seed+10000*index,
      seed_eval=seed+10000*index+1,stringsAsFactors=FALSE)
  }
  list(groups=groups,points=do.call(rbind,points),values=values)
}

# Change values to explore. Correlation/amplitude sweeps use the fixed d and n.
# Dimension sweeps use the values as d, and n=5d-1 (main) or 3d-1 (original).
# dimension_n optionally supplies one common degree or five degrees IN THAT ORDER.
# Main means are a*(1,tstar,0,...), with default rho=.95 and a=3.5.
# Original rho/amplitude means are (a,0,...), default a=2 and active rho=.9.
# Original dimension defaults: first-block rho=.7 (overridable), a=3,
# and means (a,.5,...,.5,0,...) with nine .5 entries.
# Five amplitudes share ONE covariance, null bank and evaluation draw per replicate.
# By default all results stay in memory: no directory, files or caches are written.
# Supply outdir="some_folder" explicitly only if saved results are wanted.
# In saved mode, editing parameters automatically selects a new hashed directory.
run_custom_sweep <- function(sweep=c("rho","d","amplitude"),values,
    study=c("main","original"),d=50,n=249,rho=NULL,amplitude=NULL,
    dimension_n=NULL,budget=c("pilot","full","quick"),B=NULL,R=NULL,
    workers=4L,seed=2026101500,outdir=NULL,backend="base",resume=TRUE) {
  sweep<-match.arg(sweep);study<-match.arg(study);budget<-match.arg(budget)
  backend<-match.arg(backend,c("auto","base","cpp"))
  required<-c("run_setting","make_block_cov","settings_table","write_result_tables",
              "adjust_paired","enable_acceleration",".manuscript_signature")
  if(any(!vapply(required,exists,FALSE,envir=environment(run_custom_sweep),
                 mode="function",inherits=TRUE)))
    stop("Source EA_manuscript_simulations.R before running a custom sweep",call.=FALSE)
  defaults<-switch(budget,quick=c(199L,200L),pilot=c(1999L,2000L),full=c(19999L,10000L))
  if(is.null(B))B<-defaults[1L];if(is.null(R))R<-defaults[2L]
  valid_integer<-function(x,minimum)is.numeric(x)&&!is.complex(x)&&length(x)==1L&&
    is.finite(x)&&x>=minimum&&x==floor(x)&&x<=.Machine$integer.max
  if(!valid_integer(B,1)||!valid_integer(R,2)||!valid_integer(workers,1))
    stop("Require integer B>=1, R>=2 and workers>=1",call.=FALSE)
  workers<-.ea_check_workers(workers)
  persist<-!is.null(outdir)
  if(persist && (!is.character(outdir)||length(outdir)!=1L||is.na(outdir)||!nzchar(outdir)))
    stop("outdir must be NULL (memory only) or one nonempty directory name",call.=FALSE)
  if(!is.logical(resume)||length(resume)!=1L||is.na(resume))stop("resume must be TRUE or FALSE",call.=FALSE)
  B<-as.integer(B);R<-as.integer(R)
  design<-.custom_sweep_design(sweep,values,study,d,n,rho,amplitude,dimension_n,seed)
  settings<-design$groups;enable_acceleration(backend)
  cluster<-NULL
  if(workers>1L) {
    message('Using ',workers,' simulation workers plus one coordinating R process; no nested pools.')
    cluster<-.ea_make_cluster(workers,backend)
    on.exit(parallel::stopCluster(cluster),add=TRUE)
  } else message('Using serial simulation in the coordinating R process.')
  config<-list(version="custom-five-point-v1",study=study,sweep=sweep,
    values=design$values,groups=settings,budget=budget,B=B,R=R,alpha=.05,
    seed=as.integer(seed),backend_requested=backend,backend_effective=.ea_runtime$backend,
    score_signature=if(persist).manuscript_signature() else NA_character_,family_size=75L,
    family_scope="Five requested grid-point scenarios x (8 known + 7 estimated) comparisons; includes a=0 if chosen")
  directory<-run_hash<-NULL
  if(persist) {
    temporary<-tempfile();on.exit(unlink(temporary),add=TRUE)
    saveRDS(config,temporary,compress=FALSE,version=3)
    run_hash<-unname(tools::md5sum(temporary))
    directory<-file.path(outdir,paste(study,sweep,run_hash,sep="-"))
    dir.create(directory,recursive=TRUE,showWarnings=FALSE)
    directory<-normalizePath(directory,mustWork=TRUE)
    config_file<-file.path(directory,"custom_config.rds")
    if(file.exists(config_file)&&!identical(readRDS(config_file),config))
      stop("Existing hashed directory has a different configuration",call.=FALSE)
    saveRDS(config,config_file)
  }
  workers_used<-if(is.null(cluster))1L else length(cluster)
  # One reusable pool scores replicates within each setting. No nested parallelism.
  answer<-lapply(seq_along(settings),function(j)
    run_setting(settings[[j]],names(settings)[j],B=B,R=R,outdir=directory,
                backend=backend,resume=resume,workers=1L,cluster=cluster))
  if(any(vapply(answer,inherits,FALSE,"try-error")))
    stop(if(persist)"A custom group failed; completed groups remain available to resume" else
      "A custom group failed; no results were saved to disk",call.=FALSE)
  tab<-settings_table(settings);tab$figure_set<-study
  power<-do.call(rbind,lapply(answer,`[[`,"power"))
  paired<-adjust_paired(do.call(rbind,lapply(answer,`[[`,"paired")),75L)
  paired$in_simultaneous_family<-paired$scenario!="null"
  result<-list(power=power,paired=paired,
    analytic_EA=do.call(rbind,lapply(answer,`[[`,"analytic_EA")),settings=tab,
    directory=directory,config=config,run_hash=run_hash,workers_used=workers_used,
    parallel=list(workers=if(is.null(cluster))0L else length(cluster),coordinator=1L,
      maximum_R_processes=if(is.null(cluster))1L else length(cluster)+1L,
      native_threads_per_worker=1L))
  if(persist) {
    write_result_tables(result,directory)
    if(exists("generate_tables",mode="function"))generate_tables(result,directory)
  }
  manifest<-design$points;manifest$study<-study;manifest$sweep<-sweep
  manifest$B<-B;manifest$R<-R;manifest$budget<-budget;manifest$alpha<-.05
  manifest$family_size<-75L;manifest$backend<-config$backend_effective
  manifest$score_signature<-config$score_signature
  manifest$run_directory<-if(persist)directory else NA_character_
  if(persist)write.csv(manifest,file.path(directory,"custom_manifest.csv"),row.names=FALSE)
  result$manifest<-manifest
  if(persist)message("Custom ",study," ",sweep," results: ",directory) else
    message("Custom ",study," ",sweep," results kept in R memory; no files saved.")
  invisible(result)
}

if(sys.nframe()==0L)main()
