// Optional acceleration for EA_methods.R. The moment and Lasso kernels retain the previous
// manuscript engine's formulas. The readable R implementations are the reference.
#include <Rcpp.h>
#ifdef __aarch64__
#include <arm_neon.h>
#endif
// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::export]]
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
}

// [[Rcpp::export]]
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
}

// Only the first b subdiagonals of the triangular inverse are needed for IHC.
// [[Rcpp::export]]
Rcpp::NumericMatrix banded_inverse_cpp(Rcpp::NumericMatrix L, int b) {
  int d = L.nrow();
  Rcpp::NumericMatrix U(d, d);
  for (int j = 0; j < d; ++j) {
    for (int i = j; i < std::min(d, j + b); ++i) {
      double value = (i == j ? 1.0 : 0.0);
      for (int k = j; k < i; ++k) value -= L(i, k) * U(k, j);
      U(i, j) = value / L(i, i);
    }
  }
  return U;
}
