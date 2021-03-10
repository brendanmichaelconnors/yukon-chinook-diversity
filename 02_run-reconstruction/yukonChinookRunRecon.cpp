//-----------------------------------------------------------------------------//
// yukonChinookRunRecon.cpp                                                    //
// Objective function for Yukon River Chinook run reconstruction               //
//                                                                             //
// Copyright 2019 by Landmark Fisheries Research, Ltd.                         //
//                                                                             //
// This software is provided to Essa Technologies in the hope that it will be  //
// useful, but WITHOUT ANY WARRANTY; without even the implied warranty of      //
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                        //
//                                                                             //
// ALL INTELLECTUAL PROPERTY REMAINS WITH LANDMARK FISHERIES RESEARCH, LTD.    //
// THIS SOFTWARE MAY NOT BE REDISTRIBUTED, SUBLICENCED, COPIED, OR SHARED      //
// OUTSIDE OF ESSA TECHNOLOGIES WITHOUT THE EXPRESS WRITTEN CONSENT OF         //
// LANDMARK FISHERIES RESEARCH, LTD.                                           //
//                                                                             //
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" //
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   //
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  //
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE    //
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         //
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        //
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    //
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     //
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     //
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  //
// POSSIBILITY OF SUCH DAMAGE.                                                 //
//-----------------------------------------------------------------------------//

#include <TMB.hpp>
#include <iostream>

template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}
template <class Type> 
Type square(Type x){
  return pow(x,2);
}
template <class Type> 
Type binvlogit(Type x)
{
  Type lb = -1;
  Type ub = 1;
  return invlogit(x)*(ub-lb) + lb;
}
template<class Type>
Type posfun(Type x, Type eps, Type &pen){
  pen += CppAD::CondExpLt(x, eps, Type(0.01) * pow(x-eps,2), Type(0.));
  return CppAD::CondExpGe(x, eps, x, eps/(Type(2.0)-x/eps));
}
template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace density;
  DATA_ARRAY(n_sdtg);
  DATA_ARRAY(E_dtg);
  DATA_VECTOR(I_t);
  DATA_VECTOR(day_d);
  int nS = n_sdtg.dim(0);
  int nD = n_sdtg.dim(1);
  int nT = n_sdtg.dim(2);
  int nG = n_sdtg.dim(3);
  PARAMETER_ARRAY(lnRunSize_st);
  PARAMETER_VECTOR(lnArrivMu_s);
  PARAMETER_VECTOR(lnArrivSD_s);
  PARAMETER_ARRAY(arrivErr_st);
  PARAMETER_VECTOR(lnErrSD_s);
  PARAMETER_ARRAY(logitCor_ss);
  PARAMETER_ARRAY(lnqE_sg);
  PARAMETER_VECTOR(lnqI_s);
  PARAMETER(lnWeightI);
  PARAMETER_ARRAY(lnDisp_tg);
  array<Type> runSize_st(nS,nT);
  vector<Type> arrivSD_s = exp(lnArrivSD_s);
  vector<Type> errSD_s = exp(lnErrSD_s);
  vector<Type> runSize_t(nT);
  array<Type> mu_st(nS,nT);
  array<Type> N_dt(nD,nT);
  array<Type> N_dst(nD,nS,nT);
  array<Type> rho_dst(nD,nS,nT);
  matrix<Type> cor_ss(nS,nS);
  array<Type> Ihat_dtg(nD,nT,nG);
  vector<Type> mrIhat_t(nT);
  array<Type> Phat_sdtg(nS,nD,nT,nG);
  array<Type> nllI_tg(nT,nG);
  vector<Type> nllP_g(nG);
  Type nllMR = 0;
  Type nlp = 0;
  Type varPen = 0;
  N_dt.fill(0);
  nllI_tg.fill(0);
  nllP_g.fill(0);
  mu_st.col(0) = exp(lnArrivMu_s);
  for( int t=0; t<nT; t++ )
  {
    runSize_st.col(t) = exp( lnRunSize_st.col(t) );
    runSize_t(t) = runSize_st.col(t).sum();
    for( int s=0; s<nS; s++ )
    {
      for( int s2=0; s2<nS; s2++ )
        cor_ss(s,s2) = binvlogit(logitCor_ss(s,s2));
      if( t>0 )
        mu_st(s,t) = mu_st(s,t-1)*exp(arrivErr_st(s,t-1));
      for( int d=0; d<nD; d++ )
        rho_dst(d,s,t) = exp( -square(day_d(d)-mu_st(s,t))/(2*square(arrivSD_s(s))) );
      rho_dst.col(t).col(s) /= rho_dst.col(t).col(s).sum();
      N_dst.col(t).col(s) = runSize_st(s,t)*rho_dst.col(t).col(s);
      N_dt.col(t) += N_dst.col(t).col(s);
      for( int g=0; g<nG; g++ )
        Ihat_dtg.col(g).col(t) += exp(lnqE_sg(s,g))*N_dst.col(t).col(s);
    }
    mrIhat_t(t) = (exp(lnqI_s)*runSize_st.col(t)).sum();
  }
  for( int t=0; t<nT; t++ )
  {
    for( int g=0; g<nG; g++ )
    {
      for( int d=0; d<nD; d++ )
      {
        if( !isNA(E_dtg(d,t,g)) )
        {
          Type mu = posfun( Ihat_dtg(d,t,g), Type(1e-8), varPen );
          Type var = mu+mu*mu*exp(lnDisp_tg(t,g));
          if( t==16 & g==1 & d<60 ) Type var = mu + Type(1e-4);
          nllI_tg(t,g) -= dnbinom_robust(E_dtg(d,t,g),log(mu),log(var-mu),TRUE);
        }
        for( int s=0; s<nS; s++ )
          Phat_sdtg(s,d,t,g) = exp(lnqE_sg(s,g))*N_dst(d,s,t)/Ihat_dtg(d,t,g);
        if( !isNA(n_sdtg(0,d,t,g)) )
        {
          vector<Type> N_s = n_sdtg.col(g).col(t).col(d);
          vector<Type> P_s = Phat_sdtg.col(g).col(t).col(d);
          nllP_g(g) -= dmultinom( N_s, P_s, TRUE );
        }
      }
    }
    if( !isNA(I_t(t)) )
      nllMR -= dnorm( log(I_t(t)),
                      log(mrIhat_t(t)),
                      log(I_t(t))*0.06,
                      TRUE );
  }
  matrix<Type> D(nS,nS);
  D.fill(0);
  D.diagonal() = exp(lnErrSD_s);
  matrix<Type> cov_ss = D*cor_ss*D;
  MVNORM_t<Type> errDens(cov_ss);
  for( int t=0; t<(nT-1); t++ )
    nlp += errDens(arrivErr_st.col(t));
  Type objFun = nllI_tg.sum() + nllP_g.sum() + exp(lnWeightI)*nllMR + nlp + 1e3*varPen;
  ADREPORT(runSize_st);
  ADREPORT(runSize_t);
  ADREPORT(errSD_s);
  ADREPORT(mu_st);
  REPORT(n_sdtg);
  REPORT(E_dtg);
  REPORT(I_t);
  REPORT(day_d);
  REPORT(nS);
  REPORT(nG);
  REPORT(nT);
  REPORT(nD);
  REPORT(lnRunSize_st);
  REPORT(lnArrivMu_s);
  REPORT(lnArrivSD_s);
  REPORT(arrivErr_st);
  REPORT(logitCor_ss);
  REPORT(cor_ss);
  REPORT(lnErrSD_s);
  REPORT(lnqI_s);
  REPORT(lnWeightI);
  REPORT(lnDisp_tg);
  REPORT(runSize_t);
  REPORT(runSize_st);
  REPORT(N_dt);
  REPORT(N_dst);
  REPORT(mu_st);
  REPORT(rho_dst);
  REPORT(lnqE_sg);
  REPORT(Ihat_dtg);
  REPORT(mrIhat_t);
  REPORT(Phat_sdtg);
  REPORT(objFun);
  REPORT(nllI_tg);
  REPORT(nllP_g);
  REPORT(nllMR);
  REPORT(nlp);
  REPORT(cov_ss);
  REPORT(varPen);
  return objFun;
}




