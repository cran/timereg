#include <stdio.h>
#include <math.h>
#include "matrix.h"
                 
/* ====================================================== */
void twostage(times,Ntimes,
designX,nx,px, designG,ng,pg,
antpers,start,stop, betaS,Nit,cu,
vcu,Iinv,Vbeta, detail,Rvcu,RVbeta,
id,status,ratesim, score,robust, clusters,
antclust, betafixed,theta,vartheta, thetascore,inverse,
clustsize)
double *designX,*designG,*times,*betaS,*start,*stop,*cu,*Vbeta,*RVbeta,*vcu,*Rvcu,*Iinv,*score,*theta,*vartheta,*thetascore;
int *nx,*px,*ng,*pg,*antpers,*Ntimes,*Nit,*detail,*id,*status,
*robust,*clusters,*antclust,*betafixed,*inverse,*clustsize,*ratesim;
{
  matrix *ldesignX,*ldesG0,*ldesignG,*cdesX,*cdesX2,*cdesX3,*CtVUCt,*A,*AI;
  matrix *Vcov,*dYI,*Ct,*dM1M2,*M1M2t,*COV,*ZX,*ZP,*ZPX; 
  matrix *tmp1,*tmp2,*tmp3,*dS,*S1,*SI,*S2,*M1,*VU,*ZXAI,*VUI; 
  matrix *RobVbeta;
  matrix *St[*Ntimes],*M1M2[*Ntimes],*C[*Ntimes],*ZXAIs[*Ntimes],*dYIt[*Ntimes]; 
  matrix *W3t[*antclust],*W4t[*antclust],*W2t[*antclust],*AIxit[*antpers]; 
  vector *dA,*VdA,*dN,*MdA,*delta,*zav,*lamt,*lamtt;
  vector *xi,*zi,*U,*beta,*xtilde; 
  vector *Gbeta,*zcol,*one,*difzzav; 
  vector *offset,*weight,*ZXdA[*Ntimes],*Uprofile;
  vector *ahatt,*phit[*Ntimes],*dbetaNH,*dANH; 
  vector *tmpv1,*tmpv2,*rowX,*rowZ,*difX,*VdB,*atrisk[*antpers]; 
  vector *W2[*antclust],*W3[*antclust],*reszpbeta,*res1dim,*dAt[*Ntimes]; 
  int c,pers=0,i,j,k,s,it,count,sing,pmax,
      *cluster=calloc(*antpers,sizeof(int));
  double Nt[*antclust],dtime,time,dummy,ll,lle,llo;
  double tau,hati,scale,sumscore,d2Utheta=0;
  double *HeHi=calloc(*antpers,sizeof(double)),
         *Nti=calloc(*antpers,sizeof(double)),
	 *H2eHi=calloc(*antpers,sizeof(double)),
	 *Rthetai=calloc(*antpers,sizeof(double)),
	 dtheta,*Hik=calloc(*antpers,sizeof(double)),
	 *NH=calloc(*antclust,sizeof(double)),
	 *HeH=calloc(*antclust,sizeof(double)),
	 *H2eH=calloc(*antclust,sizeof(double)),
	 *Rtheta=calloc(*antclust,sizeof(double)),
	 *thetaiid=calloc(*antclust,sizeof(double)),
	 Dthetanu=0,DDthetanu=0,nu,
	 *dAiid=calloc(*antclust,sizeof(double)); 
  int *ipers=calloc(*Ntimes,sizeof(int));




  for (j=0;j<*antclust;j++) { Nt[j]=0; NH[j]=0; dAiid[j]=0; }
  for (i=0;i<*antpers;i++) { HeHi[j]=0; H2eHi[j]=0; Rthetai[j]=0; 
    Nti[i]=0; malloc_vec(*Ntimes,atrisk[i]); }

  if (*robust==1) {
    for (j=0;j<*antclust;j++) { malloc_mat(*Ntimes,*px,W3t[j]); 
      malloc_mat(*Ntimes,*px,W4t[j]); malloc_mat(*Ntimes,*pg,W2t[j]); 
      malloc_vec(*pg,W2[j]); malloc_vec(*px,W3[j]); }
    for (j=0;j<*antpers;j++) {malloc_mat(*Ntimes,*px,AIxit[j]);}
  }
  malloc_vec(1,reszpbeta); malloc_vec(1,res1dim); 
  malloc_vec(*pg,dbetaNH); malloc_vec(*px,dANH); 

  malloc_mats(*antpers,*px,&ldesignX,&cdesX,&cdesX2,&cdesX3,NULL); 
  malloc_mats(*antpers,*pg,&ZP,&ldesignG,&ldesG0,NULL); 
  malloc_mats(*px,*px,&Vcov,&COV,&A,&AI,&M1,&CtVUCt,NULL); 
  malloc_mats(*pg,*pg,&RobVbeta,&tmp1,&tmp2,&dS,&S1,&S2,&SI,&VU,&VUI,NULL); 
  malloc_mats(*pg,*px,&ZXAI,&ZX,&dM1M2,&M1M2t,NULL); 
  malloc_mats(*px,*pg,&tmp3,&ZPX,&dYI,&Ct,NULL); 

  malloc_vecs(*antpers,&weight,&lamtt,&lamt,&dN,&zcol,&Gbeta,&one,&offset,NULL); 
  malloc_vecs(*px,&ahatt,&tmpv1,&difX,&VdB,&rowX,&xi,&dA,&VdA,&MdA,NULL); 
  malloc_vecs(*px,&xtilde,NULL); 
  malloc_vecs(*pg,&tmpv2,&rowZ,&zi,&U,&beta,&delta,&zav,&difzzav,&Uprofile,NULL); 

  for(j=0;j<*Ntimes;j++) { malloc_mat(*px,*pg,C[j]); malloc_mat(*pg,*px,M1M2[j]); 
    malloc_mat(*pg,*px,ZXAIs[j]); malloc_mat(*px,*pg,dYIt[j]); malloc_vec(*px,dAt[j]);
    malloc_vec(*pg,ZXdA[j]); malloc_mat(*pg,*pg,St[j]); malloc_vec(*px,phit[j]);  } 

  if (*px>=*pg) pmax=*px; else pmax=*pg; ll=0; 
  for(j=0;j<*pg;j++) VE(beta,j)=betaS[j]; 
  for(j=0;j<*antpers;j++) {Hik[j]=0; VE(one,j)=1; VE(weight,j)=1; 
	  VE(offset,j)=1;} 

  for (it=0;it<*Nit;it++)
    {
      vec_zeros(U); mat_zeros(S1);  sumscore=0;   
      for (s=1;s<*Ntimes;s++)
	{
	  time=times[s]; vec_zeros(dN); sing=0; mat_zeros(ldesignX); mat_zeros(ldesignG); 

	  for (c=0,count=0;((c<*nx) && (count!=*antpers));c++) 
	    {
	      if ((start[c]<time) && (stop[c]>=time))  {
		cluster[id[c]]=clusters[c];
		for(j=0;j<pmax;j++) {
		  if (j<*px) ME(ldesignX,id[c],j)=designX[j*(*nx)+c];
		  if (j<*pg) ME(ldesignG,id[c],j)=designG[j*(*ng)+c]; } 
		if (time==stop[c] && status[c]==1) {VE(dN,id[c])=1; pers=id[c];} 
		if (it==*Nit-1) VE(atrisk[id[c]],s)=1; 
		count=count+1; }
	    }
	  ipers[s]=pers;

	  Mv(ldesignG,beta,Gbeta); 

	  for (j=0;j<*antpers;j++)
	  {extract_row(ldesignX,j,xi); 
	   dummy=exp(VE(Gbeta,j));   /* W D(exp(z*beta))X */ 
	   scl_vec_mult(VE(weight,j)*dummy,xi,xtilde); 
	   replace_row(cdesX,j,xtilde); 
	  }
	  if (it==(*Nit-1) && s==1) { 
		  ldesG0=mat_copy(ldesignG,ldesG0); 
	          cdesX2=mat_copy(cdesX,cdesX2);}

	  scale=VE(weight,pers); 
  
	  MtA(cdesX,ldesignX,A); invert(A,AI); 
	  if (ME(AI,0,0)==0) {
		  printf("X'X not invertible at time %lf \n",time);
	  }

	  extract_row(ldesignX,pers,xi); scl_vec_mult(scale,xi,xi); 
	  Mv(AI,xi,dA); MtA(ldesignG,cdesX,ZX); 
	  MxA(ZX,AI,ZXAIs[s]); 
	  Mv(ZX, dA, ZXdA[s]);  scl_vec_mult(1,dA,dAt[s]); 

	  if (s<1) {print_mat(A); print_mat(AI); print_mat(ZX); print_mat(ZXAIs[s]); }

	  /* First derivative U and Second derivative S  */ 
  
	  extract_row(ldesignG,pers,zi); scl_vec_mult(scale,zi,zi); 
	  Mv(ZX, dA, zav);  vec_subtr(zi,zav,difzzav); vec_add(difzzav,U,U); 

	  if (s<1) {print_vec(zi); print_vec(zav); print_vec(difzzav);}

	  Mv(cdesX,dA,lamt);  
	  for (j=0;j<*antpers;j++)
	    {extract_row(ldesignG,j,zi); 
	      scl_vec_mult(VE(lamt,j),zi,zi); replace_row(ZP,j,zi);} 

	  MtA(ldesignX,ZP,ZPX); MxA(ZXAIs[s],ZPX,tmp2); 

	  MtA(ZP,ldesignG, tmp1); mat_subtr( tmp1,tmp2, dS); 
	  mat_add(dS,S1,S1);  St[s]=mat_copy(S1,St[s]);

	  /* varians beregninger */ 
	  if (it==((*Nit)-1)) { 

	    for (i=0;i<*px;i++) for (j=0;j<*pg;j++) ME(dM1M2,j,i)=VE(dA,i)*VE(difzzav,j);
	    if (*betafixed==0) 
	      for (i=0;i<*pg;i++) { 
		for (j=0;j<*pg;j++) ME(VU,i,j)=ME(VU,i,j)+VE(difzzav,i)*VE(difzzav,j); }

	    MxA(AI,ZPX,dYIt[s]); mat_subtr(Ct,dYIt[s],Ct); C[s]=mat_copy(Ct,C[s]); 

	    vec_star(dA,dA,VdA); mat_add(dM1M2,M1M2t,M1M2t); M1M2[s]=mat_copy(M1M2t,M1M2[s]); 

	    for (k=1;k<=*px;k++) {cu[k*(*Ntimes)+s]=cu[k*(*Ntimes)+s-1]+VE(dA,k-1); 
	      vcu[k*(*Ntimes)+s]=VE(VdA,k-1)+vcu[k*(*Ntimes)+s-1];}
	    if (*robust==1) 
	      for (j=0;j<*antpers;j++)
		{extract_row(ldesignX,j,xi); Mv(AI,xi,rowX); replace_row(AIxit[j],s,rowX);}


	    for (i=0;i<*antpers;i++) {
	      extract_row(cdesX,i,rowX); scl_vec_mult(1/VE(weight,i),rowX,rowX); 
	      vec_star(rowX,dAt[s],tmpv1); hati=vec_sum(tmpv1); 
	      /*  printf(" %ld %ld %lf %lf \n",i,s,Hik[i],hati); */
	      Hik[i]=Hik[i]+hati; 
	      extract_row(ldesignG,i,rowZ); 
	    }
	    Nti[pers]=Nti[pers]+1; 

	    for (j=0;j<*antclust;j++) {
	      for (i=0;i<*antpers;i++)  if (cluster[i]==j) Nt[j]=Nt[j]+(i==pers); }

	  } /* it==*Nit-1 */ 
	} /* Ntimes */ 


      /* for (k=0;k<*pg;k++) ME(S1,k,k)=ME(S1,k,k)+*ridge;  */
      invert(S1,SI); 

      Mv(SI,U,delta); if (*betafixed==0) vec_add(beta,delta,beta); 
      MxA(SI,VU,S2); MxA(S2,SI,VU); 

      for (k=0;k<*pg;k++) sumscore= sumscore+VE(U,k); 

      if ((fabs(sumscore)<0.0000000001) & (it<*Nit-2)) it=*Nit-2; 
    } /* it */
  for (k=0;k<*pg;k++) score[k]=VE(U,k); 

  for (j=0;j<*antclust;j++) 
    {
      for (i=0;i<*antpers;i++) if (cluster[i]==j)  {
	NH[j]=NH[j]+Nti[i]*Hik[i]; 
	extract_row(ldesG0,i,zi); extract_row(cdesX2,i,xi);
	vec_add_mult(dbetaNH,zi,Nti[i]*Hik[i],dbetaNH);  
	vec_add_mult(dANH,xi,Nti[i],dANH);  }
      if (*betafixed==1) vec_zeros(dbetaNH); 
    }


  lle=0; llo=0;
  /* terms for robust variances ============================ */
  for (s=1;s<*Ntimes;s++) 
    {
      time=times[s]; vec_zeros(dN);dtime=time-times[s-1]; 
      cu[s]=times[s]; vcu[s]=times[s]; Rvcu[s]=times[s]; 

      mat_zeros(ldesignX); mat_zeros(ldesignG); 
      for (c=0,count=0;((c<*nx) && (count!=*antpers));c++) 
	{
	  if ((start[c]<time) && (stop[c]>=time))  {
	    cluster[id[c]]=clusters[c];
	    for(j=0;j<pmax;j++) {
	      if (j<*px) ME(ldesignX,id[c],j)=designX[j*(*nx)+c];
	      if (j<*pg) ME(ldesignG,id[c],j)=designG[j*(*ng)+c]; } 
	    if (time==stop[c] && status[c]==1) {pers=id[c];} 
	    count=count+1; }
	}
      Mv(ldesignG,beta,Gbeta); 

      for (j=0;j<*antpers;j++)
	{extract_row(ldesignX,j,xi); extract_row(ldesignG,j,zi);
	  dummy=exp(VE(Gbeta,j))*VE(weight,j)*VE(offset,j); 
	  scl_vec_mult(dummy,xi,xtilde); replace_row(cdesX,j,xtilde); }

      /* terms for robust variance   */ 
      if (*robust==1) {
	for (j=0;j<*antclust;j++) 
	  {
	    for (i=0;i<*antpers;i++) if (cluster[i]==j)  {

	      extract_row(cdesX,i,rowX); scl_vec_mult(1/VE(weight,i),rowX,rowX); 
	      extract_row(ldesignG,i,zi); extract_row(ldesignX,i,xi); 
	      vec_star(rowX,dAt[s],tmpv1); hati=vec_sum(tmpv1); 

	      Mv(ZXAIs[s],xi,tmpv2);  vec_subtr(zi,tmpv2,tmpv2); 
	      scl_vec_mult(VE(weight,i),tmpv2,tmpv2); 

	      if (*betafixed==0) {
		if (i==pers) vec_add(tmpv2,W2[j],W2[j]);
		if (*ratesim==1) {scl_vec_mult(hati,tmpv2,rowZ); vec_subtr(W2[j],rowZ,W2[j]); }
	      }

	      extract_row(AIxit[i],s,rowX);
	      scl_vec_mult(VE(weight,i),rowX,rowX); 

	      if (i==pers) {vec_add(rowX,W3[j],W3[j]); if (hati>0) lle=lle+log(hati);}
	      llo=llo+hati;

	      if (*ratesim==1) {scl_vec_mult(hati,rowX,rowX); vec_subtr(W3[j],rowX,W3[j]);}

	    } /* i if cluste==j */
	    replace_row(W2t[j],s,W2[j]); 
	    replace_row(W3t[j],s,W3[j]);  

	  } /* j and i=1..antclust */ 
      }

      /* MG baseret varians beregning */
      MxA(C[s],VU,tmp3); MAt(tmp3,C[s],CtVUCt);
      MxA(C[s],SI,tmp3); MxA(tmp3,M1M2[s],COV); 

      for (k=1;k<=*px;k++) {
	vcu[k*(*Ntimes)+s]=vcu[k*(*Ntimes)+s]+ME(CtVUCt,k-1,k-1)
	  +2*ME(COV,k-1,k-1); }

      /* */
    } /* s=1 ..Ntimes */ 

  ll=lle-llo; /* likelihood beregnes */


  /* ROBUST VARIANCES   */
  if (*robust==1)
    {
      for (s=1;s<*Ntimes;s++) {
	vec_zeros(VdB); mat_zeros(Vcov);vec_zeros(zi); 

	for (j=0;j<*antclust;j++) {

	  if (*betafixed==0)  {
	    if (s==1 ) {Mv(SI,W2[j],tmpv2); scl_vec_mult(1,tmpv2,W2[j]);  }
	    Mv(C[s],tmpv2,rowX);
	  }
	  if (*betafixed==1) vec_zeros(rowX); 
	  extract_row(W3t[j],s,tmpv1); 
	  vec_add(tmpv1,rowX,difX); 
	  replace_row(W4t[j],s,difX);
	  vec_star(difX,difX,tmpv1); vec_add(tmpv1,VdB,VdB);

	  if (s==1) { for (c=0;c<*pg;c++) for (k=0;k<*pg;k++)
			ME(RobVbeta,c,k)=ME(RobVbeta,c,k)+VE(W2[j],c)*VE(W2[j],k);}
	} /* j in clusters  */

	for (k=1;k<*px+1;k++) { Rvcu[k*(*Ntimes)+s]=VE(VdB,k-1); }
      }  /*  s=1 ..Ntimes */ 

    } /* if robust==1 */ 

  for(j=0;j<*pg;j++) { betaS[j]= VE(beta,j); 
    for (k=0;k<*pg;k++){ Iinv[k*(*pg)+j]=ME(SI,j,k);
      Vbeta[k*(*pg)+j]=-ME(VU,j,k); 
      RVbeta[k*(*pg)+j]=-ME(RobVbeta,j,k); } } 


  /*===================Estimates theta, two stage approach of glidden ==== */

  if (*inverse==1) nu=1/theta[0]; else nu=theta[0]; 

  for (it=0;it<*Nit;it++)
    {
      for (j=0;j<*antclust;j++) {Rtheta[j]=1; HeH[j]=0;H2eH[j]=0;}
      thetascore[0]=0; d2Utheta=0; 

      if (*inverse==1) {
	theta[0]=1/nu; Dthetanu=-1/pow(nu,2);  DDthetanu=2/pow(nu,3); 
      } else {DDthetanu=0; Dthetanu=1; theta[0]=nu;}

      for (j=0;j<*antclust;j++) 
	if (clustsize[j]>=2) {
	  for (i=0;i<*antpers;i++) if (cluster[i]==j)  {
	    Rtheta[j]=Rtheta[j]+exp(theta[0]*Hik[i])-1; 
	    HeH[j]=HeH[j]+Hik[i]*exp(theta[0]*Hik[i]); 
	    H2eH[j]=H2eH[j]+pow(Hik[i],2)*exp(theta[0]*Hik[i]); } }
		

      for (j=0;j<*antclust;j++)  
	if (clustsize[j]>=2) {

	  sumscore=0;  ll=0; 
	  if (Nt[j]>=2) 
	    for (k=2;k<=Nt[j];k++) {
	      tau=(Nt[j]-1)/(1+theta[0]*(Nt[j]-1));
	      lle=-pow((Nt[j]-1),2)/pow((1+theta[0]*(Nt[j]-1)),2);
	      sumscore=sumscore+tau; ll=ll+lle; }

	  thetaiid[j]=sumscore+log(Rtheta[j])/(theta[0]*theta[0])
	    -(1/theta[0]+Nt[j])*HeH[j]/Rtheta[j]+NH[j]; 

	  thetascore[0]=thetascore[0]+Dthetanu*thetaiid[j]; 

	  d2Utheta=d2Utheta+(ll+(2/pow(theta[0],2))*HeH[j]/Rtheta[j]
			     -(2/pow(theta[0],3))*log(Rtheta[j])-(1/theta[0]+Nt[j])*
			     (H2eH[j]*Rtheta[j]-HeH[j]*HeH[j])/pow(Rtheta[j],2));
	}
      if (*inverse==1) {thetascore[0]=thetascore[0]; 
	d2Utheta=
	  DDthetanu/Dthetanu*thetascore[0]+pow(Dthetanu,2)*d2Utheta;}

      dtheta=nu-thetascore[0]/d2Utheta; 
      nu=dtheta; 

      if (*detail==1) { 
	printf("====================Iteration %d ==================== \n",it);
	printf("Estimate theta \n"); printf(" %lf \n",theta[0]); 
	printf("Score D l\n"); printf(" %lf \n",thetascore[0]); 
	printf("Information -D^2 l\n"); 
	printf(" %lf \n",-1/d2Utheta); 
      }

      if ((fabs(thetascore[0])<0.0000000001) & (it<*Nit-2)) it=*Nit-2; 
    } /* it theta Newton-Raphson */ 


  /* beregning af varians bidrag via iid decomp */
  if (*robust==1) {
    vec_zeros(rowZ);   
    for (k=0;k<*antclust;k++) if (clustsize[k]>=2) {

      for (i=0;i<*antpers;i++) if (cluster[i]==k)  {
	dummy=(1/(theta[0]*Rtheta[k]))*exp(theta[0]*Hik[i])-
	  (1/theta[0]+Nt[k])*(1+theta[0]*Hik[i])*exp(theta[0]*Hik[i])
	  /Rtheta[k]+Nti[i]+
	  (1+theta[0]*Nt[k])*exp(theta[0]*Hik[i])*HeH[k]/pow(Rtheta[k],2);
	extract_row(ldesG0,i,zi); 
	vec_add_mult(rowZ,zi,Hik[i]*dummy,rowZ);}


      for (s=1;s<*Ntimes;s++) {
	if (k==0)  {
	  vec_zeros(rowX); 
	  for (j=0;j<*antclust;j++) if (clustsize[j]>=2) {
	    for (i=0;i<*antpers;i++) if (cluster[i]==j && VE(atrisk[i],s)!=0)  {
	      extract_row(cdesX2,i,xi); scl_vec_mult(VE(atrisk[i],s),xi,xi); 
	      dummy=(1/(theta[0]*Rtheta[j]))*exp(theta[0]*Hik[i])-
		(1/theta[0]+Nt[j])*(1+theta[0]*Hik[i])*exp(theta[0]*Hik[i])
		/Rtheta[j]+Nti[i]+
		(1+theta[0]*Nt[j])*exp(theta[0]*Hik[i])*HeH[j]/pow(Rtheta[j],2);
	      vec_add_mult(rowX,xi,dummy,rowX); 
	    }
	  }
	}
	extract_row(W4t[k],s,tmpv1); extract_row(W4t[k],s-1,xi); 
	vec_subtr(tmpv1,xi,xi); vec_star(rowX,xi,tmpv1); 
	if (k==-1) printf("  s er %d \n",s); 
	if (k==-1) print_vec(rowX); 
	dAiid[k]=dAiid[k]+vec_sum(tmpv1); 
      } /* s=1..Ntimes */ 
    } /* k=1..antclust */ 

    /* printf(" Gtilde \n"); print_vec(rowZ);  */

    for (j=0;j<*antclust;j++) if (clustsize[j]>=2) {
      if (*betafixed==1) vec_zeros(W2[j]); 

      /* beta component to iid variation begin added */ 
      vec_star(rowZ,W2[j],zi); 
      /*
	print_vec(W2[j]); 
      */
      thetaiid[j]=thetaiid[j]+Dthetanu*vec_sum(zi); 

      /* A(tau) component to iid variation begin added */ 
      thetaiid[j]=thetaiid[j]+Dthetanu*dAiid[j]; 

      vartheta[0]=vartheta[0]+pow(thetaiid[j],2);
    }
    vartheta[0]=vartheta[0]/pow(d2Utheta,2); 
    /*
      printf("theta and variance is %lf %lf \n",theta[0],vartheta[0]); 
    */
  }

  theta[0]=nu; 

  if (*robust==1) {
    for (j=0;j<*antclust;j++) {
      free_mat(W3t[j]); free_mat(W4t[j]); free_mat(W2t[j]);free_vec(W2[j]); free_vec(W3[j]);
    }
    for (j=0;j<*antpers;j++) free_mat(AIxit[j]); 
  }
  for (j=0;j<*antpers;j++) free_vec(atrisk[j]); 
  for (j=0;j<*Ntimes;j++) {
    free_mat(dYIt[j]); free_vec(dAt[j]); free_mat(C[j]);free_mat(M1M2[j]);free_mat(ZXAIs[j]);
    free_vec(ZXdA[j]);free_mat(St[j]); 
    free_vec(phit[j]);  } 
  free_vec(dbetaNH); free_vec(dANH); 

  free_mats( &ldesignX,&cdesX,&cdesX2,&cdesX3,&ZP,&ldesignG,&ldesG0, 
	       &Vcov,&COV,&A,&AI,&M1,&CtVUCt,&RobVbeta,&tmp1,&tmp2,&dS,&S1,&S2,&SI,&VU,&VUI, 
	       &ZXAI,&ZX,&dM1M2,&M1M2t,&tmp3,&ZPX,&dYI,&Ct,NULL); 
  free_vecs(&weight,&lamtt,&lamt,&dN,&zcol,&Gbeta,&one,&offset, 
	      &ahatt,&tmpv1,&difX,&VdB,&rowX,&xi,&dA,&VdA,&MdA,&xtilde, 
	      &tmpv2,&rowZ,&zi,&U,&beta,&delta,&zav,&difzzav,&Uprofile,
	      &reszpbeta,&res1dim,NULL); 

  free(ipers); free(cluster);  free(Nti); 
  free(HeHi); free(H2eHi); free(Rthetai); free(Hik); free(NH); free(HeH); 
  free(H2eH); free(Rtheta); free(thetaiid); free(dAiid); 
}
