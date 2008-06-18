#include <stdio.h>
#include <math.h>
#include "matrix.h"

void aalendes(alltimes,Nalltimes,Ntimes,designX,nx,px,designG,ng,pg,
antpers,start,stop,status,id,clusters,antclust,deltaweight,desret)
double *designX,*alltimes,*start,*stop,*designG,*desret; 
int *nx,*px,*antpers,*Nalltimes,*Ntimes,*ng,*pg,*status,*id,*clusters,*antclust,*deltaweight;
{
  matrix *X,*Z;
  int i,j,l,c,s,cluster[*antpers],count0;
  int stat,maxtime,pmax,count,pers=0;
  double time,dtime,fabs(),sqrt();
  double times[*Ntimes];

  malloc_mat(*nx,*px,X); if (*pg>0) malloc_mat(*nx,*pg,Z); if (*pg==0) malloc_mat(1,1,Z); 

  if (*px>=*pg) pmax=*px; else pmax=*pg; 
  times[0]=alltimes[0]; l=0; count0=0; 
  maxtime=alltimes[*Nalltimes]; 

  for (s=1;s<*Nalltimes;s++)
    {
      time=alltimes[s]; dtime=time-alltimes[s-1]; 
      mat_zeros(X); if (*pg>0) mat_zeros(Z); stat=0;  
      for (c=0,count=0;((c<*nx) && (count!=*antpers));c++) 
	{
	  if ((start[c]<time) && (stop[c]>=time)) {
	    cluster[id[c]]=clusters[c];
	    for(j=0;j<pmax;j++) {
	      if (j<*px) ME(X,id[c],j)=designX[j*(*nx)+c];
	      if (j<*pg) ME(Z,id[c],j)=designG[j*(*ng)+c]; 
	    }
	    if (time==stop[c] && status[c]==1) {pers=id[c];stat=1;l=l+1;}
	    count=count+1; }
	}


      for (i=0;i<(*px);i++)
	for (c=0;c<(*nx);c++)
	  {
	    l=(s-1)*(*nx)+c; 
	    desret[i*(*nx)*(*Nalltimes-1)+l]=ME(X,c,i);
	    count0=count0+1; 
	  }

      if (*pg>0) 
	for (i=*px;i<(*px)+(*pg);i++)
	  for (c=0;c<(*nx);c++)
	    {
	      l=(s-1)*(*nx)+c; 
	      desret[i*(*nx)*(*Nalltimes-1)+l]=ME(Z,c,i-(*px))*dtime;
	      count0=count0+1; 
	    }


      i=(*pg)+(*px);
      for (c=0;c<(*nx);c++)
	{ 
	  l=(s-1)*(*nx)+c; 
	  desret[i*(*nx)*(*Nalltimes-1)+l]=(c==pers)*stat;
	  count0=count0+1; 
	}

      i=(*pg)+(*px)+1; 
      for (c=0;c<(*nx);c++)
	{ 
	  l=(s-1)*(*nx)+c; 
	  desret[i*(*nx)*(*Nalltimes-1)+l]=1/dtime;
	  count0=count0+1; 
	}


      i=(*pg)+(*px)+2; 
      for (c=0;c<(*nx);c++)
	{ 
	  l=(s-1)*(*nx)+c; 
	  desret[i*(*nx)*(*Nalltimes-1)+l]=cluster[c];
	  count0=count0+1; 
	}


      i=(*pg)+(*px)+3; 
      for (c=0;c<(*nx);c++)
	{ 
	  l=(s-1)*(*nx)+c; 
	  desret[i*(*nx)*(*Nalltimes-1)+l]=time;
	  count0=count0+1; 
	}


    } /* s in Nalltimes */

  free_mat(X); 
  if (*pg>0) free_mat(Z); 
}

void aalendesL(alltimes,Nalltimes,Ntimes,designX,nx,px,designG,ng,pg,
antpers,start,stop,status,id,clusters,antclust,deltaweight,desret)
double *designX,*alltimes,*start,*stop,*designG,*desret; 
int *nx,*px,*antpers,*Nalltimes,*Ntimes,*ng,*pg,*status,*id,*clusters,
*antclust,*deltaweight;
{
  matrix *X,*Z,*A,*AI,*XZ,*AIXZ,*avZ;
  int i,j,l,c,s,count0;
  int stat,pmax,count,pers=0; 
  double time,dtime,fabs(),sqrt();

  pmax=max(*px,*pg); 
  malloc_mat(*nx,*px,X); 
  malloc_mat(*nx,*pg,Z); malloc_mat(*nx,*pg,avZ); 
  malloc_mat(*px,*px,A); malloc_mat(*px,*px,AI); 
  malloc_mat(*px,*pg,XZ); malloc_mat(*px,*pg,AIXZ); 
  pmax=max(*px,*pg); 

  for (s=1;s<*Nalltimes;s++)
    {
      time=alltimes[s]; dtime=time-alltimes[s-1]; 
      mat_zeros(X); mat_zeros(Z); stat=0;  
      /* printf(" %ld \n",s);  */ 

      for (c=0,count=0;((c<*nx) && (count!=*antpers));c++) 
	{
	  if ((start[c]<time) && (stop[c]>=time)) {
	    for(j=0;j<pmax;j++) {
	      if (j<*px) ME(X,id[c],j)=designX[j*(*nx)+c];
	      if (j<*pg) ME(Z,id[c],j)=designG[j*(*ng)+c]; 
	      if (time==stop[c] && status[c]==1) {pers=id[c];
		stat=1;l=l+1;
		/* ipers[l]=pers; * ls[l]=s;*/}
	    }
	    count=count+1;
	  }
	}

      MtA(X,X,A); invert(A,AI); 
      if (ME(AI,0,0)==0.0) printf(" X'X not invertible at time %lf \n",time);

      /* printf("0a %lf  %ld \n",time,s); */
      MtA(X,Z,XZ); MxA(AI,XZ,AIXZ); MxA(X,AIXZ,avZ); 

      mat_subtr(Z,avZ,avZ); 

      for (i=0;i<(*pg);i++)
	for (c=0;c<(*nx);c++)
	  {
	    l=(s-1)*(*nx)+c; 
	    desret[i*(*nx)*(*Nalltimes-1)+l]=ME(avZ,c,i)*sqrt(dtime);
	    count0=count0+1; 
	  }

      for (c=0;c<(*nx);c++)
	{ 
	  l=(s-1)*(*nx)+c; 
	  desret[i*(*nx)*(*Nalltimes-1)+l]=(c==pers)*stat/sqrt(dtime);
	  desret[(i+1)*(*nx)*(*Nalltimes-1)+l]=dtime;
	  desret[(i+2)*(*nx)*(*Nalltimes-1)+l]=c; 
	  desret[(i+3)*(*nx)*(*Nalltimes-1)+l]=time;
	  count0=count0+1; 
	}
      /* printf("%lf  %ld \n",time,s);  */
    } /* s =1...Ntimes */ 

  free_mats(&avZ,&X,&Z,&A,&AI,&XZ,&AIXZ,NULL); 
}
