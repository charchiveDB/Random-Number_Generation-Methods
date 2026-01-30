/* takes around 30 nanosecs  (~300 million/second) */
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

static unsigned long jxr=521288629;
static  int *K,*P;
static double *V;
static int J[256];
static FILE *fp;
static int size,offset,last;

void PoissP(double lam)  /* Creates Poisson Probabilites */
{int i,j=-1,nlam;        /* P's are 30-bit integers, assumed denominator 2^30 */
double p=1,c,t=1.;
    /* generate  P's from 0 if lam<21.4 */
if(lam<21.4){p=t=exp(-lam); for(i=1;t*2147483648.>1;i++) t*=(lam/i);
          size=i-1; last=i-2;
          /* Given size, malloc and fill P array, (30-bit integers) */
          P=malloc(size*sizeof(int)); P[0]=exp(-lam)*(1<<30)+.5;
          for(i=1;i<size;i++) {p*=(lam/i); P[i]=p*(1<<30)+.5; }
            }
    /* If lam>=21.4, generate from largest P up,then largest down */
if(lam>=21.4)
{nlam=lam;      /*first find size */
c=lam*exp(-lam/nlam);
for(i=1;i<=nlam;i++) t*=(c/i);
p=t;
for(i=nlam+1;t*2147483648.>1;i++) t*=(lam/i);
last=i-2;
t=p; j=-1;
for(i=nlam-1;i>=0;i--){t*=((i+1)/lam); if(t*2147483648.<1){j=i;break;} }
offset=j+1;  size=last-offset+1;
   /* malloc and fill P array as 30-bit integers */
P=malloc(size*sizeof(int));
t=p; P[nlam-offset]=p*(1<<30)+0.5;
for(i=nlam+1;i<=last;i++){t*=(lam/i); P[i-offset]=t*(64<<24)+.5;}
t=p;
for(i=nlam-1;i>=offset;i--){t*=((i+1)/lam);P[i-offset]=t*(1<<30)+0.5;}
} // end lam>=21.4
} //end PoissP

void BinomP(int n, double p)  /*Creates Binomial Probabilities */
         /* Note: if n large and p>.5, generate j=Binom(n,1-p), return n-j*/
{double h,t,p0;
 int i,kb=0,ke,k=0,s;
/* first find size of P array */
p0=t=exp(n*log(1.-p));
h=p/(1.-p);
ke=n;
if(t*2147483648.>1) {k=1;kb=0;}
for(i=1;i<=n;i++)
  { t*=(n+1-i)*h/i;
  if(k==0 && t*2147483648.>1) {k=1;kb=i;}
  if(k==1 && t*2147483648.<1) {k=2;ke=i-1;}
  }
  size=ke-kb+1; offset=kb;
 /* Then malloc and assign P values as 30-bit integers */
P=malloc(size*sizeof(int));
t=p0; for(i=1;i<=kb;i++) t*=(n+1-i)*h/i;
        s=t*1073741824.+.5; P[0]=s;
for(i=kb+1;i<=ke;i++) {t*=(n+1-i)*h/i;
                       P[i-kb]=t*1073741824.+.5;
                       s+=P[i-kb];}
i=n*p-offset; P[i]+=(s-1073741824);
} //end BinomP








void DSQset(double *p,int Hsize)    //sets table-lookup J[0],...J[255]
 {  int k,h,L=0,i=0,j=0;       // then V and K tables for square histogram
    double t,min,max,ss=0,a=1./Hsize;
    V=malloc(Hsize*sizeof(double));
    K=malloc(Hsize*sizeof(int));

    for(i=0;i<Hsize;i++)  // set up 256 table
        {  k=256*(t=p[i]); p[i]=256*t-k; ss+=p[i];
           for(j=0;j<k;j++) J[L+j]=i;  L+=k;
		}
         for(i=L;i<256;i++) J[i]=-1;    //fill empty cells of 256 table with -1's
         for(i=0;i<Hsize;i++) p[i]/=ss;     // make new P's sum to 1
        printf("J table is %5.2f percent full\n",100.*L/256.);

         for (j=0;j<Hsize;++j){K[j]=j; V[j]=(j+1)*a;}  //Initialize V[],K[]
         // Apply Robin Hood rule n-1 times:
            for(L=1;L<=Hsize;L++)
        { min=a;         //Find the minimum and maximum columns
            max=a;
             for (h=0;h<Hsize; h++)
               {  t=p[h];
                  if(t<min) {min=t;i=h;}
                  else if(t>max) {max=t;j=h;}
               }
                 //Take from maximum to bring minimum to average.
                 //Store donating index in K[] and division point in V[]
                 V[i]=min + i*a;
                 K[i]=j;
                 p[j]=max + min - a;
                 p[i]=a;
         }  //end L loop
}

int Dran()  // Table+SquareHisto discrete variate generator.
{int d; float U;
  jxr^=(jxr<<13); jxr^=(jxr>>17); jxr^=(jxr<<5);  //get random (xorshift) 32-bit integer
     d=J[jxr&255];                 //use rightmost byte to access table J[]
     if(d>=0) return d+offset;     //return if tabled value >=0, incrementing with offset.
       U=jxr*2.328306437e-10;
       d=(int)(size*U);              // use square histogram
       if (U<V[d]) return d+offset;
       return K[d]+offset;
}

double Phi(double x)
 {long double s=x,t=0,b=x,x2=x*x,i=1;
    while(s!=t) s=(t=s)+(b*=x2/(i+=2));
    return  .5+s*exp(-.5*x2-.91893853320467274178L);
 }

double chisq(double z,int n)
{double s=0.,t=1.,q,h;
 int i;
 if(z<=0.) return (0.);
 if(n>3000) return Phi((exp(log(z/n)/3.)-1.+2./(9*n))/sqrt(2./(9*n)));
 h=.5*z;
 if(n&1){q=sqrt(z); t=2*exp(-.5*z-.918938533204673)/q;
    for(i=1;i<=(n-2);i+=2){ t=t*z/i; s+=t;}
    return(2*Phi(q)-1-s);
        }
 for(i=1;i<n/2;i++) { t=t*h/i; s+=t;}
 return (1.-(1+s)*exp(-h));
}


void Dtest(n)
/* requires static 'size', static int array P[size] */
/* generates n Dran()'s, tests output */
{ double x=0,y,s=0,*E;
  int kb=0,ke=1000,i,j=0,*M;
 E=malloc(size*sizeof(double));
 M=malloc(size*sizeof(int));
 for(i=0;i<size;i++) {E[i]=(n+0.)*P[i]/1073741824.;M[i]=0;}
 s=0; for(i=0;i<size;i++) {s+=E[i]; if(s>10){kb=i;E[kb]=s;break;} }
 s=0; for(i=size-1;i>0;i--) {s+=E[i]; if(s>10){ke=i;E[ke]=s;break;} }
 j=0; for(i=0;i<=kb;i++) j+=M[i]; M[kb]=j;
 j=0; for(i=ke;i<size;i++) j+=M[i]; M[ke]=j;
 for(i=0;i<size;i++) M[i]=0; s=0; x=0;
 for(i=0;i<n;i++) {j=Dran(); if(j<kb+offset) j=kb+offset;
                           if(j>ke+offset) j=ke+offset;
                        M[j-offset]++;}
printf("\n   D     Observed     Expected    (O-E)^2/E   sum\n");
fprintf(fp,"\n   D     Observed     Expected    (O-E)^2/E   sum\n");
 for(i=kb;i<=ke;i++){ y=M[i]-E[i]; y=y*y/E[i]; s+=y;

 printf("%4d %10d  %12.2f   %7.2f   %7.2f\n",i+offset,M[i],E[i],y,s);
 fprintf(fp,"%4d %10d   %12.2f   %7.2f   %7.2f\n",i+offset,M[i],E[i],y,s);
                  }
 printf("    chisquare for %d d.f.=%7.2f, p=%7.5f\n",ke-kb,s,chisq(s,ke-kb));
 fprintf(fp,"     chisquare for %d d.f.=%7.2f, p=%7.5f\n",ke-kb,s,chisq(s,ke-kb));
}




int main(){
int i,j=0,n,nsmpls=100000000;
double *p;
double lam,bp;
fp=fopen("tests.out","w");
printf("  Enter lambda:\n");
scanf("%lf",&lam);
PoissP(lam);
p=malloc(size*sizeof(double));
for(i=0;i<size;i++) p[i]=P[i]*9.313225746e-10;
DSQset(p,size);
//printf("start"); for(n=0;n<1000000000;n++) j+=Dran(); printf("END\n"); //15 nanos
Dtest(nsmpls);
fprintf(fp," Above results for sample of %d from Poisson, lambda=%3.2f\n\n",nsmpls,lam);
free(P); free(p);
printf("\n Enter n and p for Binomial:\n");
scanf("%d %lf",&n,&bp);
BinomP(n,bp);
p=malloc(size*sizeof(double));
for(i=0;i<size;i++) p[i]=P[i]*9.313225746e-10;
DSQset(p,size);
Dtest(nsmpls);
fprintf(fp," Above result for sample of %d from Binomial(%d,%3.3f)\n\n",nsmpls,n,bp);
free(P); free(p);
printf(" Test results sent to file tests.out\n");
return 0;
}
