#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "inc\\myStatFunc.h"
#define FpRand  0
#define FpRandn  1
#define FpGibbs 2

#define ARRAY_SIZE 1000000
#define Ns_SIZE 5000000

double rand_u[ARRAY_SIZE];
double rand_n[ARRAY_SIZE];
double rand_g[Ns_SIZE-ARRAY_SIZE][2];

typedef struct NORMAL_2D_PARAMETER
{
    double u[2][2];
    double cov[2][2][2];
    double rou[2];
}NORM_2D_PARA;
typedef struct NORMAL_PARAMETER
{
    double u;
    double sgm;
    double Amp;
}NORM_PARA;
typedef struct MIX_NORMAL_PARA
{
    NORM_PARA normPara[2];
}MIX_NORM_PARA;

double fmarginv(NORM_2D_PARA* pa,double t,char marginvar)
{
    double tmp;
    if(marginvar!='x'&&marginvar!='y')
        return 0.0;
    int i = marginvar - 'x';

    tmp = 0.5*normpdf(t,pa->u[0][i],sqrt(pa->cov[0][i][i]));
    tmp += 0.5*normpdf(t,pa->u[1][i],sqrt(pa->cov[1][i][i]));
    return tmp;
}
double mixPdf(MIX_NORM_PARA* mixPa,double t)
{
    double tmp;
    tmp = mixPa->normPara[0].Amp*normpdf(t,mixPa->normPara[0].u,mixPa->normPara[0].sgm);
    tmp +=mixPa->normPara[1].Amp*normpdf(t,mixPa->normPara[1].u,mixPa->normPara[1].sgm);
    return tmp;
}
double randU(FILE* fp)
{
    double tmp;
    int flag = fscanf(fp,"%lf",&tmp);
    if(flag==EOF)
    {
        if(fseek(fp,0L,SEEK_SET))
            printf("Error! randU return to file begin position failure!");
        if(EOF == fscanf(fp,"%lf",&tmp))
        {
            printf("Error! ranN is empty!");
            exit(EXIT_FAILURE);
        }
    }
    return tmp;
}
double randN(FILE* fp)
{
    double tmp = 0.0;
    int flag = fscanf(fp,"%lf",&tmp);
    if(flag==EOF)
    {
        if(fseek(fp,0L,SEEK_SET))
            printf("Error! ranN return to file begin position failure!");
        if(EOF == fscanf(fp,"%lf",&tmp))
        {
            printf("Error! ranN is empty!");
            exit(EXIT_FAILURE);
        }
    }
    return tmp;
}
double GeneRandByHMcontinus(double (*p)(MIX_NORM_PARA*,double),MIX_NORM_PARA* mixPa,long Nth)
{
    double xt=0.0;
    ++Nth;
    double y,tmp,alp;
    for(int i=0;i<Nth;++i)
    {
        y = rand_n[i]+xt;
        tmp = rand_u[i];
        alp = p(mixPa,y)*normpdf(xt,y,1.0)/p(mixPa,xt)/normpdf(y,xt,1.0);
        alp = alp<1.0?alp:1.0;
        if(tmp<alp)
            xt = y;
    }
    return xt;
}
double RandSampOnConditionPdf(NORM_2D_PARA* pa,double t,char condVar)
{
    if(condVar!='x'&&condVar!='y')
        return 0.0;
    int j = 1- (condVar - 'x');
    MIX_NORM_PARA mixPa;
    for(int i=0;i<2;++i)
    {
        mixPa.normPara[i].u = pa->u[i][j]+pa->rou[j]*sqrt(pa->cov[i][j][j]/pa->cov[i][1-j][1-j])
            *(t-pa->u[i][1-j]);
        mixPa.normPara[i].sgm = sqrt(pa->cov[i][j][j]*(1-pow(pa->rou[j],2)));
        mixPa.normPara[i].Amp = normpdf(t,pa->u[i][1-j],sqrt(pa->cov[i][1-j][1-j]))
            /fmarginv(pa,t,'x'+(1-j))/2;
    }
    return GeneRandByHMcontinus(mixPdf,&mixPa,30000);
}
int main(int argc,char** argv)
{
    //srand((unsigned)time(NULL));
    /*printf("Hello world!\n");
    char* fileRand[]="Rand.dat";
    char* fileRandn[]="Randn.dat";
    char* fileGibbs[]="Gibbs.dat"*/
    if(argc<4)
    {
        printf("Useage: GibbsSampling.exe randU.dat randN.dat randGibbs.dat\n");
        exit(EXIT_FAILURE);
    }
    FILE** fp = (FILE**)malloc(sizeof(FILE*)*(argc-1));
    for(int i=1;i<argc;++i)
    {
        if(i<argc-1)
            fp[i-1] = fopen(argv[i],"r");
        else
            fp[i-1] = fopen(argv[i],"w");
        if(!fp[i-1])
            printf("%s can not open!\n",argv[i]);
        else
            printf("%s is open!\n",argv[i]);
    }
	printf("begin to read rand file, loading data...\n");
	for(long i=0;i<ARRAY_SIZE;++i)
	{
		rand_u[i] = randU(fp[FpRand]);
		rand_n[i] = randN(fp[FpRandn]);
	}
    long start,end;
    start = clock();
    printf("Program start running, Please wait...\n");

    NORM_2D_PARA pa={{{0,0},{2,3}},
        {{{1.5,0.2},{0.2,1.4}},{{1.1,0.4},{0.4,1.3}}},{0,0}};
    pa.rou[0] = pa.cov[0][0][1]/(sqrt(pa.cov[0][0][0]*pa.cov[0][1][1]));
    pa.rou[1] = pa.cov[1][0][1]/(sqrt(pa.cov[1][0][0]*pa.cov[1][1][1]));
    long Ns = Ns_SIZE;long Nth = ARRAY_SIZE;
    double Z[2]={0.0,0.0};
    long runtime;double lasttime = 0.0;int barnum = 0;double totaltime = 0.0;
    for(long k=0,q=0;k<Ns;++k)
    {
        Z[0]=RandSampOnConditionPdf(&pa,Z[1],'y');
        Z[1]=RandSampOnConditionPdf(&pa,Z[0],'x');
        if(k>Nth)
        {
			rand_g[q][0] = Z[0];
			rand_g[q++][1] = Z[1];
        }
		if((k+1)%100==0)
		{
			runtime = clock();
			totaltime =(double)(runtime-start)*Ns/(k+1)/CLOCKS_PER_SEC;
			lasttime =(double)(runtime-start)*(Ns-k+1)/(k+1)/CLOCKS_PER_SEC;
			barnum = (int)((totaltime-lasttime)/totaltime*50.0);
			printf("\rProgress: [");
			for(int j=1;j<=50;++j)
			if(j<=barnum)
				printf("#");
			else
				printf(" ");
			if(lasttime>3600)
				printf("] : last time maybe %.3f h",lasttime/3600.0);
			else if(lasttime>60)
				printf("] : last time maybe %.3f min",lasttime/60);
			else
				printf("] : last time maybe %.3f",lasttime);
			fflush(stdout);
		}
    }
	for(long i=0;i<Ns_SIZE-ARRAY_SIZE;++i)
	{
		fprintf(fp[FpGibbs],"%lf\t%lf\n",rand_g[i][0],rand_g[i][1]);
	}
    end = clock();
    printf("\nAll done!\n");
    printf("\nProgram runtime = %.3f s.\n",(double)(end-start)/CLOCKS_PER_SEC);
    for(int i=1;i<argc;++i)
    {
        fclose(fp[i-1]);
    }
    printf("All file closed.\n");
    return 0;
}
