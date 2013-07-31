
#include"TFIM.h"
#include"Calculate.h"
#include"Regulation.h"
#include"EasytoDebug.h"
#include"Sequence.h"
#include"GRN.h"




void Calculate::RandMatrix(double a[][DIMENS],double b[][DIMENS],const int n)
//
{
    int i1,j1,k1,m1;
    //srand((unsigned)time(0));
    for(i1=0;i1!=n+1;++i1)
        for(j1=0;j1!=n+1;++j1)
    {
        if(i1<n&&j1<n)                    //original network coefficient
        {
            m1=(int)a[i1][j1];
            switch(m1)
            {
                case 1:b[i1][j1]=(double)(rand()%101)/25.0+2.0;break;//positive regulate
                case -1:b[i1][j1]=-(double)(rand()%101)/25.0-2.0;break;//negtive regulate
                case 2:b[i1][j1]=0;break;               //no regulation
                default :m1=rand()%2;                
                if(m1)b[i1][j1]=(double)(rand()%101)/25.0+2.0;
                else b[i1][j1]=-(double)(rand()%101)/25.0-2.0;
            }
        }
        else                              
        {
            m1=(int)(a[i1][j1]*100);
            if(rand()%100<m1)k1=m1>0?rand()%6+5:-rand()%6-5;      
            else if(rand()%5)k1=0;
            else k1=-rand()%6-5;
            b[i1][j1]=k1;
        }
    }
}

double Calculate::FaNexVal(double Matr[][DIMENS],double a[],const int n,const int i,    double p[],double q[],double nn[],double r[])
{
    int j;
    double m[3];m[0]=0;                
    m[0]=0-r[i]*a[i];   
    m[2]=0;m[1]=0;
    for(j=0;j<n;++j)
    {
        if(Matr[j][i]>0)m[1]+=pow(a[j],Matr[j][i]);    
        else if(Matr[j][i]<0)m[2]+=pow(a[j],-Matr[j][i]); 
    }
    if(m[1]==0)
    m[0]+=(p[i]/(nn[i])+q[i]/(nn[i]+m[2]));                  
    else m[0]+=(p[i]*m[1]/(nn[i]+m[1])+q[i]/(nn[i]+m[2]));
    return m[0];
}

double Calculate::NextValue(double Matr[][DIMENS],double a[],const int n,const int i,double step ,   double p[],double q[],double nn[],double r[])
{
    int j;double b[DIMENS];
    for(j=0;j<n;++j)b[j]=a[j]+step*FaNexVal(Matr,a,n,j,p,q,nn,r)/2.0;
    return FaNexVal(Matr,b,n,i,p,q,nn,r);
}

void Calculate::Network_1(double ReguMatrix[][DIMENS],int n)
{
    double a[DIMENS],b[DIMENS],c[DIMENS],d[DIMENS],e[DIMENS],f[DIMENS],sum(0),AbsValue(1),MaxScore(-1),Score(0);
    double FenShu[101];
    int i(0),j=0,cou(0),k(0);double TempMatrix[DIMENS][DIMENS];double step=STEP;int pets=PETS;
    for(i=0;i<101;++i)FenShu[i]=0;
    for(i=0;i<n;++i)a[i]=b[i]=e[i]=INITIALVALUE;           
    while(k<NN)                                
    {
        ++k;
        RandMatrix(ReguMatrix,TempMatrix,n); 
        AbsValue=10;j=0;step=STEP;pets=PETS;
        while(AbsValue>0.000001&&cou<MAXTIME) 
        {
            ++j;
            for(i=0;i<n;++i)
            {
                a[i]+=NextValue(TempMatrix,b,n,i,step,p,q,nn,r)*step;
                if(a[i]<0.000001){cou=MAXTIME+2;break;}
            }
            if(cou>=MAXTIME)break;                
            for(i=0;i<n;++i)b[i]=a[i];
            if(j%(pets/8)==0)
            {
                AbsValue=0;
                for(i=0;i<n+1;++i){AbsValue+=fabs(a[i]-e[i]);}
                for(i=0;i<n+1;++i)e[i]=a[i];
                if(AbsValue<5&&pets==PETS&&j==PETS){step=STEP*4.0;pets=PETS/4;++cou;j=0;continue;}
            }
            if(j%pets==0){++cou;j=0;}
        }
        step=STEP;pets=PETS;
        if(cou>=MAXTIME)break;
        for(i=0;i<n;++i){c[i]=d[i]=a[i];c[n]+=a[i];}
        c[n]/=n;d[n]=c[n];cou=0;AbsValue=10;j=0;
        while(AbsValue>0.000001&&cou<MAXTIME)
        {
            ++j;
            for(i=0;i<n+1;++i)
            {
                c[i]+=NextValue(TempMatrix,d,n+1,i,step,p,q,nn,r)*step;
                if(c[i]<0.000001){cou=MAXTIME+2;break;}
            }
            for(i=0;i<n+1;++i)d[i]=c[i];
            if(j%(pets/8)==0)
            {
                AbsValue=0;
                for(i=0;i<n+1;++i){AbsValue+=fabs(c[i]-f[i]);}
                for(i=0;i<n+1;++i)f[i]=c[i];
                if(AbsValue<5&&pets==PETS&&j==PETS){step=STEP*4.0;pets=PETS/4;++cou;j=0;continue;}
            }
            if(j%pets==0){j=0;++cou;}
        }
        if(cou>=MAXTIME){sum+=1;FenShu[0]+=1;break;}                        
        AbsValue=0;
        sum+=1;
        for(i=0;i<n;++i)AbsValue+=(fabs(c[i]-a[i]))/(c[i]>a[i]?a[i]:c[i]);             
        Score=1-AbsValue*AbsValue/(n*n/9.0+AbsValue*AbsValue); 
        Score*=100;
        Score*=(1-(pow(((double)cou)/MAXTIME,3)));
        if(Score>MaxScore)                                 
        {
            for(i=0;i<n+1;++i)
            for(j=0;j<n+1;++j)
            MaxMa[i][j]=TempMatrix[i][j];
            MaxScore=Score;
        }
        FenShu[(int)Score]+=1;               
    }
    for(i=0;i<101;++i)FenShu[i]/=sum;                 
    ofstream fi1;
    fi1.open("Score");
    for(i=0;i<101;++i)
        fi1<<i+0.5<<' '<<FenShu[i]<<endl;
    fi1.close();
    for(i=0;i<n;++i)a[i]=b[i]=e[i]=INITIALVALUE;                              
    AbsValue=10;
    j=0;cou=0;step=STEP;pets=PETS;
    while(AbsValue>0.000001&&cou<MAXTIME)                     
    {
        ++j;
        for(i=0;i<n;++i){a[i]+=NextValue(MaxMa,b,n,i,step,p,q,nn,r)*step;}
        for(i=0;i<n;++i)b[i]=a[i];
        if(j%(pets/8)==0)
        {
            AbsValue=0;
            for(i=0;i<n+1;++i){AbsValue+=fabs(a[i]-e[i]);}
            for(i=0;i<n+1;++i)e[i]=a[i];
            if(AbsValue<5&&pets==PETS&&j==PETS){step=STEP*4.0;pets=PETS/4;++cou;j=0;continue;}
        }
        if(j%pets==0){++cou;j=0;}
    }
    for(i=0;i<n;++i){c[i]=d[i]=f[i]=a[i];c[n]+=a[i];}
    c[n]/=n;d[n]=f[n]=c[n];
    ofstream igemSfw;
    igemSfw.open("ustcsoftware.txt");
    if(!igemSfw)exit(0);
    j=0;cou=0;AbsValue=10;step=STEP;pets=PETS;
    while(AbsValue>0.00000001&&cou<MAXTIME)             
    {
        ++j;
        for(i=0;i<n+1;++i){c[i]+=NextValue(MaxMa,d,n+1,i,step,p,q,nn,r)*step;}
        for(i=0;i<n+1;++i)d[i]=c[i];
        if(j%(pets/8)==0)
        {
            igemSfw<<cou+j*step<<' '<<flush;
            for(i=0;i<n;++i)igemSfw<<c[i]<<' '<<flush;
            igemSfw<<c[n]<<endl;
            AbsValue=0;
            for(i=0;i<n+1;++i){AbsValue+=fabs(c[i]-f[i]);}
            for(i=0;i<n+1;++i)f[i]=c[i];
            if(AbsValue<5&&pets==PETS&&j==PETS){step=STEP*4.0;pets=PETS/4;++cou;j=0;continue;}
        }
        if(j%pets==0)
        {
            ++cou;j=0;
        }
    }
}



void Calculate::Network_2(double Matr[][DIMENS],int n)
{
    int i,j=0,cou=0;
    double b[DIMENS],AbsValue=10, c[DIMENS];
    for(i=0;i<n;++i){nong[i]=INITIALVALUE;b[i]=nong[i];c[i]=b[i];}
    int pets=PETS;double step=STEP;
    while(AbsValue>0.0000001&&cou<MAXTIME)  
    {
        ++j;
        for(i=0;i<n;++i){nong[i]+=NextValue(Matr,b,n,i,step,p,q,nn,r)*step;if(nong[i]<0.0000001){cou=MAXTIME+2;break;}}
        if(cou>=MAXTIME)break;                
        for(i=0;i<n;++i)b[i]=nong[i];
        if(j%(pets/8)==0)
        {
            AbsValue=0;
            for(i=0;i<n;++i){AbsValue+=fabs(nong[i]-c[i]);}
            for(i=0;i<n;++i)c[i]=nong[i];
            if((AbsValue<5)&&(pets==PETS)&&(j==PETS)){step=STEP*4.0;pets=PETS/4;++cou;j=0;continue;}
        }
        if(j%pets==0){++cou;j=0;}
    }
    if(cou>=MAXTIME)                                          
    for(i=0;i<n+1;++i)nong[i]=-1;
}