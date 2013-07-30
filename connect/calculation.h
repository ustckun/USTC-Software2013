#define NN 100                                 
#define PETS 256                         
#define STEP (1.0/PETS)                          
#define MAXTIME 1000                           
#define INITIALVALUE 2                         
#define DIMENS 170                              


using namespace std;


inline void RandMatrix(double a[][DIMENS],double b[][DIMENS],const int n)

{
    int i1,j1,k1,m1;
    srand((unsigned)time(0));
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
            if(rand()%100<(50+abs(m1)/2))k1=m1>0?rand()%6+5:-rand()%6-5;    
            else if(rand()%2)k1=0;
            else k1=m1<0?rand()%6+5:-rand()%6-5;
            b[i1][j1]=k1;
        }
    }
}

double FaNexVal(double Matr[][DIMENS],double a[],const int n,const int i,    double p[],double q[],double nn[],double r[])
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

double NextValue(double Matr[][DIMENS],double a[],const int n,const int i,    double p[],double q[],double nn[],double r[])
{
    int j;double b[DIMENS];
    for(j=0;j<n;++j)b[j]=a[j]+STEP*FaNexVal(Matr,a,n,j,p,q,nn,r)/2.0;
    return FaNexVal(Matr,b,n,i,p,q,nn,r);
}

void Network_1(double ReguMatrix[][DIMENS],double MaxMa[][DIMENS],int n)
{
    double a[DIMENS],b[DIMENS],c[DIMENS],d[DIMENS],e[DIMENS],f[DIMENS],sum(0),AbsValue(1),MaxScore(0),Score(0);
    double p[DIMENS],q[DIMENS],r[DIMENS],nn[DIMENS],FenShu[100];
    int i(0),j=0,cou(0),k(0);double TempMatrix[DIMENS][DIMENS];
    for(i=0;i<DIMENS;++i)                      
    {
        p[i]=6;q[i]=3;r[i]=1.5;nn[i]=3;
    }
    for(i=0;i<100;++i)FenShu[i]=0;
    for(i=0;i<n;++i)a[i]=b[i]=e[i]=INITIALVALUE;               
    while(k<NN)                                 
    {
        ++k;
        RandMatrix(ReguMatrix,TempMatrix,n+1);  
        AbsValue=1;j=0;
        while(AbsValue>0.0000001&&cou<MAXTIME)  
        {
            ++j;
            for(i=0;i<n;++i){a[i]+=NextValue(TempMatrix,b,n,i,p,q,nn,r)*STEP;if(a[i]<0.0000001){cou=MAXTIME+2;break;}}
            if(cou>=MAXTIME)break;                  
            for(i=0;i<n;++i)b[i]=a[i];
            if(j%(PETS/8)==0)
            {
                AbsValue=0;
                for(i=0;i<n+1;++i){AbsValue+=fabs(a[i]-e[i]);}
                for(i=0;i<n+1;++i)e[i]=a[i];
            }
            if(j%PETS==0){++cou;j=0;}
        }
        if(cou>=MAXTIME)break;
        for(i=0;i<n;++i){c[i]=d[i]=a[i];c[n]+=a[i];}
        c[n]/=n;d[n]=c[n];cou=0;AbsValue=1;j=0;
        while(AbsValue>0.0000001&&cou<MAXTIME)
        {
            ++j;
            for(i=0;i<n;++i){c[DIMENS]+=NextValue(TempMatrix,d,n,i,p,q,nn,r)*STEP;if(c[i]<0.0000001){cou=MAXTIME+2;break;}}
            if(cou>=MAXTIME)break;
            for(i=0;i<n;++i)d[i]=c[i];
            if(j%(PETS/8)==0)
            {
                AbsValue=0;
                for(i=0;i<n+1;++i){AbsValue+=fabs(c[i]-f[i]);}
                for(i=0;i<n+1;++i)f[i]=c[i];
            }
            if(j%PETS==0){j=0;++cou;}
        }
        if(cou>=MAXTIME){sum+=1;FenShu[0]+=1;break;}                        
        AbsValue=0;
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
    for(i=0;i<100;++i)FenShu[i]/=sum;                
    ofstream fi1;
    fi1.open("fenshu.txt");
    for(i=0;i<100;++i)
        fi1<<i+0.5<<' '<<FenShu[i]<<endl;
    fi1.close();
    for(i=0;i<n;++i)a[i]=b[i]=e[i]=INITIALVALUE;                             
    AbsValue=10;
    j=0;cou=0;
    while(AbsValue>0.0000001&&cou<MAXTIME)                   
    {
        ++j;
        for(i=0;i<n;++i){a[i]+=NextValue(MaxMa,b,n,i,p,q,nn,r)*STEP;}
        for(i=0;i<n;++i)b[i]=a[i];
        if(j%(PETS/8)==0)
        {
            AbsValue=0;
            for(i=0;i<n+1;++i){AbsValue+=fabs(a[i]-e[i]);}
            for(i=0;i<n+1;++i)e[i]=a[i];
        }
        if(j%PETS==0){++cou;j=0;}
    }
    for(i=0;i<n;++i){c[i]=d[i]=f[i]=a[i];c[n]+=a[i];}
    c[n]/=n;d[n]=f[n]=c[n];
    ofstream igemSfw;
    igemSfw.open("ustcsoftware.txt");
    if(!igemSfw)exit(0);
    j=0;cou=0;AbsValue=1;
    while(AbsValue>0.0000001&&cou<MAXTIME)            
    {
        ++j;
        for(i=0;i<n+1;++i){c[i]+=NextValue(TempMatrix,d,n+1,i,p,q,nn,r)*STEP;}
        for(i=0;i<n+1;++i)d[i]=c[i];
        igemSfw<<cou+j*STEP<<' ';
        for(i=0;i<n;++i)igemSfw<<c[i]<<' '<<flush;
        igemSfw<<c[n]<<endl;
        if(j%100==0)
        {
            AbsValue=0;
            for(i=0;i<n+1;++i){AbsValue+=fabs(c[i]-f[i]);}
            for(i=0;i<n+1;++i)f[i]=c[i];
        }
        if(j%(PETS)==0)
        {
            ++cou;j=0;
        }
    }
}



void Network_2(double Matr[][DIMENS],int n,double a[])
{
    int i,j=0,cou=0;double p[DIMENS], q[DIMENS], nn[DIMENS], r[DIMENS];
    double b[DIMENS],AbsValue=10, c[DIMENS];
    for(i=0;i<n;++i){a[i]=INITIALVALUE;b[i]=a[i];c[i]=b[i];}
    for(i=0;i<DIMENS;++i)                       
    {
        p[i]=6;q[i]=3;r[i]=1.5;nn[i]=3;
    }
    while(AbsValue>0.0000001&&cou<MAXTIME)  
    {
        ++j;
        for(i=0;i<n;++i){a[i]+=NextValue(Matr,b,n,i,p,q,nn,r)*STEP;if(a[i]<0.0000001){cou=MAXTIME+2;break;}}
        if(cou>=MAXTIME)break;                  
        for(i=0;i<n;++i)b[i]=a[i];
        if(j%(PETS/8)==0)
        {
            AbsValue=0;
            for(i=0;i<n;++i){AbsValue+=fabs(a[i]-c[i]);}
            for(i=0;i<n;++i)c[i]=a[i];
        }
        if(j%PETS==0){++cou;j=0;}
    }
    if(cou>=MAXTIME)                                            
    for(i=0;i<n+1;++i)a[i]=-1;
}