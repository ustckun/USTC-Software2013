#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <ctime>
#define NN 1                                   //统计次数
#define PETS 128                                //步长的倒数
#define STEP (1.0/PETS)                           //步长
#define MAXTIME 100                            //最大时间
#define INITIALVALUE 2.5                          //初值
#define DIMENS 170                              //网络最大规模

//说明:主函数为Network_1,调用NextValue来实现迭代,调用RandMatrix来实现把a[][]变成b[][];NextValue调用FaNexVal来实现梯形积分;
//结构体S_FList存放统计表;0——100分分为20个区间，score为3,8,13,18.....
using namespace std;

class mainFunc
{
public:
    mainFunc()
    {
        int ik;
        for(ik=0;ik<DIMENS;++ik)
        {
            p[ik]=18;q[ik]=9;r[ik]=2;nn[ik]=6;
        }
    }
    double MaxMa[DIMENS][DIMENS];
    double nong[DIMENS];
    void Network_1(double ReguMatrix[][DIMENS],int n);
    void Network_2(double Matr[][DIMENS],int n);
private:
    double p[DIMENS],q[DIMENS],r[DIMENS],nn[DIMENS];
};




inline void RandMatrix(double a[][DIMENS],double b[][DIMENS],const int n)
//拿一个调控矩阵产生系数矩阵，系数取+-6 7 8 9 10或者0;新基因调控关系可以取到+-3 4 5
//a[][]为输入矩阵,b[][]为输出矩阵,用数组形参输出,旧网络的规模为n.
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
                default :m1=rand()%2;                //其他,即正负都有
                if(m1)b[i1][j1]=(double)(rand()%101)/25.0+2.0;
                else b[i1][j1]=-(double)(rand()%101)/25.0-2.0;
            }
        }
        else                               //新基因调控关系的系数
        {
            m1=(int)(a[i1][j1]*100);
            if(rand()%100<m1)k1=m1>0?rand()%6+5:-rand()%6-5;       //调控关系和预测关系相同
            else if(rand()%5)k1=0;
            //调控关系和预测关系不同,相反关系和没有关系对半分
            else k1=-rand()%6-5;
            b[i1][j1]=k1;
        }
    }
}

double FaNexVal(double Matr[][DIMENS],double a[],const int n,const int i,    double p[],double q[],double nn[],double r[])
//矩形积分的下一个增量值,Matr[][]给出系数,a[]给出自己的浓度,n为网络规模,i为要增加的分量,p[]nn[]r[]为系数数组
{
    int j;
    double m[3];m[0]=0;                 //m[1]为正调控,m[2]为负调控
    m[0]=0-r[i]*a[i];    //降解量
    m[2]=0;m[1]=0;
    for(j=0;j<n;++j)
    {
        if(Matr[j][i]>0)m[1]+=pow(a[j],Matr[j][i]);       //正调控
        else if(Matr[j][i]<0)m[2]+=pow(a[j],-Matr[j][i]);  //负调控
    }
    if(m[1]==0)
    m[0]+=(p[i]/(nn[i])+q[i]/(nn[i]+m[2]));                   //总公式,即增量总和
    else m[0]+=(p[i]*m[1]/(nn[i]+m[1])+q[i]/(nn[i]+m[2]));
    return m[0];
}

double NextValue(double Matr[][DIMENS],double a[],const int n,const int i,double step ,   double p[],double q[],double nn[],double r[])
//梯形积分，形参意义同上
{
    int j;double b[DIMENS];//作为暂时空间
    for(j=0;j<n;++j)b[j]=a[j]+step*FaNexVal(Matr,a,n,j,p,q,nn,r)/2.0;
    //a[j]作为下一步高度的估计值与a[j]作平均,然后用这个估计值来取下一步增量
    return FaNexVal(Matr,b,n,i,p,q,nn,r);
}

void mainFunc::Network_1(double ReguMatrix[][DIMENS],int n)
//参数意义:ReguMatrix[][]作为调控矩阵,MaxMa[][]输出最大分数的系数矩阵
//S_Freq[]输出分数-频率表,结构体见前;n为旧网络规模.其余数组意义同上
{
    double a[DIMENS],b[DIMENS],c[DIMENS],d[DIMENS],e[DIMENS],f[DIMENS],sum(0),AbsValue(1),MaxScore(-1),Score(0);
    double FenShu[101];
    //a[]b[]c[]d[]临时记录表达产物的浓度值,交替迭代
    //AbsValue用来记录反应过程中相邻两个时刻的差值,检查平衡
    //MaxScore记录最大分数,Score记录分数
    int i(0),j=0,cou(0),k(0);double TempMatrix[DIMENS][DIMENS];double step=STEP;int pets=PETS;
    //cou记录反应时间,j记录迭代次数,k记录分数-频率统计次数,TempMatrix[][]记录临时产生的(n+1)*(n+1)的系数矩阵


    for(i=0;i<101;++i)FenShu[i]=0;
    for(i=0;i<n;++i)a[i]=b[i]=e[i]=INITIALVALUE;               //初始化浓度
    while(k<NN)                                  //统计次数小于预定值
    {
        ++k;
        RandMatrix(ReguMatrix,TempMatrix,n);  //产生临时随机矩阵放到TempMatrix里面
        AbsValue=10;j=0;step=STEP;pets=PETS;
        while(AbsValue>0.000001&&cou<MAXTIME)  // 没有稳定且反应时间没有达到预定时间
        {
            ++j;
            for(i=0;i<n;++i)
            {
                a[i]+=NextValue(TempMatrix,b,n,i,step,p,q,nn,r)*step;
                if(a[i]<0.000001){cou=MAXTIME+2;break;}
            }
            if(cou>=MAXTIME)break;                  //某一物质浓度太低直接终止
            for(i=0;i<n;++i)b[i]=a[i];
            if(j%(pets/8)==0)
            {
                AbsValue=0;
                for(i=0;i<n+1;++i){AbsValue+=fabs(a[i]-e[i]);}//检查平衡
                for(i=0;i<n+1;++i)e[i]=a[i];
                if(AbsValue<5&&pets==PETS&&j==PETS){step=STEP*4.0;pets=PETS/4;++cou;j=0;continue;}
            }
            if(j%pets==0){++cou;j=0;}
        }
        step=STEP;pets=PETS;
        if(cou>=MAXTIME)break;//没有稳定直接跳过,没有打分之类的
        for(i=0;i<n;++i){c[i]=d[i]=a[i];c[n]+=a[i];}//赋初值，上次稳定的值.ps:一般来说老网络稳定值都是一样的
        c[n]/=n;d[n]=c[n];cou=0;AbsValue=10;j=0;//最新的物质浓度设为其他的平均值
        //跑新网络
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
                for(i=0;i<n+1;++i){AbsValue+=fabs(c[i]-f[i]);}//检查平衡
                for(i=0;i<n+1;++i)f[i]=c[i];
                if(AbsValue<5&&pets==PETS&&j==PETS){step=STEP*4.0;pets=PETS/4;++cou;j=0;continue;}
            }
            if(j%pets==0){j=0;++cou;}
        }
        if(cou>=MAXTIME){sum+=1;FenShu[0]+=1;break;}                         //如果新网络没有稳定,直接跳过
        AbsValue=0;
        sum+=1;
        for(i=0;i<n;++i)AbsValue+=(fabs(c[i]-a[i]))/(c[i]>a[i]?a[i]:c[i]);              //计算新旧网络的浓度差
        Score=1-AbsValue*AbsValue/(n*n/9.0+AbsValue*AbsValue);  //每一种物质都差一倍基本上就废了
        Score*=100;
        Score*=(1-(pow(((double)cou)/MAXTIME,3)));
        if(Score>MaxScore)                                      //记录最大分数的矩阵
        {
            for(i=0;i<n+1;++i)
            for(j=0;j<n+1;++j)
            MaxMa[i][j]=TempMatrix[i][j];
            MaxScore=Score;
        }
        FenShu[(int)Score]+=1;                   //记录分数
    }
    for(i=0;i<101;++i)FenShu[i]/=sum;                  //sum为有结果的次数
    ofstream fi1;
    fi1.open("fenshu.txt");
    for(i=0;i<101;++i)
        fi1<<i+0.5<<' '<<FenShu[i]<<endl;//输出频率分数表
    fi1.close();
    //由最大分数的系数矩阵跑新网络的时间-浓度放到txt里面
    for(i=0;i<n;++i)a[i]=b[i]=e[i]=INITIALVALUE;                               //赋初值
    AbsValue=10;
    j=0;cou=0;step=STEP;pets=PETS;
    while(AbsValue>0.000001&&cou<MAXTIME)                      //跑旧网络
    {
        ++j;
        for(i=0;i<n;++i){a[i]+=NextValue(MaxMa,b,n,i,step,p,q,nn,r)*step;}
        for(i=0;i<n;++i)b[i]=a[i];
        if(j%(pets/8)==0)
        {
            AbsValue=0;
            for(i=0;i<n+1;++i){AbsValue+=fabs(a[i]-e[i]);}//检查平衡
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
    while(AbsValue>0.00000001&&cou<MAXTIME)              //跑新网络
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
            for(i=0;i<n+1;++i){AbsValue+=fabs(c[i]-f[i]);}//检查平衡
            for(i=0;i<n+1;++i)f[i]=c[i];
            if(AbsValue<5&&pets==PETS&&j==PETS){step=STEP*4.0;pets=PETS/4;++cou;j=0;continue;}
        }
        if(j%pets==0)
        {
            ++cou;j=0;
        }
    }
}



void mainFunc::Network_2(double Matr[][DIMENS],int n)
//输入系数矩阵和网络规模,输出稳定物质.若不能稳定或者稳定时间太长,则a[]所有值置为-1,
//p[] q[] nn[] r[]意义同上
{
    int i,j=0,cou=0;
    double b[DIMENS],AbsValue=10, c[DIMENS];
    for(i=0;i<n;++i){nong[i]=INITIALVALUE;b[i]=nong[i];c[i]=b[i];}
    int pets=PETS;double step=STEP;
    while(AbsValue>0.0000001&&cou<MAXTIME)  // 没有稳定且反应时间没有达到预定时间
    {
        ++j;
        for(i=0;i<n;++i){nong[i]+=NextValue(Matr,b,n,i,step,p,q,nn,r)*step;if(nong[i]<0.0000001){cou=MAXTIME+2;break;}}
        if(cou>=MAXTIME)break;                  //某一物质浓度太低直接终止
        for(i=0;i<n;++i)b[i]=nong[i];
        if(j%(pets/8)==0)
        {
            AbsValue=0;
            for(i=0;i<n;++i){AbsValue+=fabs(nong[i]-c[i]);}//检查平衡
            for(i=0;i<n;++i)c[i]=nong[i];
            if((AbsValue<5)&&(pets==PETS)&&(j==PETS)){step=STEP*4.0;pets=PETS/4;++cou;j=0;continue;}
        }
        if(j%pets==0){++cou;j=0;}
    }
    if(cou>=MAXTIME)                                            //如果没有稳定,值全为-1
    for(i=0;i<n+1;++i)nong[i]=-1;
}