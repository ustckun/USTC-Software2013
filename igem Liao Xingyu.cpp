#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <ctime>
#define NN 1                                   //ͳ�ƴ���
#define PETS 128                                //�����ĵ���
#define STEP (1.0/PETS)                           //����
#define MAXTIME 100                            //���ʱ��
#define INITIALVALUE 2.5                          //��ֵ
#define DIMENS 170                              //��������ģ

//˵��:������ΪNetwork_1,����NextValue��ʵ�ֵ���,����RandMatrix��ʵ�ְ�a[][]���b[][];NextValue����FaNexVal��ʵ�����λ���;
//�ṹ��S_FList���ͳ�Ʊ�;0����100�ַ�Ϊ20�����䣬scoreΪ3,8,13,18.....
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
//��һ�����ؾ������ϵ������ϵ��ȡ+-6 7 8 9 10����0;�»�����ع�ϵ����ȡ��+-3 4 5
//a[][]Ϊ�������,b[][]Ϊ�������,�������β����,������Ĺ�ģΪn.
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
                default :m1=rand()%2;                //����,����������
                if(m1)b[i1][j1]=(double)(rand()%101)/25.0+2.0;
                else b[i1][j1]=-(double)(rand()%101)/25.0-2.0;
            }
        }
        else                               //�»�����ع�ϵ��ϵ��
        {
            m1=(int)(a[i1][j1]*100);
            if(rand()%100<m1)k1=m1>0?rand()%6+5:-rand()%6-5;       //���ع�ϵ��Ԥ���ϵ��ͬ
            else if(rand()%5)k1=0;
            //���ع�ϵ��Ԥ���ϵ��ͬ,�෴��ϵ��û�й�ϵ�԰��
            else k1=-rand()%6-5;
            b[i1][j1]=k1;
        }
    }
}

double FaNexVal(double Matr[][DIMENS],double a[],const int n,const int i,    double p[],double q[],double nn[],double r[])
//���λ��ֵ���һ������ֵ,Matr[][]����ϵ��,a[]�����Լ���Ũ��,nΪ�����ģ,iΪҪ���ӵķ���,p[]nn[]r[]Ϊϵ������
{
    int j;
    double m[3];m[0]=0;                 //m[1]Ϊ������,m[2]Ϊ������
    m[0]=0-r[i]*a[i];    //������
    m[2]=0;m[1]=0;
    for(j=0;j<n;++j)
    {
        if(Matr[j][i]>0)m[1]+=pow(a[j],Matr[j][i]);       //������
        else if(Matr[j][i]<0)m[2]+=pow(a[j],-Matr[j][i]);  //������
    }
    if(m[1]==0)
    m[0]+=(p[i]/(nn[i])+q[i]/(nn[i]+m[2]));                   //�ܹ�ʽ,�������ܺ�
    else m[0]+=(p[i]*m[1]/(nn[i]+m[1])+q[i]/(nn[i]+m[2]));
    return m[0];
}

double NextValue(double Matr[][DIMENS],double a[],const int n,const int i,double step ,   double p[],double q[],double nn[],double r[])
//���λ��֣��β�����ͬ��
{
    int j;double b[DIMENS];//��Ϊ��ʱ�ռ�
    for(j=0;j<n;++j)b[j]=a[j]+step*FaNexVal(Matr,a,n,j,p,q,nn,r)/2.0;
    //a[j]��Ϊ��һ���߶ȵĹ���ֵ��a[j]��ƽ��,Ȼ�����������ֵ��ȡ��һ������
    return FaNexVal(Matr,b,n,i,p,q,nn,r);
}

void mainFunc::Network_1(double ReguMatrix[][DIMENS],int n)
//��������:ReguMatrix[][]��Ϊ���ؾ���,MaxMa[][]�����������ϵ������
//S_Freq[]�������-Ƶ�ʱ�,�ṹ���ǰ;nΪ�������ģ.������������ͬ��
{
    double a[DIMENS],b[DIMENS],c[DIMENS],d[DIMENS],e[DIMENS],f[DIMENS],sum(0),AbsValue(1),MaxScore(-1),Score(0);
    double FenShu[101];
    //a[]b[]c[]d[]��ʱ��¼�������Ũ��ֵ,�������
    //AbsValue������¼��Ӧ��������������ʱ�̵Ĳ�ֵ,���ƽ��
    //MaxScore��¼������,Score��¼����
    int i(0),j=0,cou(0),k(0);double TempMatrix[DIMENS][DIMENS];double step=STEP;int pets=PETS;
    //cou��¼��Ӧʱ��,j��¼��������,k��¼����-Ƶ��ͳ�ƴ���,TempMatrix[][]��¼��ʱ������(n+1)*(n+1)��ϵ������


    for(i=0;i<101;++i)FenShu[i]=0;
    for(i=0;i<n;++i)a[i]=b[i]=e[i]=INITIALVALUE;               //��ʼ��Ũ��
    while(k<NN)                                  //ͳ�ƴ���С��Ԥ��ֵ
    {
        ++k;
        RandMatrix(ReguMatrix,TempMatrix,n);  //������ʱ�������ŵ�TempMatrix����
        AbsValue=10;j=0;step=STEP;pets=PETS;
        while(AbsValue>0.000001&&cou<MAXTIME)  // û���ȶ��ҷ�Ӧʱ��û�дﵽԤ��ʱ��
        {
            ++j;
            for(i=0;i<n;++i)
            {
                a[i]+=NextValue(TempMatrix,b,n,i,step,p,q,nn,r)*step;
                if(a[i]<0.000001){cou=MAXTIME+2;break;}
            }
            if(cou>=MAXTIME)break;                  //ĳһ����Ũ��̫��ֱ����ֹ
            for(i=0;i<n;++i)b[i]=a[i];
            if(j%(pets/8)==0)
            {
                AbsValue=0;
                for(i=0;i<n+1;++i){AbsValue+=fabs(a[i]-e[i]);}//���ƽ��
                for(i=0;i<n+1;++i)e[i]=a[i];
                if(AbsValue<5&&pets==PETS&&j==PETS){step=STEP*4.0;pets=PETS/4;++cou;j=0;continue;}
            }
            if(j%pets==0){++cou;j=0;}
        }
        step=STEP;pets=PETS;
        if(cou>=MAXTIME)break;//û���ȶ�ֱ������,û�д��֮���
        for(i=0;i<n;++i){c[i]=d[i]=a[i];c[n]+=a[i];}//����ֵ���ϴ��ȶ���ֵ.ps:һ����˵�������ȶ�ֵ����һ����
        c[n]/=n;d[n]=c[n];cou=0;AbsValue=10;j=0;//���µ�����Ũ����Ϊ������ƽ��ֵ
        //��������
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
                for(i=0;i<n+1;++i){AbsValue+=fabs(c[i]-f[i]);}//���ƽ��
                for(i=0;i<n+1;++i)f[i]=c[i];
                if(AbsValue<5&&pets==PETS&&j==PETS){step=STEP*4.0;pets=PETS/4;++cou;j=0;continue;}
            }
            if(j%pets==0){j=0;++cou;}
        }
        if(cou>=MAXTIME){sum+=1;FenShu[0]+=1;break;}                         //���������û���ȶ�,ֱ������
        AbsValue=0;
        sum+=1;
        for(i=0;i<n;++i)AbsValue+=(fabs(c[i]-a[i]))/(c[i]>a[i]?a[i]:c[i]);              //�����¾������Ũ�Ȳ�
        Score=1-AbsValue*AbsValue/(n*n/9.0+AbsValue*AbsValue);  //ÿһ�����ʶ���һ�������Ͼͷ���
        Score*=100;
        Score*=(1-(pow(((double)cou)/MAXTIME,3)));
        if(Score>MaxScore)                                      //��¼�������ľ���
        {
            for(i=0;i<n+1;++i)
            for(j=0;j<n+1;++j)
            MaxMa[i][j]=TempMatrix[i][j];
            MaxScore=Score;
        }
        FenShu[(int)Score]+=1;                   //��¼����
    }
    for(i=0;i<101;++i)FenShu[i]/=sum;                  //sumΪ�н���Ĵ���
    ofstream fi1;
    fi1.open("fenshu.txt");
    for(i=0;i<101;++i)
        fi1<<i+0.5<<' '<<FenShu[i]<<endl;//���Ƶ�ʷ�����
    fi1.close();
    //����������ϵ���������������ʱ��-Ũ�ȷŵ�txt����
    for(i=0;i<n;++i)a[i]=b[i]=e[i]=INITIALVALUE;                               //����ֵ
    AbsValue=10;
    j=0;cou=0;step=STEP;pets=PETS;
    while(AbsValue>0.000001&&cou<MAXTIME)                      //�ܾ�����
    {
        ++j;
        for(i=0;i<n;++i){a[i]+=NextValue(MaxMa,b,n,i,step,p,q,nn,r)*step;}
        for(i=0;i<n;++i)b[i]=a[i];
        if(j%(pets/8)==0)
        {
            AbsValue=0;
            for(i=0;i<n+1;++i){AbsValue+=fabs(a[i]-e[i]);}//���ƽ��
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
    while(AbsValue>0.00000001&&cou<MAXTIME)              //��������
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
            for(i=0;i<n+1;++i){AbsValue+=fabs(c[i]-f[i]);}//���ƽ��
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
//����ϵ������������ģ,����ȶ�����.�������ȶ������ȶ�ʱ��̫��,��a[]����ֵ��Ϊ-1,
//p[] q[] nn[] r[]����ͬ��
{
    int i,j=0,cou=0;
    double b[DIMENS],AbsValue=10, c[DIMENS];
    for(i=0;i<n;++i){nong[i]=INITIALVALUE;b[i]=nong[i];c[i]=b[i];}
    int pets=PETS;double step=STEP;
    while(AbsValue>0.0000001&&cou<MAXTIME)  // û���ȶ��ҷ�Ӧʱ��û�дﵽԤ��ʱ��
    {
        ++j;
        for(i=0;i<n;++i){nong[i]+=NextValue(Matr,b,n,i,step,p,q,nn,r)*step;if(nong[i]<0.0000001){cou=MAXTIME+2;break;}}
        if(cou>=MAXTIME)break;                  //ĳһ����Ũ��̫��ֱ����ֹ
        for(i=0;i<n;++i)b[i]=nong[i];
        if(j%(pets/8)==0)
        {
            AbsValue=0;
            for(i=0;i<n;++i){AbsValue+=fabs(nong[i]-c[i]);}//���ƽ��
            for(i=0;i<n;++i)c[i]=nong[i];
            if((AbsValue<5)&&(pets==PETS)&&(j==PETS)){step=STEP*4.0;pets=PETS/4;++cou;j=0;continue;}
        }
        if(j%pets==0){++cou;j=0;}
    }
    if(cou>=MAXTIME)                                            //���û���ȶ�,ֵȫΪ-1
    for(i=0;i<n+1;++i)nong[i]=-1;
}