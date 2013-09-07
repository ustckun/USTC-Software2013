#include "widget.h"
#include "ui_widget.h"
#include "console.h"
#include "ui_console.h"

//delete margin
Widget::Widget(QWidget *parent,Qt::WindowFlags flags) :
    QWidget(parent,flags),
    ui(new Ui::Widget)
{
    ui->setupUi(this);
    setWindowFlags(Qt::FramelessWindowHint);
    ui->Do->setText("Ready?");
    ui->doing->setText("Software for pioneer...\n");
    setAttribute(Qt::WA_TranslucentBackground, true);
    setWindowTitle("gNAP");
    ui->gene->close();
    ui->Gene_done->close();
    ui->promoter->close();
    ui->Promoter_done->close();
    ui->progress_bar_start->close();
    ui->Take_a_nap->close();
    ui->Promoter_Info->close();
    ui->TU_Info->close();
    ui->Gene_Info->close();
    ui->TF_Gene->close();
    console_ui=new console;
    QAction::connect(console_ui->ui->back,SIGNAL(clicked()),this,SLOT(back_func()));
}

Widget::~Widget()
{
    delete ui;
}

void Widget::back_func()
{
    this->show();
    console_ui->hide();
    ui->gene->hide();
    ui->Gene_done->hide();
}

void Widget::on_TF_TF_clicked()
{
    QString TF_TF_name = QFileDialog::getOpenFileName(this);
    TF_TF_address=TF_TF_name.toStdString();
    int length=(int)TF_TF_address.size();
    string temp;
    int i;
    for(i=0;TF_TF_address[length-i]!='/';i++)
    {
    }
    char *p=&TF_TF_address[length-i+1];
    temp=p;
    if(temp=="TF-TF")
    {
        ui->TF_TF->close();
        ui->TF_Gene->show();
        ui->Do->setText("Loading TF-Gene...");
        ui->doing->setText(ui->doing->text()+"TF-TF file: Done!\n");
    }
    else
    {
        ui->Do->setText("Please try again...");
        ui->doing->setText(ui->doing->text()+"TF-TF file: Failed!\n");
    }
}

void Widget::on_TF_Gene_clicked()
{
    QString TF_Gene_name = QFileDialog::getOpenFileName(this);
    TF_Gene_address=TF_Gene_name.toStdString();
    int length=(int)TF_Gene_address.size();
    string temp;
    int i;
    for(i=0;TF_Gene_address[length-i]!='/';i++)
    {
    }
    char *p=&TF_Gene_address[length-i+1];
    temp=p;
    if(temp=="TF-Gene")
    {
        ui->TF_Gene->close();
        ui->Gene_Info->show();
        ui->Do->setText("Loading Gene Info...");
        ui->doing->setText(ui->doing->text()+"TF-Gene file: Done!\n");
    }
    else
    {
        ui->Do->setText("Please try again...");
        ui->doing->setText(ui->doing->text()+"TF-Gene file: Failed!\n");
    }
}

void Widget::on_Gene_Info_clicked()
{
    QString Gene_Info_name = QFileDialog::getOpenFileName(this);
    gene_info_address=Gene_Info_name.toStdString();
    int length=(int)gene_info_address.size();
    string temp;
    int i;
    for(i=0;gene_info_address[length-i]!='/';i++)
    {
        //temp[i]=gene_info_address[length-i];
    }
    char *p=&gene_info_address[length-i+1];
    temp=p;
    if(temp=="Gene Info")
    {
        ui->Gene_Info->close();
        ui->TU_Info->show();
        ui->Do->setText("Loading TU Info...");
        ui->doing->setText(ui->doing->text()+"Gene Info file: Done!\n");
    }
    else
    {
        ui->Do->setText("Please try again...");
        ui->doing->setText(ui->doing->text()+"Gene Info file: Failed!\n");
    }
}

void Widget::on_TU_Info_clicked()
{
    QString TU_Info_name = QFileDialog::getOpenFileName(this);
    TU_info_address=TU_Info_name.toStdString();
    int length=(int)TU_info_address.size();
    string temp;
    int i;
    for(i=0;TU_info_address[length-i]!='/';i++)
    {
        //temp[i]=TU_info_address[length-i];
    }
    char *p=&TU_info_address[length-i+1];
    temp=p;
    if(temp=="TU Info")
    {
        ui->TU_Info->close();
        ui->Promoter_Info->show();
        ui->Do->setText("Loading Promoter Info...");
        ui->doing->setText(ui->doing->text()+"TU Info file: Done!\n");
    }
    else
    {
        ui->Do->setText("Please try again...");
        ui->doing->setText(ui->doing->text()+"TU Info file: Failed!\n");
    }
}

void Widget::on_Promoter_Info_clicked()
{
    QString Promoter_Info_name = QFileDialog::getOpenFileName(this);
    promoter_info_address=Promoter_Info_name.toStdString();
    int length=(int)promoter_info_address.size();
    string temp;
    int i;
    for(i=0;promoter_info_address[length-i]!='/';i++)
    {
        //temp[i]=promoter_info_address[length-i];
    }
    char *p=&promoter_info_address[length-i+1];
    temp=p;
    if(temp=="Promoter Info")
    {
        ui->Promoter_Info->close();
        ui->Take_a_nap->show();
        ui->Do->setText("Ready to read...");
        ui->doing->setText(ui->doing->text()+"Promoter Info file: Done!\n");
    }
    else
    {
        ui->Do->setText("Please try again...");
        ui->doing->setText(ui->doing->text()+"Promoter Info file: Failed!\n");
    }
}

void Widget::mousePressEvent(QMouseEvent *event)
{
    if (event->button() == Qt::LeftButton) {
        m_DragPosition = event->globalPos() - this->pos();
        event->accept();
    }
}

void Widget::mouseMoveEvent(QMouseEvent *event)
{
    if (event->buttons() && Qt::LeftButton) {
        move(event->globalPos() - m_DragPosition);
        event->accept();
    }
}

void Widget::on_close_clicked()
{
    this->close();
}

void Widget::on_Take_a_nap_clicked()
{
    GetReady setup;
    ui->progress_bar_start->show();
    ui->Do->setText("Get Regu...");
    setup.getRegulationMatrix(ecoli,TF_TF_address,TF_Gene_address);
    setup.inputUncertainGene();
    int G_A=setup.getGeneAmount();
    ui->progress_bar_start->setMaximum(G_A*3+GENEAM);
    map<string,string> dict;
    ui->doing->setText("Get Regu: Done!\n");
    ui->Do->setText("Map Gene Info...");
    dict=setup.mapTFIM(gene_info_address);
    for(int i=0;i<setup.getGeneAmount();i++)
    {
        ui->progress_bar_start->setValue(i);
        ecoli[i].getGeneInformation(dict);
    }
    setup.readTUPosition(TU_info_address);
    setup.getGenePromoter(ecoli);
    ui->doing->setText(ui->doing->text()+"Map Gene Info: Done!\n");
    ui->Do->setText("Map Promoter Info...");
    dict=setup.mapPromoter(promoter_info_address);
    for(int i=0;i<setup.getGeneAmount();i++)
    {
        ui->progress_bar_start->setValue(i+G_A);
        ecoli[i].getPromoterIF(dict);
    }
    /*for(int i=0;i<50000;i++)
    {
        for(int j=0;j<2000;j++)
            bar->ui->DNA_info->setValue(i);
    }*/
    ui->doing->setText(ui->doing->text()+"Map Promoter Info: Done!\n");
    //ui->Do->setText("Write out Info...");
    ofstream all_info("all_info");
    for(int i=0;i<setup.getGeneAmount();i++)
    {
        ui->progress_bar_start->setValue(i+G_A*2);
        all_info<<ecoli[i].gene_number<<'\t';
        all_info<<ecoli[i].getPromoterSequence()<<'\t';
        all_info<<ecoli[i].getGeneSequence()<<'\t';
        all_info<<ecoli[i].getGeneTrueName()<<'\t';
        all_info<<ecoli[i].getID()<<'\t';
        all_info<<ecoli[i].getLeftPosition()<<'\t';
        all_info<<ecoli[i].getRightPosition()<<'\t';
        all_info<<ecoli[i].getPromoterName()<<'\t';
        all_info<<ecoli[i].getGeneDescription()<<endl;
    }
    all_info.close();
    ofstream Regu_Matrix("old_GRN");
    for(int i=0;i<GENEAM;i++)
    {
        ui->progress_bar_start->setValue(i+G_A*3);
        for(int j=0;j<TFScale;j++)
        {
            Regu_Matrix<<setup.originalGRN[i][j]<<" ";
        }
        Regu_Matrix<<endl;
    }
    ui->doing->setText(ui->doing->text()+"Write out Info: Done!\n");
    ui->Do->setText("Read Sequence...");
    flag=1;
    ui->promoter->show();
    ui->Promoter_done->show();
    QAction::connect(this,SIGNAL(pass_amount(int,int)),console_ui,SLOT(getAmount(int,int)));
    emit pass_amount(setup.getGeneAmount(),setup.getTFAmount());
}

void Widget::on_Promoter_done_clicked()
{
    QString temp=ui->promoter->text();
    promoter_sequence=temp.toStdString();
    ui->gene->show();
    ui->Gene_done->show();
}

void Widget::on_Gene_done_clicked()
{
    QString temp=ui->gene->text();
    gene_sequence=temp.toStdString();
    QAction::connect(this,SIGNAL(pass_sequence(string,string)),console_ui,SLOT(getSequence(string,string)));
    emit pass_sequence(promoter_sequence,gene_sequence);
    for(int i=0;i<GENEAM;i++)
    {
        QAction::connect(this,SIGNAL(pass_information(GeneIM,int)),console_ui,SLOT(getInformation(GeneIM,int)));
        emit pass_information(ecoli[i],i);
    }
    this->hide();
    console_ui->show();
    if((int)promoter_sequence.size()<50||(int)gene_sequence.size()<50)
    {
        console_ui->ui->analyze->setStyleSheet("#analyze{border-image: url(:/picture/picture/gray1.png);color:white;font: 75 18pt \"Arial\";}#analyze:hover{border-image: url(:/picture/picture/gray2.png);}");
        console_ui->analyze_flag=1;
        console_ui->ui->tips->setText("Sleepy?");
    }
    else
    {
        console_ui->ui->analyze->setStyleSheet("#analyze{border-image: url(:/picture/picture/blue1.png);color:white;font: 75 18pt \"Arial\";}#analyze:hover{border-image: url(:/picture/picture/blue2.png);}");
        console_ui->analyze_flag=0;
        console_ui->ui->tips->setText("Sleepy?");
    }
}
