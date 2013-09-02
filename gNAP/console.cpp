#include "console.h"
#include "ui_console.h"
#include "result.h"
#include "ui_result.h"

console::console(QWidget *parent,Qt::WindowFlags flags) :
    QWidget(parent,flags),
    ui(new Ui::console)
{
    ui->setupUi(this);
    setWindowFlags(Qt::FramelessWindowHint);
    setAttribute(Qt::WA_TranslucentBackground, true);
    setWindowTitle("gNAP");
    n=-1;
    flag=0;
    analyze_flag=0;
    if_analyze=0;
    if_predict=0;
    ui->cout->insertPlainText("$>Hello iGEM!\n$>_");
    ui->tips->setText("Sleepy?");
    result_ui=new result;
    result_ui->ui->choose->setMaxVisibleItems(10);
    result_ui->ui->choose->setEditable(true);
    QValidator *validator = new QRegExpValidator;
    result_ui->ui->choose->lineEdit()->setValidator(validator);
    result_ui->hide();
    QAction::connect(result_ui->ui->back,SIGNAL(clicked()),this,SLOT(show_out()));
    QAction::connect(result_ui->ui->done,SIGNAL(clicked()),this,SLOT(done_clicked()));
    QAction::connect(result_ui->ui->sbol,SIGNAL(clicked()),this,SLOT(sbol_clicked()));
}

console::~console()
{
    delete ui;
}

void console::show_out()
{
    this->show();
    result_ui->hide();
}

void console::mousePressEvent(QMouseEvent *event)
{
    if (event->button() == Qt::LeftButton) {
        m_DragPosition = event->globalPos() - this->pos();
        event->accept();
    }
}

void console::mouseMoveEvent(QMouseEvent *event)
{
    if (event->buttons() && Qt::LeftButton) {
        move(event->globalPos() - m_DragPosition);
        event->accept();
    }
}

void console::getAmount(int row,int column)
{
    row_number=row;
    column_number=column;
}

void console::getSequence(string promoter,string gene)
{
    new_promoter_sequence=promoter;
    new_gene_sequence=gene;
}

void console::getInformation(GeneIM ecoli,int position)
{
    ecoli_k12[position]=ecoli;
}

void console::set_high()
{
    low_high[n]=1;
    choose_high[n]->setStyleSheet("#choose_high{border-image: url(:/picture/picture/up1.png);}#choose_high:hover{border-image:url(:/picture/picture/up2.png);}");
    choose_low[n]->setStyleSheet("#choose_low{border-image: url(:/picture/picture/down_hold1.png);}#choose_low:hover{border-image:url(:/picture/picture/down_hold2.png);}");
}

void console::set_low()
{
    low_high[n]=0;
    choose_high[n]->setStyleSheet("#choose_high{border-image: url(:/picture/picture/up_hold1.png);}#choose_high:hover{border-image:url(:/picture/picture/up_hold2.png);}");
    choose_low[n]->setStyleSheet("#choose_low{border-image: url(:/picture/picture/down1.png);}#choose_low:hover{border-image:url(:/picture/picture/down2.png);}");
}

void console::on_analyze_clicked()
{
    if(analyze_flag==1)
    {
        ui->tips->setText("No sequences input. Please enter sequences in backword page.");
    }
    else{
    if_analyze=1;
    ui->progress_bar->setMaximum(6500);
    QString output_temp;
    //QObject::connect(,SIGNAL());
    double **old_GRN;
    old_GRN = new double*[GENEAM];
    for (int i = 0; i != GENEAM; ++i) {
        old_GRN[i] = new double[TFScale];
    }
    ifstream infile;
    ui->tips->setText("Reading gene regulation network...");
    infile.open("old_GRN");
    for (int i = 0; i != GENEAM; ++i) {
        for (int j = 0; j != TFScale; ++j) {
            infile >> old_GRN[i][j];
            //ui->cout->setPlainText(ui->cout->toPlainText()+QString::number(old_GRN[i][j])+"\n$>_");

        }
        output_temp="Get No."+QString::number(i)+"'s regulation...\n$>";
        ui->cout->insertPlainText(output_temp);
        ui->cout->moveCursor(QTextCursor::End);
        ui->progress_bar->setValue(i);
    }
    infile.close();
    GRN network;
    ui->cout->insertPlainText("Initialize_GRN...\n$>");
    ui->cout->moveCursor(QTextCursor::End);
    ui->progress_bar->setValue(ui->progress_bar->value()+20);
    ui->tips->setText("Initializing gene sequence...");
    network.initialize_GRN(old_GRN, row_number, column_number);
    for (int i = 0; i != row_number; ++i) {
        sequence_array[i].initialize_Sequence(i,
                                              ecoli_k12[i].getPromoterSequence(),
                                              (int)ecoli_k12[i].getPromoterSequence().size(),
                                              ecoli_k12[i].getGeneSequence(),
                                              (int)ecoli_k12[i].getGeneSequence().size());
        //ui->cout->setPlainText(ui->cout->toPlainText()+QString::fromStdString(ecoli_k12[i].getPromoterSequence())+"\n$>_");
        //ui->cout->setPlainText(ui->cout->toPlainText()+QString::fromStdString(ecoli_k12[i].getGeneSequence())+"\n$>_");
        output_temp="Initialize sequence: No."+QString::number(i)+"\n$>";
        ui->cout->insertPlainText(output_temp);
        ui->cout->moveCursor(QTextCursor::End);
        ui->progress_bar->setValue(1820+i);
    }
    int new_num = row_number;
    sequence_array[row_number].initialize_Sequence(new_num,
                                             new_promoter_sequence,
                                             (int)new_promoter_sequence.size(),
                                             new_gene_sequence,
                                             (int)new_gene_sequence.size());
    ui->tips->setText("Constructing New GRN...");
    ui->cout->insertPlainText("Construct New GRN...\n$>_");
    ui->cout->moveCursor(QTextCursor::End);
    ui->progress_bar->setValue(ui->progress_bar->value()+20);
    network.construct_new_GRN(sequence_array);
    ofstream new_Network("new_GRN");
    for(int i=0;i<GENEAM;i++)
    {
        for(int j=0;j<TFScale;j++)
        {
            new_Network<<network.new_GRN[i][j]<<'\t';
            //ui->cout->setPlainText(ui->cout->toPlainText()+QString::number(network.new_GRN[i][j])+"\n$>");

        }
        new_Network<<endl;
        output_temp="Output No."+QString::number(i)+"'s regulation...\n$>";
        ui->cout->insertPlainText(output_temp);
        ui->cout->moveCursor(QTextCursor::End);
        ui->progress_bar->setValue(3640+i);
    }
    ModleNetwork calculate;
    ui->tips->setText("Calculating GRN...");
    ui->cout->insertPlainText("Calculate Network...\n$>_");
    ui->cout->moveCursor(QTextCursor::End);
    ui->progress_bar->setValue(ui->progress_bar->value()+20);
    calculate.Network_1(network.new_GRN,column_number,row_number);
    ui->cout->insertPlainText("Score:\n$>");
    ui->cout->moveCursor(QTextCursor::End);
    ui->progress_bar->setValue(ui->progress_bar->value()+20);
    ifstream score("Score");
    double temp_score1,temp_score2;
    for(int i=0;i<101;i++)
    {
        score>>temp_score1;
        ui->cout->insertPlainText(QString::number(temp_score1)+"\t");
        ui->cout->moveCursor(QTextCursor::End);
        score>>temp_score2;
        if(most_score<temp_score2)
            best_score=temp_score1;
        ui->cout->insertPlainText(QString::number(temp_score2)+"\n$>");
        ui->cout->moveCursor(QTextCursor::End);
        ui->progress_bar->setValue(5485+i*10);
    }
    ui->tips->setText("More information could be get in \"Result\"!");
    ui->cout->insertPlainText("Analyze Done!\n$>");
    ui->cout->moveCursor(QTextCursor::End);
    ui->cout->insertPlainText("Ready to work...\n$>_");
    ui->cout->moveCursor(QTextCursor::End);
    ui->progress_bar->setValue(6500);
    ui->progress_bar->hide();
    }
}

void console::on_add_clicked()
{
    n++;
    int a,b,c,d;
    a=ui->add->pos().rx();
    c=a+270;
    d=a+310;
    b=ui->add->pos().ry();
    pick[n]=new QComboBox(this);
    pick[n]->setGeometry(a,b,270,40);
    pick[n]->setStyleSheet("font: 87 18pt \"Arial\";");
    for(int i=0;i<row_number;i++)
    {
        pick[n]->addItem(QString::fromStdString(ecoli_k12[i].getGeneName()));
    }
    pick[n]->setMaxVisibleItems(10);
    pick[n]->setEditable(true);
    QValidator *validator = new QRegExpValidator;
    pick[n]->lineEdit()->setValidator(validator);
    pick[n]->show();
    hold_low[n]=new QFrame(this);
    hold_low[n]->setGeometry(c,b,40,40);
    hold_low[n]->setStyleSheet("border-image: url(:/picture/picture/down3.png);");
    hold_low[n]->show();
    choose_low[n]=new QPushButton(this);
    choose_low[n]->setGeometry(c,b,40,40);
    choose_low[n]->setObjectName("choose_low");
    choose_low[n]->setStyleSheet("#choose_low{border-image: url(:/picture/picture/down1.png);}#choose_low:hover{border-image:url(:/picture/picture/down2.png);}");
    choose_low[n]->setCursor(Qt::PointingHandCursor);
    choose_low[n]->show();
    hold_high[n]=new QFrame(this);
    hold_high[n]->setGeometry(d,b,40,40);
    hold_high[n]->setStyleSheet("border-image: url(:/picture/picture/up3.png);");
    hold_high[n]->show();
    choose_high[n]=new QPushButton(this);
    choose_high[n]->setGeometry(d,b,40,40);
    choose_high[n]->setObjectName("choose_high");
    choose_high[n]->setStyleSheet("#choose_high{border-image: url(:/picture/picture/up1.png);}#choose_high:hover{border-image:url(:/picture/picture/up2.png);}");
    choose_high[n]->setCursor(Qt::PointingHandCursor);
    choose_high[n]->show();
    ui->add->move(a,b+40);
    QAction::connect(choose_low[n],SIGNAL(clicked()),this,SLOT(set_low()));
    QAction::connect(choose_high[n],SIGNAL(clicked()),this,SLOT(set_high()));
    if(n==6)
        ui->add->close();
}

void console::on_predict_clicked()
{
    if_predict=1;
    ui->progress_bar->show();
    ui->progress_bar->setValue(0);
    ui->progress_bar->setMaximum(1950);
    QString temp_string;
    ModleNetwork calculate;
    ui->tips->setText("Initializing PSO...");
    ui->cout->setText("");
    ui->cout->insertPlainText("$>Initialize PSO method:_\n");
    ui->progress_bar->setValue(10);
    PSO predict_regu(calculate,row_number,column_number);
    ui->progress_bar->setValue(ui->progress_bar->value()+10);
    ui->cout->insertPlainText("$>Get gene expression range...\n");
    predict_regu.getRange(row_number,column_number,calculate);
    for(int i=0;i<n+1;i++)
    {
        if(low_high[i]==0)
            predict_regu.target[pick[i]->currentIndex()]=predict_regu.random_max[pick[i]->currentIndex()];
        else
            predict_regu.target[pick[i]->currentIndex()]=predict_regu.random_min[pick[i]->currentIndex()];
        temp_string="$>No."+QString::number(pick[i]->currentIndex())+"need change\n";
        ui->cout->insertPlainText(temp_string);
        ui->progress_bar->setValue(ui->progress_bar->value()+10);
    }
    ui->tips->setText("Runing PSO...");
    temp_string="$>Using PSO to predict...\n$>";
    ui->cout->insertPlainText(temp_string);
    ui->progress_bar->setValue(ui->progress_bar->value()+10);
    predict_regu.getPrediction(calculate,row_number,column_number);
    ui->tips->setText("Putting out predicted interaction...");
    for(int i=0;i<column_number;i++)
    {
        if(predict_regu.edPick[i]>0.9)
        {
            text=text+QString::fromStdString(ecoli_k12[i].getGeneName())+"->new"+'\t'
                    +"+"+'\n';
        }
        if(predict_regu.edPick[i]<-0.9)
        {
            text=text+QString::fromStdString(ecoli_k12[i].getGeneName())+"->new"+'\t'
                    +"-"+'\n';
        }
        temp_string=QString::fromStdString(ecoli_k12[i].getGeneName())+"->new"+'\t'
                +QString::number(predict_regu.edPick[i])+"\n$>";
        ui->cout->insertPlainText(temp_string);
        ui->cout->moveCursor(QTextCursor::End);
        ui->progress_bar->setValue(ui->progress_bar->value()+1);
    }
    for(int i=0;i<row_number;i++)
    {
        if(predict_regu.toPick[i]>0.9)
        {
            text=text+"new->"+QString::fromStdString(ecoli_k12[i].getGeneName())+'\t'
                    +"+"+'\n';
        }
        if(predict_regu.toPick[i]<-0.9)
        {
            text=text+"new->"+QString::fromStdString(ecoli_k12[i].getGeneName())+'\t'
                    +"-"+'\n';
        }
        temp_string="new->"+QString::fromStdString(ecoli_k12[i].getGeneName())+'\t'
                +QString::number(predict_regu.toPick[i])+"\n$>";
        ui->cout->insertPlainText(temp_string);
        ui->cout->moveCursor(QTextCursor::End);
        ui->progress_bar->setValue(ui->progress_bar->value()+1);
    }
    ui->cout->insertPlainText("Predict Done!\n$>");
    ui->cout->moveCursor(QTextCursor::End);
    ui->cout->insertPlainText("Ready to work...\n$>_");
    ui->cout->moveCursor(QTextCursor::End);
    ui->tips->setText("More information could be get in \"Result\"!");
    ui->progress_bar->hide();
}

void console::on_result_clicked()
{
    result_ui->show();
    this->hide();
    if(if_analyze==1)
    {
    result_ui->ui->analyze->show();
    for(int i=0;i<row_number;i++)
    {
        result_ui->ui->choose->addItem(QString::fromStdString(ecoli_k12[i].getGeneName()));
    }
    int star_number=0;
    if(best_score>0&&best_score<66)
    {
        star_number=1;
        result_ui->ui->star2->hide();
        result_ui->ui->star3->hide();
        result_ui->ui->star4->hide();
        result_ui->ui->star5->hide();
    }
    if(best_score>66&&best_score<72)
    {
        star_number=2;
        result_ui->ui->star5->hide();
        result_ui->ui->star4->hide();
        result_ui->ui->star3->hide();
    }
    if(best_score>72&&best_score<78)
    {
        star_number=3;
        result_ui->ui->star5->hide();
        result_ui->ui->star4->hide();
    }
    if(best_score>78&&best_score<84)
    {
        star_number=4;
        result_ui->ui->star5->hide();
    }
    if(best_score>84&&best_score<101)
    {
        star_number=5;
    }
    for(int i=0;i<star_number;i++)
    {
        star=new QFrame(result_ui->ui->star);
        star->setGeometry(i*50,0,50,50);
        star->setStyleSheet("border-image: url(:/picture/picture/star.png);");
    }
    ifstream gene_value("EachGeneChange");
    string nouse;
    double first,last;
    int choose_number=result_ui->ui->choose->currentIndex();
    for(int i=0;i<choose_number+1;i++)
        getline(gene_value,nouse);
    gene_value>>first;
    char *p=&nouse[(int)nouse.size()-7];
    nouse=p;
    last=atof(nouse.c_str());
    if(last-first>0.0001)
        result_ui->ui->up_down->setStyleSheet("border-image: url(:/picture/picture/to_high.png);");
    else if(last-first<-0.0001)
        result_ui->ui->up_down->setStyleSheet("border-image: url(:/picture/picture/to_low.png);");
    else
        result_ui->ui->up_down->setStyleSheet("border-image: url(:/picture/picture/to_eq.png);");
    QString temp_text,link;
    temp_text="<br/>:<br/>:<br/>:<br/>:<br/>Gene Name:"+QString::fromStdString(ecoli_k12[choose_number].getGeneName())
            +"<br/>RegulonDB ID:E"+QString::fromStdString(ecoli_k12[choose_number].getID())
            +"<br/>Promoter Name:"+QString::fromStdString(ecoli_k12[choose_number].getPromoterName())
            +"<br/>Left Position:"+QString::number(ecoli_k12[choose_number].getLeftPosition())
            +"<br/>Right Position:"+QString::number(ecoli_k12[choose_number].getRightPosition());
    link=" http://regulondb.ccg.unam.mx/gene?term=E"
            +QString::fromStdString(ecoli_k12[choose_number].getID())
            +"&organism=ECK12&format=jsp&type=gene";
    QString html;
    html="<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\"><html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">p, li { white-space: pre-wrap; }</style></head><body style=\" font-family:'Courier'; font-size:16pt; font-weight:72; font-style:normal;\"><p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><br />"
            +temp_text+"</p>"
            +"<a href=\""+link+"\" target=\"_blank\"><span style=\" text-decoration: underline; color:#ffffff;\">Click here for Info on RegulonDB</a>"
            +"</body></html>";
    result_ui->ui->Information->setHtml(html);
    }
    else
    {
        result_ui->ui->analyze->hide();
    }
    if(if_predict==1)
    {
        result_ui->ui->predict->show();
        QString temp;
        for(int i=0;i<n+1;i++)
        {
            if(low_high[i]==0)
                temp=temp+QString::fromStdString(ecoli_k12[pick[i]->currentIndex()].getGeneName())+"\tLow\n";
            else
                temp=temp+QString::fromStdString(ecoli_k12[pick[i]->currentIndex()].getGeneName())+"\tHigh\n";
        }
        result_ui->ui->target_gene->setText(temp);
        result_ui->ui->prediction->setText(text);
    }
    else
    {
        result_ui->ui->predict->hide();
    }
}

void console::done_clicked()
{
    ifstream gene_value("EachGeneChange");
    string nouse;
    double first,last;
    int choose_number=result_ui->ui->choose->currentIndex();
    for(int i=0;i<choose_number+1;i++)
        getline(gene_value,nouse);
    gene_value>>first;
    char *p=&nouse[(int)nouse.size()-7];
    nouse=p;
    last=atof(nouse.c_str());
    if(last-first>0.001)
        result_ui->ui->up_down->setStyleSheet("border-image: url(:/picture/picture/to_high.png);");
    else if(last-first<-0.001)
        result_ui->ui->up_down->setStyleSheet("border-image: url(:/picture/picture/to_low.png);");
    else
        result_ui->ui->up_down->setStyleSheet("border-image: url(:/picture/picture/to_eq.png);");
    QString temp_text,link;
    temp_text="<br/>:<br/>:<br/>:<br/>:<br/>Gene Name:"+QString::fromStdString(ecoli_k12[choose_number].getGeneName())
            +"<br/>RegulonDB ID:E"+QString::fromStdString(ecoli_k12[choose_number].getID())
            +"<br/>Promoter Name:"+QString::fromStdString(ecoli_k12[choose_number].getPromoterName())
            +"<br/>Left Position:"+QString::number(ecoli_k12[choose_number].getLeftPosition())
            +"<br/>Right Position:"+QString::number(ecoli_k12[choose_number].getRightPosition());
    link="http://regulondb.ccg.unam.mx/gene?term=E"
            +QString::fromStdString(ecoli_k12[choose_number].getID())
            +"&organism=ECK12&format=jsp&type=gene";
    QString html;
    html="<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\"><html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">p, li { white-space: pre-wrap; }</style></head><body style=\" font-family:'Courier'; font-size:16pt; font-weight:72; font-style:normal;\"><p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><br />"
            +temp_text+"</p>"
            +"<a href=\""+link+"\" target=\"_blank\"><span style=\" text-decoration: underline; color:#ffffff;\">Click here for Info on RegulonDB</a>"
            +"</body></html>";
    result_ui->ui->Information->setHtml(html);
}

void console::sbol_clicked()
{
    SBOL creat_SBOL;
    int m=result_ui->ui->choose->currentIndex();
    string left=to_string(ecoli_k12[m].getLeftPosition());
    string right=to_string(ecoli_k12[m].getRightPosition());
    creat_SBOL.CreatSBOL(ecoli_k12[m].getGeneName(),
                         ecoli_k12[m].getID(),
                         left,
                         right,
                         ecoli_k12[m].getGeneDescription(),
                         ecoli_k12[m].getGeneSequence());
}

void console::on_close_clicked()
{
    this->close();
}
