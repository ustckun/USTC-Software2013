#include "result.h"
#include "ui_result.h"

result::result(QWidget *parent,Qt::WindowFlags flags) :
    QWidget(parent,flags),
    ui(new Ui::result)
{
    ui->setupUi(this);
    setWindowFlags(Qt::FramelessWindowHint);
    setAttribute(Qt::WA_TranslucentBackground, true);
    setWindowTitle("gNAP");
    ui->star1->show();
    ui->star2->show();
    ui->star3->show();
    ui->star4->show();
    ui->star5->show();
}

result::~result()
{
    delete ui;
}

void result::mousePressEvent(QMouseEvent *event)
{
    if (event->button() == Qt::LeftButton) {
        m_DragPosition = event->globalPos() - this->pos();
        event->accept();
    }
}

void result::mouseMoveEvent(QMouseEvent *event)
{
    if (event->buttons() && Qt::LeftButton) {
        move(event->globalPos() - m_DragPosition);
        event->accept();
    }
}

void result::on_close_clicked()
{
    this->close();
}
