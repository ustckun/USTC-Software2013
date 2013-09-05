#ifndef WIDGET_H
#define WIDGET_H

#include <QWidget>
#include <QFileDialog>
#include <QMouseEvent>
#include <QApplication>

#include "Code/define.h"
#include "Code/GeneIM.h"
#include "Code/GetReady.h"

#include "console.h"


namespace Ui {
class Widget;
}

class Widget : public QWidget
{
    Q_OBJECT

public:
    explicit Widget(QWidget *parent = 0,Qt::WindowFlags flags=0);
    ~Widget();
    Ui::Widget *ui;

signals:
    void pass_amount(int row,int column);

    void pass_sequence(string promoter,string gene);

    void pass_information(GeneIM ecoli,int position);

private slots:
    void on_TF_TF_clicked();

    void on_TF_Gene_clicked();

    void on_Gene_Info_clicked();

    void on_TU_Info_clicked();

    void on_Promoter_Info_clicked();

    void on_close_clicked();

    void on_Take_a_nap_clicked();

    void on_Promoter_done_clicked();

    void on_Gene_done_clicked();

    void back_func();

private:
    GeneIM ecoli[GENEAM];
    string TF_TF_address;
    string TF_Gene_address;
    string gene_info_address;
    string TU_info_address;
    string promoter_info_address;
    QPoint m_DragPosition;
    void mousePressEvent(QMouseEvent *event);
    void mouseMoveEvent(QMouseEvent *event);
    string promoter_sequence;
    string gene_sequence;
    int flag;
    console *console_ui;
};

#endif // WIDGET_H
