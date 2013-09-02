#ifndef CONSOLE_H
#define CONSOLE_H

#include <QWidget>
#include <QFileDialog>
#include <QMouseEvent>
#include <QApplication>
#include <QLineEdit>
#include <QComboBox>
#include <QFrame>
#include <QString>
#include <stdlib.h>
#include <assert.h>
#include <sstream>

#include "Code/define.h"
#include "Code/GeneIM.h"
#include "Code/GRN.h"
#include "Code/Sequence.h"
#include "Code/RandSeq.h"
#include "Code/ModleNetwork.h"
#include "Code/PSO.h"
#include "Code/SBOL.h"

#include "result.h"

namespace Ui {
class console;
}

class console : public QWidget
{
    Q_OBJECT

public:
    Ui::console *ui;
    explicit console(QWidget *parent = 0,Qt::WindowFlags flags=0);
    ~console();
    int analyze_flag;

private slots:

    void getAmount(int row,int column);

    void getSequence(string promoter,string gene);

    void on_analyze_clicked();

    void on_add_clicked();

    void set_low();

    void set_high();

    void getInformation(GeneIM ecoli,int position);

    void on_predict_clicked();

    void on_close_clicked();

    void show_out();

    void on_result_clicked();

    void done_clicked();

    void sbol_clicked();

private:
    int row_number,column_number;
    string new_promoter_sequence;
    string new_gene_sequence;
    QPoint m_DragPosition;
    void mousePressEvent(QMouseEvent *event);
    void mouseMoveEvent(QMouseEvent *event);
    Sequence sequence_array[GENEAM];
    QComboBox *pick[7];
    QPushButton *choose_low[7];
    QPushButton *choose_high[7];
    QFrame *hold_low[7];
    QFrame *hold_high[7];
    int n;
    string gene_name[GENEAM];
    int low_high[7];
    int flag;
    GeneIM ecoli_k12[GENEAM];
    result *result_ui;
    double best_score,most_score;
    QFrame *star;
    int if_analyze,if_predict;
    QString text;
};

#endif // CONSOLE_H
