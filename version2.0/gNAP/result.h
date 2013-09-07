#ifndef RESULT_H
#define RESULT_H

#include <QWidget>
#include <QMouseEvent>

namespace Ui {
class result;
}

class result : public QWidget
{
    Q_OBJECT

public:
    explicit result(QWidget *parent = 0,Qt::WindowFlags flags=0);
    ~result();

    Ui::result *ui;

private slots:

    void on_close_clicked();

private:
    QPoint m_DragPosition;
    void mousePressEvent(QMouseEvent *event);
    void mouseMoveEvent(QMouseEvent *event);
};

#endif // RESULT_H
