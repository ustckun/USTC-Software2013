#include "widget.h"
#include <QApplication>

int main(int argc, char *argv[])
{
    GeneIM gene_information[GENEAM];
    GetReady prepare;
    QApplication a(argc, argv);
    Widget w;
    w.show();
    return a.exec();
}
