#ifndef PLOTTER_H
#define PLOTTER_H

#include <QWidget>
#include <QFile>
#include <vector>
#include <iostream>
#include <QPainter>
#include <QFileDialog>

using namespace std;

class plotter : public QWidget
{
    Q_OBJECT

public:
    plotter(QWidget *parent = 0);
    ~plotter();
    vector<double> step, rmsd, energy;
    enum { Margin = 50 };

private:
    void paintEvent(QPaintEvent* event);
};

#endif // PLOTTER_H
