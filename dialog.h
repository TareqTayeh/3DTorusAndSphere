/* Name: Tareq Tayeh
 * Student Number: 250725776
 * Course: CS3388A Computer Graphics
 * Assignment Number: Assignment #2
 * Date: 17 October 2017
 * Program purpose: The purpose of this program is to create and render 3D wiremesh objects defined with parametric functions.
 *                  The program will display spheres and tori with the use of Bresenham's algorithm for 2D line segments,
 *                  the synthetic camera, and the parametric functions for the sphere and the torus.
 * File: dialog.h
 */

#ifndef DIALOG_H
#define DIALOG_H

#include <QDialog>
#include <QKeyEvent>
#include <QtGui>
#include <QtCore>
#include <iostream>
#include <QVector3D>

namespace Ui {
class Dialog;
}

class Dialog : public QDialog
{
    Q_OBJECT

public:
    explicit Dialog(QWidget *parent = 0);
    ~Dialog();
    void keyPressEvent(QKeyEvent *event);

protected:
    void paintEvent(QPaintEvent *e);
    void Bresenham( int x1, int y1, int x2, int y2);

private:
    Ui::Dialog *ui;
};

#endif // DIALOG_H
