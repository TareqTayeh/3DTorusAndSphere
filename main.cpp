/* Name: Tareq Tayeh
 * Student Number: 250725776
 * Course: CS3388A Computer Graphics
 * Assignment Number: Assignment #2
 * Date: 17 October 2017
 * Program purpose: The purpose of this program is to create and render 3D wiremesh objects defined with parametric functions.
 *                  The program will display spheres and tori with the use of Bresenham's algorithm for 2D line segments,
 *                  the synthetic camera, and the parametric functions for the sphere and the torus.
 * File: main.cpp
 */

#include "dialog.h"
#include <QApplication>

/* Author: Tareq Tayeh
 * Date of creation: 12 October 2017
 *
 * Function: main()
 * Purpose: Instantiates an instance of Dialog, and sets its size, background and shows it.
 */
int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    Dialog w;

    w.setFixedSize(512, 512); //Opening a window of size 512 by 512 pixels
    w.setStyleSheet("QWidget { background-color: white; }"); //Setting the window background color to white
    w.show(); //Show dialog

    return a.exec();
}
