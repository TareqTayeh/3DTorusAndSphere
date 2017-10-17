/* Name: Tareq Tayeh
 * Student Number: 250725776
 * Course: CS3388A Computer Graphics
 * Assignment Number: Assignment #2
 * Date: 17 October 2017
 * Program purpose: The purpose of this program is to create and render 3D wiremesh objects defined with parametric functions.
 *                  The program will display spheres and tori with the use of Bresenham's algorithm for 2D line segments,
 *                  the synthetic camera, and the parametric functions for the sphere and the torus.
 * File: dialog.cpp
 */

#include "dialog.h"
#include "ui_dialog.h"
#include "matrix.h"
#include <vector>

#define Ex 15.0
#define Ey 15.0
#define Ez 15.0

#define Gx 0.0
#define Gy 0.0
#define Gz 0.0

#define UPx 0.0
#define UPy 0.0
#define UPz 1.0

#define NP 5.0
#define FP 50.0

#define THETA 90.0
#define ASPECT 1.0

#define W  512
#define H  512

Dialog::Dialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::Dialog)
{
    ui->setupUi(this);
}

Dialog::~Dialog()
{
    delete ui;
}

/* Author: Dr. Steven S. Beauchemin
 * Date of creation: 12 October 2017
 *
 * Function:build_camera_matrix()
 * Parameter:
 *  In: dmatrix_t *E - Camera's "Eye" matrix
 *      dmatrix_t *G - Camera's "Gaze" matrix
 *  Out: None
 * Returns: dmatrix_t - A matrix of the computed built camera
 * Desc: Computes the camera matrix
 */
dmatrix_t *build_camera_matrix(dmatrix_t *E, dmatrix_t *G) {

    dmatrix_t N ; /* Viewing axis */

    N = *dmat_normalize(dmat_sub(E,G)) ;
    N.l = 3 ;

    dmatrix_t UP ;
    dmat_alloc(&UP,4,1) ;
    UP.l = 3 ;

    UP.m[1][1] = UPx ;
    UP.m[2][1] = UPy ;
    UP.m[3][1] = UPz ;
    UP.m[4][1] = 1.0 ;

    dmatrix_t U ;

    U = *dmat_normalize(dcross_product(&UP,&N)) ;

    dmatrix_t V ;
    V = *dcross_product(&N,&U) ;

    dmatrix_t Mv ; /* Build matrix M_v */
    dmat_alloc(&Mv,4,4) ;

    Mv.m[1][1] = U.m[1][1] ;
    Mv.m[1][2] = U.m[2][1] ;
    Mv.m[1][3] = U.m[3][1] ;
    Mv.m[1][4] = -1.0*((*E).m[1][1]*U.m[1][1] + (*E).m[2][1]*U.m[2][1] + (*E).m[3][1]*U.m[3][1]) ;

    Mv.m[2][1] = V.m[1][1] ;
    Mv.m[2][2] = V.m[2][1] ;
    Mv.m[2][3] = V.m[3][1] ;
    Mv.m[2][4] = -1.0*((*E).m[1][1]*V.m[1][1] + (*E).m[2][1]*V.m[2][1] + (*E).m[3][1]*V.m[3][1]) ;

    Mv.m[3][1] = N.m[1][1] ;
    Mv.m[3][2] = N.m[2][1] ;
    Mv.m[3][3] = N.m[3][1] ;
    Mv.m[3][4] = -1.0*((*E).m[1][1]*N.m[1][1] + (*E).m[2][1]*N.m[2][1] + (*E).m[3][1]*N.m[3][1]) ;

    Mv.m[4][1] = 0.0 ;
    Mv.m[4][2] = 0.0 ;
    Mv.m[4][3] = 0.0 ;
    Mv.m[4][4] = 1.0 ;

    dmatrix_t Mp ; /* Build matrix Mp */
    dmat_alloc(&Mp,4,4) ;
    Mp = *dmat_identity(&Mp) ;

    float a = -1.0*(FP + NP)/(FP - NP) ;
    float b = -2.0*(FP*NP)/(FP - NP) ;

    Mp.m[1][1] = NP ;
    Mp.m[2][2] = NP ;
    Mp.m[3][3] = a ;
    Mp.m[3][4] = b ;
    Mp.m[4][3] = -1.0 ;
    Mp.m[4][4] = 0.0 ;

    /* Build matrices T_1 and S_1 */

    /* Work out coordinates of near plane corners */

    float top = NP*tan(M_PI/180.0*THETA/2.0) ;
    float right = ASPECT*top ;
    float bottom = -top ;
    float left = -right ;

    dmatrix_t T1 ;
    dmat_alloc(&T1,4,4) ;

    T1 = *dmat_identity(&T1) ;
    T1.m[1][4] = -(right + left)/2.0 ;
    T1.m[2][4] = -(top + bottom)/2.0 ;

    dmatrix_t S1 ;
    dmat_alloc(&S1,4,4) ;

    S1 = *dmat_identity(&S1) ;
    S1.m[1][1] = 2.0/(right - left) ;
    S1.m[2][2] = 2.0/(top - bottom) ;

    /* Build matrices T2, S2, and W2 */

    dmatrix_t T2 ;
    dmatrix_t S2 ;
    dmatrix_t W2 ;

    dmat_alloc(&T2,4,4) ;
    dmat_alloc(&S2,4,4) ;
    dmat_alloc(&W2,4,4) ;

    T2 = *dmat_identity(&T2) ;
    S2 = *dmat_identity(&S2) ;
    W2 = *dmat_identity(&W2) ;

    T2.m[1][4] = 1.0 ;
    T2.m[2][4] = 1.0 ;

    S2.m[1][1] = W/2.0 ;
    S2.m[2][2] = H/2.0 ;

    W2.m[2][2] = -1.0 ;
    W2.m[2][4] = (double)H ;

    return dmat_mult(&W2,dmat_mult(&S2,dmat_mult(&T2,dmat_mult(&S1,dmat_mult(&T1,dmat_mult(&Mp,&Mv)))))) ;
}

/* Author: Dr. Steven S. Beauchemin
 * Date of creation: 12 October 2017
 *
 * Function:perspective_projection()
 * Parameter:
 *  In: dmatrix_t *P - A matrix pointer
 *  Out: None
 * Returns: dmatrix_t - A matrix of the computed perspective projection
 * Desc: Computes the perspective projection of a matrix and returns it
 */
dmatrix_t *perspective_projection(dmatrix_t *P) {

    (*P).m[1][1] /= (*P).m[4][1] ;
    (*P).m[2][1] /= (*P).m[4][1] ;
    (*P).m[3][1] /= (*P).m[4][1] ;
    (*P).m[4][1] /= (*P).m[4][1] ;

    return P ;
}

/* Author: Dr. Steven S. Beauchemin
 * Date of creation: 12 October 2017
 *
 * Function: camera_setup()
 * Parameter:
 *  In: None
 *  Out: None
 * Returns: dmatrix_t - A matrix of the camera
 * Desc: Sets up the camera matrix by using its E (eye) and G (gaze), and calling the build_camera_matrix() function
 */
dmatrix_t camera_setup(){
    dmatrix_t E ; /* The centre of projection for the camera */

    dmat_alloc(&E,4,1) ;

    E.m[1][1] = Ex ;
    E.m[2][1] = Ey ;
    E.m[3][1] = Ez ;
    E.m[4][1] = 1.0 ;

    dmatrix_t G ; /* Point gazed at by camera */

    dmat_alloc(&G,4,1) ;

    G.m[1][1] = Gx ;
    G.m[2][1] = Gy ;
    G.m[3][1] = Gz ;
    G.m[4][1] = 1.0 ;

    dmatrix_t C ; /* The camera matrix */

    dmat_alloc(&C,4,4) ;
    C = *build_camera_matrix(&E,&G) ;

    return C;
}

/* Author: Tareq Tayeh
 * Date of creation: 16 October 2017
 *
 * Function: DrawTorus()
 * Parameter:
 *  In: None
 *  Out: None
 * Returns: std::vector<dmatrix_t> - A matrix with all the torus points from the parametric equation
 * Desc: Algorithm to draw a torus. Computes points from the parametric equations and stores them in a matrix
 */
std::vector<dmatrix_t> DrawTorus(){
    //Variables and Initializations for the Torus' Algorithm
    double theta = 0, phi = 0; //Angles
    float DTOR = 0.01745329252; //(2 * PI / 36)
    dmatrix_t P1,P2,P3,P4; //Points
    std::vector<dmatrix_t> torusPVector;

    //Algorithm to Draw a Torus
    int du = 10, dv = 10; //u and v increments
    double outerRadius = 220.0, innerRadius = 55.0; //Tori's Radii

    for (int u = 0;u < 360;u += du) {
        for (int v = 0;v < 360;v += dv) {

            //Segment 1
            theta = (u) * DTOR;
            phi   = (v) * DTOR;
            dmat_alloc(&P1,4,1) ;
            P1.m[1][1] = 256 + cos(theta) * ( outerRadius + innerRadius * cos(phi) ); //x1
            P1.m[2][1] = 256 + sin(theta) * ( outerRadius + innerRadius * cos(phi) ); //y1
            P1.m[3][1] = 256 + innerRadius * sin(phi); //z1
            P1.m[4][1] = 1.0 ;
            torusPVector.push_back(P1);

            //Segment 2
            theta = (u+du) * DTOR;
            phi   = (v) * DTOR;
            dmat_alloc(&P2,4,1) ;
            P2.m[1][1] = 256 + cos(theta) * ( outerRadius + innerRadius * cos(phi) ); //x2
            P2.m[2][1] = 256 + sin(theta) * ( outerRadius + innerRadius * cos(phi) ); //y2
            P2.m[3][1] = 256 + innerRadius * sin(phi); //z2
            P2.m[4][1] = 1.0 ;
            torusPVector.push_back(P2);

            //Segment 3
            theta = (u+du) * DTOR;
            phi   = (v+dv) * DTOR;
            dmat_alloc(&P3,4,1) ;
            P3.m[1][1] = 256 + cos(theta) * ( outerRadius + innerRadius * cos(phi) ); //x3
            P3.m[2][1] = 256 + sin(theta) * ( outerRadius + innerRadius * cos(phi) ); //y3
            P3.m[3][1] = 256 + innerRadius * sin(phi); //z3
            P3.m[4][1] = 1.0 ;
            torusPVector.push_back(P3);

            //Segment 4
            theta = (u) * DTOR;
            phi   = (v+dv) * DTOR;
            dmat_alloc(&P4,4,1) ;
            P4.m[1][1] = 256 + cos(theta) * ( outerRadius + innerRadius * cos(phi) ); //x4
            P4.m[2][1] = 256 + sin(theta) * ( outerRadius + innerRadius * cos(phi) ); //y4
            P4.m[3][1] = 256 + innerRadius * sin(phi); //z4
            P4.m[4][1] = 1.0 ;
            torusPVector.push_back(P4);
        }
    }

    return torusPVector;
}

/* Author: Tareq Tayeh
 * Date of creation: 16 October 2017
 *
 * Function: DrawSphere()
 * Parameter:
 *  In: None
 *  Out: None
 * Returns: std::vector<dmatrix_t> - A matrix with all the sphere points from the parametric equation
 * Desc: Algorithm to draw a sphere. Computes points from the parametric equations and stores them in a matrix
 */
std::vector<dmatrix_t> DrawSphere(){
    //Variables and Initializations for the Sphere's Algorithms
    double theta = 0, phi = 0; //Angles
    float DTOR = 0.01745329252; //(2 * PI / 36)
    dmatrix_t P1,P2,P3,P4; //Points
    std::vector<dmatrix_t> spherePVector;

    //Algorithm to Draw a Sphere
    int di = 10, dj = 10; //i and j increments
    int sphereRadius = 100; //Sphere Radius

    for (int i = 0; i < 360; i += di){
        for (int j = 0; j < 180; j += dj){

            //Segment 1
            theta = i * DTOR;
            phi = j * DTOR;
            dmat_alloc(&P1,4,1) ;
            P1.m[1][1] = 256 + cos(theta) * sin(phi) * sphereRadius; //x1
            P1.m[2][1] = 256 + sin(theta) * sin(phi) * sphereRadius; //y1
            P1.m[3][1] = 256 + cos(phi) * sphereRadius; //z1
            P1.m[4][1] = 1.0 ;
            spherePVector.push_back(P1);

            //Segment 2
            theta = (i + di) * DTOR;
            phi = j * DTOR;
            dmat_alloc(&P2,4,1) ;
            P2.m[1][1] = 256 + cos(theta) * sin(phi) * sphereRadius; //x2
            P2.m[2][1] = 256 + sin(theta) * sin(phi) * sphereRadius; //y2
            P2.m[3][1] = 256 + cos(phi) * sphereRadius; //z2
            P2.m[4][1] = 1.0 ;
            spherePVector.push_back(P2);

            //Segment 3
            theta = (i + di) * DTOR;
            phi = (j + dj) * DTOR;
            dmat_alloc(&P3,4,1) ;
            P3.m[1][1] = 256 + cos(theta) * sin(phi) * sphereRadius; //x3
            P3.m[2][1] = 256 + sin(theta) * sin(phi) * sphereRadius; //y3
            P3.m[3][1] = 256 + cos(phi) * sphereRadius; //z3
            P3.m[4][1] = 1.0 ;
            spherePVector.push_back(P3);

            //Segment 4
            theta = i * DTOR;
            phi = (j + dj) * DTOR;
            dmat_alloc(&P4,4,1) ;
            P4.m[1][1] = 256 + cos(theta) * sin(phi) * sphereRadius; //x4
            P4.m[2][1] = 256 + sin(theta) * sin(phi) * sphereRadius; //y4
            P4.m[3][1] =256 + cos(phi) * sphereRadius; //z4
            P4.m[4][1] = 1.0 ;
            spherePVector.push_back(P4);
        }
    }

    return spherePVector;
}

/* Author: Tareq Tayeh
 * Date of creation: 12 October 2017
 *
 * Function: keyPressEvent()
 * Parameter:
 *  In: QKeyEvent *event - Any key event when the dialog is open
 *  Out: None
 * Returns: None
 * Desc: Quits the program when the key "q" is pressed
 */
void Dialog::keyPressEvent(QKeyEvent *event)
{
    switch(event->key())
    {
    case Qt::Key_Q:
        close();
        break;
    default:
        QDialog::keyPressEvent(event);
    }
}

/* Author: Tareq Tayeh
 * Date of creation: 16 October 2017
 *
 * Function: paintEvent()
 * Parameter:
 *  In: QPaintEvent - A paint event when the dialog is open
 *  Out: None
 * Returns: None
 * Desc: Dialog's painter function where the shape algorithms are called, perspective projections are calculated,
 *       then the Bresenham function is called where the dots are "painted"
 */
void Dialog::paintEvent(QPaintEvent *e)
{
    //Camera's Matrix and Setup
    printf("Camera Matrix:\n") ;
    dmatrix_t C = camera_setup();
    write_dmatrix(&C) ;

    //Variables and Initializations for the Shapes
    dmatrix_t *PP1,*PP2; //Perspective Projections
    std::vector<dmatrix_t> torusPVector, spherePVector; //Shape's points vector

    //Draws the Torus
    torusPVector = DrawTorus();
    for(int a = 0; a < torusPVector.size() - 1; a++){
        PP1 = perspective_projection(dmat_mult(&C,&(torusPVector.at(a))));
        PP2 = perspective_projection(dmat_mult(&C,&(torusPVector.at(a+1))));
        Bresenham(PP1->m[1][1],PP1->m[2][1],PP2->m[1][1],PP2->m[2][1]);
    }

    //Draws the Sphere
    spherePVector = DrawSphere();
    for(int a = 0; a < spherePVector.size() - 1; a++){
        PP1 = perspective_projection(dmat_mult(&C,&(spherePVector.at(a))));
        PP2 = perspective_projection(dmat_mult(&C,&(spherePVector.at(a+1))));
        Bresenham(PP1->m[1][1],PP1->m[2][1],PP2->m[1][1],PP2->m[2][1]);
    }

    //Freeing Memory
    delete (PP1->m);
    delete (PP2->m);
    delete PP1;
    delete PP2;
}

/* Author: Tareq Tayeh
 * Date of creation: 21 September 2017
 *
 * Function: Bresenham()
 * Parameter:
 *  In: int x1 - coordinate x1 as an integer
 *      int y1 - coordinate y1 as an integer
 *      int x2 - coordinate x2 as an integer
 *      int y2 - coordinate y2 as an integer
 *  Out: None
 * Returns: None
 * Desc: Bresenham's line algorithm, and where the dots are "painted" on the dialog
 */
void Dialog::Bresenham(int x1, int y1, int x2, int y2)
{
  QPainter painter(this);
  int dx, dy, error, ystep, y, maxX;

  const bool steep = (fabs(y2 - y1) > fabs(x2 - x1));
  if(steep)
  {
    std::swap(x1, y1);
    std::swap(x2, y2);
  }

  if(x1 > x2)
  {
    std::swap(x1, x2);
    std::swap(y1, y2);
  }

  dx = x2 - x1;
  dy = fabs(y2 - y1);

  error = dx / 2.0f;
  ystep = (y1 < y2) ? 1 : -1;
  y = y1;

  maxX = x2;

  for(int x=x1; x<maxX; x++)
  {
    if(steep)
    {
        painter.drawPoint(y,x);
    }
    else
    {
        painter.drawPoint(x,y);
    }

    error -= dy;
    if(error < 0)
    {
        y += ystep;
        error += dx;
    }
  }
}
