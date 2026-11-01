/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   MarchingCube.h
 * Author: obublik
 *
 * Created on 9. ledna 2019, 10:41
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifndef MARCHINGCUBE_H
#define MARCHINGCUBE_H

typedef struct {
   double x,y,z;
} XYZ;

typedef struct {
   XYZ p[3];
} TRIANGLE;

typedef struct {
   XYZ p[8];
   double val[8];
} GRIDCELL;

int Polygonise(GRIDCELL grid,double isolevel,TRIANGLE *triangles);
XYZ VertexInterp(double isolevel, XYZ p1, XYZ p2, double valp1, double valp2);


#endif /* MARCHINGCUBE_H */

