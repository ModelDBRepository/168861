//
//
// File author(s):  <Julian Andres Garcia Grajales>, (C) 2014
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//

#include "GL/glut.h"
#include "RgbImage.h"
#include "neurite.h"
#include "time.h"

#ifndef _graphics_H_
#define _graphics_H_

/* local functions */
void runGL();
void display();
void captureScreen();
void wait(int microseconds);
void renderString(double x, double y, char *string, char orientation);
void renderAxes(char *xlabel, char *ylabel);
void renderCurves(double (*color)[3]);
void renderMaxInds(char *magName, double (*color)[3], int curvOne, int curvTwo);
void renderNeurites(double (*color)[3], int curvOne, int curvTwo, char *oneName, char *twoName);

#endif
