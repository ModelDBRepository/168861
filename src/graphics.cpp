//
//
// File author(s):  <Julian Andres Garcia Grajales>, (C) 2014
//
// Copyright: this software is licenced under the terms stipulated in the license.txt file located in its root directory
//
//

/*!\file graphics.cpp
  \brief This file contains the functions related to the visualization of the electrical part of Neurite by means of OpenGL. It is related to the "g" option of Neurite
*/

#include "graphics.h"
#include "configuration.h"
   #include <sys/stat.h>
   #include <sys/types.h>


#define WIN_X 400
#define WIN_Y 400
#define WIN_PX_H 800
#define WIN_PX_W 800
#define MAX_V 0.08
#define MIN_V -0.01
#define UNIT_V 0.01

/* global variables */
int img=0, i=0;
extern int numStepsX, numStepsT, open_GL_steps;
extern TYPE xspan, tspan, unit_x;
extern TYPE *neuriteE, **neuriteR, ***neurite,* open_GL_graphE, **open_GL_graphR,*** open_GL_graph;

void runGL()
{
	char winName[128] = "Discrete Cable Equation with Hodgkin-Huxley Channels";
	
	/* Initiate Glut */
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
	glutInitWindowSize(WIN_PX_W, WIN_PX_H);
	glutInitWindowPosition(100, 100);
	glutCreateWindow(winName);

	// Clear the background of our window to white
	glClearColor(1,1,1,0);

	// Define viewing area
	gluOrtho2D(-5,WIN_X+5,-5,WIN_Y+5);

	//Glut drawing functions
	glutDisplayFunc(display);
	glutIdleFunc(display);

	//The main loop
	glutMainLoop();

}

void display()
{
	double color[][3] = {{0,0.5,0}, {0.7,0,0.7}, {0,0,0.7}, {0.7,0,0}};
	int sleep = 33333;

	//Clear the color buffer 
	glClear(GL_COLOR_BUFFER_BIT);
	glColor3f(0,0,0); //Set color to black

	// create graph-space and put the names at the axis
	renderAxes("Distance from Voltage Source (meters)", "VOLTS");

	// render curves
	renderCurves(color);

	// render neurite
	renderNeurites(color, 0, 1, "Control", "Stretched");

	// Flush the OpenGL buffers to the window
	glFlush();
	glutPostRedisplay();

	// capture the first set of frames
	if(img < open_GL_steps)
	{
		captureScreen();
	}

	// properly iterate i
	if(i<open_GL_steps-2) {
		i++;
	} else {
		i=0;
	}

	sleep = (int)(pow(10,23)*tspan/(9*open_GL_steps));

	wait(sleep);
}

extern void wait(int microseconds)
{
  clock_t endwait;
  endwait = clock () + microseconds * CLOCKS_PER_SEC / 1000000;
  while (clock() < endwait) {}
}

extern void renderAxes(char *xlabel, char *ylabel)
{
	char tmp[8];
	int j;

	// render graph axes
	glBegin(GL_LINE_STRIP);
		glVertex2f(5*0.025*WIN_X, 39*0.025*WIN_Y);
		glVertex2f(5*0.025*WIN_X, 12*0.025*WIN_Y);
		glVertex2f(35*0.025*WIN_X, 12*0.025*WIN_Y);
	glEnd();

	// render y-axis tick marks and labels
	for(j=-1; j<(MAX_V/UNIT_V+1); j++)
	{
		// tick mark
		glBegin(GL_LINES);
			glVertex2f(5*0.025*WIN_X, (27*0.025*WIN_Y/(MAX_V - MIN_V))*UNIT_V*j + 15*0.025*WIN_Y);
			glVertex2f(4.5*0.025*WIN_X, (27*0.025*WIN_Y/(MAX_V - MIN_V))*UNIT_V*j +15*0.025*WIN_Y);
		glEnd();

		// label
		if(j < 0)
		{
			snprintf(tmp, 7, "%f", UNIT_V*j-0.065);

		       	renderString(2*0.025*WIN_X, (27*0.025*WIN_Y/(MAX_V - MIN_V))*UNIT_V*j + 15*0.025*WIN_Y, tmp, 'h');
		} else {
			snprintf(tmp, 7, "%f", UNIT_V*j-0.065);
			//fprintf(stdout,"yo=%f\n",UNIT_V*j-0.065);
			 // if(j==1)exit(0);
			renderString(2.4*0.025*WIN_X, (27*0.025*WIN_Y/(MAX_V - MIN_V))*UNIT_V*j + 15*0.025*WIN_Y, tmp, 'h');
		}
	}

	// render x-axis tick marks and labels
	for(j=0; j<(xspan/unit_x+1); j++)
	{
		// tick mark
		glBegin(GL_LINES);
			glVertex2f((30*0.025*WIN_X/xspan)*unit_x*j + 5*0.025*WIN_X, 12*0.025*WIN_Y);
			glVertex2f((30*0.025*WIN_X/xspan)*unit_x*j + 5*0.025*WIN_X, 11.5*0.025*WIN_Y);
		glEnd();

		// label
		snprintf(tmp, 4, "%f", unit_x*1E3*j);
		renderString((30*0.025*WIN_X/xspan)*unit_x*j + 4.5*0.025*WIN_X, 10.5*0.025*WIN_Y, tmp, 'h');
	}
		strcpy(tmp, "E-3");
//	snprintf(tmp, 8, "%.1e", xspan);  //  ****FIX****
	renderString(35.5*0.025*WIN_X, 10.5*0.025*WIN_Y, tmp, 'h');

	// render x and y axis labels
	renderString(0.6*0.025*WIN_X, 29.0*0.025*WIN_Y, ylabel, 'v');  // y-axis
	renderString(14*0.025*WIN_X, 9*0.025*WIN_Y, xlabel, 'h');  // x-axis
}

extern void renderString(double x, double y, char *string, char orientation)
{
	// orientation == (h)orizontal or (v)ertical

	int len, i;
	if(orientation == 'h')
		glRasterPos2f(x, y);
	len = (int) strlen(string);
	for (i = 0; i < len; i++)
	{
		if(orientation == 'v')
			glRasterPos2f(x, y-i*5);
		glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, string[i]);
	}
}

extern void renderCurves(double (*color)[3])
{
	int j, curv;
	int numCurves=1;

	glLineWidth(2);
	for(curv=0; curv<numCurves; curv++)
	{
		if (curv < 4)
			glColor3f(color[curv][0],color[curv][1],color[curv][2]);
		else
			glColor3f((rand()%10)/10,(rand()%10)/10,(rand()%10)/10);

		glBegin(GL_LINE_STRIP);
			for(j=1; j<numStepsX-1; j++)
			{

			  glVertex2f(j*(30*0.025*WIN_X/numStepsX) + 5*0.025*WIN_X, (open_GL_graph[curv][i][j]+0.065)*(27*0.025*WIN_Y/(MAX_V - MIN_V)) + 15*0.025*WIN_Y);
			}

		glEnd();
	}
	glLineWidth(1);
}


extern void renderNeurites(double (*color)[3], int curvOne, int curvTwo, char *oneName, char *twoName)
{
	int c=0, j, curv;
	double Ncolor;
	int numCurves=1;

	for(curv=0; curv<numCurves; curv++)
	{
		if (curv < 4)
			glColor3f(color[curv][0],color[curv][1],color[curv][2]);
		else
			glColor3f((rand()%10)/10,(rand()%10)/10,(rand()%10)/10);

		if((curv == curvOne) || (curv == curvTwo))
		{
			if(curv == curvOne)
				renderString(0.7*0.025*WIN_X, (-c*2 + 3.2)*0.025*WIN_Y, oneName, 'h');
			if(curv == curvTwo)
				renderString(0.7*0.025*WIN_X, (-c*2 + 3.2)*0.025*WIN_Y, twoName, 'h');
			renderString(0.7*0.025*WIN_X, (-c*2 + 2.4)*0.025*WIN_Y, "Neurite", 'h');

			for(j=0; j<numStepsX; j++)
			{
				Ncolor = open_GL_graph[curv][i][j]*200; // should be 12.5, but boosting for effect
				glColor3f(Ncolor,0,1-Ncolor);

				glBegin(GL_QUADS);
					glVertex2f(j*(0.75*WIN_X/numStepsX) + 0.125*WIN_X, (-c*2 + 2.5)*0.025*WIN_Y);
					glVertex2f(j*(0.75*WIN_X/numStepsX) + 0.125*WIN_X, (-c*2 + 3.5)*0.025*WIN_Y);
					glVertex2f((j+1)*(0.75*WIN_X/numStepsX) + 0.125*WIN_X, (-c*2 + 3.5)*0.025*WIN_Y);
					glVertex2f((j+1)*(0.75*WIN_X/numStepsX) + 0.125*WIN_X, (-c*2 + 2.5)*0.025*WIN_Y);
				glEnd();
			}

			c++;
		}
	}
}

extern void captureScreen()
{
	char fname[250];
	char tmp[8];

	char *execution_folder = getenv("working_directory");

	if (execution_folder != NULL){

	  sprintf(tmp, "%03d", img);
	  strcpy (fname, execution_folder);
	  strcat (fname, "/outputs/GUI_captures/capture_");
	  strcat (fname, tmp);
	  strcat (fname, ".bmp");
	  
	  RgbImage theImage(WIN_PX_H, WIN_PX_W);
	  theImage.LoadFromOpenglBuffer();
	  theImage.WriteBmpFile(fname);
	  img++;
	}else{
	  fprintf(stdout,"ERROR. You did not define the folder to save the pictures for the video\n");
	  exit(0);
	}
	if(img==open_GL_steps)
	  exit(0);

}
