#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>

#include "Lennar_Jones.h"

#define RED 1.0, 0.0, 0.0
#define BLACK 0.0, 0.0, 0.0

GLint windowWidth  = 1000;
GLint windowHeight = 1000;

GLint FPS = 24;
GLfloat left = 0.0;
GLfloat right = 1.0;
GLfloat bottom = 0.0;
GLfloat top = 1.0;
GLfloat systemWidth = 27.386127875258;
GLfloat systemHeight = 27.386127875258;


world *pW;

void display(){

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glLoadIdentity();

	GLfloat xSize = (right - left) / systemWidth;
	GLfloat ySize = (top - bottom) / systemHeight;

	glPointSize(10.0f);
	glBegin(GL_POINTS);
		int i; for(i = 0; i < N; i++) {

			GLfloat x = pW -> X[i][0];
			GLfloat y = pW -> X[i][1];
		
			glColor3f(RED);
			glVertex2f(x * xSize, y * ySize);

			
		}
	glEnd();
	
	glFlush();
	glutSwapBuffers();	
}


void reshape(int w, int h){

	windowWidth = w;
	windowHeight = h;

	glViewport(0, 0, windowWidth, windowHeight);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(left, right, bottom, top);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	glutPostRedisplay();
}


void update(int value){

	stepMC(pW);

	glutPostRedisplay();
	glutTimerFunc(1000 / FPS, update, 0);
}


int main(int argc, char **argv){

	glutInit(&argc, argv);
	glutInitWindowSize(windowWidth, windowHeight);
	glutInitWindowPosition(0, 0);
	glutCreateWindow("Aplication: Lennar-Jones (OFF-Lattice MC)");
	glClearColor(1, 1, 1, 1);
	
	glutReshapeFunc(reshape);
	glutDisplayFunc(display);

	pW = makeSystem();
		
	update(0);
	glutMainLoop();

	return 0;
}