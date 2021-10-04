#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>


#include "glauber_system.c"

#define RED 1.0, 0.0, 0.0
#define BLACK 0.0, 0.0, 0.0

GLint windowWidth  = 1000;
GLint windowHeight = 1000;

GLint FPS = 24;
GLfloat left = 0.0;
GLfloat right = 1.0;
GLfloat bottom = 0.0;
GLfloat top = 1.0;
GLint systemWidth = 500;
GLint systemHeight = 500;


GlauberModel pG;

void display(){

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glLoadIdentity();

	GLfloat xSize = (right - left) / systemWidth;
	GLfloat ySize = (top - bottom) / systemHeight;

	glBegin(GL_QUADS);
		for(GLint x = 0; x < systemWidth; ++x) {
			for(GLint y = 0; y < systemHeight; ++y){

				(pG -> spins[getId(x,y)] == UP)?glColor3f(BLACK):glColor3f(RED);

				glVertex2f(    x*xSize+left,    y*ySize+bottom);
				glVertex2f((x+1)*xSize+left,    y*ySize+bottom);
				glVertex2f((x+1)*xSize+left,(y+1)*ySize+bottom);
				glVertex2f(    x*xSize+left,(y+1)*ySize+bottom);
			}
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

	stepMC(pG);
	printf("t=%lf\tH=%lf\n", pG->step - 0.2*N_SWEEPS, pG -> H *cos((2 * M_PI * (pG -> step - 0.2*N_SWEEPS)) / pG -> T)) ;
	glutPostRedisplay();
	glutTimerFunc(1000 / FPS, update, 0);
}


int main(int argc, char **argv){

	glutInit(&argc, argv);
	glutInitWindowSize(windowWidth, windowHeight);
	glutInitWindowPosition(0, 0);
	glutCreateWindow("Aplication: Glauber Model");
	glClearColor(1, 1, 1, 1);
	
	glutReshapeFunc(reshape);
	glutDisplayFunc(display);

	pG = makeSystem(0.48, 1000, 0.03);
		
	update(0);
	glutMainLoop();

	return 0;
}