#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>


#include "ising_cpu.h"

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


IsingModel pS;


static void keyFunc(unsigned char key, int x, int y ) { 
  switch (key) { 
    case 'j':
        pS -> J -= 0.01;
        printf("J = %lf\n", pS -> J);
        break;    
    case 'J':
        pS -> J += 0.01;
        printf("J = %lf\n", pS -> J);
		break;
	case 'h':
		pS -> H -= 0.01;
        printf("H = %lf\n", pS -> H);
        break;
    case 'H':
		pS -> H += 0.01;
        printf("H = %lf\n", pS -> H);
        break;        
    default:
        break;
  } 
} 

// Mouse event handlers
int mouseOldX, mouseOldY;
int mouseButtons = 0;

float rotateX = 0.0;
float rotateY = 0.0;

float traslateZ = -3.0;

int sqX, sqY;
int sign = 1;

void mouse(int button, int state, int x, int y){

	if (state == GLUT_DOWN){
		mouseButtons |= 1 << button;
	}
	else if (state == GLUT_UP){
		mouseButtons = 0;
	}

	mouseOldX = x;
	mouseOldY = y;

	sqX=(x>0 && x<L)?(x):(sqX); sqY=(y>0 && y<L)?(y):(sqY);
	printf("(cliks) ---> (x,y)=(%d, %d)\n", sqX, sqY);
}


void motion(int x, int y){
    
    float dx, dy;
    dx = (float)(x - mouseOldX);
    dy = (float)(y - mouseOldY);

    if (mouseButtons & 1){
        
        rotateX += dy * 0.2f;
        rotateY += dx * 0.2f;
    } else if (mouseButtons & 4){
		traslateZ += dy * 0.01f;
	}

    mouseOldX = x;
    mouseOldY = y;

    sqX=(x>0 && x<L)?(x):(sqX); sqY=(y>0 && y<L)?(y):(sqY);
    addSquare(sqX ,L - sqY, 10, sign, pS);
}


void display(){

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glLoadIdentity();

	GLfloat xSize = (right - left) / systemWidth;
	GLfloat ySize = (top - bottom) / systemHeight;

	glBegin(GL_QUADS);
		for(GLint x = 0; x < systemWidth; ++x) {
			for(GLint y = 0; y < systemHeight; ++y){

				(pS -> spins[getId(x,y)] == UP)?glColor3f(BLACK):glColor3f(RED);

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

	stepMC(pS);

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

	glutKeyboardFunc(keyFunc);
	glutMouseFunc(mouse);
 	glutMotionFunc(motion);

	pS = makeSystem(0.48, 1000, 0.03);
		
	update(0);
	glutMainLoop();

	return 0;
}