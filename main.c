/*
* Tsunami simulation - main file with graphics
*
* We use finite elements methods with Navier-Stokes equations for the model
*
* By Arbai Faycal and Sofyan Mimouni.
*
* With the contribution and work of Vincent Legat, Université Catholique de Louvain (UCL)
*/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <string.h>

#define GLFW_INCLUDE_GLU
#include "GL/glfw.h"

int tsunamiReadFile(const char *baseResultName, int iter, double *U, double *V, double *E, int nelem);
void tsunamiAnimate(double dt, int nmax, int sub, const char *meshFileName, const char *baseResultName);
void tsunamiCompute(double dt, int nmax, int sub, const char *meshFileName, const char *baseResultName);

int main( void )
{
    int nElem,nNode,i,j,index,trash,*elem;
    double *X,*Y,*H,*E,*U,*V;
    int width,height,mouse;
    tsunamiCompute(4.0,10000,100, "PacificTriangleFine.txt","data/tsunami");

    int nout = 100;
    double t;
    double R = 6371220.0;
    double BathMax = 9368.0;
    GLfloat colors[9], coord[9];

    FILE* file = fopen("PacificTriangleFine.txt","r");
    if (file == NULL) {
    	printf("Error : cannot open mesh file :-) \n");
        exit(0); }
 	fscanf(file, "Number of nodes %d \n", &nNode);
  	X = malloc(sizeof(double)*nNode);
  	Y = malloc(sizeof(double)*nNode);
  	H = malloc(sizeof(double)*nNode);
	for (i = 0; i < nNode; i++)
    	fscanf(file,"%d : %le %le %le  \n",&trash,&X[i],&Y[i],&H[i]);
    fscanf(file, "Number of triangles %d \n", &nElem);
  	elem = malloc(sizeof(int)*3*nElem);
 	U = malloc(sizeof(double)*3*nElem);
 	V = malloc(sizeof(double)*3*nElem);
 	E = malloc(sizeof(double)*3*nElem);
  	for (i = 0; i < nElem; i++)
    	fscanf(file,"%d : %d %d %d \n", &trash,&elem[i*3],&elem[i*3+1],&elem[i*3+2]);
  	fclose(file);


   	glfwInit();
   	glfwOpenWindow(640,480,0,0,0,0,1,0,GLFW_WINDOW );
	glfwSetWindowTitle( "MECA1120 Tsunami" );

    glfwEnable( GLFW_STICKY_KEYS );
    glfwSwapInterval( 1 );

    GLfloat mat_specular[] = { 1.0, 1.0, 1.0, 0.0 };
    GLfloat mat_shininess[] = { 50.0 };
    GLfloat light_position[] = { 8.0, 8.0, 8.0, 0.0 };

	glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
    glMaterialfv(GL_FRONT, GL_SHININESS, mat_shininess);
    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
    glLightfv(GL_LIGHT0, GL_POSITION, light_position);
    GLfloat light_radiance[] = {1., 1., 1., 1.};

    glLightfv(GL_LIGHT0, GL_DIFFUSE, light_radiance);
    glLightfv(GL_LIGHT0, GL_SPECULAR, light_radiance);
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glDepthFunc(GL_LEQUAL);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_COLOR_MATERIAL);
    glEnable(GL_NORMALIZE);

    double t0 = 0;
    int frame0,frame = -1;

    do {
        t = glfwGetTime();
        frame0 = frame;
    	frame = (int) ((t-t0) * 2);

        if (frame0 != frame) {

            glfwGetMousePos( &mouse, NULL );  mouse = 389;

            char filename[256];
            const char *basename = "data/tsunami-%08d.txt";
            sprintf(filename, basename, frame * nout);
    //A decommenter lorsque vous avez des fichiers de resultats
            if (access(filename, F_OK)) {
               glfwTerminate(); exit( EXIT_SUCCESS ); }

            printf("===  Reading local file %s %d %f \n",filename,frame,t);
            tsunamiReadFile("data/tsunami",frame*nout,U,V,E,nElem);

            glfwGetWindowSize( &width, &height );
            height = height > 0 ? height : 1;
            glViewport( 0, 0, width, height );

            glClearColor( 0.9f, 0.9f, 0.8f, 0.0f );
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

            glMatrixMode(GL_PROJECTION);
            glLoadIdentity();
            gluPerspective(65.0f,(GLfloat)width/(GLfloat)height,1.0f,100.0f);

            glMatrixMode(GL_MODELVIEW);
            glLoadIdentity();
            gluLookAt(0.0f,1.0f,0.0f,0.0f, 20.0f, 0.0f,0.0f,0.0f,1.0f);
            glTranslatef(0.0f,14.0f,0.0f);
            double tt = 0;
            glRotatef(0.3f*(GLfloat)mouse + (GLfloat)tt*10.0f,0.0f,0.0f,1.0f);

            GLUquadricObj *quadratic = gluNewQuadric();
            gluQuadricNormals(quadratic, GLU_SMOOTH);
            glColor3f(1.0,1.0,1.0);
            gluSphere(quadratic,5.95,400,200);
            // A commenter pour supprimer la sphere interieure
            // Conseille pour les maillages grossiers :-)

            for (i=0; i < nElem; ++i) {
                for (j=0; j < 3; ++j) {
                    index = elem[3*i+j];
                    double value = H[index]/BathMax;
                    value = E[3*i+j]*10;
                    if (value < 0) value = 0;
                    if (value > 1) value = 1;
                    colors[j*3+0] = 3.5*(value)*(value);
                    colors[j*3+1] = (1-value)*(value)*3.5;
                    colors[j*3+2] = (1-value)*(1-value);
                    double x = X[index];
                    double y = Y[index];
                    double Factor = (4*R*R + x*x + y*y)*(R/6);
                    coord[j*3+0] = 4*R*R * x / Factor;
                    coord[j*3+1] = 4*R*R * y / Factor;
                    coord[j*3+2] = (4*R*R - x*x - y*y)*R / Factor;  }

                glEnableClientState(GL_VERTEX_ARRAY);
                glEnableClientState(GL_COLOR_ARRAY);
                glEnableClientState(GL_NORMAL_ARRAY);
                glVertexPointer(3, GL_FLOAT, 0, coord);
                glNormalPointer(GL_FLOAT, 0, coord);
                glColorPointer(3, GL_FLOAT, 0, colors);
                glDrawArrays(GL_TRIANGLES, 0, 3);
                glDisableClientState(GL_NORMAL_ARRAY);
                glDisableClientState(GL_COLOR_ARRAY);
                glDisableClientState(GL_VERTEX_ARRAY);

                glColor3f(0.0, 0.0, 0.0);
                glEnableClientState(GL_VERTEX_ARRAY);
                for (j=0; j < 9; ++j)
                     coord[j] = coord[j] * 1.001;
                glVertexPointer(3, GL_FLOAT, 0, coord);
                glDrawArrays(GL_LINE_LOOP, 0, 3);
                glDisableClientState(GL_VERTEX_ARRAY); }
            glfwSwapBuffers();

    }} while(glfwGetKey( GLFW_KEY_ESC ) != GLFW_PRESS && glfwGetWindowParam( GLFW_OPENED ) );

    glfwTerminate();
    exit( EXIT_SUCCESS );

}
