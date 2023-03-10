// ABCSketch1.cpp : Defines the entry point for the console application.
//

#include <iostream>
#include <stdlib.h>

#include <time.h>
#include <conio.h>
#include <math.h>
#include <io.h>
#include <fcntl.h>
#include <cstdio>
#include <fstream>
#include <crtdbg.h>
#include <chrono>

#include "SketchManager.h"
#include "HapticsEventManager.h"
#include "opencv2\opencv.hpp"

#define GL_WIN_SIZE_X 256	
#define GL_WIN_SIZE_Y 256

using namespace std;
using namespace midl;

int drillStatus;

// Event Manager
static HapticsEventManagerGeneric *hapticsEvent = 0;
static ABCSketchManager *skEvent = 0;

float eye[] = { 0.0,0.0,5.0 };
float center[] = { 0.0,0.0,0.0 };
float head[] = { 0.0, 1.0, 0.0 };

//---- Open GL View
PerspectiveView view;

//---- Open GL Lighting
Light light1, light2;
GLfloat lightPos[] = { 0.0, 0.0, 10.0, 1.0 };
GLfloat diffuse1[] = { 0.54, 0.0, 0.0, 1.0 };
GLfloat diffuse2[] = { 0.0, 0.8, 0.0, 0.5 };
GLfloat ambient[] = { 0.5, 0.5, 0.5, 1.0 };
GLfloat specular[] = { 1.0, 1.0, 1.0, 1.0 };

float scale[] = { 1.0 ,1.0, 1.0 };

void CleanupAndExit()
{
	if (hapticsEvent)
	{
		hapticsEvent->Cleanup();
		HapticsEventManagerGeneric::Delete(hapticsEvent);
	}
	exit(0);
}

void glutDisplay(void)
{
	hapticsEvent->UpdateState();

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glMatrixMode(GL_MODELVIEW);
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_NORMALIZE);

	view.Bind();
	
	glPolygonMode(GL_FRONT, GL_FILL);
	glEnable(GL_POLYGON_OFFSET_FILL);
	glPolygonOffset(2.0f, 2.0f);		
	
	skEvent->Update();

	glDisable(GL_POLYGON_OFFSET_LINE);

	view.Unbind();

	glutSwapBuffers();
}

void glutReshape(int w, int h)
{
	view.Reshape(w, h);
}

void glutIdle(void)
{
	// We just set the backgournd to White!
	// and do nothing.
	glClearColor(1.0, 1.0, 1.0, 1.0);
	glutPostRedisplay();
}

void mouseClickCallback(int button, int state, int x, int y)
{
	// Example usage:
	if (button == GLUT_LEFT_BUTTON)
	{
		if (state == GLUT_DOWN)
		{
			// When the left button is pressed
			// down, do something
		}
		else if (state == GLUT_UP)
		{
			// When the left button is lifted
			// up, do something else
			//skEvent->SaveLog();
		}
	}
}

void mouseActiveMotionCallback(int x, int y)
{
	// You can define anythging here
	// Right now, we do not want to 
	// do anything when the mouse moves!
	// So, let us leave this function empty!
}

void mousePassiveMotionCallback(int x, int y)
{
}

void glutKeyboard(unsigned char key, int x, int y)
{
	switch (key)
	{
	case 27: // 27 is ASCII code for the "Esc" key
		CleanupAndExit();
		break;

	case 's':
		skEvent->SaveLog();
		break;

	case 'S':
		skEvent->SaveLog();
		break;	
	case 'r':
		skEvent->RecordData();
		drillStatus = -1;
		skEvent->SetDrillingStatus(drillStatus);
		break;

	case 'R':
		skEvent->RecordData();
		drillStatus = -1;
		skEvent->SetDrillingStatus(drillStatus);
		break;

	case 32:
		drillStatus++;
		if (drillStatus == 0) skEvent->SetDrillingStatus(drillStatus);
		else if (drillStatus == 1) skEvent->SetDrillingStatus(drillStatus);
		else if (drillStatus == 2) skEvent->SetDrillingStatus(drillStatus);
		else if (drillStatus == 3) skEvent->SetDrillingStatus(drillStatus);
		else if (drillStatus == 4) skEvent->SetDrillingStatus(drillStatus);
		else if (drillStatus == 5) skEvent->SaveLog();
		else if (drillStatus > 5) { cerr << "Data already saved" << endl; }
		break;

	case 'x':
		drillStatus = -1;
		skEvent->Reset();
		break;

	case 'X':
		drillStatus = -1;
		skEvent->Reset();
		break;
	}
}

void specialKeyCallback(unsigned char key, int x, int y)
{
	switch (key)
	{
	case GLUT_KEY_LEFT:
		break;
	case GLUT_KEY_RIGHT:
		break;
	}
}

void glInit(int * pargc, char ** argv)
{
	glutInit(pargc, argv);
	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH | GLUT_MULTISAMPLE);
	
	//glBlendFunc(GL_ONE_MINUS_DST_ALPHA, GL_DST_ALPHA);
	// Create a window of size GL_WIN_SIZE_X, GL_WIN_SIZE_Y
	// See where GL_WIN_SIZE_X and GL_WIN_SIZE_Y are defined!!
	glutInitWindowSize(GL_WIN_SIZE_X, GL_WIN_SIZE_Y);

	// Create a window with some name
	glutCreateWindow("THIS WINDOWS NEEDS TO BE ACTIVE");

	// GLEW is another library that
	// will be used as a standard for
	// initializing Shaders. 
	// Do not worry about this now!
	GLenum err = glewInit();
	if (GLEW_OK != err) // If GLEW was not properly initialized, we create a pop-up window and tell the user!
	{
		printf("Error: %s\n", glewGetErrorString(err));
		MessageBox(NULL, L"GLEW unable to initialize.", L"ERROR", MB_OK);
		return;
	}

	//---- Initialize Haptics device
	hapticsEvent = HapticsEventManager::Initialize();
	hapticsEvent->Setup(skEvent);
	
	//---- Initialize the openGL camera 
	//---- and projection parameters
	view.SetParameters(35, 0.1, 10000.0);
	view.SetCameraCenter(0.0, 0.0, 0.0);
	view.SetCameraEye(0.0, 0.2, 2.0);
	view.SetCameraHead(0.0, 1.0, 0.0);

	// Shader Initialization:
	skEvent->InitShaders();

	light1.SetPosition(lightPos);
	light1.SetDiffuseColor(diffuse1);
	light1.SetAmbientColor(ambient);
	light1.SetSpecularColor(specular);

	light2.SetPosition(lightPos);
	light2.SetDiffuseColor(diffuse2);
	light2.SetAmbientColor(ambient);
	light2.SetSpecularColor(specular);

	glutReshapeFunc(glutReshape);
	glutDisplayFunc(glutDisplay);
	glutKeyboardFunc(glutKeyboard);
	glutMouseFunc(mouseClickCallback);
	//glutPassiveMotionFunc(mousePassiveMotionCallback);
	glutIdleFunc(glutIdle);

}

int main(int argc, char** argv)
{
	skEvent = new ABCSketchManager();


	cerr << "Please enter the User ID ?" << endl;
	cin >> skEvent->userID;

	cerr << "Please enter the bone variant ? (OB/YB)" << endl;
	cin >> skEvent->boneType;

	//cerr << "Please enter the trial number ?" << endl;
	//cin >> skEvent->trial;	
	
	atexit(CleanupAndExit);

	glInit(&argc, argv);
	glutMainLoop();

	return 0;
}

