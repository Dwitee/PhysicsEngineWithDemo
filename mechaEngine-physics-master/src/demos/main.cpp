/*
 * The main entry point for all demos.
 *
 * Part of the mechaEngine physics system.
 *
 * Copyright (c) Dwitee Krishna Panda. All Rights Reserved.
 *
 * This software is distributed under licence. Use of this software
 * implies agreement with all terms and conditions of the accompanying
 * software licence.
 */

#include <gl/glut.h>

// Include the general application structure.
#include "app.h"

// Include the timing functions
#include "timing.h"

// Forward declaration of the function that will return the
// application object for this particular demo. This should be
// implemented in the demo's .cpp file.
extern Application* getApplication();

const int windowWidth = 1280;
const int windowHeight = 640;
// Store the global application object.
//The Game project cpp should return this pointer when getApplication is called.
Application* app;


//Called each frame to update the 3D scene. Delegates to
//the application.
void update()
{
    // Update the timing.
    TimingData::get().update();
    
    // Delegate to the application.
    app->update();
}

//Creates a window in which to display the scene.
void createWindow(const char* title)
{
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize(windowWidth,windowHeight);
    glutInitWindowPosition(0,0);
    glutCreateWindow(title);
}


//Called when a mouse button is pressed. Delegates to the
//application.
void mouse(int button, int state, int x, int y)
{
    app->mouse(button, state, x, y);
}


void display()
{
    app->display();

    // Update the displayed content.
    glFlush();
    glutSwapBuffers();
}

//Called when the mouse is dragged.
void motion(int x, int y)
{
    app->mouseDrag(x, y);
}

//Called when the display window changes size.
void reshape(int width, int height)
{
    app->resize(width, height);
}

//Called when a key is pressed.
void keyboard(unsigned char key, int x, int y)
{
    // Note we omit passing on the x and y: they are rarely needed.
    app->key(key);
}

//The main entry point. We pass arguments onto GLUT.

int main(int argc, char** argv)
{
    // Sets up GLUT and the timers
    glutInit(&argc, argv);
    TimingData::init();

    // Creates the application and its window
    app = getApplication();
    createWindow(app->getTitle());

    // Sets up the appropriate handler functions
    glutReshapeFunc(reshape);
    glutKeyboardFunc(keyboard);
    glutDisplayFunc(display);
    glutIdleFunc(update);
    glutMouseFunc(mouse);
    glutMotionFunc(motion);

    // Runs the application
    app->initGraphics();
    glutMainLoop();

    // Clean up the application
    app->deinit();
    delete app;
    TimingData::deinit();
}
