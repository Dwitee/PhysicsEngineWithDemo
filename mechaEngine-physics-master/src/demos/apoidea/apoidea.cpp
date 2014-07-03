/*
 * The honey bee behaviour.
 */

#include <gl/glut.h>
#include <mechaEngine/mechaEngine.h>
#include "../app.h"
#include "../timing.h"

#include <stdio.h>
#include <cassert>

#define BEE_COUNT 5
#define PLATFORM_COUNT 100
#define BEE_RADIUS 0.4f

static const mechaEngine::real cameraSpeedX = 0.5;
static bool isGameOver = false;
static bool isGameStarted = false;

/**
 * The Platforms are actually two dimensional: lines on which the
 * particles can rest. Platforms are also contact Creators for the physics collision
 */
class Platform : public mechaEngine::ParticleContactGenerator
{
public:
    mechaEngine::Vector3 startPoint;
    mechaEngine::Vector3 endPoint;

    /**
     * Holds a pointer to the particles we're checking for collisions with. 
     */
    mechaEngine::Particle *particles;

    virtual unsigned addContact(
        mechaEngine::ParticleContact *contact, 
        unsigned limit
        ) const;
};

unsigned Platform::addContact(mechaEngine::ParticleContact *ParticleContact,
                              unsigned limit) const
{
    const static mechaEngine::real restitutionCoef = 0.0f;

    unsigned used = 0;
    for (unsigned i = 0; i < BEE_COUNT; i++)
    {
        if (used >= limit) break;
        
        // Check for collision detection for line and circle
        mechaEngine::Vector3 toParticleVec = particles[i].getPosition() - startPoint;
        mechaEngine::Vector3 lineDirection = endPoint - startPoint;
        mechaEngine::real projected = toParticleVec * lineDirection;
        mechaEngine::real platformSqLength = lineDirection.squareMagnitude();
        if (projected <= 0)
        {
            // The bee is nearest to the start point
            if (toParticleVec.squareMagnitude() < BEE_RADIUS*BEE_RADIUS)
            {
                // collision is detected then
                ParticleContact->contactNormal = toParticleVec.unit();
                ParticleContact->contactNormal.z = 0;
                ParticleContact->restitution = restitutionCoef;
                ParticleContact->particle[0] = particles + i;
                ParticleContact->particle[1] = 0;
                ParticleContact->penetration = BEE_RADIUS - toParticleVec.magnitude();
                used ++;
                ParticleContact ++;
            }
            
        }
        else if (projected >= platformSqLength)
        {
            // The bee is nearest to the end point
            toParticleVec = particles[i].getPosition() - endPoint;
            if (toParticleVec.squareMagnitude() < BEE_RADIUS*BEE_RADIUS)
            {
                // collision is detected then
                ParticleContact->contactNormal = toParticleVec.unit();
                ParticleContact->contactNormal.z = 0;
                ParticleContact->restitution = restitutionCoef;
                ParticleContact->particle[0] = particles + i;
                ParticleContact->particle[1] = 0;
                ParticleContact->penetration = BEE_RADIUS - toParticleVec.magnitude();
                used ++;            
                ParticleContact ++;
            }
        }
        else
        {
            // the bee is nearest to the middle.
            mechaEngine::real distanceToPlatform = 
                toParticleVec.squareMagnitude() -
                projected*projected / platformSqLength;
            if (distanceToPlatform < BEE_RADIUS*BEE_RADIUS)
            {
                // collision is detected then
                mechaEngine::Vector3 closestPoint = 
                    startPoint + lineDirection*(projected/platformSqLength);

                ParticleContact->contactNormal = (particles[i].getPosition()-closestPoint).unit();
                ParticleContact->contactNormal.z = 0;
                ParticleContact->restitution = restitutionCoef;
                ParticleContact->particle[0] = particles + i;
                ParticleContact->particle[1] = 0;
                ParticleContact->penetration = BEE_RADIUS - real_sqrt(distanceToPlatform);
                used ++;
                ParticleContact ++;
            }
        }
    }
    return used;
}

/**
 * A force creator for near field attraction.
 */
class BeeForceGenerator : public mechaEngine::ParticleForceGenerator
{
public:

    mechaEngine::Particle *particles;


    mechaEngine::real maxReplusion;


    mechaEngine::real maxAttraction;
    
 
    mechaEngine::real minNaturalDistance;
    
    mechaEngine::real maxNaturalDistance;

    mechaEngine::real floatHead;

    unsigned maxFloat;


    mechaEngine::real maxDistance;

    virtual void updateForce(
        mechaEngine::Particle *particle, 
        mechaEngine::real duration
        );
};

void BeeForceGenerator::updateForce(mechaEngine::Particle *beeParticle,
                                      mechaEngine::real duration)
{
    unsigned joinCount = 0;
    // checking for each of the bee which are particles
    for (unsigned int i = 0; i < BEE_COUNT; i++)
    {
        // should not  attract it self
        if (particles + i == beeParticle) continue;

        // Work out the separation distance between queen bee and other bee
        mechaEngine::Vector3 separationDistance =
            particles[i].getPosition() - beeParticle->getPosition();
        
        //make both the bee in same plane same distacne from camera to make simple 2D
        separationDistance.z = 0.0f;
        
        // get the 2D distance
        mechaEngine::real distance = separationDistance.magnitude();

        // the following code give the behaviour between the two bees
        if (distance < minNaturalDistance)
        {
            // Use a repulsion force.
            distance = 1.0f - distance / minNaturalDistance;
            beeParticle->addForce(
                separationDistance.unit() * (1.0f - distance) * maxReplusion * -1.0f
                );
            joinCount++;
        }
        else if (distance > maxNaturalDistance && distance < maxDistance)
        {
            // Use an attraction force.
            distance = 
                (distance - maxNaturalDistance) / 
                (maxDistance - maxNaturalDistance);
            beeParticle->addForce(
                separationDistance.unit() * distance * maxAttraction
                );
            joinCount++;
        }
    }

    // If the particle is the head, and we've got a join count, then float it.
    if (beeParticle == particles && joinCount > 0 && maxFloat > 0)
    {
        mechaEngine::real force = mechaEngine::real(joinCount / maxFloat) * floatHead;
        if (force > floatHead) force = floatHead;
        beeParticle->addForce(mechaEngine::Vector3(0, force, 0));
    }

}

/**
 * The main demo class definition.
 */
class BeeDemo : public Application
{
    mechaEngine::Particle *bees;

    Platform *platforms;

    mechaEngine::ParticleWorld world;

    BeeForceGenerator beeForceGenerator;
    
    mechaEngine::Vector3 cameraPosition;

    /* The control for the x-axis. */
    float xAxis;

    /* The control for the y-axis. */
    float yAxis;

    void reset();

public:
    /** Creates a new demo object. */
    BeeDemo();
    virtual ~BeeDemo();

    /** Returns the window title for the demo. */
    virtual const char* getTitle();

    /** Display the particles. */
    virtual void display();

    /** Update the particle positions. */
    virtual void update();

    /** Handle a key press. */
    virtual void key(unsigned char key);

};

// Method definitions
BeeDemo::BeeDemo()
:
xAxis(0), yAxis(0),
world(PLATFORM_COUNT+BEE_COUNT, PLATFORM_COUNT)
{
    // Create the bee storage
    bees = new mechaEngine::Particle[BEE_COUNT];
    mechaEngine::Random r;

    // Create the force generator
    beeForceGenerator.minNaturalDistance = BEE_RADIUS* 1.25f;
    beeForceGenerator.maxNaturalDistance = BEE_RADIUS* 2.5f;
    beeForceGenerator.maxDistance = BEE_RADIUS * 10.0f;
    beeForceGenerator.maxFloat = 9;
    beeForceGenerator.floatHead = 8.0f;
    beeForceGenerator.particles = bees;
    beeForceGenerator.maxAttraction = 10.0f;
    beeForceGenerator.maxReplusion = 2.0f;

    
    // Create the platforms
    mechaEngine::real endPointX = 0;
    mechaEngine::real roofY = 4;
    mechaEngine::real floorY = 0;
    
    platforms = new Platform[PLATFORM_COUNT];
    //create the floor
    for (unsigned i = 0; i < PLATFORM_COUNT/4; i++)
    {
        mechaEngine::real  platformLength = 10;
        platforms[i].startPoint.x = 0.5f+ endPointX;
        platforms[i].startPoint.y = floorY;
        platforms[i].endPoint.x = platforms[i].startPoint.x + platformLength;
        platforms[i].endPoint.y = platforms[i].startPoint.y;
        
        endPointX = platforms[i].endPoint.x;


        // Make sure the platform knows which particles it 
        // should collide with.
        platforms[i].particles = bees;
        world.getContactGenerators().push_back(platforms + i);
    }
    endPointX = 0;
    //creat the roof
    for ( unsigned i = PLATFORM_COUNT/4 ; i < PLATFORM_COUNT/2; i++)
    {
        
        mechaEngine::real  platformLength = 10;
        platforms[i].startPoint.x = 0.5f+endPointX;
        platforms[i].startPoint.y = roofY;
        platforms[i].endPoint.x = platforms[i].startPoint.x + platformLength;
        platforms[i].endPoint.y = platforms[i].startPoint.y;
        endPointX = platforms[i].endPoint.x;
        platforms[i].particles = bees;
        world.getContactGenerators().push_back(platforms + i);
    }
    
    //creats floor obstacles
    mechaEngine::real offsetX = 0;
    for ( unsigned i = PLATFORM_COUNT/2 ; i < PLATFORM_COUNT*3/4; i++)
    {
        
        mechaEngine::real  platformLength = 2 ;
        platforms[i].startPoint.x = offsetX;
        platforms[i].startPoint.y = floorY;
        platforms[i].endPoint.x = platforms[i].startPoint.x ;
        platforms[i].endPoint.y = platforms[i].startPoint.y + platformLength;
        offsetX += r.randomReal(5,10);
        
        platforms[i].particles = bees;
        world.getContactGenerators().push_back(platforms + i);
    }
    //creat roof obstacles
    offsetX = 2;
    for ( unsigned i = PLATFORM_COUNT*3/4 ; i < PLATFORM_COUNT; i++)
    {
        
        mechaEngine::real  platformLength = 2 ;
        platforms[i].startPoint.x = offsetX;
        platforms[i].startPoint.y = roofY;
        platforms[i].endPoint.x = platforms[i].startPoint.x ;
        platforms[i].endPoint.y = platforms[i].startPoint.y - platformLength;
        offsetX += r.randomReal(5,15);
        
        platforms[i].particles = bees;
        world.getContactGenerators().push_back(platforms + i);
    }
    
    // Create the bees.
    Platform *p = platforms;
    mechaEngine::real fraction = (mechaEngine::real)1.0 / BEE_COUNT;
    mechaEngine::Vector3 delta = p->endPoint - p->startPoint;
    for (unsigned i = 0; i < BEE_COUNT; i++)
    {
        unsigned me = (i+BEE_COUNT/2) % BEE_COUNT;
        bees[i].setPosition(
            p->startPoint + delta * (mechaEngine::real(me)*0.8f*fraction+0.1f) +
            mechaEngine::Vector3(0, 1.0f+r.randomReal(), 0));

        bees[i].setVelocity(0,0,0);
        bees[i].setDamping(0.2f);
        bees[i].setAcceleration(mechaEngine::Vector3::GRAVITY * 0.4f);
        if( i == 0)//for the queen bee
            bees[i].setMass(0.750f);
        else
            bees[i].setMass(0.10f);
        bees[i].clearAccumulator();

        world.getParticles().push_back(bees + i);
        world.getForceRegistry().add(bees + i, &beeForceGenerator);
    }
    bees[0].getPosition(&cameraPosition);
}

void BeeDemo::reset()
{
    mechaEngine::Random r;
    Platform *p = platforms + (PLATFORM_COUNT-2);
    mechaEngine::real fraction = (mechaEngine::real)1.0 / BEE_COUNT;
    mechaEngine::Vector3 delta = p->endPoint - p->startPoint;
    for (unsigned i = 0; i < BEE_COUNT; i++)
    {
        unsigned me = (i+BEE_COUNT/2) % BEE_COUNT;
        bees[i].setPosition(
            p->startPoint + delta * (mechaEngine::real(me)*0.8f*fraction+0.1f) +
            mechaEngine::Vector3(0, 1.0f+r.randomReal(), 0));
        bees[i].setVelocity(0,0,0);
        bees[i].clearAccumulator();
    }
    isGameStarted = false;
}

BeeDemo::~BeeDemo()
{
    delete bees;
}

void BeeDemo::display()
{
    mechaEngine::Vector3 pos = cameraPosition ;

    if (abs (bees[0].getPosition().x - cameraPosition.x) > 7) // if less then viewport width
    {
        isGameOver = true;
        glColor3f(0,0,0);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glClearColor(0.5, 0.5, 0.5,0.5);
        

    }
    // Clear the view port and set the camera direction
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glLoadIdentity();
    gluLookAt(pos.x, pos.y, 6.0,  pos.x, pos.y, 0.0,  0.0, 1.0, 0.0);

    glColor3f(0,0,0);
    if (isGameOver == true) {
        renderText(width/2, height/2, "Game Over: Keep the queen bee moving by pressing space ");
    }
 

    glBegin(GL_LINES);
    glColor3f(0,0,1);
    for (unsigned i = 0; i < PLATFORM_COUNT/2; i++)
    {
        const mechaEngine::Vector3 &p0 = platforms[i].startPoint;
        const mechaEngine::Vector3 &p1 = platforms[i].endPoint;

        glVertex3f(p0.x, p0.y, p0.z);
        glVertex3f(p1.x, p1.y, p1.z);

    }
    glEnd();
    glBegin(GL_LINES);
    for (unsigned i =  PLATFORM_COUNT/2; i < PLATFORM_COUNT; i++)
    {
        const mechaEngine::Vector3 &p0 = platforms[i].startPoint;
        const mechaEngine::Vector3 &p1 = platforms[i].endPoint;
        glVertex3f(p0.x, p0.y, p0.z);
        glVertex3f(p1.x, p1.y, p1.z);
    }
    glEnd();

    glColor3f(1,0,0);
    
    for (unsigned i = 0; i < BEE_COUNT; i++)
    {

        const mechaEngine::Vector3 &p = bees[i].getPosition();
        glPushMatrix();
        glTranslatef(p.x, p.y, p.z);
        glutSolidSphere(BEE_RADIUS, 12, 12);
        glPopMatrix();
        glColor3f(0,0,0);
    }
    glColor3f(1,0,0);
    
    mechaEngine::Vector3 p = bees[0].getPosition();
    mechaEngine::Vector3 v = bees[0].getVelocity() * 0.05f;
    v.trim(BEE_RADIUS*0.5f);
    p = p + v;
    glPushMatrix();
    glTranslatef(p.x-BEE_RADIUS*0.2f, p.y, BEE_RADIUS);
    glColor3f(1,1,1);
    glutSolidSphere(BEE_RADIUS*0.2f, 8, 8);
    glTranslatef(0,0,BEE_RADIUS*0.2f);
    glColor3f(0,0,0);
    glutSolidSphere(BEE_RADIUS*0.1f, 8, 8);
    glTranslatef(BEE_RADIUS*0.4f, 0, -BEE_RADIUS*0.2f);
    glColor3f(1,1,1);
    glutSolidSphere(BEE_RADIUS*0.2f, 8, 8);
    glTranslatef(0,0,BEE_RADIUS*0.2f);
    glColor3f(0,0,0);
    glutSolidSphere(BEE_RADIUS*0.1f, 8, 8);
    glPopMatrix();
    
    
    if (isGameStarted == false) {
        renderText(width/2, height/2, " Press D to move forward\n, tap space bar to fly, \nPress any key to start,  \n,");
    }
}

void BeeDemo::update()
{
    
   if ( isGameStarted == false)
   {
       return;
   }
    // Clear accumulators
    world.startFrame();

    // Find the duration of the last frame in seconds
    float duration = (float)TimingData::get().lastFrameDuration * 0.001f;
    if (duration <= 0.0f) return;
    
    
    // Recenter the axes
    xAxis *= pow(0.1f, duration);
    yAxis *= pow(0.1f, duration);

    //move cmera at a constant speed
    if (isGameOver == false && isGameStarted == true) {
        cameraPosition.x += cameraSpeedX*duration;
    }
    
    mechaEngine::Vector3 queenBeePosition;
    bees[0].getPosition(&queenBeePosition);
    // if bee is on center camera make bee speed and camera speed
    if (queenBeePosition.x > cameraPosition.x ) {
        cameraPosition.x = queenBeePosition.x;
    }

    //move the floor obstacle up and down
    for ( unsigned i = PLATFORM_COUNT/2 ; i < PLATFORM_COUNT*3/4; i++)
    {
        
        platforms[i].endPoint.y = 2*cos(cameraPosition.x) > 0? 2*cos(cameraPosition.x) : -2*cos(cameraPosition.x)  ;


    }
    
    //move the roof obstacle left and right
    for ( unsigned i = PLATFORM_COUNT*3/4 ; i < PLATFORM_COUNT; i++)
    {
        
        platforms[i].endPoint.x += 0.05*cos(cameraPosition.x);
        platforms[i].startPoint.x += 0.05* cos(cameraPosition.x);
        
        
    }
    
    // Move the queen bee constantly forward
    bees[0].addForce(mechaEngine::Vector3(1, 0.0, 0)*5.0f);
    bees[0].addForce(mechaEngine::Vector3(xAxis, yAxis*2.0, 0)*5.0f);

    // Run the simulation
    world.runPhysics(duration);

    // Bring all the particles back to 2d
    mechaEngine::Vector3 position;
    for (unsigned i = 0; i < BEE_COUNT; i++)
    {
        bees[i].getPosition(&position);
        position.z = 0.0f;
        bees[i].setPosition(position);
    }

    
    
    Application::update();
}

const char* BeeDemo::getTitle()
{
    return "mechaEngine > Bee Demo";
}

void BeeDemo::key(unsigned char key)
{
    switch(key)
    {
    case ' ':
        yAxis = 1.0;
        break;
    case 's': case 'S':
        yAxis = -1.0;
        break;
    case 'a': case 'A':
        xAxis = -1.0f;
        break;
    case 'd': case 'D':
        xAxis = 1.0f;
        break;
    case 'r': case 'R':
        reset();
            break;
    default:
            isGameStarted = true;
        break;
    }
}

/**
 * Called by the common demo framework to create an application
 * object (with new) and return a pointer.
 */
Application* getApplication()
{
    return new BeeDemo();
}