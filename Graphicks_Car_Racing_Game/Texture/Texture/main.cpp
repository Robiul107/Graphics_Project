#include<windows.h>
#include <GL/glut.h>
#include<bits/stdc++.h>
#include <stdlib.h>
#define rad (3.1416/180)
#define EN_SIZE 20
#include "BmpLoader.h"


using namespace std;


const double PI = 3.14159265389;
unsigned int ID;
int anglex=0, angley = 0, anglez = 270;
//int anglex=0, angley = 0, anglez =0;          //rotation angles
int window;
int wired=0;
int shcpt=0;
int total_point=0;
int highest_point=-1;
int life=5;
int animat = 0;
const int L=20,L1=14;
const int dgre=3;
int ncpt=L+1;
int clikd=0;
const int nt = 40;				//number of slices along x-direction
const int ntheta = 100;






float zoom=4;
int tola[5000][5000];
float tX=0,tY=0,tZ=-5,rX=0,rY=0,rZ=4;
float tZ1=-20,tZ2=-40,tZ3=-60,tZ4=-80,tZ5=-100,tZ6=-120;
float rotX=0,rotY=0,rotZ=0;
float cosX=0,cosY=1,cosZ=0;
float angle=0;
float xEye=-1.0f,yEye=10.5f,zEye=30.0f;
float cenX=0,cenY=4,cenZ=-30,roll=0;
float radius=0;
float theta=0,slope=0;
float speed =0.6;
bool rt=true;
float angleBackFrac = 0.2;
bool saheedMinarVisible = false;
float r[] = {0.1,0.4,0.0,0.9,0.2,0.5,0.0,0.7,0.5,0.0};
float g[] = {0.2,0.0,0.4,0.5,0.2,0.0,0.3,0.9,0.0,0.2};
float b[] = {0.4,0.5,0.0,0.7,0.9,0.0,0.1,0.2,0.5,0.0};
int TIME=0;
bool START = false,START1=false;
bool l1=true,l2=false,l3=true,l4=false,l5=false,df=true,am=true,sp=true;
float torusPosX[7] = {1,-2,3,-4,-2,0,2};
float torusPosY[7] = {2,3,10,6,7,4,1};

bool rot = false;



    GLfloat l_spt[] = {0,0,-1,1};
    GLfloat spt_ct[] = {10};

    GLfloat l_sptb[] = {0,-1,0,1};
    GLfloat spt_ctb[] = {50};


    GLfloat no_light[]  = { 0.0f, 0.0f, 0.0f, 1.0f };
    GLfloat light_ambient[]  = { 0.3f, 0.3f, 0.3f, 1.0f };
    GLfloat light_diffuse[]  = { 1.0f, 1.0f, 1.0f, 1.0f };
    GLfloat light_specular[] = { 1.0f, 1.0f, 1.0f, 1.0f };

    GLfloat light_position0[] = { 2.0f, 5.0f, 5.0f, 0.0f }; //sun
    GLfloat light_position1[] = { 0, 1,0, 1.0f }; //spot light
    GLfloat light_position2[] = { 0,-5,-100, 1.0f }; //normal light
    GLfloat light_position3[] = { 0,15,0, 1.0f }; //normal light


    GLfloat mat_ambient[]    = { 0.7f, 0.7f, 0.7f, 1.0f };
    GLfloat mat_diffuse[]    = { 0.8f, 0.8f, 0.8f, 1.0f };
    GLfloat mat_specular[]   = { 1.0f, 1.0f, 1.0f, 1.0f };
    GLfloat high_shininess[] = { 100.0f };




//floor
GLfloat mat_ambient2[] = { 0.5,0.5,0.5, 1.0 };
GLfloat mat_diffuse2[] = { 0.5,0.5,0.5, 1.0 };
GLfloat mat_specular2[] = { 1.0, 1.0, 1.0, 1.0 };
GLfloat mat_shininess2[] = {60};


//white
GLfloat mat_ambient3[] = { 1,1,1, 1.0 };
GLfloat mat_diffuse3[] = { 1,1,1, 1.0 };
GLfloat mat_specular3[] = { 1.0, 1.0, 1.0, 1.0 };
GLfloat mat_shininess3[] = {60};



//Green
GLfloat mat_ambient4[] = { 0,1,0, 1.0 };
GLfloat mat_diffuse4[] = { 0,1,0, 1.0 };
GLfloat mat_specular4[] = { 1.0, 1.0, 1.0, 1.0 };
GLfloat mat_shininess4[] = {60};


//Black
GLfloat mat_ambient5[] = { 0,0,0, 1.0 };
GLfloat mat_diffuse5[] = { 0,0,0, 1.0 };
GLfloat mat_specular5[] = { 1.0, 1.0, 1.0, 1.0 };
GLfloat mat_shininess5[] = {60};

//red
GLfloat mat_ambient6[] = { 1,0,0, 1.0 };
GLfloat mat_diffuse6[] = { 1,0,0, 1.0 };
GLfloat mat_specular6[] = { 1.0, 1.0, 1.0, 1.0 };
GLfloat mat_shininess6[] = {60};



static void resize(int width, int height)
{
    const float ar = (float) width / (float) height;

    glViewport(0, 0, width, height);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glFrustum(-ar, ar, -1.0, 1.0, 2.0, 1000.0);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
}
static GLfloat v_cube[8][3] =
{
    {0,0,0},
    {0,0,1},
    {0,1,0},
    {0,1,1},

    {1,0,0},
    {1,0,1},
    {1,1,0},
    {1,1,1}
};
void ownScale(double x, double y, double z)
{
    GLfloat m[16];
    for(int i=0; i<16; i++)
    {
        m[i]=0;
    }
    m[0]=x;
    m[5]=y;
    m[10]=z;
    m[15]=1;

    glMatrixMode(GL_MODELVIEW);
    glMultMatrixf(m);

}


static GLubyte c_ind[6][4] =
{
    {0,2,6,4},
    {1,5,7,3},
    {0,4,5,1},
    {2,3,7,6},
    {0,1,3,2},
    {4,6,7,5}
};

static void getNormal3p(GLfloat x1, GLfloat y1, GLfloat z1, GLfloat x2, GLfloat y2, GLfloat z2, GLfloat x3, GLfloat y3, GLfloat z3)
{
    GLfloat Ux, Uy, Uz, Vx, Vy, Vz, Nx, Ny, Nz;

    Ux = x2-x1;
    Uy = y2-y1;
    Uz = z2-z1;

    Vx = x3-x1;
    Vy = y3-y1;
    Vz = z3-z1;

    Nx = Uy*Vz - Uz*Vy;
    Ny = Uz*Vx - Ux*Vz;
    Nz = Ux*Vy - Uy*Vx;

    glNormal3f(Nx,Ny,Nz);
}

void cube(float ColR=0.5, float ColG=0.5, float ColB=0.5)
{

    GLfloat m_no[] = {0, 0, 0, 1.0};
    GLfloat m_amb[] = {ColR,ColG,ColB,1};
    GLfloat m_diff[] = {ColR,ColG,ColB,1};
    GLfloat m_spec[] = {1,1,1,1};
    GLfloat m_sh[] = {40};

    glMaterialfv(GL_FRONT, GL_AMBIENT, m_amb);
    glMaterialfv(GL_FRONT, GL_DIFFUSE, m_diff);
    glMaterialfv(GL_FRONT, GL_SPECULAR, m_spec);
    glMaterialfv(GL_FRONT, GL_SHININESS, m_sh);


    glBegin(GL_QUADS);
    for (GLint i = 0; i <6; i++)
    {
        getNormal3p(v_cube[c_ind[i][0]][0], v_cube[c_ind[i][0]][1], v_cube[c_ind[i][0]][2],
                    v_cube[c_ind[i][1]][0], v_cube[c_ind[i][1]][1], v_cube[c_ind[i][1]][2],
                    v_cube[c_ind[i][2]][0], v_cube[c_ind[i][2]][1], v_cube[c_ind[i][2]][2]);

        for (GLint j=0; j<4; j++)
        {
            glVertex3fv(&v_cube[c_ind[i][j]][0]);
        }
    }
    glEnd();
}

void cube1(double r, double g, double b)
{
    GLfloat m_no[] = {0.5, 0.5, 0.5, 1.0};
    GLfloat m_amb[] = {r,g,b,1};
    GLfloat m_diff[] = {r,g,b,1};
    GLfloat m_spec[] = {1,1,1,1};
    GLfloat m_sh[] = {100};

    glMaterialfv(GL_FRONT, GL_AMBIENT, m_amb);
    glMaterialfv(GL_FRONT, GL_DIFFUSE, m_diff);
    glMaterialfv(GL_FRONT, GL_SPECULAR, m_spec);
    glMaterialfv(GL_FRONT, GL_SHININESS, m_sh);




    glBegin(GL_QUADS);

    glVertex3f(-1, -1, 1);
    glVertex3f(1, -1, 1);
    glVertex3f(1, 1, 1);
    glVertex3f(-1, 1, 1);


    glVertex3f(1, 1, 1);
    glVertex3f(1, -1, 1);
    glVertex3f(1, -1, -1);
    glVertex3f(1, 1, -1);


    glVertex3f(-1, 1, -1);
    glVertex3f(-1, -1, -1);
    glVertex3f(-1, -1, 1);
    glVertex3f(-1, 1, 1);


    glVertex3f(-1, 1, -1);
    glVertex3f(1, 1, -1);
    glVertex3f(1, -1, -1);
    glVertex3f(-1, -1, -1);


    glVertex3f(-1, 1, -1);
    glVertex3f(-1, 1, 1);
    glVertex3f(1, 1, 1);
    glVertex3f(1, 1, -1);


    glVertex3f(-1, -1, 1);
    glVertex3f(-1, -1, -1);
    glVertex3f(1, -1, -1);
    glVertex3f(1, -1, 1);
    glEnd();
}


void Cube()
{

    glBegin(GL_QUADS);
    for (GLint i = 0; i <6; i++)
    {
        getNormal3p(v_cube[c_ind[i][0]][0], v_cube[c_ind[i][0]][1], v_cube[c_ind[i][0]][2],
                    v_cube[c_ind[i][1]][0], v_cube[c_ind[i][1]][1], v_cube[c_ind[i][1]][2],
                    v_cube[c_ind[i][2]][0], v_cube[c_ind[i][2]][1], v_cube[c_ind[i][2]][2]);

        for (GLint j=0; j<4; j++)
        {
            glVertex3fv(&v_cube[c_ind[i][0]][0]);
            glTexCoord2f(1,1);
            glVertex3fv(&v_cube[c_ind[i][1]][0]);
            glTexCoord2f(1,0);
            glVertex3fv(&v_cube[c_ind[i][2]][0]);
            glTexCoord2f(0,0);
            glVertex3fv(&v_cube[c_ind[i][3]][0]);
            glTexCoord2f(0,1);
        }
    }
    glEnd();
}

GLfloat ctrlpoints_tunnel[L1+1][3]=
{
    {73.2249, 23.684, 0},
    {73.2249, 23.684, 0},
    {-29.9166, 23.602, 0},
    {-29.9166, 23.602, 0},
    {51.3486, 23.4378, 0},
    {51.3486, 23.4378, 0},
    {58.2295, 23.2681, 0},
    {58.2295, 23.2681, 0},
    {-21.0119, 23.0602, 0},
    {-21.0119, 23.0602, 0},
    {-17.9402, 22.9124, 0},
    {-17.9402, 22.9124, 0},
    {60.5347, 22.7647, 0},
    {60.5347, 22.7647, 0},
    {48.1226, 22.6935, 0}
};




GLfloat ctrlpoints[L+1][3] =
{

    {2.15, 1.7, 0},
    {2.15, 1.7, 0},
    {2.175, 2.95, 0},
    {2.175, 2.95, 0},
    {2.2, 4.275, 0},
    {2.2, 4.275, 0},
    {0.975, 4.3, 0},
    {0.975, 4.3, 0},
    {-0.225, 4.3, 0},
    {-0.225, 4.3, 0},
    {-1.45, 4.3, 0},
    {-1.45, 4.3, 0},
    {-1.525, 2.95, 0},
    {-1.525, 2.95, 0},
    {-1.55, 1.675, 0},
    {-1.55, 1.675, 0},
    {-0.6, 1.675, 0},
    {-0.6, 1.675, 0},
    {0.85, 1.675, 0},
    {0.85, 1.675, 0},
    {2.175, 1.7, 0}
};

float wcsClkDn[3],wcsClkUp[3];
///////////////////////////////
class point1
{
public:
    point1()
    {
        x=0;
        y=0;
    }
    int x;
    int y;
} clkpt[2];
int flag=0;
GLint viewport[4]; //var to hold the viewport info
GLdouble modelview[16]; //var to hold the modelview info
GLdouble projection[16]; //var to hold the projection matrix info

//////////////////////////
void scsToWcs(float sx,float sy, float wcsv[3] );
void processMouse(int button, int state, int x, int y);
void matColor(float kdr, float kdg, float kdb,  float shiny, int frnt_Back=0, float ambFactor=1.0, float specFactor=1.0);
///////////////////////////

void scsToWcs(float sx,float sy, float wcsv[3] )
{

    GLfloat winX, winY, winZ; //variables to hold screen x,y,z coordinates
    GLdouble worldX, worldY, worldZ; //variables to hold world x,y,z coordinates

    //glGetDoublev( GL_MODELVIEW_MATRIX, modelview ); //get the modelview info
    glGetDoublev( GL_PROJECTION_MATRIX, projection ); //get the projection matrix info
    glGetIntegerv( GL_VIEWPORT, viewport ); //get the viewport info

    winX = sx;
    winY = (float)viewport[3] - (float)sy;
    winZ = 0;

    //get the world coordinates from the screen coordinates
    gluUnProject( winX, winY, winZ, modelview, projection, viewport, &worldX, &worldY, &worldZ);
    wcsv[0]=worldX;
    wcsv[1]=worldY;
    wcsv[2]=worldZ;


}
void processMouse(int button, int state, int x, int y)
{
    if(button==GLUT_LEFT_BUTTON && state==GLUT_DOWN)
    {
        if(flag!=1)
        {
            flag=1;
            clkpt[0].x=x;
            clkpt[0].y=y;
        }


        scsToWcs(clkpt[0].x,clkpt[0].y,wcsClkDn);
        //cout<<"\nD: "<<x<<" "<<y<<" wcs: "<<wcsClkDn[0]<<" "<<wcsClkDn[1];
        cout<<"{"<<wcsClkUp[0]<<", "<<wcsClkUp[1]<<", "<<0.0<<"},"<<endl;
    }
    else if(button==GLUT_LEFT_BUTTON && state==GLUT_UP)
    {
        if (flag==1)
        {
            clkpt[1].x=x;
            clkpt[1].y=y;
            flag=0;
        }
        float wcs[3];
        scsToWcs(clkpt[1].x,clkpt[1].y,wcsClkUp);
        //cout<<"\nU: "<<x<<" "<<y<<" wcs: "<<wcsClkUp[0]<<" "<<wcsClkUp[1];
        cout<<"{"<<wcsClkUp[0]<<", "<<wcsClkUp[1]<<", "<<0.0<<"},"<<endl;

        clikd=!clikd;
    }
}

//control points
long long nCr(int n, int r)
{
    if(r > n / 2)
        r = n - r; // because C(n, r) == C(n, n - r)
    long long ans = 1;
    int i;

    for(i = 1; i <= r; i++)
    {
        ans *= n - r + i;
        ans /= i;
    }

    return ans;
}

//polynomial interpretation for N points
void BezierCurve ( double t,  float xy[2])
{
    double y=0;
    double x=0;
    t=t>1.0?1.0:t;
    for(int i=0; i<=L; i++)
    {
        int ncr=nCr(L,i);
        double oneMinusTpow=pow(1-t,double(L-i));
        double tPow=pow(t,double(i));
        double coef=oneMinusTpow*tPow*ncr;
        x+=coef*ctrlpoints[i][0];
        y+=coef*ctrlpoints[i][1];

    }
    xy[0] = float(x);
    xy[1] = float(y);

    //return y;
}

///////////////////////
void setNormal(GLfloat x1, GLfloat y1,GLfloat z1, GLfloat x2, GLfloat y2,GLfloat z2, GLfloat x3, GLfloat y3,GLfloat z3)
{
    GLfloat Ux, Uy, Uz, Vx, Vy, Vz, Nx, Ny, Nz;

    Ux = x2-x1;
    Uy = y2-y1;
    Uz = z2-z1;

    Vx = x3-x1;
    Vy = y3-y1;
    Vz = z3-z1;

    Nx = Uy*Vz - Uz*Vy;
    Ny = Uz*Vx - Ux*Vz;
    Nz = Ux*Vy - Uy*Vx;

    glNormal3f(-Nx,-Ny,-Nz);
}

void bottleBezier()
{
    int i, j;
    float x, y, z, r;				//current coordinates
    float x1, y1, z1, r1;			//next coordinates
    float theta;

    const float startx = 0, endx = ctrlpoints[L][0];
    //number of angular slices
    const float dx = (endx - startx) / nt;	//x step size
    const float dtheta = 2*PI / ntheta;		//angular step size

    float t=0;
    float dt=1.0/nt;
    float xy[2];
    BezierCurve( t,  xy);
    x = xy[0];
    r = xy[1];
    //rotate about z-axis
    float p1x,p1y,p1z,p2x,p2y,p2z;
    for ( i = 0; i < nt; ++i )  			//step through x
    {
        theta = 0;
        t+=dt;
        BezierCurve( t,  xy);
        x1 = xy[0];
        r1 = xy[1];

        //draw the surface composed of quadrilaterals by sweeping theta
        glBegin( GL_QUAD_STRIP );
        for ( j = 0; j <= ntheta; ++j )
        {
            theta += dtheta;
            double cosa = cos( theta );
            double sina = sin ( theta );
            y = r * cosa;
            y1 = r1 * cosa;	//current and next y
            z = r * sina;
            z1 = r1 * sina;	//current and next z

            //edge from point at x to point at next x
            glVertex3f (x, y, z);

            if(j>0)
            {
                setNormal(p1x,p1y,p1z,p2x,p2y,p2z,x, y, z);
            }
            else
            {
                p1x=x;
                p1y=y;
                p1z=z;
                p2x=x1;
                p2y=y1;
                p2z=z1;

            }
            glVertex3f (x1, y1, z1);

            //forms quad with next pair of points with incremented theta value
        }
        glEnd();
        x = x1;
        r = r1;
    } //for i

}


void BezierCurvet ( double t,  float xy[2])
{
    double y=0;
    double x=0;
    t=t>1.0?1.0:t;
    for(int i=0; i<=L1; i++)
    {
        int ncr=nCr(L1,i);
        double oneMinusTpow=pow(1-t,double(L1-i));
        double tPow=pow(t,double(i));
        double coef=oneMinusTpow*tPow*ncr;
        x+=coef*ctrlpoints_tunnel[i][0];
        y+=coef*ctrlpoints_tunnel[i][1];

    }
    xy[0] = float(x);
    xy[1] = float(y);

    //return y;
}


void tunnelBezier()
{
    int i, j;
    float x, y, z, r;				//current coordinates
    float x1, y1, z1, r1;			//next coordinates
    float theta;

    const float startx = 0, endx = ctrlpoints_tunnel[L1][0];
    //number of angular slices
    const float dx = (endx - startx) / nt;	//x step size
    const float dtheta = 2*PI / ntheta;		//angular step size

    float t=0;
    float dt=1.0/nt;
    float xy[2];
    BezierCurvet( t,  xy);
    x = xy[0];
    r = xy[1];
    //rotate about z-axis
    float p1x,p1y,p1z,p2x,p2y,p2z;
    for ( i = 0; i < nt; ++i )  			//step through x
    {
        theta = 0;
        t+=dt;
        BezierCurvet( t,  xy);
        x1 = xy[0];
        r1 = xy[1];

        //draw the surface composed of quadrilaterals by sweeping theta
        glBegin( GL_QUAD_STRIP );
        for ( j = 0; j <= ntheta/2; ++j )
        {
            theta += dtheta;
            double cosa = cos( theta );
            double sina = sin ( theta );
            y = r * cosa;
            y1 = r1 * cosa;	//current and next y
            z = r * sina;
            z1 = r1 * sina;	//current and next z

            //edge from point at x to point at next x
            glVertex3f (x, y, z);

            if(j>0)
            {
                setNormal(p1x,p1y,p1z,p2x,p2y,p2z,x, y, z);
            }
            else
            {
                p1x=x;
                p1y=y;
                p1z=z;
                p2x=x1;
                p2y=y1;
                p2z=z1;

            }
            glVertex3f (x1, y1, z1);

            //forms quad with next pair of points with incremented theta value
        }
        glEnd();
        x = x1;
        r = r1;
    } //for i

}



void showControlPoints()
{
    glPointSize(5.0);
    glColor3f(1.0, 0.0, 1.0);
    glBegin(GL_POINTS);
    for (int i = 0; i <=L; i++)
        glVertex3fv(&ctrlpoints[i][0]);
    glEnd();
}

void matColor(float kdr, float kdg, float kdb,  float shiny, int frnt_Back, float ambFactor, float specFactor)
{
    /*
        const GLfloat mat_ambient[]    = { kdr*ambFactor, kdg*ambFactor, kdb*ambFactor, 1.0f };
        const GLfloat mat_diffuse[]    = { kdr, kdg, kdb, 1.0f };
        const GLfloat mat_specular[]   = { 1.0f*specFactor, 1.0f*specFactor, 1.0f*specFactor, 1.0f };
        const GLfloat high_shininess[] = { shiny };
        if(frnt_Back==0)
        {
            glMaterialfv(GL_FRONT, GL_AMBIENT,   mat_ambient);
            glMaterialfv(GL_FRONT, GL_DIFFUSE,   mat_diffuse);
            glMaterialfv(GL_FRONT, GL_SPECULAR,  mat_specular);
            glMaterialfv(GL_FRONT, GL_SHININESS, high_shininess);
        }
        else if(frnt_Back==1)
        {
            glMaterialfv(GL_BACK, GL_AMBIENT,   mat_ambient);
            glMaterialfv(GL_BACK, GL_DIFFUSE,   mat_diffuse);
            glMaterialfv(GL_BACK, GL_SPECULAR,  mat_specular);
            glMaterialfv(GL_BACK, GL_SHININESS, high_shininess);
        }
        else if(frnt_Back==2)
        {
            glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT,   mat_ambient);
            glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE,   mat_diffuse);
            glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR,  mat_specular);
            glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, high_shininess);
        }
    */
}


void human()
{
    //right leg
    glPushMatrix();
    glTranslated(5,0,0);
    glScaled(1,0.2,-3);
    cube(0.6,0.6,0.6);
    glPopMatrix();

    glPushMatrix();
    glTranslated(5.15,0,0);
    glScaled(0.7,4,-0.7);
    cube(0.1,0.1,0.4);
    glPopMatrix();

    glPushMatrix();
    glTranslated(7,0,0);
    glScaled(1,0.2,-3);
    cube(0.6,0.6,0.6);
    glPopMatrix();

    glPushMatrix();
    glTranslated(7.15,0,0);
    glScaled(0.7,4,-0.7);
    cube(0.1,0.1,0.4);
    glPopMatrix();

    glPushMatrix();
    glTranslated(4.5,4,1.5);
    glScaled(4,7,-3);
    cube(0.6,0.1,0.1);
    glPopMatrix();


    glPushMatrix();
    glTranslated(6,11,0.5);
    glScaled(1,1,-1);
    cube(0.7,0.6,0.6);
    glPopMatrix();


    glPushMatrix();
    glTranslated(5.5,12,.5);
    glScaled(2,2,-1);
    cube(0.8,0.6,0.6);
    glPopMatrix();

    glPushMatrix();
    glTranslated(2.8,6,0);
    glRotated(60,0,0,1);
    glPushMatrix();
    glScaled(5,0.8,-0.8);
    cube(0.5,0.2,0.2);
    glPopMatrix();
    glPopMatrix();

    glPushMatrix();
    glTranslated(7.8,10,0);
    glRotated(-60,0,0,1);
    glPushMatrix();
    glScaled(5,0.8,-0.8);
    cube(0.5,0.2,0.2);
    glPopMatrix();
    glPopMatrix();


}


void human3()
{
    glPushMatrix();
    glCullFace(GL_FRONT);
    glPushMatrix();
    human();
    glPopMatrix();

    glPushMatrix();
    glTranslated(11,0,-0);
    glRotated(60,0,1,0);
    human();
    glPopMatrix();

    glPushMatrix();
    glTranslated(-4,0,-10);
    glRotated(-60,0,1,0);
    human();
    glPopMatrix();

    glPopMatrix();
}



void wheel()
{



    const double t = glutGet(GLUT_ELAPSED_TIME) / 500.0;
    const double a = t*90.0;
    double aa=a;
    if(wired)
    {
        glPolygonMode( GL_FRONT, GL_LINE ) ;
        glPolygonMode( GL_BACK, GL_LINE ) ;

    }
    else
    {
        glPolygonMode( GL_FRONT,GL_FILL ) ;
        glPolygonMode( GL_BACK, GL_FILL ) ;
    }

    glPushMatrix();

    if(animat)
        glRotated(a,0,0,1);



    glRotatef( 90, 0.0, 0.0, 1.0);
    glTranslated(-3.5,0,0);
    glGetDoublev( GL_MODELVIEW_MATRIX, modelview ); //get the modelview info

    //matColor(1,1,1,20);   // front face color
//  matColor(1,1,1,20,1);  // back face color
    glPushMatrix();
    if(rt)
    {
       glRotated(-aa,1,0,0);
    }

    bottleBezier();
    glPopMatrix();

    if(shcpt)
    {
        matColor(0.0,0.0,0.9,20);
        showControlPoints();
    }

    glPopMatrix();

}



void singleTolaHouse(int R,int G,int B)
{

    GLfloat m_no[] = {0, 0, 0, 1.0};
    GLfloat m_amb[] = {r[R%11],g[G%11],b[B%11],1};
    GLfloat m_diff[] = {r[R%11],g[G%11],b[B%11],1};
    GLfloat m_spec[] = {0.1,0.1,0.1,1};
    GLfloat m_sh[] = {1};


    glMaterialfv(GL_FRONT, GL_AMBIENT, m_amb);
    glMaterialfv(GL_FRONT, GL_DIFFUSE, m_diff);
    glMaterialfv(GL_FRONT, GL_SPECULAR, m_spec);
    glMaterialfv(GL_FRONT, GL_SHININESS, m_sh);



    glPushMatrix();
    glTranslated(0,0,0);
    glutSolidCube(1);
    glPopMatrix();

    glMaterialfv(GL_FRONT, GL_AMBIENT,mat_ambient5);
    glMaterialfv(GL_FRONT, GL_DIFFUSE,mat_diffuse5);
    glMaterialfv(GL_FRONT, GL_SPECULAR,mat_specular5);
    glMaterialfv(GL_FRONT, GL_SHININESS,mat_shininess5);
    glPushMatrix();
    glTranslated(0.2,0,0);
    glScaled(0.3,0.3,1.001);
    glutSolidCube(1);
    glPopMatrix();

    glMaterialfv(GL_FRONT, GL_AMBIENT,mat_ambient5);
    glMaterialfv(GL_FRONT, GL_DIFFUSE,mat_diffuse5);
    glMaterialfv(GL_FRONT, GL_SPECULAR,mat_specular5);
    glMaterialfv(GL_FRONT, GL_SHININESS,mat_shininess5);
    glPushMatrix();
    glTranslated(-0.2,0,0);
    glScaled(0.3,0.3,1.001);
    glutSolidCube(1);
    glPopMatrix();

    glMaterialfv(GL_FRONT, GL_AMBIENT,mat_ambient5);
    glMaterialfv(GL_FRONT, GL_DIFFUSE,mat_diffuse5);
    glMaterialfv(GL_FRONT, GL_SPECULAR,mat_specular5);
    glMaterialfv(GL_FRONT, GL_SHININESS,mat_shininess5);
    glPushMatrix();
    glTranslated(0,0,0.2);
    glScaled(1.001,0.3,0.3);
    glutSolidCube(1);
    glPopMatrix();

    glMaterialfv(GL_FRONT, GL_AMBIENT,mat_ambient5);
    glMaterialfv(GL_FRONT, GL_DIFFUSE,mat_diffuse5);
    glMaterialfv(GL_FRONT, GL_SPECULAR,mat_specular5);
    glMaterialfv(GL_FRONT, GL_SHININESS,mat_shininess5);
    glPushMatrix();
    glTranslated(0,0,-0.2);
    glScaled(1.001,0.3,0.3);
    glutSolidCube(1);
    glPopMatrix();

}

void house(int n,int R,int G)
{
    //road();
    for(int i=0; i<n; i++)
    {
        glPushMatrix();
        glTranslated(0,0.8+i,0);
        singleTolaHouse(R,G,i);

        glPopMatrix();
    }
}

void colony()
{
    for(int i=-(EN_SIZE/2)+1; i<(EN_SIZE/2); i+=2)
    {
        for(int j=-(EN_SIZE/2)+1; j<(EN_SIZE/2); j+=4)
        {
            if(tola[i+(EN_SIZE/2)+1][j+(EN_SIZE/2)+1]!=0)
            {
                glPushMatrix();
                glTranslated(i,0,j);
                house(tola[i+(EN_SIZE/2)+1][j+(EN_SIZE/2)+1],i,j);
                glPopMatrix();
            }
            else if(i>=-5&&i<=5) {}
            else
            {
                tola[i+(EN_SIZE/2)+1][j+(EN_SIZE/2)+1]=(rand()%5)+1;
                glPushMatrix();
                glTranslated(i,0,j);

                house(tola[i+(EN_SIZE/2)+1][j+(EN_SIZE/2)+1],i,j);

                glPopMatrix();
            }
        }
    }
}


void colony2()
{
    /*
    glPushMatrix();
    glTranslated(0,0,0);
    house(5,2,3);
    glPopMatrix();

    glPushMatrix();
    glTranslated(2,0,0);
    house(6,1,2);
    glPopMatrix();

    glPushMatrix();
    glTranslated(4,0,0);
    house(6,2,3);
    glPopMatrix();*/

    glPushMatrix();
    glTranslated(0,0,2);
    house(3,1,2);
    glPopMatrix();


    glPushMatrix();
    glTranslated(2,0,2);
    house(6,2,3);
    glPopMatrix();

    glPushMatrix();
    glTranslated(4,0,2);
    house(6,1,2);
    glPopMatrix();

    glPushMatrix();
    glTranslated(0,0,4);
    house(8,2,3);
    glPopMatrix();

    glPushMatrix();
    glTranslated(2,0,4);
    house(2,1,2);
    glPopMatrix();


    glPushMatrix();
    glTranslated(4,0,4);
    house(4,2,3);
    glPopMatrix();

    glPushMatrix();
    glTranslated(0,0,6);
    house(6,1,2);
    glPopMatrix();

    glPushMatrix();
    glTranslated(2,0,6);
    house(4,2,3);
    glPopMatrix();

    glPushMatrix();
    glTranslated(4,0,6);
    house(5,1,2);
    glPopMatrix();


    glPushMatrix();
    glTranslated(0,0,8);
    house(2,7,2);
    glPopMatrix();

    glPushMatrix();
    glTranslated(2,0,8);
    house(2,2,7);
    glPopMatrix();

    glPushMatrix();
    glTranslated(4,0,8);
    house(5,8,2);
    glPopMatrix();

}

void triangle(double r, double g, double b)
{

    glPushMatrix();

    GLfloat m_no[] = {0, 0, 0, 1.0};
    GLfloat m_amb[] = {r,g,b,1};
    GLfloat m_diff[] = {r,g,b,1};
    GLfloat m_spec[] = {1,1,1,1};
    GLfloat m_sh[] = {40};

    glMaterialfv(GL_FRONT, GL_AMBIENT, m_amb);
    glMaterialfv(GL_FRONT, GL_DIFFUSE, m_diff);
    glMaterialfv(GL_FRONT, GL_SPECULAR, m_spec);
    glMaterialfv(GL_FRONT, GL_SHININESS, m_sh);


    glBegin(GL_POLYGON);

    glMaterialfv(GL_FRONT, GL_AMBIENT, m_amb);
    glMaterialfv(GL_FRONT, GL_DIFFUSE, m_diff);
    glMaterialfv(GL_FRONT, GL_SPECULAR, m_spec);
    glMaterialfv(GL_FRONT, GL_SHININESS, m_sh);
    glVertex3f(-0.3,5,-20.1);


    glMaterialfv(GL_FRONT, GL_AMBIENT, m_amb);
    glMaterialfv(GL_FRONT, GL_DIFFUSE, m_diff);
    glMaterialfv(GL_FRONT, GL_SPECULAR, m_spec);
    glMaterialfv(GL_FRONT, GL_SHININESS, m_sh);
    glVertex3f(10.3,5,-20.1);


    glMaterialfv(GL_FRONT, GL_AMBIENT, m_amb);
    glMaterialfv(GL_FRONT, GL_DIFFUSE, m_diff);
    glMaterialfv(GL_FRONT, GL_SPECULAR, m_spec);
    glMaterialfv(GL_FRONT, GL_SHININESS, m_sh);
    glVertex3f(5, 9.9, -20.1);


    glEnd();
    glPopMatrix();
}

void galarry()
{
    glPushMatrix();
    glScaled(1,0.5,0.5);

    glPushMatrix();
    glTranslated(-10,0,-10);
    glScaled(40,3,6);
    cube(0.9,0.5,0.5);
    glPopMatrix();

    glPushMatrix();
    glTranslated(-10,3,-7);
    glScaled(40,3,6);
    cube(0.5,0.9,0.5);
    glPopMatrix();

    glPushMatrix();
    glTranslated(-10,6,-4);
    glScaled(40,3,6);
    cube(0.5,0.5,0.9);
    glPopMatrix();


    glPushMatrix();
    glTranslated(-10,9,-1);
    glScaled(40,3,6);
    cube(0.9,0.5,0.5);
    glPopMatrix();

    glPushMatrix();
    glTranslated(-10,12,2);
    glScaled(40,3,6);
    cube(0.5,0.9,0.5);
    glPopMatrix();

    glPushMatrix();
    glTranslated(-10,15,5);
    glScaled(40,3,6);
    cube(0.5,0.5,0.9);
    glPopMatrix();

    glPopMatrix();

}


void hospital()
{
    glPushMatrix();
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D,6);
    glMaterialfv(GL_FRONT, GL_AMBIENT,mat_ambient3);
    glMaterialfv(GL_FRONT, GL_DIFFUSE,mat_diffuse3);
    glMaterialfv(GL_FRONT, GL_SPECULAR,mat_specular3);
    glMaterialfv(GL_FRONT, GL_SHININESS,mat_shininess3);
    glRotated(135,0,1,0);
    glTranslated(0,-5,0);
    glScaled(15,12,-15);
    Cube();
    glPopMatrix();
    glDisable(GL_TEXTURE_2D);


    glPushMatrix();
    glRotated(135,0,1,0);
    glTranslated(0,7,0);
    glScaled(15,2,-15);
    cube(0.6,0.6,0.6);
    glPopMatrix();


    glPushMatrix();
    glTranslated(-2.5,0,0);
    glScaled(0.8,1,0.8);
    glPushMatrix();
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D,6);
    glMaterialfv(GL_FRONT, GL_AMBIENT,mat_ambient3);
    glMaterialfv(GL_FRONT, GL_DIFFUSE,mat_diffuse3);
    glMaterialfv(GL_FRONT, GL_SPECULAR,mat_specular3);
    glMaterialfv(GL_FRONT, GL_SHININESS,mat_shininess3);
    glRotated(135,0,1,0);
    glTranslated(0,9,0);
    glScaled(15,12,-15);
    Cube();
    glPopMatrix();
    glDisable(GL_TEXTURE_2D);
    glPopMatrix();

}



void Play_ground()
{

    glPushMatrix();
    glTranslated(0,0.2,0);
    glPushMatrix();
    glTranslated(-18.6,0.2,-5);
    glScaled(0.08,0.08,0.08);
    glRotated(-90,0,1,0);
    human();
    glPopMatrix();

    glPushMatrix();
    glTranslated(-16,0.2,-5);
    glScaled(0.1,0.1,0.1);
    glRotated(-90,0,1,0);
    human();
    glPopMatrix();

    glPushMatrix();
    glTranslated(-17,0.2,-3);
    glScaled(0.1,0.1,0.1);
    glRotated(-90,0,1,0);
    human();
    glPopMatrix();

    glPushMatrix();
    glTranslated(-17,0.2,-7);
    glScaled(0.1,0.1,0.1);
    glRotated(-90,0,1,0);
    human();
    glPopMatrix();

    glPushMatrix();
    glTranslated(-15,0.2,-1);
    glScaled(0.1,0.1,0.1);
    glRotated(-90,0,1,0);
    human();
    glPopMatrix();

    glPushMatrix();
    glTranslated(-14,0.2,-3);
    glScaled(0.1,0.1,0.1);
    glRotated(-90,0,1,0);
    human();
    glPopMatrix();

    glPushMatrix();
    glTranslated(-14,0.2,-5);
    glScaled(0.1,0.1,0.1);
    glRotated(-90,0,1,0);
    human();
    glPopMatrix();

    glPushMatrix();
    glTranslated(-15,0.2,-7);
    glScaled(0.1,0.1,0.1);
    glRotated(-90,0,1,0);
    human();
    glPopMatrix();

    glPopMatrix();

//2nd team


    glPushMatrix();
    glRotated(-180,0,1,0);
    glTranslated(26,0.2,8);
    glPushMatrix();
    glTranslated(-18.6,0.2,-5);
    glScaled(0.08,0.08,0.08);
    glRotated(-90,0,1,0);
    human();
    glPopMatrix();

    glPushMatrix();
    glTranslated(-16,0.2,-5);
    glScaled(0.1,0.1,0.1);
    glRotated(-90,0,1,0);
    human();
    glPopMatrix();

    glPushMatrix();
    glTranslated(-17,0.2,-3);
    glScaled(0.1,0.1,0.1);
    glRotated(-90,0,1,0);
    human();
    glPopMatrix();

    glPushMatrix();
    glTranslated(-17,0.2,-7);
    glScaled(0.1,0.1,0.1);
    glRotated(-90,0,1,0);
    human();
    glPopMatrix();

    glPushMatrix();
    glTranslated(-15,0.2,-1);
    glScaled(0.1,0.1,0.1);
    glRotated(-90,0,1,0);
    human();
    glPopMatrix();

    glPushMatrix();
    glTranslated(-14,0.2,-3);
    glScaled(0.1,0.1,0.1);
    glRotated(-90,0,1,0);
    human();
    glPopMatrix();

    glPushMatrix();
    glTranslated(-14,0.2,-5);
    glScaled(0.1,0.1,0.1);
    glRotated(-90,0,1,0);
    human();
    glPopMatrix();

    glPushMatrix();
    glTranslated(-15,0.2,-7);
    glScaled(0.1,0.1,0.1);
    glRotated(-90,0,1,0);
    human();
    glPopMatrix();

    glPopMatrix();



    glPushMatrix();
    glTranslated(-7,0.2,1);
    glScaled(0.6,1,1);
    glRotated(90,0,1,0);


    glPushMatrix();
    glTranslated(0,0,0);
    glScaled(10,0.1,-20);
    cube(0.1,1,0.3);
    glPopMatrix();




    glPushMatrix();
    glTranslated(4,0,0);
    glScaled(0.05,1,0.05);
    cube(0.9,0.5,0.85);
    glPopMatrix();


    glPushMatrix();
    glTranslated(6,0,-20);
    glScaled(0.05,1,0.05);
    cube(0.9,0.5,0.85);
    glPopMatrix();

    glPushMatrix();
    glTranslated(4,1,0);
    glScaled(2,0.05,0.05);
    cube(0.9,0.5,0.85);
    glPopMatrix();


    glPushMatrix();
    glTranslated(4,0,-20);
    glScaled(0.05,1,0.05);
    cube(0.9,0.5,0.85);
    glPopMatrix();


    glPushMatrix();
    glTranslated(6,0,0);
    glScaled(0.05,1,0.05);
    cube(0.9,0.5,0.85);
    glPopMatrix();


    glPushMatrix();
    glTranslated(4,1,-20);
    glScaled(2,0.05,0.05);
    cube(0.9,0.5,0.85);
    glPopMatrix();

    //bilbord
    glPushMatrix();
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D,2);
    GLfloat mat_ambient3[] = { 1,1,1, 1.0 };
    GLfloat mat_diffuse3[] = { 1,1,1, 1.0 };
    GLfloat mat_specular3[] = { 1.0, 1.0, 1.0, 1.0 };
    glTranslated(-2,2,-1);
    glScaled(0.1,2,-3);
    Cube();
    glPopMatrix();
    glDisable(GL_TEXTURE_2D);

    glPushMatrix();
    glTranslated(-2,0,-1.5);
    glScaled(0.1,2,-0.1);
    cube(0,0,0);
    glPopMatrix();

    glPushMatrix();
    glTranslated(-2,0,-3.5);
    glScaled(0.1,2,-0.1);
    cube(0,0,0);
    glPopMatrix();


    glPopMatrix();

}



void stadiumEnv()
{
    /// Ground
    glPushMatrix();
    glMaterialfv(GL_FRONT, GL_AMBIENT,mat_ambient2);
    glMaterialfv(GL_FRONT, GL_DIFFUSE,mat_diffuse2);
    glMaterialfv(GL_FRONT, GL_SPECULAR,mat_specular2);
    glMaterialfv(GL_FRONT, GL_SHININESS,mat_shininess2);
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D,5);
    glCullFace(GL_BACK);
    glTranslated(0,0,0);
    glScaled(EN_SIZE*2,0.3,EN_SIZE*2);
    glutSolidCube(1);
    glPopMatrix();
    glDisable(GL_TEXTURE_2D);



    glPushMatrix();
    glRotated(180,0,1,0);
    glTranslated(10,0,11);
    glScaled(0.3,0.3,0.2);
    galarry();
    glPopMatrix();

    glPushMatrix();
    glRotated(0,0,1,0);
    glTranslated(-17,0,3);
    glScaled(0.2,0.3,0.2);
    galarry();
    glPopMatrix();

    glPushMatrix();
    glCullFace(GL_FRONT);
    Play_ground();
    glPopMatrix();



}

void sky()
{
    glPushMatrix();
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D,4);
    glMaterialfv(GL_FRONT, GL_AMBIENT,mat_ambient3);
    glMaterialfv(GL_FRONT, GL_DIFFUSE,mat_diffuse3);
    glMaterialfv(GL_FRONT, GL_SPECULAR,mat_specular3);
    glMaterialfv(GL_FRONT, GL_SHININESS,mat_shininess3);
    Cube();
    glPopMatrix();
    glDisable(GL_TEXTURE_2D);
}


void road()
{
    glPushMatrix();
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D,4);
    glMaterialfv(GL_FRONT, GL_AMBIENT,mat_ambient3);
    glMaterialfv(GL_FRONT, GL_DIFFUSE,mat_diffuse3);
    glMaterialfv(GL_FRONT, GL_SPECULAR,mat_specular3);
    glMaterialfv(GL_FRONT, GL_SHININESS,mat_shininess3);
    Cube();
    glPopMatrix();
    glDisable(GL_TEXTURE_2D);

}

void environment(int n)
{

    /// Ground
    glMaterialfv(GL_FRONT, GL_AMBIENT,mat_ambient2);
    glMaterialfv(GL_FRONT, GL_DIFFUSE,mat_diffuse2);
    glMaterialfv(GL_FRONT, GL_SPECULAR,mat_specular2);
    glMaterialfv(GL_FRONT, GL_SHININESS,mat_shininess2);
    glPushMatrix();
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D,5);

    glTranslated(0,0,0);
    glScaled(EN_SIZE*2,0.3,EN_SIZE*2);
    //glScaled(EN_SIZE,0.3,EN_SIZE*2);
    glutSolidCube(1);
    //cube();
    glPopMatrix();
    glDisable(GL_TEXTURE_2D);


    glPushMatrix();
    glPushMatrix();
    glTranslated(-(EN_SIZE*0.33),0.1,-EN_SIZE);
    glScaled(EN_SIZE*0.66,0.3,EN_SIZE*2);
    glCullFace(GL_BACK);
    road();
    glPopMatrix();
    glPushMatrix();
    // glColor3d(0,1,0.1);
    glTranslated(0,0,0);
    //glTranslated(torusPosX[n],1,0);
    glScaled(0.3,0.3,0.3);
    //glutSolidTorus(1,3,30,30);
    glPopMatrix();
    glPopMatrix();

    //  colony();

}


void busStand()
{

    //sadd
    glPushMatrix();
    glTranslated(-5,10,-10);
    glScaled(15,0.5,15);
    cube(0.6,0.3,0);
    glPopMatrix();


    glPushMatrix();
    glTranslated(2.5,0,-2.5);
    glScaled(.3,10,.3);
    cube(0.8,0.8,0.8);
    glPopMatrix();


    // seat
    glPushMatrix();
    glScaled(5,2,1.5);
    cube(0.7,0.6,0.6);
    glPopMatrix();


//helan
    glPushMatrix();
    glTranslated(0,0,1.5);
    glScaled(5,4,1);
    cube(.7,0.7,0.4);
    glPopMatrix();



    //side
    glPushMatrix();
    glTranslated(6,0,0);
    glRotated(90,0,1,0);
    // seat
    glPushMatrix();
    glScaled(5,2,1.5);
    cube(0.7,0.6,0.6);
    glPopMatrix();

//helan
    glPushMatrix();
    glTranslated(0,0,1.5);
    glScaled(5,4,1);
    cube(.7,0.7,0.4);
    glPopMatrix();
    glPopMatrix();


    //side
    glPushMatrix();
    glTranslated(-1,0,-5);
    glRotated(-90,0,1,0);
    // seat
    glPushMatrix();
    glScaled(5,2,1.5);
    cube(0.7,0.6,0.6);
    glPopMatrix();

//helan
    glPushMatrix();
    glTranslated(0,0,1.5);
    glScaled(5,4,1);
    cube(.7,0.7,0.4);
    glPopMatrix();
    glPopMatrix();

    glPushMatrix();
    glScaled(0.5,0.5,0.5);
    human3();
    glPopMatrix();

}



void bari()
{
    glPushMatrix();


    glPushMatrix();
    glTranslated(0,0,0);
    glScaled(10,0.1,-20);
    cube(0.6,0.5,0.4);
    glPopMatrix();


// 4 wall
    glPushMatrix();
    glTranslated(0,0,0);
    glScaled(0.1,5,-20);
    cube(1,0.6,0.6);
    glPopMatrix();


    glPushMatrix();
    glTranslated(10,0,0);
    glScaled(0.1,5,-20);
    cube(1,0.6,0.6);
    glPopMatrix();


    glPushMatrix();
    glColor3d(1,1,1);
    glTranslated(0,0,0);
    glScaled(10,5,-0.1);
    cube(1,0.6,0.6);
    glPopMatrix();



    glPushMatrix();
    glColor3d(1,1,1);
    glTranslated(0,0,0);
    glScaled(10,5,-20);
    cube(1,0.6,0.6);
    glPopMatrix();


    //tin

    glPushMatrix();
    glEnable(GL_TEXTURE_2D);
    glMaterialfv(GL_FRONT, GL_AMBIENT,mat_ambient3);
    glMaterialfv(GL_FRONT, GL_DIFFUSE,mat_diffuse3);
    glMaterialfv(GL_FRONT, GL_SPECULAR,mat_specular3);
    glMaterialfv(GL_FRONT, GL_SHININESS,mat_shininess3);
    glBindTexture(GL_TEXTURE_2D,3);

    glPushMatrix();
    glTranslated(3.1,0.8,0);
    glRotated(45,0,0,1);
    glPushMatrix();
    glTranslated(0,5.1,0);
    glScaled(7,0.1,-20);
    Cube();
    //cube(0.9,0.9,0.9);
    glPopMatrix();
    glPopMatrix();
    glPopMatrix();
    glDisable(GL_TEXTURE_2D);

    glPushMatrix();
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D,3);
    glMaterialfv(GL_FRONT, GL_AMBIENT,mat_ambient3);
    glMaterialfv(GL_FRONT, GL_DIFFUSE,mat_diffuse3);
    glMaterialfv(GL_FRONT, GL_SPECULAR,mat_specular3);
    glMaterialfv(GL_FRONT, GL_SHININESS,mat_shininess3);
    glTranslated(-5,13,0);
    glRotated(-45,0,0,1);
    glPushMatrix();
    glTranslated(10,5.1,0);
    glScaled(7,0.1,-20);
    Cube();
    glPopMatrix();
    glPopMatrix();
    glDisable(GL_TEXTURE_2D);

    glPushMatrix();
    glTranslated(4.4,9.5,0);
    glScaled(1.4,0.1,-20);
    cube(0.4,0.4,0.4);
    glPopMatrix();


    glPushMatrix();
    glTranslated(10.3,1,-16);
    glScaled(0.1,2.5,2.5);
    cube(1,1,1);
    glPopMatrix();


    glPushMatrix();
    glTranslated(10.3,1,-8);
    glScaled(0.1,2.5,2.5);
    cube(1,1,1);
    glPopMatrix();


    glPushMatrix();
    glTranslated(5,1,-20.2);
    glScaled(2.5,2.5,0.1);
    cube(1,1,1);
    glPopMatrix();


    glPushMatrix();
    glTranslated(0.2,-0.2,0.2);
    triangle(1,0.2,0.2);
    glPopMatrix();

    glPushMatrix();
    glCullFace(GL_BACK);
    glTranslated(18,0,-5);
    busStand();
    glPopMatrix();


    glPushMatrix();
    glCullFace(GL_BACK);
    glTranslated(80,23,-15);
    busStand();
    glPopMatrix();

    glPopMatrix();

}



void game_over()
{


    glPushMatrix();
    ///Environment
    if(tX>=5)
    {
        START=false;
        START1=true;
        tX=0;
    }


    if(tX<=-5.1)
    {
        START=false;
        START1=true;
        tX=0;
    }



    //car7
    if((tZ5<=-108))
    {
        if(tX>=-3.5)
        {
            if(tX<=-1.5)
            {
                START=false;
                START1=true;
                tX=0;
            }
        }

    }


    //car2
    if((tZ<=-108))
    {
        if(tX>=1)
        {
            if(tX<=3)
            {
                START=false;
                START1=true;
                tX=0;
            }
        }

    }

//car1 ok

    if((tZ6<=-108))
    {
        if(tX>=1)
        {
            if(tX<=3)
            {
                START=false;
                START1=true;
                tX=0;
            }
        }

    }


//car4 ok
    if((tZ2<=-108))
    {
        if(tX>=-1)
        {
            if(tX<=1)
            {
                START=false;
                START1=true;
                tX=3;
                // tZ1=-5;tZ1=-20;tZ2=-40;tZ3=-60;tZ4=-80;tZ5=-100;tZ6=-120;
            }
        }

    }



//car6 ok
    if((tZ4<=-108))
    {
        if(tX>=-1)
        {
            if(tX<=1)
            {
                START=false;
                START1=true;
                tX=3;
                // tZ1=-5;tZ1=-20;tZ2=-40;tZ3=-60;tZ4=-80;tZ5=-100;tZ6=-120;
            }
        }

    }



//car3 ok
    if((tZ1<=-108))
    {
        if(tX>=-3.5)
        {
            if(tX<=-1.5)
            {
                START=false;
                START1=true;
                tX=0;
            }
        }

    }


    //car5 ok
    if((tZ3<=-108))
    {
        if(tX>=-3.5)
        {
            if(tX<=-1.5)
            {
                START=false;
                START1=true;
                tX=0;
            }
        }

    }



    if(tY>0.1)
        tY= 0.1;
    if(tY<-15)
        tY= -15;
    glPopMatrix();


}


void car_tire(double r,double g,double b)
{



    GLfloat m_no[] = {0, 0, 0, 1.0};
    GLfloat m_amb[] = {r,g,b,1};
    GLfloat m_diff[] = {r,g,b,1};
    GLfloat m_spec[] = {1,1,1,1};
    GLfloat m_sh[] = {40};

    glMaterialfv(GL_FRONT, GL_AMBIENT, m_amb);
    glMaterialfv(GL_FRONT, GL_DIFFUSE, m_diff);
    glMaterialfv(GL_FRONT, GL_SPECULAR, m_spec);
    glMaterialfv(GL_FRONT, GL_SHININESS, m_sh);

    glBegin(GL_POLYGON);
    glMaterialfv(GL_FRONT, GL_AMBIENT, m_amb);
    glMaterialfv(GL_FRONT, GL_DIFFUSE, m_diff);
    glMaterialfv(GL_FRONT, GL_SPECULAR, m_spec);
    glMaterialfv(GL_FRONT, GL_SHININESS, m_sh);


    for(int i=0; i<360; i=i+10)
    {
        float x,y;
        x = 3.2 * cos(((-3.1416)*i)/180);
        y = 3.2 * sin(((-3.1416)*i)/180);
        glVertex3f(x, y, -.1);

    }

    glEnd();

    glBegin(GL_POLYGON);
    glMaterialfv(GL_FRONT, GL_AMBIENT, m_amb);
    glMaterialfv(GL_FRONT, GL_DIFFUSE, m_diff);
    glMaterialfv(GL_FRONT, GL_SPECULAR, m_spec);
    glMaterialfv(GL_FRONT, GL_SHININESS, m_sh);


    for(int i=0; i<360; i=i+10)
    {
        float p,q;
        p = 3.2 * cos(((3.1416)*i)/180);
        q = 3.2 * sin(((3.1416)*i)/180);
        glVertex3f(p, q, .1);

    }

    glEnd();

    glBegin(GL_QUAD_STRIP);
    glColor3d(r+.5,g+.5,b);


    for(int i=0; i<=360; i=i+10)
    {
        float m,n;
        m= 3.2 * cos(((-3.1416)*i)/180);
        n = 3.2 * sin(((-3.1416)*i)/180);
        glVertex3f(m,n, -.1);
        glVertex3f(m,n, .1);
        glMaterialfv(GL_FRONT, GL_AMBIENT, m_amb);
        glMaterialfv(GL_FRONT, GL_DIFFUSE, m_diff);
        glMaterialfv(GL_FRONT, GL_SPECULAR, m_spec);
        glMaterialfv(GL_FRONT, GL_SHININESS, m_sh);
    }
    glEnd();
}




void frontlight(double r, double g, double b)
{

    GLfloat m_no[] = {0, 0, 0, 1.0};
    GLfloat m_amb[] = {r,g,b,1};
    GLfloat m_diff[] = {r,g,b,1};
    GLfloat m_spec[] = {1,1,1,1};
    GLfloat m_sh[] = {40};

    glMaterialfv(GL_FRONT, GL_AMBIENT, m_amb);
    glMaterialfv(GL_FRONT, GL_DIFFUSE, m_diff);
    glMaterialfv(GL_FRONT, GL_SPECULAR, m_spec);
    glMaterialfv(GL_FRONT, GL_SHININESS, m_sh);

    glBegin(GL_POLYGON);
    glMaterialfv(GL_FRONT, GL_AMBIENT, m_amb);
    glMaterialfv(GL_FRONT, GL_DIFFUSE, m_diff);
    glMaterialfv(GL_FRONT, GL_SPECULAR, m_spec);
    glMaterialfv(GL_FRONT, GL_SHININESS, m_sh);


    for(int i=0; i<360; i=i+10)
    {
        float z,y;
        z = .6 * cos(((-3.1416)*i)/180);
        y = .6 * sin(((-3.1416)*i)/180);
        glVertex3f(-.1, y, z);

    }

    glEnd();

    glBegin(GL_POLYGON);
    glMaterialfv(GL_FRONT, GL_AMBIENT, m_amb);
    glMaterialfv(GL_FRONT, GL_DIFFUSE, m_diff);
    glMaterialfv(GL_FRONT, GL_SPECULAR, m_spec);
    glMaterialfv(GL_FRONT, GL_SHININESS, m_sh);


    for(int i=0; i<360; i=i+10)
    {
        float z,y;
        z = .4 * cos(((-3.1416)*i)/180);
        y = .4 * sin(((-3.1416)*i)/180);
        glVertex3f(.08, y, z);

    }

    glEnd();

    glBegin(GL_QUAD_STRIP);
    glMaterialfv(GL_FRONT, GL_AMBIENT, m_amb);
    glMaterialfv(GL_FRONT, GL_DIFFUSE, m_diff);
    glMaterialfv(GL_FRONT, GL_SPECULAR, m_spec);
    glMaterialfv(GL_FRONT, GL_SHININESS, m_sh);


    for(int i=0; i<=360; i=i+10)
    {
        float z,y,z1,y1;
        z = .8 * cos(((3.1416)*i)/180);
        y = .8 * sin(((3.1416)*i)/180);

        z1 = .4 * cos(((3.1416)*i)/180);
        y1 = .4 * sin(((3.1416)*i)/180);

        glVertex3f(-.1,y,z);
        glVertex3f(.08,y1,z1);

    }

    glEnd();
}

void wheel_4()
{


    glPushMatrix();
    glTranslated(0,0,0);
    glRotated(90,1,0,0);
    glScaled(0.5,0.5,0.5);
    wheel();
    glPopMatrix();




    glPushMatrix();
    glTranslated(0,0,-7);
    glRotated(90,1,0,0);
    glScaled(0.5,0.5,0.5);
    wheel();
    glPopMatrix();


    glPushMatrix();
    glTranslated(0,7,0);
    glRotated(90,1,0,0);
    glScaled(0.5,0.5,0.5);
    wheel();
    glPopMatrix();




    glPushMatrix();
    glTranslated(0,7,-7);
    glRotated(90,1,0,0);
    glScaled(0.5,0.5,0.5);
    wheel();
    glPopMatrix();


}

void body()

{

    const double t = glutGet(GLUT_ELAPSED_TIME) / 1000.0;
    const double a = t*180.0;

    glPushMatrix();
    glTranslated(0,2.5,0);

    glPushMatrix();

    cube1(0, 1, 0);

    glTranslated(0, -2, 0);
    ownScale(2.3, 1, 1.1);
    cube1(0, 1, 0);
// Four Glass
    glTranslated(0,2,.95);
    ownScale(.4, .8, 0.05);
    cube1(1,1,0);

    glTranslated(0,0,-38);
    ownScale(1, 1, 1);
    cube1(1,1,0);

    glTranslated(1.15,0,19);
    ownScale(.07, 1, 16);
    cube1(1,1,0);

    glTranslated(-33,0,0);
    ownScale(1,1,1);
    cube1(1,1,0);


//Backside Bumper

    glTranslated(-20,-1.8,0);
    ownScale(1,.2,1.1);
    cube1(1,0,0);

    glTranslated(0,-7,0);
    ownScale(1,1.4,1);
    cube1(0.5,.5,0.5);


//Front side Bumper

    glTranslated(73,0,0);
    ownScale(1,1,1);
    cube1(0.5,.5,0.5);


//Tires




    /*

        glPushMatrix();
        glTranslated(-56,-2,1.32);
        ownScale(3,1.5,2);
        glRotated(a,0,0,1);
        car_tire(0, 0, 0);
        glPopMatrix();


        glPushMatrix();
        glTranslated(-18,-2,1.32);
        ownScale(3,1.5,2);
        glRotated(a,0,0,1);
        car_tire(0, 0, 0);
        glPopMatrix();


        glPushMatrix();
        glTranslated(-56,-2,-1.35);
        ownScale(3,1.5,2);
        glRotated(a,0,0,1);
        car_tire(0, 0, 0);
        glPopMatrix();
        glPushMatrix();
        glTranslated(-18,-2,-1.35);
        ownScale(3,1.5,2);
        glRotated(a,0,0,1);
        car_tire(0, 0, 0);
        glPopMatrix();

    */
//Front side Light

    glTranslated(.8,4.6,.77);
    ownScale(15,2,.2);
    frontlight(0.9,0.8,0);

    glTranslated(0,0,-7.7);
    ownScale(1,1,1);
    frontlight(0.9,0.8,0.7);


    glPopMatrix();

    glPopMatrix();


}

void car1()
{


    wheel_4();


    glPushMatrix();
    glTranslated(-1.1,4,-5);
    glRotated(90,0,0,1);
    glScaled(2.5,2.5,2.5);
    body();
    glPopMatrix();

}

void shopingmall()
{


    glPushMatrix();
    glScaled(25,0.1,-25);
    cube(0.4,0.4,0.4);
    glPopMatrix();
    glPushMatrix();

    glTranslated(0,10,0);
    glScaled(25,0.5,-25);
    cube(0.7,0.7,1);
    glPopMatrix();


    glPushMatrix();
    glTranslated(16,0,0);
    glScaled(0.1,10,-16);
    cube(0,1,0);
    glPopMatrix();

    glPushMatrix();
    glTranslated(0,0,-16);
    glScaled(16,10,-0.1);
    cube(0,1,0);
    glPopMatrix();


    glPushMatrix();
    glTranslated(0,0,0);
    glScaled(0.1,10,-16);
    cube(0,1,0);
    glPopMatrix();

    glPushMatrix();
    glTranslated(0,0,0);
    glScaled(16,10,-0.1);
    cube(0,1,0);
    glPopMatrix();



    //khuti

    glPushMatrix();
    glTranslated(3,0,-22);
    glScaled(1,10,-1);
    cube(0,1,1);
    glPopMatrix();


    glPushMatrix();
    glTranslated(22,0,-3);
    glScaled(1,10,-1);
    cube(0,1,1);
    glPopMatrix();

    glPushMatrix();
    glTranslated(22,0,-22);
    glScaled(1,10,-1);
    cube(0,1,1);
    glPopMatrix();



    //left side janala

    glPushMatrix();
    glTranslated(16.1,3,-2);
    glScaled(0.1,5,-3);
    cube(0,0,0);
    glPopMatrix();


    glPushMatrix();
    glTranslated(16.1,3,-6);
    glScaled(0.1,5,-3);
    cube(0,0,0);
    glPopMatrix();

    glPushMatrix();
    glTranslated(16.1,3,-10);
    glScaled(0.1,5,-3);
    cube(0,0,0);
    glPopMatrix();

    //glass
    glPushMatrix();
    glTranslated(16.2,3.25,-2.25);
    glScaled(0.1,4.5,-2.5);
    cube(1,1,1);
    glPopMatrix();


    glPushMatrix();
    glTranslated(16.2,3.25,-6.25);
    glScaled(0.1,4.5,-2.5);
    cube(1,1,1);
    glPopMatrix();

    glPushMatrix();
    glTranslated(16.2,3.25,-10.25);
    glScaled(0.1,4.5,-2.5);
    cube(1,1,1);
    glPopMatrix();



    //front side

    glPushMatrix();
    glTranslated(6,0.5,-16.1);
    glScaled(4,8,-0.1);
    cube(0,0,0);
    glPopMatrix();

    glPushMatrix();
    glTranslated(6.5,1,-16.2);
    glScaled(3,7,-0.1);
    cube(1,1,1);
    glPopMatrix();



    glPushMatrix();
    glTranslated(11,3,-16.1);
    glScaled(3,3,-0.1);
    cube(0,0,0);
    glPopMatrix();

    glPushMatrix();
    glTranslated(11.2,3.2,-16.2);
    glScaled(2.6,2.6,-0.1);
    cube(1,1,1);
    glPopMatrix();



    glPushMatrix();
    glTranslated(2,3,-16.1);
    glScaled(3,3,-0.1);
    cube(0,0,0);
    glPopMatrix();

    glPushMatrix();
    glTranslated(2.2,3.2,-16.2);
    glScaled(2.6,2.6,-0.1);
    cube(1,1,1);
    glPopMatrix();



//2ndflor


    glPushMatrix();
    glTranslated(0,10.5,0);
    glScaled(16,12,-16);
    cube(0,1,0);
    glPopMatrix();


    //left side janala

    glPushMatrix();
    glTranslated(16.1,15,-2);
    glScaled(0.1,5,-3);
    cube(0,0,0);
    glPopMatrix();


    glPushMatrix();
    glTranslated(16.1,15,-6);
    glScaled(0.1,5,-3);
    cube(0,0,0);
    glPopMatrix();

    glPushMatrix();
    glTranslated(16.1,15,-10);
    glScaled(0.1,5,-3);
    cube(0,0,0);
    glPopMatrix();

    //glass
    glPushMatrix();
    glTranslated(16.2,15.25,-2.25);
    glScaled(0.1,4.5,-2.5);
    cube(1,1,1);
    glPopMatrix();


    glPushMatrix();
    glTranslated(16.2,15.25,-6.25);
    glScaled(0.1,4.5,-2.5);
    cube(1,1,1);
    glPopMatrix();

    glPushMatrix();
    glTranslated(16.2,15.25,-10.25);
    glScaled(0.1,4.5,-2.5);
    cube(1,1,1);
    glPopMatrix();



    //front side

    glPushMatrix();
    glTranslated(6,12.5,-16.1);
    glScaled(4,8,-0.1);
    cube(0,0,0);
    glPopMatrix();

    glPushMatrix();
    glTranslated(6.5,13,-16.2);
    glScaled(3,7,-0.1);
    cube(1,1,1);
    glPopMatrix();



    glPushMatrix();
    glTranslated(11,15,-16.1);
    glScaled(3,3,-0.1);
    cube(0,0,0);
    glPopMatrix();

    glPushMatrix();
    glTranslated(11.2,15.2,-16.2);
    glScaled(2.6,2.6,-0.1);
    cube(1,1,1);
    glPopMatrix();



    glPushMatrix();
    glTranslated(2,15,-16.1);
    glScaled(3,3,-0.1);
    cube(0,0,0);
    glPopMatrix();

    glPushMatrix();
    glTranslated(2.2,15.2,-16.2);
    glScaled(2.6,2.6,-0.1);
    cube(1,1,1);
    glPopMatrix();



    //sader khuti

    glPushMatrix();
    glTranslated(3,10.5,-22);
    glScaled(1,4,-1);
    cube(0,1,1);
    glPopMatrix();


    glPushMatrix();
    glTranslated(22,10.5,-3);
    glScaled(1,4,-1);
    cube(0,1,1);
    glPopMatrix();

    glPushMatrix();
    glTranslated(22,10.5,-22);
    glScaled(1,4,-1);
    cube(0,1,1);
    glPopMatrix();



    //bera

    glPushMatrix();
    glTranslated(4,11.5,-22.5);
    glScaled(18,0.2,-0.2);
    cube(1,0,1);
    glPopMatrix();

    glPushMatrix();
    glTranslated(4,12.5,-22.5);
    glScaled(18,0.2,-0.2);
    cube(1,0,1);
    glPopMatrix();

    glPushMatrix();
    glTranslated(4,13.5,-22.5);
    glScaled(18,0.2,-0.2);
    cube(1,0,1);
    glPopMatrix();



    glPushMatrix();
    glTranslated(45,0,0);
    glRotated(90,0,1,0);
    glPushMatrix();
    glTranslated(4,11.5,-22.5);
    glScaled(18,0.2,-0.2);
    cube(1,0,1);
    glPopMatrix();

    glPushMatrix();
    glTranslated(4,12.5,-22.5);
    glScaled(18,0.2,-0.2);
    cube(1,0,1);
    glPopMatrix();

    glPushMatrix();
    glTranslated(4,13.5,-22.5);
    glScaled(18,0.2,-0.2);
    cube(1,0,1);
    glPopMatrix();

    glPopMatrix();


    glPushMatrix();
    glTranslated(3.5,11.5,-16);
    glScaled(0.2,0.2,-6);
    cube(1,0,1);
    glPopMatrix();

    glPushMatrix();
    glTranslated(3.5,12.5,-16);
    glScaled(0.2,0.2,-6);
    cube(1,0,1);
    glPopMatrix();

    glPushMatrix();
    glTranslated(3.5,13.5,-16);
    glScaled(0.2,0.2,-6);
    cube(1,0,1);
    glPopMatrix();


    glPushMatrix();
    glTranslated(16,11.5,-3.5);
    glScaled(6,0.2,-0.2);
    cube(1,0,1);
    glPopMatrix();

    glPushMatrix();
    glTranslated(16,12.5,-3.5);
    glScaled(6,0.2,-0.2);
    cube(1,0,1);
    glPopMatrix();

    glPushMatrix();
    glTranslated(16,13.5,-3.5);
    glScaled(6,0.2,-0.2);
    cube(1,0,1);
    glPopMatrix();



    ///human

    glPushMatrix();
    glTranslated(10,0.2,-18);
    glRotated(20,0,1,0);
    glScaled(0.4,0.4,0.4);
    human3();
    glPopMatrix();


    glPushMatrix();
    glTranslated(5,10.6,-18);
    glRotated(-20,0,1,0);
    glScaled(0.4,0.4,0.4);
    human3();
    glPopMatrix();

    glPushMatrix();
    glTranslated(17,10.6,-10);
    glRotated(-90,0,1,0);
    glScaled(0.4,0.4,0.4);
    human3();
    glPopMatrix();



}




void point_object()
{

    glPushMatrix();
    glRotated(-90,1,0,0);
    glMaterialfv(GL_BACK, GL_AMBIENT,mat_ambient6);
    glMaterialfv(GL_BACK, GL_DIFFUSE,mat_diffuse6);
    glMaterialfv(GL_BACK, GL_SPECULAR,mat_specular6);
    glMaterialfv(GL_BACK, GL_SHININESS,mat_shininess6);

    glPushMatrix();
    glTranslated(0,0,0);
    wheel();
    glPopMatrix();

    glPushMatrix();
    glTranslated(0,7,5);
    wheel();
    glPopMatrix();

    glPushMatrix();
    glTranslated(-2,14,0);
    wheel();
    glPopMatrix();

    glPushMatrix();
    glTranslated(-5,21,5);
    wheel();
    glPopMatrix();


    glPopMatrix();
}

void point_collection()
{
    glPushMatrix();




    //car7
    if((tZ5<=-108))
    {
        if(tX>=-1)
        {
            if(tX<=1)
            {
                total_point+=1;
            }
        }

    }


    //car2
    if((tZ<=-108))
    {
        if(tX>=3)
        {
            if(tX<=5)
            {
                total_point+=1;
            }
        }

    }

//car1 ok

    if((tZ6<=-108))
    {
        if(tX>=-3)
        {
            if(tX<=-1)
            {
                total_point+=1;
            }
        }

    }


//car4 ok
    if((tZ2<=-108))
    {
        if(tX>=-5)
        {
            if(tX<=-3)
            {
                total_point+=1;
            }
        }

    }



//car6 ok
    if((tZ4<=-108))
    {
        if(tX>=3)
        {
            if(tX<=5)
            {
                total_point+=1;

            }
        }

    }



//car3 ok
    if((tZ1<=-108))
    {
        if(tX>=-1)
        {
            if(tX<=1)
            {
                total_point+=1;
            }
        }

    }


    //car5 ok
    if((tZ3<=-108))
    {
        if(tX>=1)
        {
            if(tX<=3)
            {
                total_point+=1;
            }
        }

    }

    glPopMatrix();

}

void truck(double r, double g, double b)
{

    const double t = glutGet(GLUT_ELAPSED_TIME) / 1000.0;
    const double a = t*180.0;

    //Tires

    glPushMatrix();
    // glCullFace(GL_FRONT);


    glPushMatrix();
    glTranslated(-2,0,-2);
    ownScale(0.2,.2,1);
    glRotated(a,0,0,1);
    car_tire(0, 0, 0);
    glPopMatrix();

    glPushMatrix();
    glTranslated(-2,0,0);
    ownScale(0.2,.2,1);
    glRotated(a,0,0,1);
    car_tire(0, 0, 0);
    glPopMatrix();



    glPushMatrix();
    glTranslated(-0.5,0,-2);
    ownScale(0.2,.2,1);
    glRotated(a,0,0,1);
    car_tire(0, 0, 0);
    glPopMatrix();



    glPushMatrix();
    glTranslated(-0.5,0,0);
    ownScale(0.2,.2,1);
    glRotated(a,0,0,1);
    car_tire(0, 0, 0);
    glPopMatrix();



    glPushMatrix();
    glTranslated(2.5,0,-2);
    ownScale(0.2,.2,1);
    glRotated(a,0,0,1);
    car_tire(0, 0, 0);
    glPopMatrix();

    glPushMatrix();
    glTranslated(2.5,0,0);
    ownScale(0.2,.2,1);
    glRotated(a,0,0,1);
    car_tire(0, 0, 0);
    glPopMatrix();



    //main body

    glPushMatrix();
    glTranslated(-0.6,1.65,-1);
    ownScale(2,1,-1.3);
    cube1(r,g,b);
    glPopMatrix();


    glPushMatrix();
    glTranslated(2.4,1.35,-1);
    ownScale(0.8,0.7,-1.3);
    cube1(0.5,0.4,1);
    glPopMatrix();


    //connector

    glPushMatrix();
    glTranslated(1.5,1.75,-0.5);
    ownScale(0.1,0.2,-0.2);
    cube1(0,1,1);
    glPopMatrix();


    glPushMatrix();
    glTranslated(1.5,1.75,-2);
    ownScale(0.1,0.2,-0.2);
    cube1(0,1,1);
    glPopMatrix();


    glPopMatrix();




}

void obstacle()
{


    glPushMatrix();
    glTranslated(tX-2,tY,tZ);
    glPushMatrix();
    glTranslated(0,1,0);
    glRotated(90,0,1,0);
    glRotated(5,0,0,1);
    glRotated(rotX,1,0,0);
    glRotated(rotY,0,1,0);
    glRotated(rotZ,0,0,1);
    glScaled(0.4,0.4,0.4);
    truck(1,1,0);
    glPopMatrix();
    glPopMatrix();

    glPushMatrix();
    glTranslated(tX+2,tY+1,tZ);
    glPushMatrix();
    glScaled(0.1,0.1,0.1);
    point_object();
    glPopMatrix();
    glPopMatrix();




    glPushMatrix();
    glTranslated(tX-2,tY,tZ1);
    glPushMatrix();
    glTranslated(0,1,0);
    glRotated(90,0,1,0);
    glRotated(5,0,0,1);
    glRotated(rotX,1,0,0);
    glRotated(rotY,0,1,0);
    glRotated(rotZ,0,0,1);
    glScaled(0.4,0.4,0.4);
    truck(1,0.5,1);
    glPopMatrix();
    glPopMatrix();


    glPushMatrix();
    glTranslated(tX-4,tY+1,tZ1);
    glPushMatrix();
    glScaled(0.1,0.1,0.1);
    point_object();
    glPopMatrix();
    glPopMatrix();


    glPushMatrix();
    glTranslated(tX+2,tY,tZ2);
    glPushMatrix();
    glTranslated(0,1,0);
    glRotated(90,0,1,0);
    glRotated(5,0,0,1);
    glRotated(rotX,1,0,0);
    glRotated(rotY,0,1,0);
    glRotated(rotZ,0,0,1);
    glScaled(0.4,0.4,0.4);
    truck(0.2,0.5,1);
    glPopMatrix();
    glPopMatrix();

    glPushMatrix();
    glTranslated(tX,tY+1,tZ2);
    glPushMatrix();
    glScaled(0.1,0.1,0.1);
    point_object();
    glPopMatrix();
    glPopMatrix();




    glPushMatrix();
    glTranslated(tX,tY,tZ3);
    glPushMatrix();
    glTranslated(0,1,0);
    glRotated(90,0,1,0);
    glRotated(5,0,0,1);
    glRotated(rotX,1,0,0);
    glRotated(rotY,0,1,0);
    glRotated(rotZ,0,0,1);
    glScaled(0.4,0.4,0.4);
    truck(0.7,0.5,0.3);
    glPopMatrix();
    glPopMatrix();

    glPushMatrix();
    glTranslated(tX+4,tY+1,tZ3);
    glPushMatrix();
    glScaled(0.1,0.1,0.1);
    point_object();
    glPopMatrix();
    glPopMatrix();



    glPushMatrix();
    glTranslated(tX+2,tY,tZ4);
    glPushMatrix();
    glTranslated(0,1,0);
    glRotated(90,0,1,0);
    glRotated(5,0,0,1);
    glRotated(rotX,1,0,0);
    glRotated(rotY,0,1,0);
    glRotated(rotZ,0,0,1);
    glScaled(0.4,0.4,0.4);
    truck(0.3,0.5,0);
    glPopMatrix();
    glPopMatrix();

    glPushMatrix();
    glTranslated(tX-2,tY+1,tZ4);
    glPushMatrix();
    glScaled(0.1,0.1,0.1);
    point_object();
    glPopMatrix();
    glPopMatrix();



    glPushMatrix();
    glTranslated(tX,tY,tZ5);
    glPushMatrix();
    glTranslated(0,1,0);
    glRotated(90,0,1,0);
    glRotated(5,0,0,1);
    glRotated(rotX,1,0,0);
    glRotated(rotY,0,1,0);
    glRotated(rotZ,0,0,1);
    glScaled(0.4,0.4,0.4);
    truck(1,0.5,0.8);
    glPopMatrix();
    glPopMatrix();

    glPushMatrix();
    glTranslated(tX-4,tY+1,tZ5);
    glPushMatrix();
    glScaled(0.1,0.1,0.1);
    point_object();
    glPopMatrix();
    glPopMatrix();




    glPushMatrix();
    glTranslated(tX+2,tY,tZ6);
    //glRotated(180,0,1,0);
    glPushMatrix();
    glTranslated(0,1,0);
    glRotated(90,0,1,0);
    glRotated(5,0,0,1);
    glRotated(rotX,1,0,0);
    glRotated(rotY,0,1,0);
    glRotated(rotZ,0,0,1);
    glScaled(0.4,0.4,0.4);
    truck(1,0.5,0);
    glPopMatrix();
    glPopMatrix();


    glPushMatrix();
    glTranslated(tX,tY+1,tZ6);
    glPushMatrix();
    glScaled(0.1,0.1,0.1);
    point_object();
    glPopMatrix();
    glPopMatrix();


}


void speedUp()
{


        if(0<=total_point<20)
        {
            speed=0.6;
        }

                if(20<=total_point<40)
        {
            speed=1;
        }

                if(40<=total_point<90)
        {
            speed=2;
        }

                if(90<=total_point<120)
        {
            speed=3;
        }

                if(120<=total_point<200)
        {
            speed=5;
        }

                if(200<=total_point<300)
        {
            speed=10;
        }
}


void draw()
{
    double t = glutGet(GLUT_ELAPSED_TIME) / 1000.0;
    double a = t*90.0;

    TIME = t;
    if(rotX>11)
        rotX=11;
    if(rotX<-11)
        rotX=-11;
    if(rotZ>10)
        rotZ=10;
    if(rotZ<-15)
        rotZ=-15;


    game_over();
    point_collection();
   // speedUp();
    glPushMatrix();
    obstacle();
    glPopMatrix();

    glPushMatrix();
    glTranslated(0,1,0);
    glRotated(90,0,1,0);
    glRotated(5,0,0,1);
    glRotated(anglex,1,0,0);
    glRotated(angley,0,1,0);
    glRotated(anglez,0,0,1);
    glScaled(0.12,0.12,0.12);
    car1();
    glPopMatrix();


    glPushMatrix();
    glTranslated(tX,tY,tZ3);
    glPushMatrix();
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D,3);
    glMaterialfv(GL_FRONT, GL_AMBIENT,mat_ambient3);
    glMaterialfv(GL_FRONT, GL_DIFFUSE,mat_diffuse3);
    glMaterialfv(GL_FRONT, GL_SPECULAR,mat_specular3);
    glMaterialfv(GL_FRONT, GL_SHININESS,mat_shininess3);

    glTranslated(-17.5,0,-21);
    glScaled(10,4,37);
    Cube();
    glPopMatrix();
    glDisable(GL_TEXTURE_2D);

    glPushMatrix();
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D,3);
    glMaterialfv(GL_FRONT, GL_AMBIENT,mat_ambient3);
    glMaterialfv(GL_FRONT, GL_DIFFUSE,mat_diffuse3);
    glMaterialfv(GL_FRONT, GL_SPECULAR,mat_specular3);
    glMaterialfv(GL_FRONT, GL_SHININESS,mat_shininess3);

    glTranslated(7.8,0,-21);
    glScaled(10,4,37);
    Cube();
    glPopMatrix();
    glDisable(GL_TEXTURE_2D);


    glPushMatrix();
    glMaterialfv(GL_FRONT, GL_AMBIENT,mat_ambient4);
    glMaterialfv(GL_FRONT, GL_DIFFUSE,mat_diffuse4);
    glMaterialfv(GL_FRONT, GL_SPECULAR,mat_specular4);
    glMaterialfv(GL_FRONT, GL_SHININESS,mat_shininess4);
    glTranslated(-0.3,0,-20);
    glRotated(-90,1,0,0);
    glRotated(0,0,1,0);
    glRotated(270,0,0,1);
    glScaled(0.7,0.3,0.3);
    tunnelBezier();
    glPopMatrix();
    environment(1);
    glPopMatrix();



    glPushMatrix();
    glTranslated(tX,tY,tZ);
    environment(2);
    glPushMatrix();
    glRotated(-180,0,1,0);
    glTranslated(-18,0.1,0);
    glScaled(0.4,0.4,0.4);
    glCullFace(GL_FRONT);
    bari();
    glPopMatrix();
    glPopMatrix();



    glPushMatrix();
    glTranslated(tX,tY,tZ1);
    stadiumEnv();
    glPopMatrix();




    glPushMatrix();
    glCullFace(GL_BACK);
    glTranslated(tX,tY,tZ2);
    environment(3);

    glPushMatrix();
    glCullFace(GL_FRONT);
    glTranslated(10,0,25);
    glScaled(0.2,0.2,0.2);
    human3();
    glPopMatrix();

    glPushMatrix();
    glRotated(315,0,1,0);
    glTranslated(25,2.5,0);
    glCullFace(GL_FRONT);
    glScaled(0.5,0.5,0.5);
    hospital();
    glPopMatrix();
    glPopMatrix();






    glPushMatrix();
    glTranslated(tX,tY,tZ4);

    environment(5);
    glPopMatrix();

    glPushMatrix();
    glTranslated(tX,tY,tZ5);
    //colony();
    glPushMatrix();
    glTranslated(10,0,-5);
    glScaled(2,2,2);
    colony2();
    glPopMatrix();

    glPushMatrix();
    glTranslated(-15,0,0);
    glScaled(1.5,1.5,1.5);
    colony2();
    glPopMatrix();
    environment(4);
    glPopMatrix();

    glPushMatrix();
    glTranslated(tX,tY,tZ6);
    glPushMatrix();
    glRotated(-90,0,1,0);
    glCullFace(GL_FRONT);
    glTranslated(-10,0.3,18);
    glScaled(0.4,0.4,0.4);
    shopingmall();
    //colony2();
    glPopMatrix();
    environment(2);
    glPopMatrix();

    tZ+=speed;
    tZ1+=speed;
    tZ2+=speed;
    tZ3+=speed;
    tZ4+=speed;
    tZ5+=speed;
    tZ6+=speed;

    if(tZ>=20)
        tZ=-110;
    if(tZ1>=20)
        tZ1=-110;
    if(tZ2>=20)
        tZ2=-110;
    if(tZ3>=20)
        tZ3=-110;
    if(tZ4>=20)
        tZ4=-110;
    if(tZ5>=20)
        tZ5=-110;
    if(tZ6>=20)
        tZ6=-110;


    if(rotX>0)
        rotX-=angleBackFrac;
    if(rotX<0)
        rotX+=angleBackFrac;
    if(rotY>0)
        rotY-=angleBackFrac;
    if(rotY<0)
        rotY+=angleBackFrac;
    if(rotZ>0)
        rotZ-=angleBackFrac;
    if(rotZ<0)
        rotZ+=angleBackFrac;

    //cout<<tX<<" "<<tY<<" "<<tZ<<endl;
    //cout<<rotX<<" "<<rotY<<" "<<rotZ<<endl;

    speed += 0.001;
}


void drawBitmapText(char *str,float x,float y,float z)
{
    char *c;
    glRasterPos3f(x,y+8,z);

    for (c=str; *c != '\0'; c++)
    {
        //glutBitmapCharacter(GLUT_BITMAP_TIMES_ROMAN_10, *c);
        glutBitmapCharacter(GLUT_BITMAP_TIMES_ROMAN_24, *c);
    }
}

void drawStrokeText(char* str,int x,int y,int z)
{
    char *c;
    glPushMatrix();
    glMaterialfv(GL_FRONT, GL_AMBIENT,mat_ambient6);
    glMaterialfv(GL_FRONT, GL_DIFFUSE,mat_diffuse6);
    glMaterialfv(GL_FRONT, GL_SPECULAR,mat_specular6);
    glMaterialfv(GL_FRONT, GL_SHININESS,mat_shininess6);
    glTranslatef(x, y,z);
    glScalef(0.02f,0.02f,z);

    for (c=str; *c != '\0'; c++)
    {
        //glutStrokeCharacter(GLUT_STROKE_ROMAN, *c);
        glutStrokeCharacter(GLUT_STROKE_MONO_ROMAN, *c);
    }
    glPopMatrix();
}

void drawStrokeText2(char* str,int x,int y,int z)
{
    char *c;
    glPushMatrix();
    glTranslatef(x, y+8,z);
    glScalef(0.005f,0.005f,z);

    for (c=str; *c != '\0'; c++)
    {
        glutStrokeCharacter(GLUT_STROKE_ROMAN, *c);
    }
    glPopMatrix();
}
void drawStrokeChar(char c,float x,float y,float z)
{
    glPushMatrix();
    glMaterialfv(GL_FRONT, GL_AMBIENT,mat_ambient6);
    glMaterialfv(GL_FRONT, GL_DIFFUSE,mat_diffuse6);
    glMaterialfv(GL_FRONT, GL_SPECULAR,mat_specular6);
    glMaterialfv(GL_FRONT, GL_SHININESS,mat_shininess6);
    glTranslatef(x, y,z);
    glScalef(0.02f,0.02f,z);
    glutStrokeCharacter(GLUT_STROKE_ROMAN, c);
    glPopMatrix();
}

void light()
{
/*
    GLfloat l_spt[] = {0,0,-1,1};
    GLfloat spt_ct[] = {10};


    GLfloat no_light[]  = { 0.0f, 0.0f, 0.0f, 1.0f };
    GLfloat light_ambient[]  = { 0.3f, 0.3f, 0.3f, 1.0f };
    GLfloat light_diffuse[]  = { 1.0f, 1.0f, 1.0f, 1.0f };
    GLfloat light_specular[] = { 1.0f, 1.0f, 1.0f, 1.0f };

    GLfloat light_position0[] = { 2.0f, 5.0f, 5.0f, 0.0f }; //sun
    GLfloat light_position1[] = { 0, 2,150, 1.0f }; //spot light
   // GLfloat light_position2[] = { 0, 50,2, 1.0f }; //normal light
    GLfloat light_position2[] = { 0, 0.5,4.5, 1.0f }; //normal light


    GLfloat mat_ambient[]    = { 0.7f, 0.7f, 0.7f, 1.0f };
    GLfloat mat_diffuse[]    = { 0.8f, 0.8f, 0.8f, 1.0f };
    GLfloat mat_specular[]   = { 1.0f, 1.0f, 1.0f, 1.0f };
    GLfloat high_shininess[] = { 100.0f };
*/


    if(l1)
    {
        glEnable(GL_LIGHT0);
        if(am)
        {
            glLightfv(GL_LIGHT0, GL_AMBIENT,  light_ambient);
        }
        else
        {
            glLightfv(GL_LIGHT0, GL_AMBIENT,no_light);
        }

        if(df)
        {
            glLightfv(GL_LIGHT0, GL_DIFFUSE,  light_diffuse);
        }
        else
        {
            glLightfv(GL_LIGHT0, GL_DIFFUSE,no_light);
        }

        if(sp)
        {
            glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular);
        }
        else
        {
            glLightfv(GL_LIGHT0, GL_SPECULAR,no_light);
        }
        glLightfv(GL_LIGHT0, GL_POSITION, light_position0);

    }
    else
    {
        glDisable(GL_LIGHT0);
    }


    if(l4)
    {
        glEnable(GL_LIGHT1);
        glLightfv(GL_LIGHT1, GL_SPOT_DIRECTION, l_spt);
        glLightfv(GL_LIGHT1, GL_SPOT_CUTOFF, spt_ct);
        glEnable(GL_LIGHT1);
        if(am)
        {
            glLightfv(GL_LIGHT1, GL_AMBIENT,  light_ambient);
        }
        else
        {
            glLightfv(GL_LIGHT1, GL_AMBIENT,no_light);
        }

        if(df)
        {
            glLightfv(GL_LIGHT1, GL_DIFFUSE,  light_diffuse);
        }
        else
        {
            glLightfv(GL_LIGHT1, GL_DIFFUSE,no_light);
        }

        if(sp)
        {
            glLightfv(GL_LIGHT1, GL_SPECULAR, light_specular);
        }
        else
        {
            glLightfv(GL_LIGHT1, GL_SPECULAR,no_light);
        }
        glLightfv(GL_LIGHT1, GL_POSITION, light_position1);
    }
    else
    {
        glDisable(GL_LIGHT1);
    }


    if(l2)
    {
        glEnable(GL_LIGHT2);
        if(am)
        {
            glLightfv(GL_LIGHT2, GL_AMBIENT,  light_ambient);
        }
        else
        {
            glLightfv(GL_LIGHT2, GL_AMBIENT,no_light);
        }

        if(df)
        {
            glLightfv(GL_LIGHT2, GL_DIFFUSE,  light_diffuse);
        }
        else
        {
            glLightfv(GL_LIGHT2, GL_DIFFUSE,no_light);
        }

        if(sp)
        {
            glLightfv(GL_LIGHT2, GL_SPECULAR, light_specular);
        }
        else
        {
            glLightfv(GL_LIGHT2, GL_SPECULAR,no_light);
        }
        glLightfv(GL_LIGHT2, GL_POSITION, light_position2);
    }
    else
    {
        glDisable(GL_LIGHT2);
    }



}




static void display(void)
{
    const double t = glutGet(GLUT_ELAPSED_TIME) / 1000.0;
    double a = t*90.0;
    double aa=a;

    if(!rot)
    {
        a=0;
    }

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glLoadIdentity();

    gluLookAt(xEye,yEye,zEye, cenX,cenY,cenZ, 0,1,0);


    light();
    if(START)
    {
        glClearColor(1,1,1,1);
        glPushMatrix();
        glTranslated(0,0,0);
        glScaled(zoom,zoom,zoom);
        glRotated(a,0,1,0);
        draw();
        glPopMatrix();


       // drawStrokeText("UP: W, DOWN: S, LEFT: A, RIGHT: D, MAIN MENU: M",-8,1.9,0);
       /* drawStrokeText("Point : ",-20,18,2);
        int temp=total_point, mod,number=0;

        float tmp=0;

        while(temp){
            mod=temp%10;
            drawStrokeChar(mod+48,-10.6+tmp,18,2);
            temp/=10;
            tmp+=2;
        }
*/

    int temp=total_point, mod,number=0;
        while(TIME)
        {
            mod=TIME%10;
            number=number*10+mod;
            TIME/=10;
        }
        float tmp=0;
        while(number)
        {
            mod=number%10;
            drawStrokeChar(mod+48,10+tmp,18,0);
            number/=10;
            tmp+=2;
        }
    }
    else if(START1)
    {

          if(total_point>highest_point)
        {
            highest_point=total_point;
        }
        glClearColor(0,1,1,1);

        int temp=total_point, mod,number=0;

        float tmp=0;

        while(temp){
            mod=temp%10;
            drawStrokeChar(mod+48,10.6-tmp,14,2);
            temp/=10;
            tmp+=2;
        }

        drawStrokeText("Your Score is:",-20,14,2);
        drawStrokeText("Press s for Start Again",-25,0,2);

        drawStrokeText("Hightest Score is:",-25,9.5,2);
       temp=highest_point, mod,number=0;
       tmp=0;

        while(temp){
            mod=temp%10;
            drawStrokeChar(mod+48,15-tmp,8.5,2);
            temp/=10;
            tmp+=2;
        }

    }

    else
    {
        glClearColor(0,0,0,1);
        glPushMatrix();

        // glMaterialfv(GL_FRONT, GL_AMBIENT,mat_ambient4);
        // glMaterialfv(GL_FRONT, GL_DIFFUSE,mat_diffuse4);
        // glMaterialfv(GL_FRONT, GL_SPECULAR,mat_specular4);
        // glMaterialfv(GL_FRONT, GL_SHININESS,mat_shininess4);
        // glRotated(aa,0,1,0);
        glPushMatrix();
        glScaled(0.5,0.5,0.5);
        glTranslated(40,2,0);
        glRotated(anglex,1,0,0);
        //glRotated(angley,0,1,0);
        glRotated(aa,0,1,0);
        glRotated(anglez,0,0,1);

        car1();

        glPopMatrix();
        glPopMatrix();


        glPushMatrix();
        glScaled(0.7,0.7,0.7);
        drawStrokeText("Press s for Start the game",-40,0,0);
        glPopMatrix();

        glPushMatrix();
        glScaled(0.7,0.7,0.7);
        drawStrokeText("INSTRUCTION:",-37,25,0);
        glPopMatrix();



        glPushMatrix();
        glTranslated(-10,7.5,0);
        glScaled(0.6,0.6,0.6);
        glPushMatrix();
        glScaled(0.4,0.4,0.4);
        drawStrokeText("For eye position Movement press < i,t,r,o,p,u >",-55,35,0);
        drawStrokeText("For Car Movement  press < c , m >",-55,30,0);
        drawStrokeText("For Rotated Car  press < x , y , z >",-55,25,0);
        drawStrokeText("For Speed Control  press < h , j , k >",-55,20,0);
        drawStrokeText("For Zoom in Zoom Out  press < + , - >",-55,15,0);
        drawStrokeText("For show only conneted point  press < w >",-55,10,0);
        drawStrokeText("For Lighting on press Respectively < 1 , 2 , 4 >",-55,5,0);
        drawStrokeText("For Ambiuent ,defuse and Specular light Press < 7 , 8 , 9 >",-55,0,0);
        glPopMatrix();
        glPopMatrix();



        glPushMatrix();
        glTranslated(4.5,0.5,0);
        glScaled(0.7,0.7,0.7);
        drawStrokeText("Devoloped By:",0,25,0);
        glPopMatrix();

        glPushMatrix();
        glScaled(0.4,0.4,0.4);
        drawStrokeText("Md Robiul Islam",25,40,0);
        drawStrokeText("Roll:1607107",25,35,0);
        drawStrokeText("KUET,CSE",25,30,0);

        glPopMatrix();


    }

    glutSwapBuffers();
}


static void key(unsigned char key, int x, int y)
{
    float frac = 0.2;
    float rotFrac = 1;
    switch (key)
    {
    case 27 :
    case 'q':
        exit(0);
        break;
    case '+':
        zoom+=0.05;
        break;
    case '-':
        zoom-=0.05;

    case 'c':
        tX+=frac;
        //rotX-=rotFrac*3;
        //rotY+=rotFrac/2;
        break;
    case 'm':
        tX-=frac;
       // rotX+=rotFrac*3;
       // rotY-=rotFrac/2;
        break;

    case 's':
        START=true;
        START1=false;
        total_point=0;
        speed=0.6;
        break;
    case 'S':
        START=false;
        START1=false;
        total_point=0;
        break;


    //eye position
    case 't':
        xEye+=0.1;
        break;
    case 'i':
        xEye-=0.1;
        break;

    case 'r':
        yEye+=0.1;
        break;

    case '0':
        yEye-=0.1;
        break;

    case 'p':
        zEye-=.1;
        break;
    case 'u':
        zEye+=.1;
        break;

    //spedd control
    case 'h':
        speed +=0.2;
        rt=true;
        break;
    case 'j':
        speed -=0.2;
        rt=true;
        break;
    case 'k':
        speed =0;
        rt=false;
        break;

    //light operation
    case '1':
        l1=l1-1;
        break;
    case '2':
        l2=l2-1;
        break;
    case '3':
        l3=l3-1;
        break;

    case '4':
        l4=l4-1;
        break;


    case '7':
        am =am-1;
        break;
    case '8':
        df=df-1;
        break;
    case '9':
        sp=sp-1;
        break;




    case 'A':
        animat=!animat;
        break;

    case 'g':
        shcpt=!shcpt;
        break;

    case 'w':
        wired=!wired;
        break;

    case 'x':
        anglex = ( anglex + 3 ) % 360;
        break;
    case 'X':
        anglex = ( anglex - 3 ) % 360;
        break;

    case 'y':
        angley = ( angley + 3 ) % 360;
        break;
    case 'Y':
        angley = ( angley - 3 ) % 360;
        break;

    case 'z':
        anglez = ( anglez + 3 ) % 360;
        break;
    case 'Z':
        anglez = ( anglez - 3 ) % 360;
        break;

    case 'V':


        break;

    }

    glutPostRedisplay();
}


static void idle(void)
{
    glutPostRedisplay();
}

void LoadTexture(const char*filename)
{
    glGenTextures(1, &ID);
    glBindTexture(GL_TEXTURE_2D, ID);
    glPixelStorei(GL_UNPACK_ALIGNMENT, ID);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    BmpLoader bl(filename);
    gluBuild2DMipmaps(GL_TEXTURE_2D, GL_RGB, bl.iWidth, bl.iHeight, GL_RGB, GL_UNSIGNED_BYTE, bl.textureData );
}

int main(int argc, char *argv[])
{
    glutInit(&argc, argv);
    glutInitWindowPosition(0,0);
    glutInitWindowSize(1366,720);
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH | GLUT_RGBA);

    glutCreateWindow("GLUT Shapes");
    LoadTexture("C:\\Users\\robiu\\Desktop\\Thesis\\1607107_Car_Racing\\field.bmp");//field
    LoadTexture("C:\\Users\\robiu\\Desktop\\Thesis\\1607107_Car_Racing\\lio.bmp");
    LoadTexture("C:\\Users\\robiu\\Desktop\\Thesis\\1607107_Car_Racing\\wall1.bmp");
    LoadTexture("C:\\Users\\robiu\\Desktop\\Thesis\\1607107_Car_Racing\\road11.bmp");
    LoadTexture("C:\\Users\\robiu\\Desktop\\Thesis\\1607107_Car_Racing\\gras.bmp");
    LoadTexture("C:\\Users\\robiu\\Desktop\\Thesis\\1607107_Car_Racing\\hos2.bmp");

    glutReshapeFunc(resize);
    glutDisplayFunc(display);
    glutKeyboardFunc(key);
    glutMouseFunc(processMouse);
    glutIdleFunc(idle);

    //PlaySound("starwars.wav", NULL, SND_ASYNC|SND_FILENAME|SND_LOOP);

    glClearColor(1,1,1,1);
    glEnable(GL_CULL_FACE);
    glCullFace(GL_BACK);

    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LESS);

    // glEnable(GL_LIGHT0);
    glEnable(GL_NORMALIZE);
    glEnable(GL_LIGHTING);


    glMaterialfv(GL_FRONT, GL_AMBIENT,   mat_ambient);
    glMaterialfv(GL_FRONT, GL_DIFFUSE,   mat_diffuse);
    glMaterialfv(GL_FRONT, GL_SPECULAR,  mat_specular);
    glMaterialfv(GL_FRONT, GL_SHININESS, high_shininess);
    glutMainLoop();

    return EXIT_SUCCESS;
}
