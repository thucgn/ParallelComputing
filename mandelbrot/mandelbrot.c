
#include <stdlib.h>
#include <GL/glut.h>
//#include <openGL/gl.h>
//#include <openGL/glu.h>
//#include <GLUT/glut.h>

#ifdef USE_OPENMP
#include <omp.h>
#endif

#ifdef PRINT_TIME
#include <stdio.h>
#include <stdarg.h>
#include <sys/time.h>
#define LOG(format, ...) printf("[%d:%s] " format "\n", \
        __LINE__, __FUNCTION__, ##__VA_ARGS__)
#define MARK_TIME(t) gettimeofday(&t, NULL)
#define DIFF_TIME(tt, ts) (((tt).tv_sec-(ts).tv_sec) + ((tt).tv_usec*1e-6 - (ts).tv_usec*1e-6))
typedef struct timeval TIME_T;
#endif


/* Defaut data via command line */
/* Can enter other values via command line arguments */

#define CENTERX -0.5
#define CENTERY 0.5
#define HEIGHT 0.5
#define WIDTH 0.5
#define MAX_ITER 100

/* N x M array to be generated */

#define N 500
#define M 500

float height = HEIGHT; /* size of window in complex plane */
float width = WIDTH;
float cx = CENTERX; /* center of window in complex plane */
float cy = CENTERY; 
int max = MAX_ITER; /* number of interations per point */

int n=N;
int m=M;

/* Use unsigned bytes for image */

GLubyte image[N][M];

/* Complex data type and complex add, mult, and magnitude functions */
/* Probably not worth overhead */

typedef float complex[2];

void add(complex a, complex b, complex p)
{
    p[0]=a[0]+b[0];
    p[1]=a[1]+b[1];
}

void mult(complex a, complex b, complex p)
{
    p[0]=a[0]*b[0]-a[1]*b[1];
    p[1]=a[0]*b[1]+a[1]*b[0];
}

float mag2(complex a)
{
    return(a[0]*a[0]+a[1]*a[1]);
}

void form(float a, float b, complex p)
{
    p[0]=a;
    p[1]=b;
}

void display()
{
    glClear(GL_COLOR_BUFFER_BIT);
    glDrawPixels(n,m,GL_COLOR_INDEX, GL_UNSIGNED_BYTE, image);
    glFlush();
}


void myReshape(int w, int h)
{
    glViewport(0, 0, w, h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    if (w <= h)
    gluOrtho2D(0.0, 0.0, (GLfloat) n, (GLfloat) m* (GLfloat) h / (GLfloat) w);
    else
    gluOrtho2D(0.0, 0.0, (GLfloat) n * (GLfloat) w / (GLfloat) h,(GLfloat) m);
    glMatrixMode(GL_MODELVIEW);
    display();
}

void myinit()
{
    float redmap[256], greenmap[256],bluemap[256];
    int i;

    glClearColor (1.0, 1.0, 1.0, 1.0);
    gluOrtho2D(0.0, 0.0, (GLfloat) n, (GLfloat) m);

/* Define pseudocolor maps, ramps for red and blue,
   random for green */

/*#ifdef PRINT_TIME
    TIME_T ts, tt;
    MARK_TIME(ts);
#endif

#ifdef USE_OPENMP
#pragma omp parallel for schedule(static) num_threads(THREADS)
#endif*/
    for(i=0;i<256;i++) 
    {
         redmap[i]=i/255.;
         greenmap[i]=drand48();
         bluemap[i]=1.0-i/255.;
    }

/*#ifdef PRINT_TIME
    MARK_TIME(tt);
    LOG("init for time: %.5f", DIFF_TIME(tt, ts));
#endif*/

    glPixelMapfv(GL_PIXEL_MAP_I_TO_R, 256, redmap);
    glPixelMapfv(GL_PIXEL_MAP_I_TO_G, 256, greenmap);
    glPixelMapfv(GL_PIXEL_MAP_I_TO_B, 256, bluemap); 
}


main(int argc, char *argv[])
{
    int i, j, k;
    float x, y, v;
    complex c0, c, d;

    if(argc>1) cx = atof(argv[1]); /* center x */
    if(argc>2) cy = atof(argv[2]);  /* center y */
    if(argc>3) height=width=atof(argv[3]); /* rectangle height and width */
    if(argc>4) max=atoi(argv[4]); /* maximum iterations */

#ifdef PRINT_TIME
    TIME_T t1, t2, t3;
    MARK_TIME(t1);
#endif


#ifdef USE_OPENMP
#pragma omp parallel for private(j,k,x,y,v,c0,c,d) schedule(static) num_threads(THREADS)
#endif
    for (i=0; i<n; i++) 
        for(j=0; j<m; j++) 
        {
    
    /* starting point */
    
            x= i *(width/(n-1)) + cx -width/2;
            y= j *(height/(m-1)) + cy -height/2;
    
            form(0,0,c);
            form(x,y,c0);
    
    /* complex iteration */
    
            for(k=0; k<max; k++)
            {
                mult(c,c,d);
                add(d,c0,c);
                v=mag2(c);
                if(v>4.0) break; /* assume not in set if mag > 4 */
            }
    
    /* assign gray level to point based on its magnitude */
            if(v>1.0) v=1.0; /* clamp if > 1 */
                image[i][j]=255*v;
        }

#ifdef PRINT_TIME
    MARK_TIME(t2);
#endif


    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB );
    glutInitWindowSize(N, M);
    glutCreateWindow("mandlebrot");

    myinit();

#ifdef PRINT_TIME
    MARK_TIME(t3);
    LOG("time of for(main) %.5lf, time total: %.5f", DIFF_TIME(t2, t1), DIFF_TIME(t3, t1));
#endif

    glutReshapeFunc(myReshape);
    glutDisplayFunc(display);

    glutMainLoop();


}
