#include "myFace.h"
#include "myvector3d.h"
#include "myHalfedge.h"
#include "myVertex.h"
#include <GL/glew.h>

myFace::myFace(void)
{
	adjacent_halfedge = NULL;
	normal = new myVector3D(1.0, 1.0, 1.0);
}

myFace::~myFace(void)
{
	if (normal) delete normal;
}

void myFace::computeNormal()
{
    if (!this->adjacent_halfedge) return;

    myVector3D n(0.0, 0.0, 0.0);
    myHalfedge* start = this->adjacent_halfedge;
    myHalfedge* he = start;

    do {
        myPoint3D* p0 = he->source->point;
        myPoint3D* p1 = he->next->source->point;
        n.dX += (p0->Y - p1->Y) * (p0->Z + p1->Z);
        n.dY += (p0->Z - p1->Z) * (p0->X + p1->X);
        n.dZ += (p0->X - p1->X) * (p0->Y + p1->Y);
        he = he->next;
    } while (he != start);

    *this->normal = n;
    this->normal->normalize();
}
