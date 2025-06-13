#include "myVertex.h"
#include "myvector3d.h"
#include "myHalfedge.h"
#include "myFace.h"

myVertex::myVertex(void)
{
	point = NULL;
	originof = NULL;
	normal = new myVector3D(1.0,1.0,1.0);
}

myVertex::~myVertex(void)
{
	if (normal) delete normal;
}

void myVertex::computeNormal()
{
    myVector3D n(0.0, 0.0, 0.0);
    myHalfedge* h = originof;
    myHalfedge* h0 = h;
    int count = 0;

    do {
        if (h->adjacent_face && h->adjacent_face->normal) {
            n.dX += h->adjacent_face->normal->dX;
            n.dY += h->adjacent_face->normal->dY;
            n.dZ += h->adjacent_face->normal->dZ;
            count++;
        }
        if (!h->twin || !h->twin->next) {
            break;
        }
        h = h->twin->next;
        count++;

    } while (h != h0);

    if (count > 0) {
        normal->dX = n.dX / count;
        normal->dY = n.dY / count;
        normal->dZ = n.dZ / count;
        normal->normalize();
    }
}
