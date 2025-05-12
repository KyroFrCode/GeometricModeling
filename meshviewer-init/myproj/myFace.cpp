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
	//Check if structure is not null
	if (!this->adjacent_halfedge || !this->adjacent_halfedge->next || !this->adjacent_halfedge->next->next){
		return;
	}

	//Retreived the three edges of the triangle
	myHalfedge* he = this->adjacent_halfedge;
	myHalfedge* he1 = he->next;
	myHalfedge* he2 = he1->next;

	//creating the vector
	myVector3D v1(he1->source->point->X - he->source->point->X, he1->source->point->Y - he->source->point->Y, he1->source->point->Z - he->source->point->Z);
	myVector3D v2(he2->source->point->X - he1->source->point->X, he2->source->point->Y - he1->source->point->Y, he2->source->point->Z - he1->source->point->Z);

	*this->normal = v1.crossproduct(v2);
	normal->normalize();
}
