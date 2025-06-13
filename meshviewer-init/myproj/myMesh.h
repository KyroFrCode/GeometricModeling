#pragma once
#include "myFace.h"
#include "myHalfedge.h"
#include "myVertex.h"
#include <vector>
#include <string>
#include <map>

class myMesh
{
public:
	std::vector<myVertex *> vertices;
	std::vector<myHalfedge *> halfedges;
	std::vector<myFace *> faces;
	std::string name;

	void myMesh::checkMesh();
	bool myMesh::check_mesh_vertex(myVertex* v);
	bool myMesh::check_halfedges(myHalfedge* e);
	bool myMesh::check_face(myHalfedge* e);
	bool myMesh::check_twin(myHalfedge* e);

	bool readFile(std::string filename);
	void computeNormals();
	void normalize();

	void subdivisionCatmullClark();

	void splitFaceTRIS(myFace *, myPoint3D *);

	void splitEdge(myHalfedge *, myPoint3D *);
	void splitFaceQUADS(myFace *, myPoint3D *);

	void triangulate();
	bool triangulate(myFace *);

	void simplification();

	void clear();

	myMesh(void);
	~myMesh(void);
};

