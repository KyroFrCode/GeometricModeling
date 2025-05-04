#include "myMesh.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <map>
#include <utility>
#include <GL/glew.h>
#include <stdexcept>
#include "myvector3d.h"

using namespace std;

myMesh::myMesh(void)
{
	/**** TODO ****/
}


myMesh::~myMesh(void)
{
	/**** TODO ****/
}

void myMesh::clear()
{
	for (unsigned int i = 0; i < vertices.size(); i++) if (vertices[i]) delete vertices[i];
	for (unsigned int i = 0; i < halfedges.size(); i++) if (halfedges[i]) delete halfedges[i];
	for (unsigned int i = 0; i < faces.size(); i++) if (faces[i]) delete faces[i];

	vector<myVertex *> empty_vertices;    vertices.swap(empty_vertices);
	vector<myHalfedge *> empty_halfedges; halfedges.swap(empty_halfedges);
	vector<myFace *> empty_faces;         faces.swap(empty_faces);
}

void myMesh::checkMesh()
{
	vector<myHalfedge*>::iterator halfe;

	//Counting the number of errors present on the object analysed
	size_t error_count = 0;

	//Looping through all halfedges to check if the mesh structure is correct
	for (halfe = halfedges.begin(); halfe != halfedges.end(); halfe++)
	{
		try {
			check_mesh_vertex((*halfe)->source); //Check for vertex structure
			check_halfedges((*halfe)); //Check for halfdeges structure
			check_face((*halfe));//Check if the adjacent face to all halfedges are the same
			check_twin((*halfe));//Check if the twin are correct for all halfedges

		}
		catch (const std::exception& e) {
			std::cerr << "Mesh check error: " << e.what() << std::endl;
			++error_count;
		}
	}

	//Display in the terminal if passed the mesh structure test 
	if (error_count == 0) {
		std::cout << "Mesh integrity check passed: All edges and connections are valid`\n";
	}
	else {
		std::cout << "Mesh integrity check failed: " << error_count << " errors found\n";
	}

}

//Check for vertex mesh structure
bool myMesh::check_mesh_vertex(myVertex* v){
	
	//Check if the vertex v is not NULL
	if (v == nullptr) {
		throw std::runtime_error("Error: vertex is nullptr");
		return false;
	}

	//Check if the origin of the vertex v isn't empty
	if (v->originof == nullptr){
		throw std::runtime_error("Error: vertex " + std::to_string(v->index) + "originof is nullptr");
		return false;
	}

	//Check if the source of the halfedge of the vertex is not NULL
	if (v->originof->source == nullptr){
		throw std::runtime_error("Error: vertex" + std::to_string(v->index) + "halfedge source is nullptr");
		return false;
	}

	//Check if the current halfedge and the source halfedge from the originof of vertex is the same or not
	if (v->originof->source != v){
		throw std::runtime_error("Error: vertex" + std::to_string(v->index) + "and the halfedge source are not the same");
		return false;
	}

	return true;
}

//Check for halfedges connection 
bool myMesh::check_halfedges(myHalfedge* e){

	//Check if the halfedge is nullptr
	if (e == nullptr) {
		throw std::runtime_error("Error: halfedge is nullptr");
		return false;
	}

	//Check if the halfedge next and next->next is nullptr
	if ((e->next == nullptr) || (e->next->next == nullptr)) {
		throw std::runtime_error("Error: halfedge next or next->next is nullptr");
		return false;
	}

	//Check if the halfedge previous is nullptr
	if (e->prev == nullptr) {
		throw std::runtime_error("Error: halfedge prev is nullptr");
		return false;
	}

	//Check if halfedge is interconnected correctly
	if (e->next->prev != e) {
		throw std::runtime_error("Error: halfedge next->prev isn't the same as original");
		return false;
	}

	//Alt: Check if the halfedge is interconnected corretly
	if (e->prev->next != e) {
		throw std::runtime_error("Error: halfedge prev->next isn't the same as original");
		return false;
	}

	return true;
}

bool myMesh::check_face(myHalfedge* e){

	//Check if the halfedge e is nullptr
	if (e == nullptr) {
		throw std::runtime_error("Error: halfedge is nullptr");
	}

	//Check if adjacent_face is NULL
	if (e->adjacent_face == nullptr) {
		throw std::runtime_error("Error: adjacent face is nullptr");
		return false;
	}

	myFace* checked_face = e->adjacent_face;
	myHalfedge* tmp = e->next;

	//Run on all over the halfedges connected to the face
	while(tmp != e){

		if (e->adjacent_face == nullptr){
			throw std::runtime_error("Error: adjacent face is nullptr for one of the halfedges link to the face");
			return false;
		}

		if (e->adjacent_face != checked_face){
			throw std::runtime_error("Error: adjacent face and checked face not the same for one of the halfedges link to the face");
			return false;
		}

		tmp = tmp->next;
	}

	return true;
	
}

bool myMesh::check_twin(myHalfedge* e){

	//Check if halfedge e is nullptr
	if (e == nullptr) {
		throw std::runtime_error("Error: halfedge is nullptr");
	}

	//Check if the twin is NULL from the halfedge given
	if (e->twin == nullptr) {
		throw std::runtime_error("Error: the first twin from the halfedge is nullptr");
		return false;
	}

	//Check if the twin are connected correctly in the structure
	if (e->twin->twin != e) {
		throw std::runtime_error("Error: twin's twin does not go back to the original halfedge e");
	}
	
	return true;
}

bool myMesh::readFile(std::string filename)
{
	string s, t, u;
	vector<int> faceids;
	myHalfedge **hedges;

	ifstream fin(filename);
	if (!fin.is_open()) {
		cout << "Unable to open file!\n";
		return false;
	}
	name = filename;

	map<pair<int, int>, myHalfedge *> twin_map;
	map<pair<int, int>, myHalfedge *>::iterator it;

	while (getline(fin, s))
	{
		stringstream myline(s);
		myline >> t;
		if (t == "g") {}
		else if (t == "v")
		{
			float x, y, z;
			myline >> x >> y >> z;
			cout << "v " << x << " " << y << " " << z << endl;

			myVertex* v = new myVertex();
			v->point = new myPoint3D(x, y, z);
			vertices.push_back(v);
		}
		else if (t == "mtllib") {}
		else if (t == "usemtl") {}
		else if (t == "s") {}
		else if (t == "f")
		{
			faceids.clear();
			while (myline >> u) // read indices of vertices from a face into a container - it helps to access them later 
				faceids.push_back(atoi((u.substr(0, u.find("/"))).c_str()) - 1);
			if (faceids.size() < 3) // ignore degenerate faces
				continue;

			hedges = new myHalfedge * [faceids.size()]; // allocate the array for storing pointers to half-edges
			for (unsigned int i = 0; i < faceids.size(); i++)
				hedges[i] = new myHalfedge(); // pre-allocate new half-edges

			myFace* f = new myFace(); // allocate the new face
			f->adjacent_halfedge = hedges[0]; // connect the face with incident edge

			for (unsigned int i = 0; i < faceids.size(); i++)
			{
				int iplusone = (i + 1) % faceids.size();
				int iminusone = (i - 1 + faceids.size()) % faceids.size();
				
				// Connect next and prev
				hedges[i]->next = hedges[iplusone];
				hedges[i]->prev = hedges[iminusone];

				// Set face
				hedges[i]->adjacent_face = f;

				// Set origin
				hedges[i]->source = vertices[faceids[i]];
				if (vertices[faceids[i]]->originof == nullptr)
					vertices[faceids[i]]->originof = hedges[i];

				// Twin edge search and setup
				pair<int, int> this_edge(faceids[i], faceids[iplusone]);
				pair<int, int> twin_edge(faceids[iplusone], faceids[i]);

				it = twin_map.find(twin_edge);
				if (it != twin_map.end())
				{
					// Twin exists, connect both
					hedges[i]->twin = it->second;
					it->second->twin = hedges[i];
					twin_map.erase(it);
				}
				else
				{
					twin_map[this_edge] = hedges[i];
				}

				// Store edge
				halfedges.push_back(hedges[i]);

			}

			// push faces to faces in myMesh
			faces.push_back(f);

			cout << "f";
			for (int idx : faceids)
				cout << " " << (idx + 1);
			cout << endl;

			delete[] hedges;

			cout << "f"; 
			while (myline >> u) cout << " " << atoi((u.substr(0, u.find("/"))).c_str());
			cout << endl;
		}
	}

	checkMesh();
	normalize();

	return true;
}


void myMesh::computeNormals(){

	//Reset all vertices normal
	for (myVertex* v : vertices){
		v->normal = new myVector3D(0, 0, 0);
	}

	/*(Normals are calculated by doing a cross product of vectors where vectors are the difference 
	between the sources points of the half-edges from the triangle represented by the face).*/
	
	//For each face, Compute face normals and added to each vertex of the face
	for (myFace* f : faces){

		//Retreive the halfedges of the face
		myHalfedge* he = f->adjacent_halfedge;
		myHalfedge* he1 = he->next;
		myHalfedge* he2 = he1->next;

		//Retreive the source point of the halfedges
		myPoint3D* p = he->source->point;
		myPoint3D* p1 = he1->source->point;
		myPoint3D* p2 = he2->source->point;

		//Computing the two vectors by using source point of the halfedges
		myVector3D vect1 = *p1 - *p;
		myVector3D vect2 = *p2 - *p;
		
		//Computing the cross product of the vectors and normalizing it
		myVector3D Normalface = cross_prod(vect1,vect2);
		Normalface.normalize();

		//Add compute normal face to each vertex normal
		*(he->source->normal) += Normalface;
		*(he1->source->normal) += Normalface;
		*(he2->source->normal) += Normalface;
		
		//Store the compute normal face in the face f
		f->normal = new myVector3D(Normalface);

		//Normalize all vertices normals
		for (myVertex* v : vertices){

			if (v->normal){
				v->normal->normalize();
			}
		}
	}
}

myVector3D myMesh::cross_prod(myVector3D vect1, myVector3D vect2){
	
	return myVector3D(vect1.dY * vect2.dZ - vect1.dZ * vect2.dY, vect1.dZ * vect2.dX - vect1.dX * vect2.dZ, vect1.dX * vect2.dY - vect1.dY * vect2.dX);
}

void myMesh::normalize()
{
	if (vertices.size() < 1) return;

	int tmpxmin = 0, tmpymin = 0, tmpzmin = 0, tmpxmax = 0, tmpymax = 0, tmpzmax = 0;

	for (unsigned int i = 0; i < vertices.size(); i++) {
		if (vertices[i]->point->X < vertices[tmpxmin]->point->X) tmpxmin = i;
		if (vertices[i]->point->X > vertices[tmpxmax]->point->X) tmpxmax = i;

		if (vertices[i]->point->Y < vertices[tmpymin]->point->Y) tmpymin = i;
		if (vertices[i]->point->Y > vertices[tmpymax]->point->Y) tmpymax = i;

		if (vertices[i]->point->Z < vertices[tmpzmin]->point->Z) tmpzmin = i;
		if (vertices[i]->point->Z > vertices[tmpzmax]->point->Z) tmpzmax = i;
	}

	double xmin = vertices[tmpxmin]->point->X, xmax = vertices[tmpxmax]->point->X,
		ymin = vertices[tmpymin]->point->Y, ymax = vertices[tmpymax]->point->Y,
		zmin = vertices[tmpzmin]->point->Z, zmax = vertices[tmpzmax]->point->Z;

	double scale = (xmax - xmin) > (ymax - ymin) ? (xmax - xmin) : (ymax - ymin);
	scale = scale > (zmax - zmin) ? scale : (zmax - zmin);

	for (unsigned int i = 0; i < vertices.size(); i++) {
		vertices[i]->point->X -= (xmax + xmin) / 2;
		vertices[i]->point->Y -= (ymax + ymin) / 2;
		vertices[i]->point->Z -= (zmax + zmin) / 2;

		vertices[i]->point->X /= scale;
		vertices[i]->point->Y /= scale;
		vertices[i]->point->Z /= scale;
	}
}


void myMesh::splitFaceTRIS(myFace *f, myPoint3D *p)
{
	/**** TODO ****/
}

void myMesh::splitEdge(myHalfedge *e1, myPoint3D *p)
{

	/**** TODO ****/
}

void myMesh::splitFaceQUADS(myFace *f, myPoint3D *p)
{
	/**** TODO ****/
}


void myMesh::subdivisionCatmullClark()
{
	/**** TODO ****/
}


void myMesh::triangulate()
{
	/**** TODO ****/
}

//return false if already triangle, true othewise.
bool myMesh::triangulate(myFace *f)
{
	/**** TODO ****/
	return false;
}

