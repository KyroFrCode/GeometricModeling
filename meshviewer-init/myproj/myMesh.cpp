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
	for (myVertex* v : vertices) {
		v->normal->clear();
	}

	/*(Normals are calculated by doing a cross product of vectors where vectors are the difference 
	between the sources points of the half-edges from the triangle represented by the face).*/
	
	//For each face, Compute face normals and added to each vertex of the face

	for (myFace* f : faces){

		//Compute the normal of the face
		f->computeNormal();

		//Retreive the halfedges of the face
		myHalfedge* he = f->adjacent_halfedge;

		//Added the compute normal to all the halfedges
		do {
			he->source->normal->dX += f->normal->dX;
			he->source->normal->dY += f->normal->dY;
			he->source->normal->dY += f->normal->dZ;
			he = he->next;

		} while (he != f->adjacent_halfedge);
	}

	//Normalize all vertices normals
	for (myVertex* v : vertices) {
		v->normal->normalize();
	}
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
	//Keep the original faces for compute of triangulation on the object
	std::vector<myFace*> original_faces = faces;
	//Clear of the faces of the object
	faces.clear();

	//For each face, do a triangulation face of the original face
	for (myFace* f : original_faces){

		//if face or adjacent_halfedge of the face is null jump to the next face
		if (!f || !f->adjacent_halfedge) {
			continue;
		}


		std::vector<myHalfedge*> edges;
		myHalfedge* current = f->adjacent_halfedge;

		//Add all the edges to the vector "edges" for computation of the triangle associate to the face later
		do {
			edges.push_back(current);
			current = current->next;
		} while (current != f->adjacent_halfedge);

		//Counter of edges present in the object
		int counter = edges.size();

		//If edges in the face are less than three return a warning and jump to the next face 
		// (lack of information for the creation of a triangle associate to the current face)
		if (counter < 3) {

			std::cerr << "Warning : Face with less than three edges." << std::endl;
			continue;
		}

		if (counter == 3) {
			faces.push_back(f);
			continue;
		}

		myVertex* v = f->adjacent_halfedge->source;
		myHalfedge* prev_link = nullptr;

		//Calculating the triangulation of the current face f and creating the triangular face
		for (int i = 1;i <= counter - 2;i++) {

			//Linking each halfedges to their relative twins
			myHalfedge* he;

			if (i == 1){
				he = edges[0];
			}
			else{
				he = prev_link;
			}

			myHalfedge* he2 = edges[i];
			myHalfedge* he3;

			if (i == counter - 2){
				he3 = edges.back();
			}
			else {
				
				he3 = new myHalfedge();
				myHalfedge* twin = new myHalfedge();

				he3->source = edges[i + 1]->source;
				twin->source = v;
				
				he3->twin = twin;
				twin->twin = he3;
				
				halfedges.push_back(he3);
				halfedges.push_back(twin);

				prev_link = twin;
			}

			//Creation of the triangle face
			myFace* face = new myFace();

			//Associating the adjacent face of halfedges to the face created previously
			he->adjacent_face = face;
			he2->adjacent_face = face;
			he3->adjacent_face = face;

			//Linking each halfedges next between them
			he->next = he2;
			he2->next = he3;
			he3->next = he;

			//Linking each halfedges prev between them
			he->prev = he3;
			he2->prev = he;
			he3->prev = he2;

			//Linking the adjacent halfedge of the face
			face->adjacent_halfedge = he;

			//Store the new face and the new halfedges to the vector "faces" & "halfedges"
			faces.push_back(face);
		}
	}
}

//return false if already triangle, true othewise.
bool myMesh::triangulate(myFace *f)
{

	//if face or adjacent_halfedge of the face is null jump to the next face
	if (!f || !f->adjacent_halfedge) {
		return true;
	}

	std::vector<myHalfedge*> edges;
	myHalfedge* current = f->adjacent_halfedge;

	//Add all the edges to the vector "edges" for computation of the triangle associate to the face later
	do {
		edges.push_back(current);
		current = current->next;
	} while (current != f->adjacent_halfedge);

	//Counter of edges present in the object
	int counter = edges.size();

	//If edges in the face are less than three return a warning and jump to the next face 
	// (lack of information for the creation of a triangle associate to the current face)
	if (counter < 3) {

		std::cerr << "Warning : Face with less than three edges." << std::endl;
		return true;
	}

	if (counter > 3) {

		std::cerr << "Warning : Face with more than three edges." << std::endl;
		return true;
	}

	return false;
}

