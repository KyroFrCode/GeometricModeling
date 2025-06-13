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

		else if (t == "f")
		{
			faceids.clear();

			while (myline >> u) // read indices of vertices from a face into a container, it helps to access them later 
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

				vertices[faceids[i]]->originof = hedges[i];

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

			delete[] hedges;
		}
	}

	checkMesh();
	normalize();

	return true;
}


void myMesh::computeNormals(){

	/*(Normals are calculated by doing a cross product of vectors where vectors are the difference
	between the sources points of the half-edges from the triangle represented by the face).*/

	//For each face, Compute face normals and add it to the normal of each vertex of the face
	for (myFace* f : faces) {

		//Compute the normal of the face
		f->computeNormal();
	}

	//Normalize all vertex normals
	for (myVertex* v : vertices) {
		v->computeNormal();
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

void myMesh::subdivisionCatmullClark(){

	//Vectors to save the new mesh of the object after the clatmull clark algorithm
	std::vector<myFace*> new_faces;
	std::vector<myHalfedge*> new_halfedges;
	std::vector<myVertex*> new_vertices;
	std::map<myFace*, myVertex*> facePoints;
	std::map<myHalfedge*, myVertex*> edgePoints;

	//Compute face points
	for (myFace* f : faces){

		myVertex* v = new myVertex();
		v->point = new myPoint3D();
		myHalfedge* e = f->adjacent_halfedge;

		int count = 0;

		do {

			*(v->point) += *(e->source->point);
			count++;
			e = e->next;

		} while (e != f->adjacent_halfedge);

		*(v->point) /= count;
		facePoints[f] = v;
	}

	//Compute edge points
	for (myHalfedge* e : halfedges) {

		if (edgePoints.find(e) != edgePoints.end()) {
			continue;
		}

		myVertex* v = new myVertex();
		v->point = new myPoint3D();

		myPoint3D p1 = *(e->source->point);
		myPoint3D p2 = *(e->twin->source->point);

		if (e->adjacent_face && e->twin->adjacent_face) {

			myPoint3D f1 = *(facePoints[e->adjacent_face]->point);
			myPoint3D f2 = *(facePoints[e->twin->adjacent_face]->point);
			*(v->point) = (p1 + p2 + f1 + f2) / 4;
		}

		edgePoints[e] = v;
		edgePoints[e->twin] = v;
	}

	// Update vertex position
	for (myVertex* v : vertices) {

		myPoint3D F(0, 0, 0), R(0, 0, 0), P = *(v->point);
		int n = 0;
		myHalfedge* e = v->originof;
		myHalfedge* start = e;

		do {

			F += *(facePoints[e->adjacent_face]->point);
			R += (*(e->source->point) + *(e->twin->source->point)) / 2;
			e = e->twin->next;
			n++;
		} while (e != start);

		F /= n;
		R /= n;
		*(v->point) = (F + R * 2 + P * (n - 3)) / n;
	}

	// Create new faces of the mesh after catmull clark algorithm
	for (myFace* f : faces) {

		myHalfedge* start = f->adjacent_halfedge;
		myHalfedge* e = start;
		do {

			myFace* new_f = new myFace();
			myVertex* v0 = e->source;
			myVertex* v1 = edgePoints[e];
			myVertex* v2 = facePoints[f];
			myVertex* v3 = edgePoints[e->prev];

			myHalfedge* he0 = new myHalfedge();
			myHalfedge* he1 = new myHalfedge();
			myHalfedge* he2 = new myHalfedge();
			myHalfedge* he3 = new myHalfedge();

			he0->source = v0; he1->source = v1; he2->source = v2; he3->source = v3;
			he0->next = he1; he1->next = he2; he2->next = he3; he3->next = he0;
			he0->prev = he3; he1->prev = he0; he2->prev = he1; he3->prev = he2;

			for (auto he : { he0, he1, he2, he3 }) {

				he->adjacent_face = new_f;
				he->source->originof = he;
			}

			new_f->adjacent_halfedge = he0;

			new_faces.push_back(new_f);
			new_halfedges.insert(new_halfedges.end(), { he0, he1, he2, he3 });

			for (auto vtx : { v0, v1, v2, v3 }) {

				if (std::find(new_vertices.begin(), new_vertices.end(), vtx) == new_vertices.end()) {

					new_vertices.push_back(vtx);
				}
			}

			e = e->next;

		} while (e != start);
	}

	// Assign twins after the creation of the new mesh
	for (myHalfedge* he1 : new_halfedges) {

		for (myHalfedge* he2 : new_halfedges) {

			if (he1 != he2 && he1->source == he2->next->source && he1->next->source == he2->source) {

				he1->twin = he2;
				break;
			}
		}
	}

	// Update the current mesh with new one created
	faces = new_faces;
	halfedges = new_halfedges;
	vertices = new_vertices;

	checkMesh();
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

		//Calculating the triangulation of the current face f and creating the triangular faces associate 
		// to the polygone face
		for (int i = 1;i <= counter - 2;i++) {

			//Linking each halfedges to their relative twins using the fan algorithm with a vertex as point of reference
			
			//First edge of the triangle
			myHalfedge* he;

			//Retreive the first edge
			if (i == 1){
				he = edges[0];
			}
			else{ //Retrevied the first edge as the twin from the last iteration
				he = prev_link;
			}

			//Second edge of the triangle
			myHalfedge* he2 = edges[i];
			//Third edge of the triangle
			myHalfedge* he3;

			//Retreived the last edge of the polygone as the third edge of the triangle
			if (i == counter - 2){
				he3 = edges.back();
			}
			else {
				
				//Creating the intern edges for triangulation and link the new edge created
				//with the previous triangle associate twin
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

			*(face->normal) = *(f->normal);

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


void myMesh::simplification() {

	double length_min = 9999999; //minimum length initialise to a number equivalent to a inifinity 
	myHalfedge* edge_to_remove = nullptr; //the small edge to collapse
	double length;
	myVertex* k;
	myVertex* n;

	for (myHalfedge* e : halfedges) {

		//If the edge or the edge next, twin is null we got the next halfedges
		if (!e || !e->next || !e->twin) {
			continue;
		}

		//If the source of the edge or next edge is empty or the point structure is not present we go to the next halfedges
		if (!e->source || !e->next->source || !e->source->point || !e->next->source->point) {
			continue;
		}

		//Taking the two end point of the current edge e
		k = e->source;
		n = e->next->source;

		//length of the edge
		length = k->point->dist(*n->point);

		//If the condition is true we catch the current edge as the small length of edge 
		if (length < length_min) {

			length_min = length;
			edge_to_remove = e;
		}
		
	}
	
	//After we get the edge with the minimal length in the mesh we start the simplification as
	//the fusion_edge to be the edge to collapse in the mesh

	//We retreive the two vertex from the fusion edge (we reuse the k and n variable)
	k = edge_to_remove->source;
	n = edge_to_remove->next->source;

	//We create the middle point of the edge
	myPoint3D* mid_point = new myPoint3D((k->point->X + n->point->X)*0.5, (k->point->Y + n->point->Y)*0.5, (k->point->Z + n->point->Z)*0.5);

	//Creating the vertex for the middle point
	myVertex* mid_vertex = new myVertex();
	mid_vertex->point = mid_point;
	mid_vertex->normal = new myVector3D();

	vertices.push_back(mid_vertex);

	//Creating the vector<myHalfedges> of out_k and out_n to redirect the halfedges for the simplification
	std::vector<myHalfedge*> out_k, out_n, remove;

	//Remove a the end the edge with the small length choose previously
	remove.push_back(edge_to_remove);
	remove.push_back(edge_to_remove->twin);

	//push back in to the vector associate with k and n
	for (myHalfedge* he : halfedges) {

		if (he->source == k) {

			out_k.push_back(he);
		}
		else if (he->source == n) {

			out_n.push_back(he);
		}
	}

	//We redirect all the halfedges pointing to the vertices of the edge to remove to the middle vertex created previously
	
	for (myHalfedge* he : out_k) {

		if (he != edge_to_remove && he != edge_to_remove->twin) {
			he->source = mid_vertex;
		}
	}

	for (myHalfedge* he : out_n) {
		
		if (he != edge_to_remove && he != edge_to_remove->twin) {
			he->source = mid_vertex;
		}

	}

	//Update of the faces
	auto updateFace = [](myHalfedge* edge){

		if (!edge || !edge->adjacent_face)
			return;

		myHalfedge* prev = edge->prev;
		myHalfedge* next = edge->next;

		if (prev && next) {

			prev->next = next;
			next->prev = prev;
			edge->adjacent_face->adjacent_halfedge = next;
		}
	};

	//Update the face link the edge to remove 
	updateFace(edge_to_remove);
	updateFace(edge_to_remove->twin);

	//Creation of vector for degenerate faces
	std::vector<myFace*> degenerateFaces;

	//Find those degenerate faces
	for (myFace* face : faces) {

		if (!face->adjacent_halfedge) {
			continue;
		}

		int count = 0;
		myHalfedge* start = face->adjacent_halfedge;
		myHalfedge* curr = start;

		do {
			count++;

			if (!curr->next) {
				break;
			}

			curr = curr->next;

		} while (curr != start);

		if (count < 3) {
			degenerateFaces.push_back(face);
			curr = start;

			do {

				if (std::find(remove.begin(), remove.end(), curr) == remove.end()) {
					remove.push_back(curr);
				}

				curr = curr->next;

			} while (curr != start);
		}
	}

	//For each degenerate face remove the face
	for (myFace* f : degenerateFaces) {

		auto it = std::find(faces.begin(), faces.end(), f);

		if (it != faces.end()) {
			faces.erase(it);
			delete f;
		}
	}

	//remove the halfedges to remove of the mesh
	for (myHalfedge* he : remove) {

		auto it = std::find(halfedges.begin(), halfedges.end(), he);

		if (it != halfedges.end()) {

			halfedges.erase(it);
			delete he;
		}
	}

	//Sub function to remove the vertex of the mesh
	auto removeVertex = [&](myVertex* v) {

		auto it = std::find(vertices.begin(), vertices.end(), v);

		if (it != vertices.end()) {
			vertices.erase(it);
			delete v;
		}
	};

	//Remove the vertex k & n
	removeVertex(k);
	removeVertex(n);

	std::map<std::pair<myVertex*, myVertex*>, myHalfedge*> edgePairs;
	
	//For each halfedge update the edges and their twins
	for (myHalfedge* he : halfedges) {

		if (!he->next) continue;

		myVertex* a = he->source;
		myVertex* b = he->next->source;

		auto key = std::make_pair(a, b);
		auto revKey = std::make_pair(b, a);

		if (edgePairs.count(revKey)) {
			he->twin = edgePairs[revKey];
			edgePairs[revKey]->twin = he;
		}
		else {
			edgePairs[key] = he;
		}
	}

}

