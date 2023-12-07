#ifndef ObjGL_h
#define ObjGL_h

#include <stdio.h>

class ObjGL
{
public:
	ObjGL();
	
	int Load( const char *filename );
	void Draw();
	void Release();
	
private:
	float* normal_data;
	float* vertex_data;
	int num_triangles;
};

#endif /* ObjGL_h */
