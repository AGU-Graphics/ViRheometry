#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <Eigen/Dense>

#ifdef __APPLE__
#include <GLUT/glut.h>
#include <OpenGL/OpenGL.h>
#else
#include <GL/gl.h>
#include <GL/glut.h>
#endif
#include "ObjGL.h"

ObjGL::ObjGL()
: num_triangles(0)
{
}

void computeNormal( const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3, float* io_Normals )
{
	const Eigen::Vector3d n = ( ( v2 - v1 ).cross( v3 - v1 ) ).normalized();

	io_Normals[0] = io_Normals[3] = io_Normals[6] = n.x();
	io_Normals[1] = io_Normals[4] = io_Normals[7] = n.y();
	io_Normals[2] = io_Normals[5] = io_Normals[8] = n.z();
}

int ObjGL::Load( const char* in_filename )
{
	FILE* f = NULL;
  const int BUF_MAX = 4096;
  char line[ BUF_MAX ];

  char material_name[ BUF_MAX ];
  
  f = fopen( in_filename, "r" );
  if( !f )
  {
    std::cout << "Could not open obj file: " << in_filename << std::endl;
    return false;
  }

  std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d>> vertices;
  std::vector<Eigen::Vector3i, Eigen::aligned_allocator<Eigen::Vector3i>> triangles;

	while( fgets( line, BUF_MAX, f ) != NULL )
  {
		if ( line[0] == 'v' )
    {
      if ( line[1] != 'n' && line[1] != 't' )
      {
        float x, y, z;
        sscanf( line, "v %f %f %f", &x, &y, &z );
        vertices.emplace_back( x, y, z );
      }
    }

    if ( line[0] == 'f' )
    {
      char* tp = strtok( &(line[2]), " " );
      std::vector<int> _vertices;
      while( tp != NULL )
      {
        int vid = -1;
        sscanf( tp, "%d", &vid );
        
        _vertices.emplace_back( vid-1 );
        
        tp = strtok( NULL, " " );
      }
      
      for( int i=0; i<_vertices.size() - 2; i++ )
      {
        triangles.emplace_back( _vertices[0], _vertices[i+1], _vertices[i+2] );
      }
    }
  }

  vertex_data = (float*) malloc( sizeof( float ) * triangles.size() * 9 );
  normal_data = (float*) malloc( sizeof( float ) * triangles.size() * 9 );

  for( int i=0; i<triangles.size(); i++ )
  {
  	for( int j=0; j<3; j++ )
  	{
	  	vertex_data[ i*9 + j*3     ] = vertices[ triangles[i](j) ].x();
			vertex_data[ i*9 + j*3 + 1 ] = vertices[ triangles[i](j) ].y();
			vertex_data[ i*9 + j*3 + 2 ] = vertices[ triangles[i](j) ].z();
		}

		computeNormal( vertices[ triangles[i](0) ], vertices[ triangles[i](1) ], vertices[ triangles[i](2) ], &normal_data[ i * 9 ] );
  }

  num_triangles = triangles.size();
}

void ObjGL::Release()
{
	free( vertex_data );
	free( normal_data );
}

void ObjGL::Draw()
{
	glEnableClientState( GL_VERTEX_ARRAY );
	glEnableClientState( GL_NORMAL_ARRAY );
	glVertexPointer( 3, GL_FLOAT,	0, vertex_data );
	glNormalPointer( GL_FLOAT, 0, normal_data );
	glDrawArrays( GL_TRIANGLES, 0, num_triangles*3 );
	glDisableClientState( GL_VERTEX_ARRAY );
	glDisableClientState( GL_NORMAL_ARRAY );
}



