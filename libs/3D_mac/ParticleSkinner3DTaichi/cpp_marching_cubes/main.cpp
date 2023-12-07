#include <Eigen/Dense>

#include <iostream>
#include <fstream>
                                                                                
#include <iomanip>
#include <vector>
#define _USE_MATH_DEFINES
#include <math.h>

#include "MarchingCubes.h"

struct PhiData
{
  Eigen::Vector3d bb_min;
  double h;
  Eigen::Vector3i resolution;
  Eigen::VectorXd data;
};

PhiData g_Phi;
float g_MaxPhi = 1.0; //excluding 1.0e10
float g_MinPhi = -1.0e9;

void showUsage( const char* argv0 )
{
  std::cout << "<" << argv0 << "> <phi dat file> <output obj file>" << std::endl;
}

template<class T>
std::string convertToString( const T& tostring )
{
	std::string out_string;
	std::stringstream ss;

	ss << tostring;
	ss >> out_string;

	return out_string;
}

void loadData( const char* in_fn )
{
  FILE* f = fopen( in_fn, "rb" );

  float data_f;
  fread( &data_f, sizeof(float), 1, f );
  g_Phi.bb_min.x() = data_f;
  fread( &data_f, sizeof(float), 1, f );
  g_Phi.bb_min.y() = data_f;
  fread( &data_f, sizeof(float), 1, f );
  g_Phi.bb_min.z() = data_f;

  fread( &data_f, sizeof(float), 1, f );
  g_Phi.h = data_f;

  int32_t data_i;
  fread( &data_i, sizeof(int32_t), 1, f );
  g_Phi.resolution.x() = data_i;
  fread( &data_i, sizeof(int32_t), 1, f );
  g_Phi.resolution.y() = data_i;
  fread( &data_i, sizeof(int32_t), 1, f );
  g_Phi.resolution.z() = data_i;

  std::cout << "resolution: " << g_Phi.resolution.transpose() << std::endl;

  const int num_elems = g_Phi.resolution.x() * g_Phi.resolution.y() * g_Phi.resolution.z();
  g_Phi.data.resize( num_elems );

  g_MaxPhi = 0.0;
  g_MinPhi = 1.0e10;

  for( int i=0; i<num_elems; i++ )
  {
    fread( &data_f, sizeof(float), 1, f );
    g_Phi.data(i) = data_f;
    if( data_f < 1.0e9 && data_f > g_MaxPhi )
      g_MaxPhi = std::max<float>( g_MaxPhi, data_f );

    g_MinPhi = std::min<float>( g_MinPhi, data_f );
  }

  std::cout << "first 5 data: " << std::endl;
  for( int i=0; i<5; i++ )
  {
    std::cout << g_Phi.data(i) << ", ";
  }
  std::cout << std::endl;

  std::cout << "data range: " << g_Phi.data.minCoeff() << " - " << g_Phi.data.maxCoeff() << std::endl;

  std::cout << "min_phi: " << g_MinPhi << std::endl;
  std::cout << "max_phi: " << g_MaxPhi << std::endl;

  fclose(f);
}

bool exportOBJMesh( const std::string& output_mesh_name, const std::vector<Eigen::Array3i>& triangles, const std::vector<Eigen::Vector3d>& vertices, std::string& save_status )
{
  std::ofstream obj_file( output_mesh_name.c_str() );
  assert( obj_file.good() );
  if( !obj_file.good() )
  {
    save_status = "Failed to open file";
    return false;
  }

  for( std::vector<Eigen::Vector3d>::size_type vrt_idx = 0; vrt_idx < vertices.size(); ++vrt_idx )
  {
    obj_file << "v " << std::setprecision(200) << vertices[vrt_idx].x() << " " << std::setprecision(200) << vertices[vrt_idx].y() << " " << std::setprecision(200) << vertices[vrt_idx].z() << std::endl;
    if( !obj_file.good() )
    {
      save_status = "Failed to save vertex " + convertToString(vrt_idx);
      return false;
    }
  }

  for( std::vector<Eigen::Vector3i>::size_type tri_idx = 0; tri_idx < triangles.size(); ++tri_idx )
  {
    obj_file << "f " << (triangles[tri_idx]+1).transpose() << std::endl;
    if( !obj_file.good() )
    {
      save_status = "Failed to save triangle " + convertToString(tri_idx);
      return false;
    }
  }

  obj_file.close();

  return true;
}

int main( int argc, char* argv[] )
{
  if( argc != 3 )
  {
    showUsage( argv[0] );
    exit(0);
  }

  loadData( argv[1] );

  std::vector<Eigen::Array3i> triangles;
  std::vector<Eigen::Vector3d> vertices;

  Eigen::VectorXd phi = Eigen::VectorXd::Zero( g_Phi.resolution.x() * g_Phi.resolution.y() * g_Phi.resolution.z() );
  for( unsigned int i = 0; i < g_Phi.resolution.x(); ++i )
  {
    for( unsigned int j = 0; j < g_Phi.resolution.y(); ++j )
    {
      for( unsigned int k = 0; k < g_Phi.resolution.z(); ++k )
      {
        phi( i * g_Phi.resolution.y() * g_Phi.resolution.z() + j * g_Phi.resolution.z() + k ) = g_Phi.data( k * g_Phi.resolution.x() * g_Phi.resolution.y() + j * g_Phi.resolution.x() + i );
      }
    }
  }

  MarchingCubes mc( g_Phi.resolution.x(), g_Phi.resolution.y(), g_Phi.resolution.z() );
  mc.set_method( true );
  mc.set_ext_data( phi.data() );
  mc.init_all();
  
  const double level_set_value = 0;
  mc.run( level_set_value );

  triangles.resize( mc.ntrigs() );
  vertices.resize( mc.nverts() );

  for( int i = 0; i < mc.nverts(); ++i )
  {
    const Vertex* cur_vert = mc.vert( i );
    assert( cur_vert != NULL );
    vertices[i] << cur_vert->x, cur_vert->y, cur_vert->z;
    vertices[i] = vertices[i].array() * g_Phi.h + Eigen::Array3d( g_Phi.bb_min.x(), g_Phi.bb_min.y(), g_Phi.bb_min.z() );
  }

  for( int i = 0; i < mc.ntrigs(); ++i )
  {
    const Triangle* cur_tri = mc.trig(i);
    assert( cur_tri != NULL );
    triangles[i] << cur_tri->v1, cur_tri->v2, cur_tri->v3;
  }

  mc.clean_all();

  std::string save_status;
  bool res_save_obj = exportOBJMesh( std::string(argv[2]), triangles, vertices, save_status );
  if( !res_save_obj )
    std::cout << save_status << std::endl;

  return 0;
}
