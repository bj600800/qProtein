/**
 *  @file       biom.h
 *  @ingroup    global_gamer
 *  @brief      Some parameter/datatype definitions
 *  @author     Zeyun Yu (zeyun.yu@gmail.com)
 *  @note       None
 *  @version    $Id: biom.h,v 1.28 2010/10/14 19:25:09 fetk Exp $
 *  
 *  @attention
 *  @verbatim
 *
 * GAMER = < Geometry-preserving Adaptive MeshER >
 * Copyright (C) 1994-- Michael Holst and Zeyun Yu
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 * 
 *  @endverbatim
 */

#ifndef _BIOM_H_
#define _BIOM_H_

#include <gamer/gamer_base.h>

#include <gamer/triangle.h>
#include <gamer/tetgen.h>

#if 0
#include <stdio.h>       /* One of 15 ISO-C headers -- get via MALOC */
#include <stdlib.h>      /* One of 15 ISO-C headers -- get via MALOC */
#include <math.h>        /* One of 15 ISO-C headers -- get via MALOC */
#include <time.h>        /* One of 15 ISO-C headers -- get via MALOC */
#include <string.h>      /* One of 15 ISO-C headers -- get via MALOC */
#include <sys/times.h>   /* System-dependent; include in .c file via MALOC */
#include <sys/stat.h>    /* System-dependent; include in .c file via MALOC */
#include <memory.h>      /* System-dependent; include in .c file via MALOC */
#endif

// critical parameters for mesh generation

/** @brief Isovalue used in the Marching Cube method */
#define IsoValue          2.5      
/** @brief Blurring blobyness used in conversion from PDB/PQR to 3D volumes */
#define BLOBBYNESS        -0.2f    
/** @brief Discretization rate of 3D volumes */
#define DIM_SCALE         1.99   
/** @brief Coarsening Rate in surface post-processing */
#define CoarsenRate       0.1666   
/** @brief The minimal volumes (in voxels) of islands to be removed */
#define MIN_VOLUME        333333   
/** @brief The size of the bounding sphere (= object size X the following rate) */
#define SphereRatio       40

#define MaxVal               999999
#define MaxAtom              10

// Other definitions and data structures

/** @brief Other definition */
#define PIE              3.14159265358979f 
/** @brief Other definition */
#define IndexVect(i,j,k) ((k)*xdim*ydim + (j)*xdim + (i))
#define IndexVect1(i,j,k) ((k)*xdim1*ydim1 + (j)*xdim1 + (i))
/** @brief Other definition */
#define max(x, y)        ((x>y) ? (x):(y))
/** @brief Other definition */
#define min(x, y)        ((x<y) ? (x):(y))

/** @brief Other definition */
#define _LITTLE_ENDIAN   1

/** @brief Other data structure FLT2VECT (float) */
typedef struct {
  float x;   /**< @brief x-coordinate */
  float y;   /**< @brief y-coordinate */
}FLT2VECT;

/** @brief Other data structure FLTVECT (float) */
typedef struct {
  float x;   /**< @brief x-coordinate */
  float y;   /**< @brief y-coordinate */
  float z;   /**< @brief z-coordinate */
}FLTVECT;

/** @brief Other data structure DBLVECT (double float) */
typedef struct {
  double x;   /**< @brief x-coordinate */
  double y;   /**< @brief y-coordinate */
  double z;   /**< @brief z-coordinate */
}DBLVECT;

/** @brief Other data structure INT3VECT (int) */
typedef struct {
  int a;   /**< @brief first integer */
  int b;   /**< @brief second integer */
  int c;   /**< @brief third integer */
}INT3VECT;

/** @brief Other data structure INT4VECT (int) */
typedef struct {
  int a;   /**< @brief first integer */
  int b;   /**< @brief second integer */
  int c;   /**< @brief third integer */
  int d;   /**< @brief fourth integer */
}INT4VECT;

/** @brief Other data structure NPNT3 */
typedef struct NeighborPoint3 NPNT3;
/** @brief Other data structure NeighborPoint3 */
struct NeighborPoint3{
  int a;   /**< @brief first integer */
  int b;   /**< @brief second integer */
  int c;   /**< @brief third integer */
  NPNT3 *next;   /**< @brief pointer to the next triangle */
};

/** @brief Other data structure Triface */
struct Triface{
  int points[3];   /**< @brief point list */
  int marker;      /**< @brief marker integer */
  Triface *next;   /**< @brief pointer to the next Triface */
};

/** @brief Other data structure SurfaceMesh (for surface meshes) */
typedef struct {
  int nv;   /**< @brief number of vertices */
  int nf;   /**< @brief number of triangles */
  int nvm;   /**< @brief number of vertex markers */
  int nfm;   /**< @brief number of face markers */
  float avglen;   /**< @brief average edge length */
  float min[3];   /**< @brief minimal coordinate of nodes */
  float max[3];   /**< @brief maximal coordinate of nodes */
  FLTVECT* vertex;   /**< @brief pointer to the vertices */
  INT3VECT* face;    /**< @brief pointer to the triangles */
  INT3VECT* neighbor;    /**< @brief pointer to the neighbors (triangles) */
  NPNT3** neighbor_list; /**< @brief pointer to the neighbor list */
  int* vertex_markers;  /**< @brief pointer to any vertex markers */
  int* face_markers;  /**< @brief pointer to any face markers */
}SurfaceMesh;

typedef struct FETK_VX{
  int id;   
  int chrt; 
  float x;  
  float y;  
  float z;  
}FETK_VX;

typedef struct FETK_SS{
  int id;  
  int grp;
  int mat;
  int fa;  
  int fb;  
  int fc;  
  int fd;  
  int na;  
  int nb;  
  int nc;  
  int nd;  
}FETK_SS;

typedef struct GemMesh{
  int dim;       
  int dimii;     
  int vertices;  
  int simplices; 
  FETK_VX *vv;  
  FETK_SS *ss;   
}GemMesh;


/** @brief Other data structure TeTraMesh (for tetrahedral meshes) */
typedef struct {
  int nv;   /**< @brief number of nodes */
  int nf;   /**< @brief numner of tetrahedra */
  FLTVECT *vertex;   /**< @brief pointer to the vertices */
  INT4VECT *face;   /**< @brief pointer to the tetrahedra */
  INT4VECT *neighbor;   /**< @brief pointer to the neighbors (tetrahedra) */
}TeTraMesh;

/** @brief Other data structure SPNT */
typedef struct SamplePoint SPNT;
/** @brief Other data structure SamplePoint */
struct SamplePoint{
  float x;   /**< @brief x-coordinate */
  float y;   /**< @brief y-coordinate */
  float z;   /**< @brief z-coordinate */
  SPNT *next;   /**< @brief pointer to next vertex */
};

/** @brief Other data structure NPNT2 */
typedef struct NeighborPoint2 NPNT2;
/** @brief Other data structure NeighborPoint2 */
struct NeighborPoint2{
  int a;   /**< @brief first integer */
  int b;   /**< @brief second integer */
  NPNT2 *next;   /**< @brief pointer to the next point */
};

/** @brief Other data structure ATOM */
typedef struct {
  float x;   /**< @brief x-coordinate */
  float y;   /**< @brief y-coordinate */
  float z;   /**< @brief z-coordinate */
  float radius;   /**< @brief radius */
}ATOM;


/** @brief Other data structure EIGENVECT */
typedef struct {
  float x1;   /**< @brief x-coordinate of first eigenvector */
  float y1;   /**< @brief y-coordinate of first eigenvector */
  float z1;   /**< @brief z-coordinate of first eigenvector */
  float x2;   /**< @brief x-coordinate of second eigenvector */
  float y2;   /**< @brief y-coordinate of second eigenvector */
  float z2;   /**< @brief z-coordinate of second eigenvector */
  float x3;   /**< @brief x-coordinate of third eigenvector */
  float y3;   /**< @brief y-coordinate of third eigenvector */
  float z3;   /**< @brief z-coordinate of third eigenvector */
}EIGENVECT;


/** @brief Other data structure MinHeapS */
typedef struct {
  unsigned short *x; /**< @brief x-coordinate */
  unsigned short *y; /**< @brief y-coordinate */
  unsigned short *z; /**< @brief z-coordinate */
  int *seed;         /**< @brief seed */
  float *dist;       /**< @brief distance */
  int size;          /**< @brief size */
}MinHeapS;


/** @brief Other data structure SEEDS */
typedef struct{
  float seedx;       /**< @brief x-coordinate */
  float seedy;       /**< @brief y-coordinate */
  float seedz;       /**< @brief z-coordinate */
  int atom[MaxAtom]; /**< @brief atom array */
}SEEDS;

// Forward declarations
struct GemMesh;
class tetgenio;

// GemMesh functions
GemMesh* GemMesh_fromTetgen(tetgenio& tetio);
GemMesh* GemMesh_fromSurfaceMesh(SurfaceMesh* surfmesh, char* tetgen_params);
GemMesh* GemMesh_fromPdb(tetgenio* out, float radius, float centerx, 
			 float centery, float centerz, char *ActiveSite, 
			 int output_flag);
void GemMesh_writeMcsf(GemMesh* out, char* filename);
void GemMesh_dtor(GemMesh*);
void GemMesh_writeOFF(GemMesh* Gem_mesh, char* filename);

// Surface mesh constructors
SurfaceMesh* SurfaceMesh_ctor(unsigned int, unsigned int);
SurfaceMesh* SurfaceMesh_readOFF(char* filename);
SurfaceMesh* SurfaceMesh_readPDB(char* filename);
SurfaceMesh* SurfaceMesh_sphere(int);
SurfaceMesh* SurfaceMesh_marchingCube(int, int, int, float*, float, SPNT **);
SurfaceMesh* SurfaceMesh_marchingCube(int, int, int, float*, float, float*, float, SPNT **);
SurfaceMesh* SurfaceMesh_readLattice(char*, float, bool);
SurfaceMesh* SurfaceMesh_readLattice(char*, char*, float, bool);

// Methods working on a SurfaceMesh
void SurfaceMesh_dtor(SurfaceMesh* surfmesh);
void SurfaceMesh_createNeighborlist(SurfaceMesh* surfmesh);
void SurfaceMesh_destroyNeighborlist(SurfaceMesh* surfmesh);
void SurfaceMesh_releaseData(SurfaceMesh* surfmesh);
void SurfaceMesh_correctNormals(SurfaceMesh* surfmesh);
void SurfaceMesh_writeOFF(SurfaceMesh* surfmesh, char* filename);
void SurfaceMesh_createVertexMarkers(SurfaceMesh* surfmesh);
void SurfaceMesh_destroyVertexMarkers(SurfaceMesh* surfmesh);
void SurfaceMesh_resetVertexMarkers(SurfaceMesh* surfmesh);
void SurfaceMesh_createFaceMarkers(SurfaceMesh* surfmesh);
void SurfaceMesh_destroyFaceMarkers(SurfaceMesh* surfmesh);
void SurfaceMesh_resetFaceMarkers(SurfaceMesh* surfmesh);
void SurfaceMesh_getMinMaxAngles(SurfaceMesh*  surfmesh, float*, float*, int*, int*, 
				 int, int);
void SurfaceMesh_refine(SurfaceMesh* surfmesh);
bool SurfaceMesh_smooth(SurfaceMesh* surfmesh, unsigned int, unsigned int, 
			unsigned int, bool);
void SurfaceMesh_normalSmooth(SurfaceMesh* surfmesh);
char SurfaceMesh_coarse(SurfaceMesh* surfmesh, float, float, float, float);
void SurfaceMesh_assignActiveSites(SurfaceMesh* surfmesh,
				   ATOM*, unsigned int, unsigned int*);
SurfaceMesh* SurfaceMesh_merge(SurfaceMesh*, SurfaceMesh*);
ATOM SurfaceMesh_getCenterRadius(SurfaceMesh* surfmesh);
void SurfaceMesh_translate(SurfaceMesh* surfmesh, float dx, float dy, float dz);
void SurfaceMesh_scale(SurfaceMesh* surfmesh, float scale_x, 
		       float scale_y, float scale_z);
void SurfaceMesh_scale(SurfaceMesh* surfmesh, float scale);
void SurfaceMesh_centeralize(SurfaceMesh* surfmesh);

//void ReadActiveSiteFile(char*, unsigned int&, ATOM*&, unsigned int*&, float*&);
void ReadActiveSiteFile(char*, unsigned int&, ATOM*&, unsigned int*&);
void ReadRawiv(int *, int *, int *, float **, char *,float *, float *);
void SurfaceExtract(TeTraMesh *, SurfaceMesh *);
float PDB2Volume(char *, float **, int *, int *, int *,float *, float *, ATOM **,int *,char);

#endif /* _BIOM_H_ */
