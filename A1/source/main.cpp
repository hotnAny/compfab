//Computational Fabrication Assignment #1
// By David Levin 2014

// Modified by Xiang 'Anthony' Chen
// Last updated: Feb-02-2015

#include <iostream>
#include <vector>
#include <ctime>
#include "../include/CompFab.h"
#include "../include/Mesh.h"

// debugging mode
#define D true

/* ========================================================================================================================
    
    program parameters    

  ========================================================================================================================*/

#define PARAMSLICING    "s"     // specifying slicing layer
#define PARAMCONTOUR    "c"     // outputing contour only
#define PARAMMULTIDIR   "m"     // using multiple rays
#define PARAMPRESCREEN  "p"     // prescreening triangles
#define PARAMDIIM       "d"     // voxel grid dimension, defualt 32
#define PARAMLOG        "l"     // log performance result, default false

int slicedLayer = -1;
bool contourOnly = false;
bool multiDirections = false;
bool prescreenTriangles = false;
unsigned int dim = 32; //dimension of voxel grid (e.g. 32x32x32)
bool logPerf = false;

void sendUsageMessageAndExit() {
    if(!logPerf) {
        std::cout << "Usage: ./voxelizer [ -d <int> | -s <int> | -cmpl ] <input_mesh_file_name> <output_mesh_file_name> \n"
                << "\t -d \t dimension of the voxel grid (default 32).\n"
                << "\t -s \t only output voxels of s specific layer (0<= s < d).\n"
                << "\t -c \t only output contour.\n"
                << "\t -m \t casting multiple axis-aligned rays when voxelizing (default using only (1, 0, 0)).\n"
                << "\t -p \t using triangle prescreening to accelerate voxlization.\n"
                << "\t -l \t logging performance result (default false)"
                ;
        exit(0);
    }
}

/* ========================================================================================================================

    a series of useful custom functions

  ======================================================================================================================== */

double min3(double x, double y, double z) {
    return fmin(x, fmin(y, z));
}

double max3(double x, double y, double z) {
    return fmax(x, fmax(y, z));
}

CompFab::Vec3 mulScalar(CompFab::Vec3 v0, float s) {
    CompFab::Vec3 v1(v0.m_x * s, v0.m_y * s, v0.m_z * s);
    return v1;
}

double norm(CompFab::Vec3 v) {
    return sqrt(v.m_x * v.m_x + v.m_y * v.m_y + v.m_z * v.m_z);
}

void printVector(CompFab::Vec3 v, std::string description = "") {
    std::cout << description << "(" << v.m_x << ", " << v.m_y << ", " << v.m_z << ")\n";
}

void printTriangle(CompFab::Triangle tri, std::string description = "") {
    std::cout << description << "\n(" << tri.m_v1.m_x << ", " << tri.m_v1.m_y << ", " << tri.m_v1.m_z << ")\n"
                             << "(" << tri.m_v2.m_x << ", " << tri.m_v2.m_y << ", " << tri.m_v2.m_z << ")\n"
                             << "(" << tri.m_v3.m_x << ", " << tri.m_v3.m_y << ", " << tri.m_v3.m_z << ")\n\n";
}

// compute if points p1 and p2 are on the same side of the line formed by points a and b
// ref: http://www.blackpawn.com/texts/pointinpoly/
bool onSameSide(CompFab::Vec3 &p1, CompFab::Vec3 &p2, CompFab::Vec3 &a, CompFab::Vec3 &b)
{
    return ((b - a) % (p1 - a)) * ((b-a) % (p2 - a)) >= 0;
}

// compute if a point p is inside a triangle
// a point p is in a triangle if for each of the edges,
// p is on the same side as the vertex that's not on that edge.
bool isInTriangle(CompFab::Vec3 &p, CompFab::Triangle &triangle) {
    return  onSameSide(p, triangle.m_v1, triangle.m_v2, triangle.m_v3) &&
            onSameSide(p, triangle.m_v2, triangle.m_v3, triangle.m_v1) &&
            onSameSide(p, triangle.m_v3, triangle.m_v1, triangle.m_v2);
}

// compute if a point p is identical to any point in an array
bool isIn(CompFab::Vec3 &p, std::vector<CompFab::Vec3> pnts) {
    for(std::vector<CompFab::Vec3>::iterator it = pnts.begin(); it != pnts.end(); ++it) {
        if(fabs(it->m_x - p.m_x) < EPSILON && fabs(it->m_y - p.m_y) < EPSILON && fabs(it->m_z - p.m_z) < EPSILON) {
            return  true;
        }
    }
    return false;
}



/* ========================================================================================================================

    assignment-provided code below (with student modifications)

  ======================================================================================================================== */

//Ray-Triangle Intersection
//Returns 1 if triangle and ray intersect, 0 otherwise
//
// function modified to provide the point of intersection in pnt    
// ref: http://geomalgorithms.com/a06-_intersect-2.html
int rayTriangleIntersection(CompFab::Ray &ray, CompFab::Triangle &triangle, CompFab::Vec3 &pnt)
{
    /********* ASSIGNMENT *********/
    /* Ray-Triangle intersection test: Return 1 if ray intersects triangle, 
     * 0 otherwise */

    /* 
        find the point where the ray intersects with the plane of the triangle
    */

    // calculate the normal
    CompFab::Vec3 nml = (triangle.m_v2 - triangle.m_v1) % (triangle.m_v3 - triangle.m_v1);

    // calculate the intersecting point
    float nmlDotDir = nml * ray.m_direction;

    // parallel to the triangle
    if(fabs(nmlDotDir) < EPSILON) {
        return 0;
    }

    float r = (nml * (triangle.m_v1 - ray.m_origin)) / nmlDotDir;

    // not on the direction of the ray
    if(r < 0) {
        return 0;
    }

    // the point intersecting the triangle's plane
    CompFab::Vec3 p = ray.m_origin + mulScalar(ray.m_direction, r);

    /*
        find if the intersecting point is inside the triangle:
        a point P is in a triangle if for the edge between the triangle's each pair of vertices, 
        P is on the same side as the 3rd vertex
    */
    
    if(isInTriangle(p, triangle)) {
        pnt.m_x = p.m_x;
        pnt.m_y = p.m_y;
        pnt.m_z = p.m_z;
        return 1;
    }

    return 0;
}

//Triangle list (global)
typedef std::vector<CompFab::Triangle> TriangleList;

TriangleList g_triangleList;
CompFab::VoxelGrid *g_voxelGrid;


//Number of intersections with surface made by a ray originating at voxel and cast in direction.
int numSurfaceIntersections(CompFab::Vec3 &voxelPos, CompFab::Vec3 &dir)
{
    int numHits = 0;
    
    /********* ASSIGNMENT *********/
    /* Check and return the number of times a ray cast in direction dir, 
     * from voxel center voxelPos intersects the surface */

    CompFab::Ray ray(voxelPos, dir);

    // this vector filter out directions that are not axis-aligned
    CompFab::Vec3 mask( (fabs(dir.m_x) == 1 && fabs(dir.m_y) == 0 && fabs(dir.m_z) == 0),
                        (fabs(dir.m_x) == 0 && fabs(dir.m_y) == 1 && fabs(dir.m_z) == 0),
                        (fabs(dir.m_x) == 0 && fabs(dir.m_y) == 0 && fabs(dir.m_z) == 1)
                    );

    // an array to collect intersecting points
    std::vector<CompFab::Vec3> intersects;
    for(std::vector<CompFab::Triangle>::iterator it = g_triangleList.begin(); it != g_triangleList.end(); ++it) {

        // process every triangle if prescreening is not turned on
        if( !prescreenTriangles ||  

            // if prescreening is on, the following conditions filter out triangles 
            // whose bounding boxes do not contain the ray cast form voxelPos
            ( ( (it->m_min.m_x <= voxelPos.m_x && voxelPos.m_x <= it->m_max.m_x) || mask.m_x) &&
              ( (it->m_min.m_y <= voxelPos.m_y && voxelPos.m_y <= it->m_max.m_y) || mask.m_y) &&
              ( (it->m_min.m_z <= voxelPos.m_z && voxelPos.m_z <= it->m_max.m_z) || mask.m_z) )
          ) {
            CompFab::Vec3 pnt;
            int intersectValue = rayTriangleIntersection(ray, *it, pnt);
            if(intersectValue == 1) {
                // printVector(pnt, "intersection:\n");
                // printTriangle(*it, "triangle: ");

                // in cases where the ray intersects two triangles at the same point
                // the conclusion is marked ambiguous (-1)
                if(isIn(pnt, intersects)) {
                    return -1;
                }
                else {
                    intersects.push_back(pnt);
                }
            }
        }
    }

    // the number of hits = # of unique intersection points
    numHits = intersects.size();

    return numHits;
}


int loadMesh(char *filename, unsigned int dim)
{
    g_triangleList.clear();
    
    Mesh *tempMesh = new Mesh(filename, true);
    
    CompFab::Vec3 v1, v2, v3;

    //copy triangles to global list
    for(unsigned int tri =0; tri<tempMesh->t.size(); ++tri)
    {
        v1 = tempMesh->v[tempMesh->t[tri][0]];
        v2 = tempMesh->v[tempMesh->t[tri][1]];
        v3 = tempMesh->v[tempMesh->t[tri][2]];
        CompFab::Triangle triangle(v1,v2,v3);

        // pre-compute each triangle's bounding box
        if(prescreenTriangles) {
            triangle.m_min = CompFab::Vec3(min3(v1.m_x, v2.m_x, v3.m_x), min3(v1.m_y, v2.m_y, v3.m_y), min3(v1.m_z, v2.m_z, v3.m_z));
            triangle.m_max = CompFab::Vec3(max3(v1.m_x, v2.m_x, v3.m_x), max3(v1.m_y, v2.m_y, v3.m_y), max3(v1.m_z, v2.m_z, v3.m_z));
        }

        g_triangleList.push_back(triangle);
    }

    //Create Voxel Grid
    CompFab::Vec3 bbMax, bbMin;
    BBox(*tempMesh, bbMin, bbMax);
    
    //Build Voxel Grid
    double bbX = bbMax[0] - bbMin[0];
    double bbY = bbMax[1] - bbMin[1];
    double bbZ = bbMax[2] - bbMin[2];
    double spacing;
    
    if(bbX > bbY && bbX > bbZ)
    {
        spacing = bbX/(double)(dim-2);
    } else if(bbY > bbX && bbY > bbZ) {
        spacing = bbY/(double)(dim-2);
    } else {
        spacing = bbZ/(double)(dim-2);
    }
    
    CompFab::Vec3 hspacing(0.5*spacing, 0.5*spacing, 0.5*spacing);
    
    g_voxelGrid = new CompFab::VoxelGrid(bbMin-hspacing, dim, dim, dim, spacing);

    int numTriangles = tempMesh->t.size();

    delete tempMesh;
    
    return numTriangles;
   
}

void saveVoxelsToObj(const char * outfile)
{
 
    Mesh box;
    Mesh mout;
    int nx = g_voxelGrid->m_dimX;
    int ny = g_voxelGrid->m_dimY;
    int nz = g_voxelGrid->m_dimZ;
    double spacing = g_voxelGrid->m_spacing;
    
    CompFab::Vec3 hspacing(0.5*spacing, 0.5*spacing, 0.5*spacing);
    
    for (int ii = 0; ii < nx; ii++) {
        for (int jj = 0; jj < ny; jj++) {
            for (int kk = 0; kk < nz; kk++) {
                if(!g_voxelGrid->isInside(ii,jj,kk)){
                    continue;
                }
                CompFab::Vec3 coord(0.5f + ((double)ii)*spacing, 0.5f + ((double)jj)*spacing, 0.5f+((double)kk)*spacing);
                CompFab::Vec3 box0 = coord - hspacing;
                CompFab::Vec3 box1 = coord + hspacing;
                makeCube(box, box0, box1);
                mout.append(box);
            }
        }
    }

    mout.save_obj(outfile);
}

int main(int argc, char **argv)
{

    //Load OBJ
    // and parsing program parameters
    if(argc < 3)
    {
        // std::cout<<"Usage: Voxelizer InputMeshFilename OutputMeshFilename \n";
        sendUsageMessageAndExit();
        exit(0);
    } else if(argc > 3) {
        try {
            for(int i=1; i<argc-2; i++) {
                if(argv[i][0] != '-') {
                    throw 0;
                }

                std::string strOpt(argv[i]);

                if(strOpt.find(PARAMDIIM) != std::string::npos) {
                    dim = std::stoi(argv[++i]);
                } 

                if(strOpt.find(PARAMSLICING) != std::string::npos) {
                    slicedLayer = std::stoi(argv[++i]);
                    if(0 > slicedLayer || slicedLayer >= dim) {
                        if(!logPerf) {
                            std::cout << "layer number out of range (should be 0 to " << dim << ")\n";
                        }
                        exit(0);
                    }
                } 

                if(strOpt.find(PARAMCONTOUR) != std::string::npos) {
                    contourOnly = true;
                } 

                if(strOpt.find(PARAMMULTIDIR) != std::string::npos) {
                    multiDirections = true;
                } 

                if(strOpt.find(PARAMPRESCREEN) != std::string::npos) {
                    prescreenTriangles = true;
                } 

                if(strOpt.find(PARAMLOG) != std::string::npos) {
                    logPerf = true;
                }

            }
        } catch(int e) {
            sendUsageMessageAndExit();
        } catch(...) {
            sendUsageMessageAndExit();
        }
    }
    
    if(!logPerf) {
        std::cout<<"Load Mesh : "<<argv[argc - 2]<<"\n";
        std::cout << "Voxel grid dimension: " << dim << "x" << dim << "x" << dim << "\n";

        if(multiDirections) {
            std::cout << "Casting multi-directional rays ...\n";
        }

        if(prescreenTriangles) {
            std::cout << "Using triangle prescreening ...\n";
        }
    }

    clock_t begin = clock();
    int numTriangles = loadMesh(argv[argc - 2], dim);

    if(!logPerf) {
        std::cout<< "Number of triangles: " << numTriangles << "\n";
    }
    
    //Cast ray, check if voxel is inside or outside
    //even number of surface intersections = outside (OUT then IN then OUT)
    // odd number = inside (IN then OUT)

    CompFab::Vec3 voxelPos;
    CompFab::Vec3 direction(1.0,0.0,0.0);
    std::vector<CompFab::Vec3> directions;

    if(multiDirections) {
        directions.push_back(CompFab::Vec3(1.0, 0.0, 0.0));
        directions.push_back(CompFab::Vec3(-1.0, 0.0, 0.0));
        directions.push_back(CompFab::Vec3(0.0, 1.0, 0.0));
        directions.push_back(CompFab::Vec3(0.0, -1.0, 0.0));
        directions.push_back(CompFab::Vec3(0.0, 0.0, 1.0));
        directions.push_back(CompFab::Vec3(0.0, 0.0, -1.0));
    }
    
    /********* ASSIGNMENT *********/
    /* Iterate over all voxels in g_voxelGrid and test whether they are inside our outside of the
     * surface defined by the triangles in g_triangleList */

    int nx = g_voxelGrid->m_dimX;
    int ny = g_voxelGrid->m_dimY;
    int nz = g_voxelGrid->m_dimZ;
    double spacing = g_voxelGrid->m_spacing;
    
    CompFab::Vec3 hspacing(0.5*spacing, 0.5*spacing, 0.5*spacing);
    
    int numDim = nx * ny * nz;
    int cntVoxels = 0;
    int numVoxels = 0;
    for (int ii = 0; ii < nx; ii++) {
        for (int jj = 0; jj < ny; jj++) {
            for (int kk = 0; kk < nz; kk++) {

                // show progress
                if(!logPerf) {
                    std::cout << "Voxelizing: " << ++cntVoxels*100/numDim << "%\r";
                }

                // skip layers if they aren't the specified one
                if(0 <= slicedLayer && slicedLayer < dim && kk != slicedLayer) {
                    continue;
                }

                int idx = kk * (nx * ny) + jj * ny + ii;
                CompFab::Vec3 coord(((double)ii)*spacing, ((double)jj)*spacing, ((double)kk)*spacing);

                // cast single or multiple rays
                if(!multiDirections) {
                    g_voxelGrid->m_insideArray[idx] = (numSurfaceIntersections(coord, direction) % 2 == 1);
                } else {
                    int votesInside = 0;
                    int numValidVotes = 0;
                    for(std::vector<CompFab::Vec3>::iterator it = directions.begin(); it != directions.end(); ++it) {
                        int numIntersects = numSurfaceIntersections(coord, *it);
                        if(numIntersects >= 0) {
                            votesInside += (numIntersects % 2 == 1) ? 1 : 0;
                            numValidVotes++;
                        }
                    }
                    g_voxelGrid->m_insideArray[idx] = votesInside > numValidVotes / 2;
                    
                }
                
                if(g_voxelGrid->m_insideArray[idx] == 1) {
                    numVoxels++;
                }
            }
        }
    }

    if(!logPerf) {
        std::cout << "\n";
    }

    /*
        computing contour by removing voxels that have all the four neighbours (N, E, S, W)
    */
    if(0 <= slicedLayer && slicedLayer < dim && contourOnly) {
        if(!logPerf) {
            std::cout << "Extracting contour ...\n";
        }
        std::vector<int> indicesInside;
        for (int ii = 1; ii < nx-1; ii++) {
            for (int jj = 1; jj < ny-1; jj++) {
                int kk = slicedLayer;
                int idx = kk * (nx * ny) + jj * ny + ii;

                if(g_voxelGrid->m_insideArray[idx] == 1 &&
                    g_voxelGrid->m_insideArray[idx - 1] == 1 &&
                    g_voxelGrid->m_insideArray[idx + 1] == 1 &&
                    g_voxelGrid->m_insideArray[idx - ny] == 1 &&
                    g_voxelGrid->m_insideArray[idx + ny] == 1) 
                {
                    indicesInside.push_back(idx);
                }
            }
        }

        for(std::vector<int>::iterator it = indicesInside.begin(); it != indicesInside.end(); ++it) {
            g_voxelGrid->m_insideArray[*it] = 0;
            numVoxels--;
        }
    }

    clock_t end = clock();
    double elapsed = double(end - begin) / CLOCKS_PER_SEC;

    if(!logPerf) {
        std::cout << "Writing to file ...\n";
        saveVoxelsToObj(argv[argc - 1]);
    }
    //Write out voxel data as obj

    if(!logPerf) {
        std::cout << "Done! " << numVoxels << " voxels actually used. Time elapsed: " << elapsed << " sec.\n";
    }

    if(logPerf) {
        std::cout << numTriangles << "," 
                << dim << "," 
                << (multiDirections ? "true" : "false") << "," 
                << (prescreenTriangles ? "true" : "false") << "," 
                << elapsed << ",\n";
    }
    
    delete g_voxelGrid;
}