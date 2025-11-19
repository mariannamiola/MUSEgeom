#include <iostream>
#include <ostream>
#include "helpVolume.h"


void printHelpVectortoMesh()
{
    std::cout << "\n";
    std::cout << "================ Vector Data Processing Help ================\n";
    std::cout << "This module loads, processes, and converts geospatial vector data\n";
    std::cout << "(POINT, LINESTRING, POLYGON, MULTI) into surface meshes (off, obj).\n\n";

    std::cout << "Available Flags:\n";
    std::cout << "  -V                      Activate vector data processing (.shp, .gpkg, etc.)\n";
    std::cout << "  --tri                   Generate triangular mesh (Trimesh)\n";
    std::cout << "  --grid                  Generate quadrilateral grid mesh (Quadmesh)\n";
    std::cout << "  --polygon               Generate simple polygon mesh\n";
    std::cout << "  --convex                Triangulate points using convex hull\n";
    std::cout << "  --concave               Triangulate points using concave hull\n";
    std::cout << "  --boundary <file>       Use an external boundary file for constrained triangulation\n";
    std::cout << "  --opt <params>          Optional triangulation parameters (e.g., 'q1.2a0.05')\n";
    std::cout << "  --save                  Save extracted boundary/data points before meshing\n";
    std::cout << "  --saveAttributesTable   Export attribute table to CSV (from shapefile)\n";
    std::cout << "  --tolerance <val>       Set tolerance for duplicate point removal\n";
    std::cout << "  --resX --resY           Set grid resolution for Quadmesh generation\n";
    std::cout << "  --newZ <val>            Set fixed Z coordinate value for additional vertices\n";
    std::cout << "  --setEPSG <code>        Apply coordinate transformation to EPSG:<code>\n";
    std::cout << "  --rotate <angles>       Apply data rotation (if configured)\n";
    std::cout << "\n";

    std::cout << "Notes:\n";
    //std::cout << "  * Input can be a single vector file or a directory containing multiple files.\n";
    std::cout << "  * Supported formats: GDAL-vector formats.\n";
    std::cout << "  * Coordinate transformation is applied if the input EPSG differs from project EPSG.\n";
    //std::cout << "  * Multiple subdomains (layers) are automatically detected and processed.\n";
    std::cout << "\n";

    std::cout << "Examples:\n";
    std::cout << "  ./MUSE-geometry -V -v data/shapes --tri --convex --opt q1.2a0.05\n";
    std::cout << "  ./MUSE-geometry -V -v lakes.geojson --polygon --save --saveAttributesTable\n";
    std::cout << "  ./MUSE-geometry -V -v points.shp --tri --concave --tolerance 0.001\n";
    std::cout << "  ./MUSE-geometry -V -v region_dir --grid --resX 100 --resY 80 --newZ 50\n";
    std::cout << "  ./MUSE-geometry -V -v data --tri --boundary boundary.shp --opt q1.1a0.05\n";
    std::cout << "\n";

    std::cout << "=============================================================\n";
    std::cout << std::endl;
}


void printHelpVolume()
{
    std::cout << "\n";
    std::cout << "================ Volume Mesh Generation Help ================\n";
    std::cout << "This module allows you to convert a surface mesh into a 3D volume mesh.\n\n";
    std::cout << "Available Flags:\n";
    std::cout << "  -M                      Activate volume mesh generation\n";
    std::cout << "  -m                      Input surface mesh (only one file allowed)\n";
    std::cout << "  --tet                   Generate tetrahedral mesh using TetGen\n";
    std::cout << "  --vox                   Generate voxelized hexahedral mesh\n";
    std::cout << "  --hex                   Generate structured hexahedral mesh from a surface\n";
    std::cout << "  --opt <params>          Optional TetGen parameters (e.g., 'pq1.2a0.1')\n";
    std::cout << "  --flag <params>         Optional Triangle parameters\n";
    std::cout << "  --save                  Save the translated surface mesh (before tetrahedralization)\n";
    std::cout << "  --maxVoxelPerSide <N>   Set voxel resolution per side (for voxel grid)\n";
    std::cout << "  --resX --resY --resZ    Resolution in X/Y/Z (for structured hex mesh)\n";
    std::cout << "  --plane                 Input plane to cut tetrahedral mesh\n";
    std::cout << "\n";
    std::cout << "Examples:\n";
    std::cout << "  ./MUSE-geometry -M --tet -m model.obj --opt pq1.2a0.1 --vtk\n";
    std::cout << "  ./MUSE-geometry -M --hex -m model.obj --resX 50 --resY 50 --resZ 30\n";
    std::cout << "\n";
    std::cout << "==============================================================\n";
    std::cout << std::endl;
}
