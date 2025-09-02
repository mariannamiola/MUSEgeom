/********************************************************************************
*  This file is part of MUSEgeom                                                *
*  Copyright(C) 2025: Marianna Miola                                            *
*                                                                               *
*  Author:                                                                      *
*                                                                               *
*     Marianna Miola (marianna.miola@cnr.it)                                    *
*                                                                               *
*     Italian National Research Council (CNR)                                   *
*     Institute for Applied Mathematics and Information Technologies (IMATI)    *
*     Via de Marini, 6                                                          *
*     16149 Genoa,                                                              *
*     Italy                                                                     *
*                                                                               *
*  This program is free software: you can redistribute it and/or modify it      *
*  under the terms of the GNU General Public License as published by the        *
*  Free Software Foundation, either version 3 of the License, or (at your       *
*  option) any later version.                                                   *
*                                                                               *
*  This program is distributed in the hope that it will be useful, but          *
*  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY   *
*  or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for  *
*  more details.                                                                *
*                                                                               *
*  You should have received a copy of the GNU General Public License along      *
*  with this program. If not, see <https://www.gnu.org/licenses/>.              *
*                                                                               *
*********************************************************************************/

#include <iostream>
#include <algorithm>
#include <filesystem>
#include <float.h>
#include <set>

#include <tclap/CmdLine.h>

#include <concaveman.h>
// https://adared.ch/concaveman-cpp-a-very-fast-2d-concave-hull-maybe-even-faster-with-c-and-python/
// MODIFICHE:
//1) Comment out make_unique function definition
//2) Add using std::make_unique instead

#include "muselib/utils.h"
#include "muselib/colors.h"
#include "muselib/geometry/mesh.h"
#include "muselib/geometry/tools.h"
#include "muselib/geometry/grid_mesh.h"
#include "muselib/geometry/hexmesh.h"
#include "muselib/geometry/merge_meshes.h"
#include "muselib/geometry/polygon_mesh.h"

#include "muselib/data_structures/point.h"
#include "muselib/data_structures/project.h"
#include "muselib/data_structures/geometry.h"
#include "muselib/data_structures/surface.h"
#include "muselib/data_structures/volume.h"

#include "muselib/metadata/geometry_meta.h"
#include "muselib/metadata/surface_meta.h"
#include "muselib/metadata/volume_meta.h"

#include "muselib/geometry/surface_mesh.h"
#include "muselib/geometry/volume_mesh.h"

//#include "muselib/interpolation/plane.h"

#include "muselib/input/load_vector.h"
#include "muselib/input/load_raster.h"
#include "muselib/input/load_xyz.h"

#include "muselib/reference_system/coordinate_systems.h"

//for filesystem
#ifdef __APPLE__
//#include <filesystem>
using namespace std::__fs;
#else
//#include <experimental/filesystem>
using namespace std;
#endif

#ifdef WIN32
const std::string sep = "\\";
#else
const std::string sep = "/";
#endif

using namespace TCLAP;
using namespace MUSE;

int main(int argc, char** argv)
{
    try {
    CmdLine cmd("MUSEGeom tool", ' ', "version 0.0");

    std::cout << std::endl;
    std::cout << "============================================================" << std::endl;
    std::cout << "================== STARTING MUSE-GEOMETRY ==================" << std::endl;
    std::cout << "============================================================" << std::endl;
    std::cout << std::endl;


    // Option 0. New project creation
    ValueArg<std::string> projectFolder ("p", "pdir", "Project directory", false, "Directory", "path", cmd);

    // Option 0a. Project creation - optional: setting project EPSG
    ValueArg<std::string> setEPSG       ("", "setEPSG", "Set project EPSG", false, "Unknown", "authority", cmd);


    // Option 1. Reading vector file (+ flag for triangulation)
    // Include: shape (.shp), geopackage (.gpkg)
    SwitchArg loadVector                ("V", "vector", "Load Vector file", cmd, false); //booleano
    SwitchArg setSave                   ("", "save", "Saving data content of geospatial files", cmd, false); //booleano
    SwitchArg setSaveAttributesTable    ("", "attribute", "Save attribute table from geospatial file", cmd, false); //booleano


    // Option 2. Reading raster file (+ flag for triangulation)
    // Include: ASCIIGRID (.ASCII)
    SwitchArg loadRaster                ("R", "raster", "Load Raster file", cmd, false); //booleano

    // Option 3. Reading xyz_file (point cloud)
    // Include: yxz, dat, txt
    SwitchArg loadPointCloud            ("P", "pcl", "Load point cloud", cmd, false); //booleano
    ValueArg<std::string> setPoints     ("", "points", "Load filename as POINTS geometry type", false, "", "filename", cmd);
    ValueArg<std::string> setPolygon    ("", "polygon", "Load filename as POLYGON geometry type", false, "", "filename", cmd);

    SwitchArg gridData                  ("G", "gridata", "Grid data - test", cmd, false); //booleano
    MultiArg<std::string> setBBPoints   ("", "bbp", "Set bounding box points", false, "string", cmd );


    ValueArg<std::string> setRotAxis    ("", "rotaxis", "Set rotation axis", false, "NO", "rot_axis", cmd);
    ValueArg<double> setRotAngle        ("", "rotangle", "Set clockwise rotation angle (in degree)", false, 0.0, "double", cmd);
    ValueArg<double> setRotCenterX      ("", "rotcx", "Set coordinte X of rotation center", false, 0.0, "double", cmd);
    ValueArg<double> setRotCenterY      ("", "rotcy", "Set coordinte Y of rotation center", false, 0.0, "double", cmd);
    ValueArg<double> setRotCenterZ      ("", "rotcz", "Set coordinte Z of rotation center", false, 0.0, "double", cmd);


    SwitchArg triFlag                   ("", "tri", "Set triangulation for 2D meshing", cmd, false); //booleano
    SwitchArg convexFlag                ("", "convex", "Set convex hull for points triangulation", cmd, false); //booleano
    SwitchArg concaveFlag               ("", "concave", "Set concave hull for points triangulation", cmd, false); //booleano
    ValueArg<std::string> setBoundary   ("", "boundary", "Set external boundary for points triangulation", false, "", "filename", cmd);
    ValueArg<std::string> optFlag       ("", "opt", "Set optimization flags", false, "", "flag", cmd);

    // Option 4. Set grid
    SwitchArg gridFlag                  ("", "grid", "Set grid for 2D meshing", cmd, false); //booleano
    ValueArg<double> setResx            ("", "resx", "Set x resolution", false, 1.0, "double", cmd);
    ValueArg<double> setResy            ("", "resy", "Set y resolution", false, 1.0, "double", cmd);
    ValueArg<double> setResz            ("", "resz", "Set z resolution", false, 1.0, "double", cmd);

    SwitchArg polygonFlag               ("", "poly", "Set generic polygon mesh for 2D meshing", cmd, false); //booleano
    //ValueArg<std::string> setFeatures   ("", "features", "Set features", false, "DEFAULT", "string" , cmd);

    ValueArg<int> subSet                ("", "subset", "Set (random) subset of points", false, 10,"int", cmd); //booleano

    // Option 5. Set z for new points (if mesh is a section in 3D space)
    // Option: compute variogram with variable/constant lag spacing
    std::vector<std::string> allowedMethod = {"MEAN","CONSTANT","NEAR","KRIGING"};
    ValuesConstraint<std::string> allowedValsMethod(allowedMethod);
    ValueArg<std::string> setMethodZ    ("", "meth", "Set method for z values", false, "CONSTANT", &allowedValsMethod, cmd);
    ValueArg<double> setNewZ            ("", "setz", "Set const z values for new points", false, 0.0, "double", cmd);


    // Option 6. Apply offset on mesh (only extrusion in z direction is enabled)
    SwitchArg setOffset                 ("O", "offset", "Load polygon mesh and apply offset", cmd, false); //booleano
    SwitchArg deltazExtrusion           ("", "delta", "Set DELTA offset", cmd, false); //booleano
    SwitchArg abszExtrusion             ("", "abs", "Set ABSOLUTE ELEVATION offset", cmd, false); //booleano
    ValueArg<double> zOffset            ("z", "zoffset", "Set offset in Z direction", false, 0.0, "double" , cmd);


    SwitchArg appendMeshes              ("A", "append", "Append meshes", cmd, false); //booleano


    // Option 4. Creating triobject: lateral closure of meshes
    SwitchArg createTriObject           ("T", "triobj", "Load trimeshes and create an object closed by surface meshes", cmd, false); //booleano
    MultiArg<std::string> meshFiles     ("m", "mesh", "Set (multi) mesh files", false, "string", cmd );

    SwitchArg createQuadObject          ("Q", "quadobj", "Load quadmeshes and create an object closed by surface meshes", cmd, false); //booleano
    SwitchArg cleanPoly                 ("", "clean", "Clean quadrilateral mesh from isolated polys", cmd, false); //booleano

    // Option 7. Creating volumetric object
    SwitchArg createVolObject           ("M", "volmesh", "Load polygonal mesh and create polyedral mesh", cmd, false); //booleano
    SwitchArg tetFlag                   ("", "tet", "Set tetrahedralization", cmd, false); //booleano
    SwitchArg voxFlag                   ("", "vox", "Set voxel as polyedralmesh", cmd, false); //booleano
    SwitchArg hexFlag                   ("", "hex", "Set hexmesh as polyedral", cmd, false); //booleano

    ValueArg<int> setMaxVoxelperSide    ("", "nmaxvox", "Set n max voxel per side", false, 1, "int", cmd);

    // Option 8. Loading surface mesh
    SwitchArg loadSurface               ("L", "trimesh", "Load trimesh file", cmd, false); //booleano
    ValueArg<std::string> splitMethod   ("", "splmet", "Set polys split method", false, "CENTROID", "string", cmd);
    SwitchArg setRemeshing              ("", "remesh", "Set remeshing", cmd, false); //booleano
    SwitchArg setMarkedEdge             ("", "mark", "Set marked boundary edges for remeshing", cmd, false); //booleano
    SwitchArg setEdgeCollpase           ("", "collapse", "Set collapse on edge to simplify mesh boundary", cmd, false);
    SwitchArg boundaryExtract           ("", "extractbp", "Set extract boundary points", cmd, false); //booleano
    ValueArg<int> setIterations         ("", "it", "Set number of iterations", false, 1.0, "int" , cmd);

    SwitchArg setScaleMesh              ("", "scale", "Set scale mesh", cmd, false);
    ValueArg<double> setScaleFactorX    ("", "sx", "Set scale factor in X direction", false, 1.0, "double" , cmd);
    ValueArg<double> setScaleFactorY    ("", "sy", "Set scale factor in Y direction", false, 1.0, "double" , cmd);
    ValueArg<double> setScaleFactorZ    ("", "sz", "Set scale factor in Z direction", false, 1.0, "double" , cmd);

    // Option 9. Loading volumetric mesh
    SwitchArg loadVolume                ("Z", "tetmesh", "Load tetmesh file", cmd, false); //booleano
    SwitchArg extractSurface            ("", "surf", "Extract surface from volume", cmd, false); //booleano


    // Format conversion for saving meshes
    SwitchArg objConversion             ("", "obj", "Saving mesh in obj format", cmd, false); //booleano
    SwitchArg vtkConversion             ("", "vtk", "Saving mesh in vtk format", cmd, false); //booleano



    // ---------------------------------------------------------------------------------------------------------
    // ADDITIONAL FUNCTIONALITIES:

    // Option 7. Merge two meshes
    SwitchArg mergeMeshes               ("U", "merge", "Merge two trimesh", cmd, false); //booleano
    ValueArg<int> proxThreshold         ("", "proxthresh", "Set proximaty threshold", false, 0, "int" , cmd);

    SwitchArg extractMeshes             ("S", "split", "Split two trimesh", cmd, false); //booleano


    SwitchArg createScalarField             ("F", "cscalar", "Create scalar field from centroids configuration and real samples", cmd, false); //booleano
    ValueArg<std::string> setSamplesMesh    ("","smesh","Set samples mesh associated to (real) values", false, "", "string", cmd);
    ValueArg<std::string> setSamplesValues  ("","sval","Set samples values associated to each vertex of samples mesh", false, "", "string", cmd);

    SwitchArg restoreScalarField            ("", "rscalar", "Restore scalar field from centroids configuration and real samples", cmd, false); //booleano


    //Option 8. Utility to extract id cells from a point cloud



    // Option 9. Multi-resolution approach (associate a scalar field to meshes with different resolutions)
    SwitchArg setMultiResolution            ("D", "res", "Set multiresolution", cmd, false); //booleano
    ValueArg<std::string> setScalarField    ("f", "file", "Set scalar field file", false, "Directory", "path", cmd);
    ValueArg<std::string> setRefModel       ("", "refmod", "Geometry model", false, "name_geometry", "string", cmd);

    ValueArg<std::string> setOutFolder      ("", "outf", "Set folder to save outputs", false, "Directory", "string", cmd);


    // Parse the argv array.
    cmd.parse(argc, argv);

    ///////////////////////////////////////////////

    // 0) Project settings
    MUSE::Project Project;
    Project.setFolder(projectFolder.getValue()); //cartella di progetto
    Project.setName(Project.folder.substr(Project.folder.find_last_of("/")+1, Project.folder.length()));

    // 0) Set folder (in/out)
    std::string in_geometry = Project.folder + "/in";
    std::string out_geometry = Project.folder + "/out";
    std::string out_surf = out_geometry +"/surf";
    std::string out_volume = out_geometry +"/volume";

    std::string command = "";

    // 0) Define file extension - surface
    std::string ext_surf = ".off";
    if(objConversion.isSet() == true)
        ext_surf = ".obj";

    // 0) Define file extension - volume
    std::string ext_vol = ".mesh";
    if(vtkConversion.isSet() == true)
        ext_vol = ".vtk";


    ///
    /// Surface modeling after loading vector file by GDAL
    ///
    if(loadVector.isSet())
    {
        if(!filesystem::exists(out_surf))
            filesystem::create_directory(out_surf);

        MUSE::SurfaceMeta geometa;
        geometa.setProject(Project);

        std::vector<std::string> depsgeom;
        std::vector<std::string> excommands = {command};
        geometa.setCommands(excommands);

        MUSE::GeospatialData Geometry;

        // Check on input files
        if(filesystem::is_empty(in_geometry))
        {
            std::cerr << "\033[0;31mInput ERROR: Insert file into: " << in_geometry << "\033[0m" << std::endl;
            exit(1);
        }

        //Configurazione rotazione
        bool rotation_active = setRotAxis.isSet();
        MUSE::Rotation dataRotation;

        if(rotation_active)
        {
            dataRotation.rotation = true;
            dataRotation.rotation_axis = setRotAxis.getValue();
            dataRotation.rotation_center_x = setRotCenterX.getValue();
            dataRotation.rotation_center_y = setRotCenterY.getValue();
            dataRotation.rotation_center_z = setRotCenterZ.getValue();
            dataRotation.rotation_angle = setRotAngle.getValue();

            geometa.setDataRotation(dataRotation);

            std::cout << "=== Rotation activated: axis=" << dataRotation.rotation_axis << std::endl;
            std::cout << "=== Rotation center: [" << dataRotation.rotation_center_x << "; " << dataRotation.rotation_center_y << "; " << dataRotation.rotation_center_z << "]" <<  std::endl;
            std::cout << "=== Rotation angle (degree): " << dataRotation.rotation_angle << std::endl;
            std::cout << std::endl;
        }

        ////////////////////////////////////////////////
        ///
        /// Lambda per applicare rotazione ai punti
        ///
        auto applyRotation = [&](std::vector<Point3D>& points) {
            if(!rotation_active) return;

            const cinolib::vec3d axis = set_rotation_axis(dataRotation.rotation_axis);
            const cinolib::vec3d center(dataRotation.rotation_center_x,
                                        dataRotation.rotation_center_y,
                                        dataRotation.rotation_center_z);

            for(auto& point : points) {
                cinolib::vec3d sample(point.x, point.y, point.z);
                sample = point_rotation(sample, axis, dataRotation.rotation_angle, center);
                point.x = sample.x();
                point.y = sample.y();
                point.z = sample.z();
            }
            std::cout << "=== Applying rotation ... COMPLETED." << std::endl;
        };

        ///
        /// Lambda per trasformazione coordinate
        ///
        auto applyCoordTransform = [&](std::vector<std::vector<Point3D>>& data_vectors) {
            if(!setEPSG.isSet() || setEPSG.getValue() == Project.authority) return;

            for(auto& vector : data_vectors) {
                for(auto& point : vector) {
                    coordinate_transformation(point.x, point.y, point.z,
                                              setEPSG.getValue(), Project.authority,
                                              point.x, point.y, point.z);
                }
            }
            std::cout << "=== Coordinate transformation ... COMPLETED." << std::endl;
        };

        // Elaborazione directories
        std::vector<std::string> dirs = get_directories(in_geometry);

        if(!dirs.empty())
        {
            for(const auto& dir : dirs)
            {
                auto file_list = get_vectorfiles(dir);
                if(file_list.empty()) continue; //vado alla dir_shape.at(i+1)

                for(const auto& file : file_list)
                {
                    std::vector<std::vector<Point3D>> boundaries, datasets;
                    std::string GDALtype;

                    std::cout << "Loading: " << file << std::endl;

                    if(load_vectorfile(file, boundaries, datasets, GDALtype) != IOSUCCESS)
                    {
                        std::cerr << "\033[0;31mERROR loading: " << file << "\033[0m" << std::endl;
                        exit(1);
                    }

                    // Export se richiesto
                    if(setSaveAttributesTable.isSet())
                    {
                        auto csv_path = out_surf + "/" + get_basename(get_filename(file)) + ".csv";
                        if(export_attributes_to_csv(file, csv_path) == IOSUCCESS)
                        {
                            Geometry.setAttributeTable(get_basename(get_filename(file)) + ".csv");
                        }

                        std::cout << "=== Save Attribute Table (shapefile) ... TO TEST!" << std::endl;
                        exit(1);
                    }

                    if(setSave.isSet())
                    {
                        const std::string basename = get_basename(get_filename(file));
                        const std::string ext = get_extensionND(get_filename(file));

                        if(!boundaries.empty())
                        {
                            for(uint id=0; id<boundaries.size(); id++)
                                export3d_xyz(out_surf + "/" + basename + "_" + std::to_string(id) + "@" + ext + ".dat", boundaries[id]);
                            std::cout << "=== Export content (boundaries) of geospatial file: " << file << std::endl;
                        }
                        if(!datasets.empty())
                        {
                            for(uint id=0; id<datasets.size(); id++)
                                export3d_xyz(out_surf + "/" + basename + "_" + std::to_string(id) + "@" + ext + ".dat", datasets[id]);
                            std::cout << "=== Export content (datasets) of geospatial file: " << file  << std::endl;
                        }
                    }

                    if(!triFlag.isSet() && !gridFlag.isSet() && !polygonFlag.isSet())
                    {
                        std::cout << "=== Meshing algorithms are not selected!" << std::endl;
                        exit(0);
                    }

                    // Applicazione trasformazioni
                    applyCoordTransform(boundaries);
                    applyCoordTransform(datasets);

                    // Setup Geometry
                    Geometry.name = dir.substr(dir.find_last_of("/")+1);
                    Geometry.setFormat(get_extension(file));
                    Geometry.setDomains(boundaries.size());

                    if(!boundaries.empty() && !datasets.empty()) Geometry.geom_type = MULTI;
                    else setGeometryType(Geometry, GDALtype);

                    if(setEPSG.isSet()) Geometry.setAuthority(setEPSG.getValue());

                    std::vector<std::string> deps = {filesystem::relative(file, Project.folder)};
                    geometa.setDependencies(deps);
                    geometa.setGeospatialData(Geometry);

                    // LOGICA ORIGINALE - nessuna lambda per meshing
                    switch(Geometry.geom_type)
                    {
                    case POLYGON:
                    case MULTI:
                    {
                        if(boundaries.empty()) {
                            std::cerr << "\033[0;31mERROR: Polygon requires boundaries!\033[0m" << std::endl;
                            continue;
                        }

                        if(polygonFlag.isSet())
                        {
                            //MUSE::SurfaceMeta::DataSummary dataSummary;

                            //per salvataggio mesh
                            std::string out_mesh = out_surf + "/"+ Geometry.name;
                            out_mesh = out_mesh + ext_surf;

                            MUSE::Surface Surface;
                            MUSE::Surface::Parameters paramSurface;
                            paramSurface.type = "POLYGONMESH";

                            MUSE::Polygonmesh<> polygonmesh (boundaries, polygonmesh);
                            polygonmesh.save(out_mesh.c_str());
                            std::cout << "\033[0;32m=== Saving polygon mesh file: " << out_mesh << "\033[0m" << std::endl;
                            std::cout << std::endl;

                            Surface.setSummary(polygonmesh);
                            Surface.setParameters(paramSurface);

                            geometa.setMeshSummary(Surface);
                            geometa.setGeospatialData(Geometry);
                            geometa.write(out_surf +"/"+ Geometry.name + ".json");
                        }
                        else
                        {
                            for(size_t k=0; k < boundaries.size(); k++)
                            {
                                if(boundaries[k].empty())
                                {
                                    std::cout << "### Domain " << k+1 << " empty - skipping" << std::endl;
                                    continue;
                                }

                                std::vector<Point3D> boundaries_unique, dataset_unique, data_combined;
                                remove_duplicates_test_opt(boundaries[k], boundaries_unique);

                                data_combined = boundaries_unique;
                                if(k < datasets.size() && !datasets[k].empty())
                                {
                                    remove_duplicates_test_opt(datasets[k], dataset_unique);
                                    data_combined.insert(data_combined.end(), dataset_unique.begin(), dataset_unique.end());
                                }

                                applyRotation(boundaries_unique);
                                applyRotation(dataset_unique);
                                applyRotation(data_combined);

                                MUSE::SurfaceMeta::DataSummary dataSummary;
                                dataSummary.setDataSummary(data_combined);
                                geometa.setDataSummary(dataSummary);

                                std::string out_mesh = out_surf + "/" + Geometry.name;
                                if(boundaries.size() > 1)
                                {
                                    Geometry.id_subdomain = std::to_string(k+1);
                                    out_mesh += "_" + Geometry.id_subdomain;
                                }
                                out_mesh += ext_surf;

                                MUSE::Surface Surface;
                                MUSE::Surface::Parameters paramSurface;

                                if(triFlag.isSet())
                                {
                                    if(concaveFlag.isSet())
                                    {
                                        std::cerr << "=== Concave flag is not able for GDALTYPE=POLYGON triangulation." << std::endl;
                                        exit(0);
                                    }
                                    if(convexFlag.isSet())
                                    {
                                        std::cerr << "=== Convex flag is not able for GDALTYPE=POLYGON triangulation." << std::endl;
                                        exit(0);
                                    }

                                    cinolib::Trimesh<> trimesh;
                                    paramSurface.type = "TRIMESH";
                                    paramSurface.opt = optFlag.isSet() ? optFlag.getValue() : "";
                                    paramSurface.boundary = "FIXED BOUNDARY";

                                    if(!dataset_unique.empty())
                                        trimesh = constrained_triangulation2(boundaries_unique, dataset_unique, paramSurface.opt);
                                    else
                                        trimesh = boundary_triangulation(boundaries_unique, paramSurface.opt);

                                    remove_isolate_vertices(trimesh);
                                    Surface.setSummary(trimesh);
                                    trimesh.save(out_mesh.c_str());

                                    std::cout << "\033[0;32m=== Saved: " << out_mesh << "\033[0m" << std::endl;
                                }
                                else if(gridFlag.isSet())
                                {
                                    paramSurface.type = "QUADMESH";
                                    paramSurface.resx = setResx.getValue();
                                    paramSurface.resy = setResy.getValue();
                                    paramSurface.resz = 0.0;

                                    MUSE::Quadmesh<> quadmesh(setResx.getValue(), setResy.getValue(), setNewZ.getValue(), boundaries_unique);
                                    quadmesh.save(out_mesh.c_str());
                                    Surface.setSummary(quadmesh);

                                    std::cout << "\033[0;32m=== Saved: " << out_mesh << "\033[0m" << std::endl;
                                }

                                Surface.setParameters(paramSurface);
                                geometa.setMeshSummary(Surface);
                                geometa.setGeospatialData(Geometry);

                                std::string json_path = boundaries.size() > 1
                                                            ? out_surf + "/" + Geometry.name + "_" + Geometry.id_subdomain + ".json"
                                                            : out_surf + "/" + Geometry.name + ".json";
                                geometa.write(json_path);
                            }

                        }
                        break;
                    }

                    case POINT:
                    case LINESTRING:
                    {
                        if(datasets.empty()) {
                            std::cerr << "\033[0;31mError: No point data\033[0m" << std::endl;
                            continue;
                        }

                        for(size_t k=0; k< datasets.size(); k++)
                        {
                            std::vector<Point3D> data;
                            remove_duplicates_test_opt(datasets[k], data);
                            applyRotation(data);

                            MUSE::SurfaceMeta::DataSummary dataSummary;
                            dataSummary.setDataSummary(data);
                            geometa.setDataSummary(dataSummary);

                            std::string out_mesh = out_surf + "/" + Geometry.name;
                            if(datasets.size() > 1)
                            {
                                Geometry.id_subdomain = std::to_string(k+1);
                                out_mesh += "_" + std::to_string(k+1);
                            }
                            out_mesh += ext_surf;

                            if(triFlag.isSet())
                            {
                                cinolib::Trimesh<> trimesh;
                                MUSE::Surface::Parameters paramSurface;
                                paramSurface.type = "TRIMESH";

                                std::cout << "=== WARNING: Triangulation on XY plane." << std::endl;

                                if (convexFlag.isSet())
                                {
                                    paramSurface.opt = "c";
                                    if(optFlag.isSet()) paramSurface.opt += optFlag.getValue();
                                    paramSurface.boundary = "CONVEX HULL";

                                    trimesh.clear();
                                    trimesh = points_triangulation(data, paramSurface.opt);
                                    remove_isolate_vertices(trimesh);
                                    std::cout << "\033[0;32m=== Convex hull ... COMPLETED.\033[0m" << std::endl;
                                }
                                // CONCAVE HULL - LOGICA ORIGINALE
                                else if (concaveFlag.isSet())
                                {
                                    // 1. Calcolo convex hull e conversione indici
                                    trimesh.clear();
                                    trimesh = points_triangulation(data, "c");
                                    std::vector<int> convexhull;
                                    std::vector<unsigned int> convex_uint = trimesh.get_ordered_boundary_vertices();
                                    for(int i: convex_uint)
                                        convexhull.push_back((short) i);

                                    std::vector<int> b_id;
                                    std::vector<Point3D> concavehull = computing_concavehull(data, convexhull, b_id);

                                    // 2. Rimozione punti boundary dal dataset
                                    std::sort(b_id.begin(), b_id.end());
                                    std::vector<Point3D> unique_data;
                                    for(size_t i=0; i< data.size(); i++)
                                    {
                                        if (!check_index(b_id, i))
                                        {
                                            Point3D unique_p;
                                            unique_p.x = data.at(i).x;
                                            unique_p.y = data.at(i).y;
                                            unique_p.z = data.at(i).z;
                                            unique_data.push_back(unique_p);
                                        }
                                    }

                                    if(optFlag.isSet())
                                        paramSurface.opt = optFlag.getValue();

                                    trimesh.clear();
                                    trimesh = concavehull_triangulation(concavehull, unique_data, paramSurface.opt);
                                    remove_isolate_vertices(trimesh);
                                    paramSurface.boundary = "CONCAVE HULL";

                                    std::cout << "\033[0;32mConcave hull ... COMPLETED.\033[0m" << std::endl;
                                }
                                else if (setBoundary.isSet())
                                {
                                    std::vector<std::vector<Point3D>> boundaries_b, datasets_b;
                                    std::string GDALtype_b;

                                    if(load_vectorfile(setBoundary.getValue(), boundaries_b, datasets_b, GDALtype_b) != IOSUCCESS) {
                                        std::cerr << "\033[0;31mERROR loading boundary: " << setBoundary.getValue() << "\033[0m" << std::endl;
                                        continue;
                                    }

                                    paramSurface.opt = optFlag.isSet() ? optFlag.getValue() : "";
                                    paramSurface.boundary = "EXTERNAL BOUNDARY";

                                    for(const auto& boundary : boundaries_b) {
                                        trimesh = constrained_triangulation2(boundary, data, paramSurface.opt);
                                        remove_isolate_vertices(trimesh);
                                    }
                                }
                                else
                                {
                                    std::cerr << "\033[0;31mERROR: Missing --convex, --concave or --boundary\033[0m" << std::endl;
                                    continue;
                                }

                                MUSE::Surface Surface;
                                Surface.setParameters(paramSurface);
                                Surface.setSummary(trimesh);
                                geometa.setMeshSummary(Surface);

                                trimesh.save(out_mesh.c_str());
                                std::cout << "\033[0;32mSaved: " << out_mesh << "\033[0m" << std::endl;
                            }
                            else if(gridFlag.isSet())
                            {
                                std::cerr << "\033[0;31mGRID FLAG NOT ACTIVE for POINT/LINESTRING!\033[0m" << std::endl;
                                exit(1);
                            }

                            if(datasets.size() > 1)
                                Geometry.id_subdomain = std::to_string(k+1);

                            geometa.setGeospatialData(Geometry);

                            std::string json_path = datasets.size() > 1
                                                        ? out_surf + "/" + Geometry.name + "_" + std::to_string(k+1) + ".json"
                                                        : out_surf + "/" + Geometry.name + ".json";
                            geometa.write(json_path);
                        }
                        break;
                    }
                    }
                }
            }
        }

        // Elaborazione file diretti
        auto file_list = get_vectorfiles(in_geometry);
        for(const auto& file : file_list)
        {
            std::vector<std::vector<Point3D>> boundaries, datasets;
            std::string GDALtype;

            std::cout << "Loading: " << file << std::endl;

            if(load_vectorfile(file, boundaries, datasets, GDALtype) != IOSUCCESS) {
                std::cerr << "\033[0;31mERROR loading: " << file << "\033[0m" << std::endl;
                continue;
            }

            // Export
            if(setSaveAttributesTable.isSet())
            {
                auto csv_path = out_surf + "/" + get_basename(get_filename(file)) + ".csv";
                if(export_attributes_to_csv(file, csv_path) == IOSUCCESS) {
                    Geometry.setAttributeTable(get_basename(get_filename(file)) + ".csv");
                }
            }

            if(setSave.isSet()) {
                const std::string basename = get_basename(get_filename(file));
                const std::string ext = get_extensionND(get_filename(file));

                for(uint id=0; id<boundaries.size(); id++)
                    export3d_xyz(out_surf + "/" + basename + "_" + std::to_string(id) + "@" + ext + ".dat", boundaries[id]);
                for(uint id=0; id<datasets.size(); id++)
                    export3d_xyz(out_surf + "/" + basename + "_" + std::to_string(id) + "@" + ext + ".dat", datasets[id]);
            }

            if(!triFlag.isSet() && !gridFlag.isSet() && !polygonFlag.isSet()) {
                std::cout << "=== No meshing algorithms selected!" << std::endl;
                continue;
            }

            // Trasformazioni
            applyCoordTransform(boundaries);
            applyCoordTransform(datasets);

            // Setup Geometry
            Geometry.setName(get_basename(get_filename(file)));
            Geometry.setFormat(get_extension(file));
            Geometry.setDomains(std::max(boundaries.size(), datasets.size()));

            if(!boundaries.empty() && !datasets.empty()) Geometry.geom_type = MULTI;
            else setGeometryType(Geometry, GDALtype);

            if(setEPSG.isSet()) Geometry.setAuthority(setEPSG.getValue());

            std::vector<std::string> deps = {filesystem::relative(file, Project.folder)};
            geometa.setDependencies(deps);
            geometa.setGeospatialData(Geometry);

            // LOGICA ORIGINALE per meshing
            switch(Geometry.geom_type)
            {
            case POLYGON:
            case MULTI:
            {
                if(polygonFlag.isSet()) {
                    MUSE::Surface::Parameters paramSurface;
                    paramSurface.type = "POLYGONMESH";

                    MUSE::Polygonmesh<> polygonmesh(boundaries, polygonmesh);
                    std::string out_mesh = out_surf + "/" + Geometry.name + ext_surf;
                    polygonmesh.save(out_mesh.c_str());

                    MUSE::Surface Surface;
                    Surface.setSummary(polygonmesh);
                    Surface.setParameters(paramSurface);
                    geometa.setMeshSummary(Surface);
                    geometa.write(out_surf + "/" + Geometry.name + ".json");
                }
                else {
                    for(size_t i=0; i < boundaries.size(); i++)
                    {
                        if(boundaries[i].empty()) continue;

                        std::vector<Point3D> boundaries_unique, dataset_unique, data_combined;
                        remove_duplicates_test_opt(boundaries[i], boundaries_unique);

                        data_combined = boundaries_unique;
                        if(i < datasets.size() && !datasets[i].empty()) {
                            remove_duplicates_test_opt(datasets[i], dataset_unique);
                            data_combined.insert(data_combined.end(), dataset_unique.begin(), dataset_unique.end());
                        }

                        applyRotation(boundaries_unique);
                        applyRotation(dataset_unique);

                        MUSE::SurfaceMeta::DataSummary dataSummary;
                        dataSummary.setDataSummary(data_combined);
                        geometa.setDataSummary(dataSummary);

                        std::string out_mesh = out_surf + "/" + Geometry.name;
                        if(boundaries.size() > 1) {
                            Geometry.id_subdomain = std::to_string(i+1);
                            out_mesh += "_" + Geometry.id_subdomain;
                        }
                        out_mesh += ext_surf;

                        if(triFlag.isSet())
                        {
                            cinolib::Trimesh<> trimesh;
                            MUSE::Surface::Parameters paramSurface;
                            paramSurface.type = "TRIMESH";
                            paramSurface.opt = optFlag.isSet() ? optFlag.getValue() : "";
                            paramSurface.boundary = "FIXED BOUNDARY";

                            if(!dataset_unique.empty())
                                trimesh = constrained_triangulation2(boundaries_unique, dataset_unique, paramSurface.opt);
                            else
                                trimesh = boundary_triangulation(boundaries_unique, paramSurface.opt);

                            remove_isolate_vertices(trimesh);

                            MUSE::Surface Surface;
                            Surface.setSummary(trimesh);
                            Surface.setParameters(paramSurface);
                            geometa.setMeshSummary(Surface);

                            trimesh.save(out_mesh.c_str());
                            std::cout << "\033[0;32mSaved: " << out_mesh << "\033[0m" << std::endl;
                        }

                        geometa.setGeospatialData(Geometry);
                        std::string json_path = boundaries.size() > 1
                                                    ? out_surf + "/" + Geometry.name + "_" + Geometry.id_subdomain + ".json"
                                                    : out_surf + "/" + Geometry.name + ".json";
                        geometa.write(json_path);
                    }
                }
                break;
            }

            case POINT:
            case LINESTRING:
            {
                if(datasets.empty()) {
                    std::cerr << "\033[0;31mError: No point data\033[0m" << std::endl;
                    continue;
                }

                for(size_t i=0; i< datasets.size(); i++)
                {
                    std::vector<Point3D> data;
                    remove_duplicates_test_opt(datasets[i], data);
                    applyRotation(data);

                    MUSE::SurfaceMeta::DataSummary dataSummary;
                    dataSummary.setDataSummary(data);
                    geometa.setDataSummary(dataSummary);

                    std::string out_mesh = out_surf + "/" + Geometry.name;
                    if(datasets.size() > 1) {
                        Geometry.id_subdomain = std::to_string(i+1);
                        out_mesh += "_" + std::to_string(i+1);
                    }
                    out_mesh += ext_surf;

                    if(triFlag.isSet())
                    {
                        cinolib::Trimesh<> trimesh;
                        MUSE::Surface::Parameters paramSurface;
                        paramSurface.type = "TRIMESH";

                        std::cout << "=== WARNING: Triangulation on XY plane." << std::endl;

                        if (convexFlag.isSet())
                        {
                            paramSurface.opt = "c";
                            if(optFlag.isSet()) paramSurface.opt += optFlag.getValue();
                            paramSurface.boundary = "CONVEX HULL";

                            trimesh = points_triangulation(data, paramSurface.opt);
                            remove_isolate_vertices(trimesh);
                            std::cout << "\033[0;32m=== Data triangulation with convex hull ... COMPLETED.\033[0m" << std::endl;
                        }
                        // CONCAVE HULL - CODICE ORIGINALE IDENTICO
                        else if (concaveFlag.isSet())
                        {
                            // 1. Calcolo convex hull e conversione indici
                            trimesh = points_triangulation(data, "c");
                            std::vector<int> convexhull;
                            std::vector<unsigned int> convex_uint = trimesh.get_ordered_boundary_vertices();
                            for(int i: convex_uint)
                                convexhull.push_back((short) i);

                            std::vector<int> b_id;
                            std::vector<Point3D> concavehull = computing_concavehull(data, convexhull, b_id);

                            // 2. Rimozione punti boundary dal dataset
                            std::sort(b_id.begin(), b_id.end());
                            std::vector<Point3D> unique_data;
                            for(size_t i=0; i< data.size(); i++)
                            {
                                if (!check_index(b_id, i))
                                {
                                    Point3D unique_p;
                                    unique_p.x = data.at(i).x;
                                    unique_p.y = data.at(i).y;
                                    unique_p.z = data.at(i).z;
                                    unique_data.push_back(unique_p);
                                }
                            }

                            if(optFlag.isSet())
                                paramSurface.opt = optFlag.getValue();

                            trimesh.clear();
                            trimesh = concavehull_triangulation(concavehull, unique_data, paramSurface.opt);
                            remove_isolate_vertices(trimesh);
                            paramSurface.boundary = "CONCAVE HULL";

                            std::cout << "\033[0;32m=== Data traingulation with concave hull ... COMPLETED.\033[0m" << std::endl;
                        }
                        else if (setBoundary.isSet())
                        {
                            std::vector<std::vector<Point3D>> boundaries_b, datasets_b;
                            std::string GDALtype_b;

                            if(load_vectorfile(setBoundary.getValue(), boundaries_b, datasets_b, GDALtype_b) != IOSUCCESS) {
                                std::cerr << "\033[0;31mERROR loading boundary: " << setBoundary.getValue() << "\033[0m" << std::endl;
                                continue;
                            }

                            paramSurface.opt = optFlag.isSet() ? optFlag.getValue() : "";
                            paramSurface.boundary = "EXTERNAL BOUNDARY";

                            trimesh = constrained_triangulation2(boundaries_b[0], data, paramSurface.opt);
                            std::cout << "\033[0;32m=== Data traingulation with external boundary ... COMPLETED.\033[0m" << std::endl;
                        }
                        else
                        {
                            std::cerr << "\033[0;31mERROR: Missing --convex, --concave or --boundary\033[0m" << std::endl;
                            continue;
                        }

                        MUSE::Surface Surface;
                        Surface.setParameters(paramSurface);
                        Surface.setSummary(trimesh);
                        geometa.setMeshSummary(Surface);

                        trimesh.save(out_mesh.c_str());
                        std::cout << "\033[0;32m=== Saved: " << out_mesh << "\033[0m" << std::endl;
                    }
                    else if(gridFlag.isSet())
                    {
                        std::cerr << "\033[0;31mGRID FLAG NOT ACTIVE for POINT/LINESTRING!\033[0m" << std::endl;
                        exit(1);
                    }

                    if(datasets.size() > 1)
                        Geometry.id_subdomain = std::to_string(i+1);

                    geometa.setGeospatialData(Geometry);

                    std::string json_path = datasets.size() > 1
                                                ? out_surf + "/" + Geometry.name + "_" + std::to_string(i+1) + ".json"
                                                : out_surf + "/" + Geometry.name + ".json";
                    geometa.write(json_path);
                }
                break;
            }
            }
        }

        std::cout << "\033[0;32m=== Vector processing ... COMPLETED.\033[0m" << std::endl;
    }


    ///
    /// Surface modeling from loading raster files
    /// To test: (1) removing convex/concave/boundary flags: the raster is always a grid!
    /// To test: (2) updating JSON information
    ///
    if(loadRaster.isSet())
    {
        MUSE::SurfaceMeta geometa;
        geometa.setProject(Project);

        std::vector<std::string> excommands;
        excommands.push_back(command);
        geometa.setCommands(excommands);

        std::string out_rast = out_geometry +"/rast";
        if(!filesystem::exists(out_rast))
            filesystem::create_directory(out_rast);

        if(filesystem::is_empty(in_geometry))
        {
            std::cerr << "\033[0;31m=== Input ERROR: Insert raster file into: " << in_geometry << "\033[0m" << std::endl;
            exit(1);
        }

        std::vector<std::string> dir_grid = get_directories(in_geometry);
        if(dir_grid.empty())
        {
            std::cout << "=== Input ERROR: no directories found!" << std::endl;
            dir_grid.push_back(in_geometry);
            std::cout << dir_grid.at(0) << std::endl;
        }

        MUSE::GeospatialData Geometry;
        MUSE::Surface Surface;
        MUSE::Surface::Parameters paramSurface;

        for(size_t i=0; i<dir_grid.size(); i++)
        {
            std::vector<std::string> list_grid = get_rasterfiles(dir_grid.at(i));

            if(list_grid.size() == 0)
                continue; //vado alla dir_shape.at(i+1)

            for(size_t j=0; j< list_grid.size(); j++)
            {
                std::vector<std::vector<float>> grid;
                float XOrigin, YOrigin;
                int nXSize, nYSize;
                float XSizePixel = 1.0;
                float YSizePixel = 1.0;

                // Read raster file
                load_rasterfile (list_grid.at(j), grid, XOrigin, YOrigin, nXSize, nYSize, XSizePixel, YSizePixel);

                std::cout << "=== Columns number (nXSize): " << nXSize << ", Rows number (nYSize): " << nYSize << std::endl;
                std::cout << std::fixed << std::setprecision(6) << "=== XOrigin: " << XOrigin << ", YOrigin: " << YOrigin << std::endl;
                std::cout << "=== Grid size: " << grid.size() << " x " << (grid.empty() ? 0 : grid[0].size()) << std::endl;
                std::cout << "=== X Pixel size: " << XSizePixel << ", Y Pixel size: " << YSizePixel << std::endl;
                std::cout << "\033[0;32m=== Import raster file: " << list_grid.at(j) << "... COMPLETED.\033[0m" << std::endl;

                // Set Geometry class
                Geometry.setName(list_grid.at(j).substr(list_grid.at(j).find_last_of("/")+1, list_grid.at(j).length()));
                Geometry.setFormat(get_extension(list_grid.at(j)));

                if(setEPSG.isSet())
                    Geometry.setAuthority(setEPSG.getValue());

                geometa.write(out_rast + "/" + Geometry.getName() + ".json");

                std::vector<Point3D> data, uniq_data;
                for(int row = 0; row < nYSize; row++)
                {
                    for(int col = 0; col < nXSize; col++)
                    {
                        Point3D p;
                        p.x = XOrigin + col * (/* pixel_size_x se disponibile, altrimenti assumere 1.0 */ XSizePixel);
                        p.y = YOrigin + row * (/* pixel_size_y se disponibile, altrimenti assumere 1.0 */ YSizePixel);
                        p.z = grid.at(row).at(col);
                        p.index = row * nXSize + col;

                        data.push_back(p);
                    }
                }
                std::cout << "=== Extract coordinates ... COMPLETED." << std::endl;

                if(triFlag.isSet())
                {
                    cinolib::Trimesh<> trimesh;
                    trimesh.clear();

                    // FOR JSON ...
                    paramSurface.type = "TRIMESH";

                    std::cout << "WARNING: Triangulation is performed on XY plane." << std::endl;

                    // Convex hull
                    if(convexFlag.isSet())
                    {
                        paramSurface.opt = "c";

                        if(optFlag.isSet())
                            paramSurface.opt = paramSurface.opt + optFlag.getValue();

                        remove_duplicates_test_opt(data, uniq_data);
                        trimesh.clear();
                        trimesh = points_triangulation(uniq_data, paramSurface.opt);
                        remove_isolate_vertices(trimesh);

                        paramSurface.boundary = "CONVEX HULL";

                        std::cout << "\033[0;32mTriangulation with convex hull ... COMPLETED.\033[0m" << std::endl;
                    }
                    else if(concaveFlag.isSet())
                    {
                        // 1. Calcolo il convex hull e lo trasformo in int da uint
                        trimesh = points_triangulation(data, "c");
                        std::vector<int> convexhull;
                        std::vector<unsigned int> convex_uint = trimesh.get_ordered_boundary_vertices();
                        for(unsigned int idx : convex_uint)
                            convexhull.push_back((int)idx);

                        std::vector<int> b_id;
                        std::vector<Point3D> concavehull = computing_concavehull(data, convexhull, b_id);

                        // 2. Removing points of concavehull (boundary) from datasets
                        std::sort(b_id.begin(), b_id.end());
                        std::vector<Point3D> unique_data;
                        for(size_t k = 0; k < data.size(); k++)
                        {
                            if(!check_index(b_id, k))
                            {
                                Point3D unique_p;
                                unique_p.x = data.at(k).x;
                                unique_p.y = data.at(k).y;
                                unique_p.z = data.at(k).z;

                                unique_data.push_back(unique_p);
                            }
                        }

                        if(optFlag.isSet())
                            paramSurface.opt = optFlag.getValue();
                        else
                            paramSurface.opt = "";

                        trimesh.clear();
                        trimesh = concavehull_triangulation(concavehull, unique_data, paramSurface.opt);
                        remove_isolate_vertices(trimesh);

                        paramSurface.boundary = "CONCAVE HULL";

                        std::cout << "\033[0;32mTriangulation with concave hull ... COMPLETED.\033[0m" << std::endl;
                    }
                    // External boundary from cmd
                    else if(setBoundary.isSet())
                    {
                        std::vector<Point3D> boundary;
                        load_xyzfile(setBoundary.getValue(), boundary);

                        if(setRotAxis.isSet())
                        {
                            for(size_t k = 0; k < boundary.size(); k++)
                            {
                                cinolib::vec3d sample(boundary.at(k).x, boundary.at(k).y, boundary.at(k).z);
                                cinolib::vec3d axis = set_rotation_axis(setRotAxis.getValue());
                                cinolib::vec3d c(setRotCenterX.getValue(), setRotCenterY.getValue(), setRotCenterZ.getValue());

                                sample = point_rotation(sample, axis, setRotAngle.getValue(), c);

                                boundary.at(k).x = sample.x();
                                boundary.at(k).y = sample.y();
                                boundary.at(k).z = sample.z();
                            }
                        }

                        paramSurface.opt = "";

                        if(optFlag.isSet())
                            paramSurface.opt = optFlag.getValue();

                        trimesh.clear();
                        trimesh = constrained_triangulation2(boundary, data, paramSurface.opt);

                        paramSurface.boundary = "CONSTRAINED";
                    }
                    else
                    {
                        std::cerr << "ERROR: Required argument missing: --convex, --concave or --boundary -m <filename>." << std::endl;
                        exit(1);
                    }

                    Surface.setParameters(paramSurface);
                    Surface.setSummary(trimesh);

                    std::string out_mesh = out_rast + "/" + get_basename(Geometry.getName()) + ext_surf;
                    trimesh.save(out_mesh.c_str());
                }
                else if(gridFlag.isSet())
                {
                    paramSurface.type = "QUADMESH";
                    paramSurface.resx = XSizePixel;
                    paramSurface.resy = YSizePixel;
                    paramSurface.resz = 0.0;

                    //MUSE::Quadmesh<> quadmesh(nYSize-1, nXSize-1, XSizePixel, YSizePixel, XOrigin, YOrigin);
                    MUSE::Quadmesh<> quadmesh(nYSize-1, nXSize-1, XSizePixel, YSizePixel, XOrigin, YOrigin, grid);

                    std::string out_mesh = out_rast + "/grid_" + get_basename(Geometry.getName()) + ext_surf;
                    quadmesh.save(out_mesh.c_str());

                    Surface.setParameters(paramSurface);
                    Surface.setSummary(quadmesh);

                    std::cout << "\033[0;32m=== Saved quadmesh: " << out_mesh << "\033[0m" << std::endl;
                }

                geometa.setMeshSummary(Surface);
                geometa.setGeospatialData(Geometry);

                // CORREZIONE: Usa Geometry.getName() invece di Geometry.name
                geometa.write(out_rast + "/" + Geometry.getName() + ".json");
            }
        }
    }


    } catch (ArgException &e)  // catch exceptions
    { std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; }
}
