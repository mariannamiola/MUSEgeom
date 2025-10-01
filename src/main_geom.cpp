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

#include <tclap/CmdLine.h>

#include <cinolib/tetgen_wrap.h>
#include <cinolib/voxelize.h>
#include <cinolib/voxel_grid_to_hexmesh.h>
#include <cinolib/remesh_BotschKobbelt2004.h>

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

#include "muselib/interpolation/plane.h"

#include "muselib/reference_system/coordinate_systems.h"

#include "muselib/help/helpVolume.h"

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
    ValueArg<std::string> projectFolder ("p", "pdir", "Set project directory", false, "Directory", "path", cmd);
    ValueArg<std::string> inputFolder   ("i", "indir", "Set input directory", false, "Directory", "path", cmd);

    ValueArg<std::string> setEPSG       ("", "setEPSG", "Set project EPSG. TO TEST.", false, "Unknown", "authority", cmd);


    // Option 1. Reading vector file (+ flag for triangulation)
    // Include: shape (.shp), geopackage (.gpkg)
    SwitchArg loadVector                ("V", "vector", "Load Vector file", cmd, false); //booleano
    SwitchArg setSave                   ("", "save", "Saving data content of geospatial files", cmd, false); //booleano
    SwitchArg setSaveAttributesTable    ("", "attribute", "Export attribute table from geospatial file", cmd, false); //booleano


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


    ValueArg<std::string> setAxis       ("", "axis", "Set rotation axis", false, "NO", "rot_axis", cmd);

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
    ValueArg<std::string> optAdditionalFlag ("", "flag", "Set additional optimization flags", false, "", "flag", cmd);

    // Option 4. Set grid
    SwitchArg gridFlag                  ("", "grid", "Set grid for 2D meshing", cmd, false); //booleano
    ValueArg<double> setResx            ("", "resx", "Set x resolution", false, 1.0, "double", cmd);
    ValueArg<double> setResy            ("", "resy", "Set y resolution", false, 1.0, "double", cmd);
    ValueArg<double> setResz            ("", "resz", "Set z resolution", false, 1.0, "double", cmd);

    SwitchArg polygonFlag               ("", "poly", "Set generic polygon mesh for 2D meshing", cmd, false); //booleano

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

    ValueArg<double> setFactor          ("", "f", "Set moltiplication factor", false, 0.1, "double" , cmd); // ad esempio: 10% di un valore

    // Option 9. Loading volumetric mesh
    SwitchArg loadVolume                ("Z", "tetmesh", "Load tetmesh file", cmd, false); //booleano
    SwitchArg extractSurface            ("", "surf", "Extract surface from volume", cmd, false); //booleano


    // Format conversion for saving meshes
    SwitchArg objConversion             ("", "obj", "Saving mesh in obj format", cmd, false); //booleano
    SwitchArg vtkConversion             ("", "vtk", "Saving mesh in vtk format", cmd, false); //booleano

    // Format conversion for saving text file (default: .dat)
    SwitchArg xyzFormat                 ("", "xyz", "Saving text file in xyz format", cmd, false); //booleano
    SwitchArg csvFormat                 ("", "csv", "Saving text file in csv format", cmd, false); //booleano


    // ---------------------------------------------------------------------------------------------------------
    // ADDITIONAL FUNCTIONALITIES:

    // Option 7. Merge two meshes
    SwitchArg mergeMeshes                   ("U", "merge", "Merge two trimesh", cmd, false); //booleano
    ValueArg<double> proxThreshold          ("", "thresh", "Set threshold", false, 0.0, "double" , cmd);

    SwitchArg extractMeshes                 ("S", "split", "Split two trimesh", cmd, false); //booleano


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
    ValueArg<int> setPrecision              ("", "prec", "Set precision", false, 6, "int" , cmd);
    ValueArg<double> setTolerance           ("", "tol", "Set tolerance", false, 1e-02, "double" , cmd);

    // Tetgenerator - integration
    ValueArg<std::string> xyzPlane          ("", "plane", "Plane", false, "plane", "string", cmd);
    ValueArg<double> planeShift             ("", "plane-shift", "Plane shift", false, 0.0, "double", cmd);


    SwitchArg setPerturbation               ("", "perturb", "Set perturbation", cmd, false); //booleano
    SwitchArg setPointManipulation          ("", "manip", "Set mnipulation", cmd, false); //booleano



    /// Help
    SwitchArg helpVol                       ("", "help-vol", "", cmd, false); //booleano

    // Parse the argv array.
    cmd.parse(argc, argv);

    /////////////////////////////////////////////
    /// HELP
    ///

    if (helpVol.isSet()) {
        printHelpVolume();
        return 0;
    }

    ///////////////////////////////////////////////

    // 0) Project settings
    MUSE::Project Project;
    Project.setFolder(projectFolder.getValue()); //cartella di progetto
    Project.setName(Project.folder.substr(Project.folder.find_last_of("/")+1, Project.folder.length()));


    // 0) Commands
    std::string command;

    filesystem::path abspath = argv[3];

    for(int i=1; i< argc; i++)
    {
        std::string string = argv[i];
        if(string.find(abspath) != std::string::npos)
        {
            //std::cout << "Path: " << argv[i] << std::endl;

            filesystem::path path = argv[i];
            filesystem::path relpath = filesystem::relative(path, abspath);
            //std::cout << "Relative path: " << relpath << std::endl;

            if(relpath.string().length() > 1)
                command += "./" + relpath.string();
            else
                command += relpath;
            command += " ";
        }
        else
        {
            command += argv[i];
            command += " ";
        }
    }
    std::cout << "=== Command line: " << command << std::endl;
    std::cout << "=== Number of command arguments: " << argc << std::endl;

    // 0) Set folder (in/out)
    std::string in_geometry = "";
    if(!filesystem::exists(Project.folder + "/in"))
        in_geometry = inputFolder.getValue();
    else
        in_geometry = Project.folder + "/in";

    std::string out_geometry = Project.folder + "/out";
    if(!filesystem::exists(out_geometry))
        filesystem::create_directory(out_geometry);

    std::cout << "=== Absolute path: " << abspath << std::endl;
    std::cout << "=== Input folder: " << in_geometry << std::endl;
    std::cout << "=== Output folder: " << out_geometry << std::endl;
    std::cout << std::endl;

    // 0) Define file extension - surface
    std::string ext_surf = ".off";
    if(objConversion.isSet() == true)
        ext_surf = ".obj";

    // 0) Define file extension - volume
    std::string ext_vol = ".mesh";
    if(vtkConversion.isSet() == true)
        ext_vol = ".vtk";

    // 0) Define file extension - text file
    std::string ext_txt = ".dat";
    if(xyzFormat.isSet() == true)
        ext_txt = ".xyz";
    else if(csvFormat.isSet())
        ext_txt = ".csv";



    std::string out_surf = out_geometry +"/surf";
    std::string out_volume = out_geometry +"/volume";


    ////////////////////////////////////////////////
    ///
    /// Lambda per applicare rotazione ai punti
    ///
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

        std::cout << "=== Rotation activated: axis=" << dataRotation.rotation_axis << std::endl;
        std::cout << "=== Rotation center: [" << dataRotation.rotation_center_x << "; " << dataRotation.rotation_center_y << "; " << dataRotation.rotation_center_z << "]" <<  std::endl;
        std::cout << "=== Rotation angle (degree): " << dataRotation.rotation_angle << std::endl;
        std::cout << std::endl;
    }

    auto applyRotation = [&](std::vector<Point3D>& points)
    {
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

        geometa.setDataRotation(dataRotation);

        MUSE::GeospatialData Geometry;

        // Check on input files
        if(filesystem::is_empty(in_geometry))
        {
            std::cerr << "\033[0;31mInput ERROR: Insert file into: " << in_geometry << "\033[0m" << std::endl;
            exit(1);
        }

        ///
        /// Lambda per trasformazione coordinate
        ///
        auto applyCoordTransform = [&](std::vector<std::vector<Point3D>>& data_vectors)
        {
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

                    std::cout << "=== Loading: " << file << std::endl;

                    printSpatialReferenceInfo(file, Project.authority);
                    if(load_vectorfile(file, boundaries, datasets, GDALtype) != IOSUCCESS)
                    {
                        std::cerr << "\033[0;31m=== ERROR loading: " << file << "\033[0m" << std::endl;
                        exit(1);
                    }

                    // Export se richiesto
                    if(setSaveAttributesTable.isSet())
                    {
                        auto csv_path = out_surf + "/" + get_basename(get_filename(file)) + ext_txt;
                        if(export_attributes_to_csv(file, csv_path) == IOSUCCESS)
                        {
                            Geometry.setAttributeTable(get_basename(get_filename(file)) + ext_txt);
                        }

                        //std::cout << "=== Save Attribute Table (shapefile) ... TO TEST!" << std::endl;
                        //exit(0);
                    }

                    if(setSave.isSet())
                    {
                        const std::string basename = get_basename(get_filename(file));
                        const std::string ext = get_extensionND(get_filename(file));

                        if(!boundaries.empty())
                        {
                            for(uint id=0; id<boundaries.size(); id++)
                                export3d_xyz(out_surf + "/" + basename + "_" + std::to_string(id) + "@" + ext + ext_txt, boundaries[id]);
                            std::cout << "=== Export boundary points from geospatial file: " << file << std::endl;
                        }
                        if(!datasets.empty())
                        {
                            for(uint id=0; id<datasets.size(); id++)
                                export3d_xyz(out_surf + "/" + basename + "_" + std::to_string(id) + "@" + ext + ext_txt, datasets[id]);
                            std::cout << "=== Export data points from geospatial file: " << file  << std::endl;
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
                                remove_duplicates_test_opt(boundaries[k], boundaries_unique, setTolerance.getValue());

                                data_combined = boundaries_unique;
                                if(k < datasets.size() && !datasets[k].empty())
                                {
                                    remove_duplicates_test_opt(datasets[k], dataset_unique, setTolerance.getValue());
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
                            remove_duplicates_test_opt(datasets[k], data, setTolerance.getValue());
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
                                    // for(int i: convex_uint)
                                    //     convexhull.push_back((short) i);
                                    for (unsigned int i : convex_uint)
                                        convexhull.push_back(static_cast<int>(i));

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

            std::cout << "=== Loading: " << file << std::endl;

            printSpatialReferenceInfo(file, Project.authority);
            if(load_vectorfile(file, boundaries, datasets, GDALtype) != IOSUCCESS) {
                std::cerr << "\033[0;31mERROR loading: " << file << "\033[0m" << std::endl;
                continue;
            }

            // if(!boundaries.empty())
            // {
            //     // std::cout << "TRASLAZIONEEEEEEEEEEEEEE" << std::endl;
            //     // double bbox_xmax = -DBL_MAX;
            //     // double bbox_ymax = -DBL_MAX;
            //     // double bbox_xmin = DBL_MAX;
            //     // double bbox_ymin = DBL_MAX;
            //     // for(size_t i=0; i< boundaries.at(0).size(); i++)
            //     // {
            //     //     if(boundaries.at(0).at(i).x > bbox_xmax)
            //     //         bbox_xmax = boundaries.at(0).at(i).x;

            //     //     if(boundaries.at(0).at(i).y > bbox_ymax)
            //     //         bbox_ymax = boundaries.at(0).at(i).y;

            //     //     if(boundaries.at(0).at(i).x < bbox_xmin)
            //     //         bbox_xmin = boundaries.at(0).at(i).x;

            //     //     if(boundaries.at(0).at(i).y < bbox_ymin)
            //     //         bbox_ymin = boundaries.at(0).at(i).y;
            //     // }
            //     // std::cout << "max_x = " << bbox_xmax << "; min_x = " << bbox_xmin << std::endl;
            //     // std::cout << "max_y = " << bbox_ymax << "; min_y = " << bbox_ymin << std::endl;

            //     // double deltax = (bbox_xmax - bbox_xmin)*0.5;
            //     // double deltay = (bbox_ymax - bbox_ymin)*0.5;

            //     // std::cout << deltax << std::endl;

            //     // for(size_t i=0; i< boundaries.at(0).size(); i++)
            //     // {
            //     //     boundaries.at(0).at(i).x -= deltax;
            //     //     boundaries.at(0).at(i).y -= deltay;

            //     //     std::cout << boundaries.at(0).at(i).x << std::endl;
            //     // }
            //     std::cout << "TRASLAZIONE VERSO IL CENTROIDE" << std::endl;

            //     double sum_x = 0.0;
            //     double sum_y = 0.0;
            //     size_t num_points = boundaries.at(0).size();

            //     // Calcolo della somma delle coordinate
            //     for (size_t i = 0; i < num_points; ++i)
            //     {
            //         sum_x += boundaries.at(0).at(i).x;
            //         sum_y += boundaries.at(0).at(i).y;
            //     }

            //     // Centroide = media delle coordinate
            //     double centroid_x = sum_x / num_points;
            //     double centroid_y = sum_y / num_points;

            //     std::cout << "Centroide: (" << centroid_x << ", " << centroid_y << ")" << std::endl;

            //     // Traslazione dei punti verso il centroide (centratura sull'origine)
            //     for (size_t i = 0; i < num_points; ++i)
            //     {
            //         boundaries.at(0).at(i).x -= centroid_x;
            //         boundaries.at(0).at(i).y -= centroid_y;
            //     }
            // }

            // if(!datasets.empty())
            // {
            //     double bbox_xmax = -DBL_MAX;
            //     double bbox_ymax = -DBL_MAX;
            //     double bbox_xmin = DBL_MAX;
            //     double bbox_ymin = DBL_MAX;
            //     for(size_t i=0; i< datasets.size(); i++)
            //     {
            //         if(datasets.at(0).at(i).x > bbox_xmax)
            //             bbox_xmax = datasets.at(0).at(i).x;

            //         if(datasets.at(0).at(i).y > bbox_ymax)
            //             bbox_ymax = datasets.at(0).at(i).y;

            //         if(datasets.at(0).at(i).x < bbox_xmin)
            //             bbox_xmin = datasets.at(0).at(i).x;

            //         if(datasets.at(0).at(i).y < bbox_ymin)
            //             bbox_ymin = datasets.at(0).at(i).y;
            //     }
            //     std::cout << "max_x = " << bbox_xmax << "; min_x = " << bbox_xmin << std::endl;
            //     std::cout << "max_y = " << bbox_ymax << "; min_y = " << bbox_ymin << std::endl;

            //     double deltax = bbox_xmax - bbox_xmin;
            //     double deltay = bbox_ymax - bbox_ymin;

            //     for(size_t i=0; i< datasets.at(0).size(); i++)
            //     {
            //         datasets.at(0).at(i).x += deltax;
            //         datasets.at(0).at(i).y += deltay;
            //     }
            // }


            // Export
            if(setSaveAttributesTable.isSet())
            {
                auto csv_path = out_surf + "/" + get_basename(get_filename(file)) + ext_txt;
                if(export_attributes_to_csv(file, csv_path) == IOSUCCESS) {
                    Geometry.setAttributeTable(get_basename(get_filename(file)) + ext_txt);
                }
            }

            if(setSave.isSet()) {
                const std::string basename = get_basename(get_filename(file));
                const std::string ext = get_extensionND(get_filename(file));

                if(!boundaries.empty())
                {
                    for(uint id=0; id<boundaries.size(); id++)
                        export3d_xyz(out_surf + "/" + basename + "_" + std::to_string(id) + "@" + ext + ext_txt, boundaries[id]);
                    std::cout << "=== Export boundary points from geospatial file: " << file << std::endl;
                }
                if(!datasets.empty())
                {
                    for(uint id=0; id<datasets.size(); id++)
                        export3d_xyz(out_surf + "/" + basename + "_" + std::to_string(id) + "@" + ext + ext_txt, datasets[id]);
                    std::cout << "=== Export data points from geospatial file: " << file  << std::endl;
                }
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
                        remove_duplicates_test_opt(boundaries[i], boundaries_unique, setTolerance.getValue());

                        data_combined = boundaries_unique;
                        if(i < datasets.size() && !datasets[i].empty()) {
                            remove_duplicates_test_opt(datasets[i], dataset_unique, setTolerance.getValue());
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
                        else if(gridFlag.isSet())
                        {
                            MUSE::Surface Surface;

                            MUSE::Surface::Parameters paramSurface;
                            paramSurface.type = "QUADMESH";
                            paramSurface.resx = setResx.getValue();
                            paramSurface.resy = setResy.getValue();
                            paramSurface.resz = 0.0;

                            MUSE::Quadmesh<> quadmesh(setResx.getValue(), setResy.getValue(), setNewZ.getValue(), boundaries_unique);
                            quadmesh.save(out_mesh.c_str());
                            Surface.setSummary(quadmesh);

                            std::cout << "\033[0;32m=== Saved: " << out_mesh << "\033[0m" << std::endl;
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
                    remove_duplicates_test_opt(datasets[i], data, setTolerance.getValue());
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
                            // for(int i: convex_uint)
                            //     convexhull.push_back((short) i);
                            for (unsigned int i : convex_uint)
                                convexhull.push_back(static_cast<int>(i));

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

        std::string out_rast = out_geometry +"/surf";
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
            //std::cout << "=== Input ERROR: no directories found!" << std::endl;
            dir_grid.push_back(in_geometry);
            //std::cout << dir_grid.at(0) << std::endl;
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
                printSpatialReferenceInfo(list_grid.at(j), Project.authority);
                load_rasterfile (list_grid.at(j), grid, XOrigin, YOrigin, nXSize, nYSize, XSizePixel, YSizePixel);

                std::cout << "=== Columns number (nXSize): " << nXSize << ", Rows number (nYSize): " << nYSize << std::endl;
                std::cout << std::fixed << std::setprecision(setPrecision.getValue()) << "=== XOrigin: " << XOrigin << ", YOrigin: " << YOrigin << std::endl;
                std::cout << "=== Grid size: " << grid.size() << " x " << (grid.empty() ? 0 : grid[0].size()) << std::endl;
                std::cout << "=== X Pixel size: " << XSizePixel << ", Y Pixel size: " << YSizePixel << std::endl;
                std::cout << "\033[0;32m=== Import raster file: " << list_grid.at(j) << "... COMPLETED.\033[0m" << std::endl;

                // Set Geometry class
                Geometry.setName(list_grid.at(j).substr(list_grid.at(j).find_last_of("/")+1, list_grid.at(j).length()));
                Geometry.setFormat(get_extension(list_grid.at(j)));

                if(setEPSG.isSet())
                    Geometry.setAuthority(setEPSG.getValue());

                std::vector<Point3D> data, uniq_data;
                for(int row = 0; row < nYSize; row++)
                {
                    for(int col = 0; col < nXSize; col++)
                    {
                        Point3D p;
                        p.x = XOrigin + (col + 0.5) * (/* pixel_size_x se disponibile, altrimenti assumere 1.0 */ XSizePixel);
                        p.y = YOrigin + (row + 0.5) * (/* pixel_size_y se disponibile, altrimenti assumere 1.0 */ YSizePixel);
                        p.z = grid.at(row).at(col);
                        p.index = row * nXSize + col;

                        data.push_back(p);
                    }
                }
                // std::cout << data.at(0).x << "; " << data.at(0).y << "; " << data.at(0).z <<std::endl;
                // std::cout << data.at(1).x << "; " << data.at(1).y << "; " << data.at(1).z <<std::endl;
                MUSE::SurfaceMeta::DataSummary dataSummary;
                dataSummary.setDataSummary(data);
                geometa.setDataSummary(dataSummary);
                std::cout << "=== Extract coordinates of pixel centroids ... COMPLETED." << std::endl;
                std::cout << std::endl;


                ///
                /// Computing bbox data and traslate at the center
                ///
                std::vector<cinolib::vec3d> data_for_bbox;
                for(size_t di=0; di < data.size(); di++)
                    data_for_bbox.push_back(cinolib::vec3d({data.at(di).x, data.at(di).y, data.at(di).z}));
                cinolib::AABB aabb (data_for_bbox);
                data_for_bbox.clear();

                data_for_bbox.push_back(cinolib::vec3d({aabb.min.x(), aabb.min.y(), aabb.min.z()}));
                data_for_bbox.push_back(cinolib::vec3d({aabb.max.x(), aabb.min.y(), aabb.min.z()}));
                data_for_bbox.push_back(cinolib::vec3d({aabb.max.x(), aabb.max.y(), aabb.min.z()}));
                data_for_bbox.push_back(cinolib::vec3d({aabb.min.x(), aabb.max.y(), aabb.min.z()}));
                data_for_bbox.push_back(cinolib::vec3d({aabb.min.x(), aabb.min.y(), aabb.max.z()}));
                data_for_bbox.push_back(cinolib::vec3d({aabb.max.x(), aabb.min.y(), aabb.max.z()}));
                data_for_bbox.push_back(cinolib::vec3d({aabb.max.x(), aabb.max.y(), aabb.max.z()}));
                data_for_bbox.push_back(cinolib::vec3d({aabb.min.x(), aabb.max.y(), aabb.max.z()}));
                std::cout << "=== Computing best plane on points bounding box | vector size: " << data_for_bbox.size() << std::endl;

                cinolib::vec3d center (aabb.min.x() + aabb.delta_x()/2.0, aabb.min.y() + aabb.delta_y()/2.0, aabb.min.z() + aabb.delta_z()/2.0);
                std::cout << "=== Center bbox data: " << center.x() << "; " << center.y() << "; " << center.z() << std::endl;

                for(auto &point:data)
                {
                    point.x -= center.x();
                    point.y -= center.y();
                    point.z -= center.z();
                }

                ///
                /// Starting meshing
                ///
                if(triFlag.isSet())
                {
                    cinolib::Trimesh<> trimesh;


                    // FOR JSON ...
                    paramSurface.type = "TRIMESH";

                    std::cout << "=== WARNING: Triangulation is performed on XY plane." << std::endl;

                    // Convex hull
                    if(convexFlag.isSet())
                    {
                        paramSurface.opt = "c";

                        if(optFlag.isSet())
                            paramSurface.opt = paramSurface.opt + optFlag.getValue();

                        remove_duplicates_test_opt(data, uniq_data, setTolerance.getValue());
                        trimesh.clear();
                        trimesh = points_triangulation(uniq_data, paramSurface.opt);
                        remove_isolate_vertices(trimesh);

                        paramSurface.boundary = "CONVEX HULL";

                        std::cout << "\033[0;32m=== Triangulation with convex hull ... COMPLETED.\033[0m" << std::endl;
                    }
                    else if(concaveFlag.isSet())
                    {
                        // 1. Calcolo il convex hull e lo trasformo in int da uint
                        trimesh.clear();
                        trimesh = points_triangulation(data, "c");

                        std::vector<int> convexhull;
                        std::vector<unsigned int> convex_uint = trimesh.get_ordered_boundary_vertices();
                        // for(unsigned int idx : convex_uint)
                        //     convexhull.push_back((int)idx);
                        for (unsigned int i : convex_uint)
                            convexhull.push_back(static_cast<int>(i));

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

                        geometa.setDataRotation(dataRotation);
                        applyRotation(boundary);

                        // if(setRotAxis.isSet())
                        // {
                        //     for(size_t k = 0; k < boundary.size(); k++)
                        //     {
                        //         cinolib::vec3d sample(boundary.at(k).x, boundary.at(k).y, boundary.at(k).z);
                        //         cinolib::vec3d axis = set_rotation_axis(setRotAxis.getValue());
                        //         cinolib::vec3d c(setRotCenterX.getValue(), setRotCenterY.getValue(), setRotCenterZ.getValue());

                        //         sample = point_rotation(sample, axis, setRotAngle.getValue(), c);

                        //         boundary.at(k).x = sample.x();
                        //         boundary.at(k).y = sample.y();
                        //         boundary.at(k).z = sample.z();
                        //     }
                        // }

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

                    trimesh.translate(center);
                    trimesh.save(out_mesh.c_str());
                }
                else if(gridFlag.isSet())
                {
                    paramSurface.type = "QUADMESH";
                    // paramSurface.resx = XSizePixel;
                    // paramSurface.resy = YSizePixel;

                    paramSurface.resx = setResx.isSet() ? setResx.getValue() : XSizePixel;
                    paramSurface.resy = setResy.isSet() ? setResy.getValue() : YSizePixel;
                    paramSurface.resz = 0.0;

                    std::string out_mesh = out_rast + "/" + get_basename(Geometry.getName()) + ext_surf;
                    std::vector<std::vector<float>> downsampled_grid;
                    if(setResx.isSet() && setResy.isSet())
                    {
                        float corrected_XOrigin = XOrigin;
                        float corrected_YOrigin = YOrigin;

                        /// DIREZIONE DA CONTROLLARE!!!!!!!! TO DO!
                        if((setResx.getValue() > XSizePixel) && (setResy.getValue() > YSizePixel))
                            downsampled_grid = resample_elevation_grid(grid, XSizePixel, YSizePixel, setResx.getValue(), setResy.getValue(), XOrigin, YOrigin, corrected_XOrigin, corrected_YOrigin);

                        uint rows = downsampled_grid.size() - 1;
                        uint cols = downsampled_grid[0].size() - 1;

                        MUSE::Quadmesh<> quadmesh(rows, cols, setResx.getValue(), setResy.getValue(), corrected_XOrigin, corrected_YOrigin, downsampled_grid);
                        quadmesh.save(out_mesh.c_str());
                        Surface.setSummary(quadmesh);
                    }
                    else
                    {
                        MUSE::Quadmesh<> quadmesh(nYSize-1, nXSize-1, XSizePixel, YSizePixel, XOrigin, YOrigin, grid);
                        quadmesh.save(out_mesh.c_str());
                        Surface.setSummary(quadmesh);
                    }

                    Surface.setParameters(paramSurface);
                    std::cout << "\033[0;32m=== Saved quadmesh: " << out_mesh << "\033[0m" << std::endl;
                }

                geometa.setMeshSummary(Surface);
                geometa.setGeospatialData(Geometry);

                geometa.write(out_rast + "/" + get_basename(Geometry.getName()) + ".json");
                std::cout << "\033[0;32m=== Saved json: " << out_rast + "/" + get_basename(Geometry.getName()) + ".json" << "\033[0m" << std::endl;
            }
        }
    }


    ///
    /// Surface modeling by loading point cloud
    ///
    if(loadPointCloud.isSet())
    {
        // // Check on input files (.txt, .dat)
        // if (!setPoints.isSet() && !setPolygon.isSet()) {
        //     if(filesystem::is_empty(in_geometry))
        //     {
        //         std::cerr << "\033[0;31mInput ERROR: Insert file into: " << in_geometry << "\033[0m" << std::endl;
        //         exit(1);
        //     }

        //     std::vector<std::string> file_list = get_xyzfiles(in_geometry);
        //     if (file_list.empty())
        //     {
        //         std::cerr << "\033[0;31mInput ERROR: NO datafile (.txt, .dat, .xyz) in the folder"<< in_geometry << "\033[0m" << std::endl;
        //         exit(1);
        //     }
        // }

        if(!setPoints.isSet() && !setPolygon.isSet())
        {
            std::cout << FRED("ERROR: Set --points <filename> or --polygon <filename>.") << std::endl;
            exit(1);
        }

        if(triFlag.isSet() || gridFlag.isSet())
        {
            if(!filesystem::exists(out_surf))
                filesystem::create_directory(out_surf);
        }
        else
        {
            std::cerr << "=== ERROR: Set --tri/--grid for creating mesh." << std::endl;
            exit(1);
        }

        // Setup metadata
        MUSE::SurfaceMeta geometa;
        geometa.setProject(Project);

        std::vector<std::string> deps;

        std::vector<std::string> excommands;
        excommands.push_back(command);
        geometa.setCommands(excommands);


        //apro i file con la specifica della geometria (NO LOADING AUTOMATICO!!)
        //perchè per le point cloud non ho la possibilità di definire da qualche parte il tipo (come in gdal)
        MUSE::GeospatialData Geometry;
        if(setPoints.isSet())
        {
            Geometry.setName(setPoints.getValue().substr(setPoints.getValue().find_last_of("/")+1, setPoints.getValue().length()));
            Geometry.format = get_extension(Geometry.name);
            Geometry.name = get_basename(Geometry.name);

            Geometry.geom_type = geomType::POINT;

            if(setPoints.getValue().find(abspath) != std::string::npos)
            {
                filesystem::path realpath = filesystem::relative(setPoints.getValue(), abspath);
                deps.push_back(realpath);
            }
        }

        if(setPolygon.isSet())
        {
            Geometry.setName(setPolygon.getValue().substr(setPolygon.getValue().find_last_of("/")+1, setPolygon.getValue().length()));
            Geometry.format = get_extension(Geometry.name);
            Geometry.name = get_basename(Geometry.name);

            Geometry.geom_type = geomType::POLYGON;

            if(setPolygon.getValue().find(abspath) != std::string::npos)
            {
                filesystem::path realpath = filesystem::relative(setPolygon.getValue(), abspath);
                deps.push_back(realpath);
            }
        }

        geometa.setGeospatialData(Geometry);
        geometa.setDependencies(deps);


        //per salvataggio mesh
        std::string out_mesh = out_surf + "/"+ Geometry.name;
        out_mesh = out_mesh + ext_surf;


        // Procedo con la triangolazione in base al tipo di geometria
        MUSE::Surface Surface;
        MUSE::Surface::Parameters paramSurface;

        //il file è di tipo DATA -> POINT
        switch (Geometry.geom_type)
        {
        case POLYGON:
        {
            std::vector<Point3D> boundary, boundary_tmp;
            load_xyzfile(setPolygon.getValue(), boundary);

            if(boundary.size() == 0)
            {
                std::cerr << "ERROR on loading points" << std::endl;
                exit(1);
            }

            ///Check plane orientation
            ///
            align_points_to_xyplane(boundary);

            MUSE::SurfaceMeta::DataSummary dataSummary;
            dataSummary.setDataSummary(boundary);
            geometa.setDataSummary(dataSummary);

            applyRotation(boundary);
            geometa.setDataRotation(dataRotation);

            if(triFlag.isSet())
            {
                cinolib::Trimesh<> trimesh;
                trimesh.clear();

                paramSurface.type = "TRIMESH";
                paramSurface.opt = "";
                paramSurface.boundary = "FIXED BOUNDARY";

                std::cout << "WARNING: Triangulation is performed on XY plane." << std::endl;
                std::cout << "WARNING: Set --rotaxis <axis>, --rotcx <double>, --rotcy <double>, --rotcz <double>, --rotangle <degree> to perform data rotation." << std::endl;
                std::cout << std::endl;

                if(optFlag.isSet())
                    paramSurface.opt = paramSurface.opt + optFlag.getValue();

                trimesh = boundary_triangulation(boundary, paramSurface.opt);

                if(trimesh.num_verts() > boundary.size())
                {
                    std::cout << std::endl;
                    std::cout << "Restore z for additional points ..." << std::endl;
                    if(setMethodZ.getValue().compare("CONSTANT") == 0)
                    {
                        std::cout << "Constant value " << setNewZ.getValue() << " is set for Z additional points." << std::endl;
                        for(uint vid=boundary.size(); vid < trimesh.num_verts(); vid++)
                            trimesh.vert(vid).z() = setNewZ.getValue();
                    }
                    else if (setMethodZ.getValue().compare("MEAN") == 0)
                    {
                        std::cout << "Interpolation method is adopted to set z for additional points ... NOT ACTIVE." << std::endl;
                    //     fittedPlane plane = fitPlane(boundary);
                    //     for(uint vid=boundary.size(); vid < trimesh.num_verts(); vid++)
                    //         trimesh.vert(vid).z() = (trimesh.vert(vid).x()-plane.meanX)*plane.meanA0+(trimesh.vert(vid).y()-plane.meanY)*plane.meanA1 + plane.meanZ;
                    //
                    }
                    else if (setMethodZ.getValue().compare("NEAR") == 0)
                    {
                        std::cerr << "NEAR NOT ACTIVE!" << std::endl;
                        exit(1);
                    }
                    else
                    {
                        std::cerr << "No valid interpolation method is set!" << std::endl;
                        exit(1);
                    }
                    std::cout << "Restore z for additional points ... COMPLETED." << std::endl;
                    std::cout << std::endl;
                }

                remove_isolate_vertices(trimesh);

                Surface.setSummary(trimesh);
                trimesh.save(out_mesh.c_str());
            }

            else if(gridFlag.isSet())
            {
                std::cout << std::endl;
                std::cout << "Resolution in X direction: " << setResx.getValue() << std::endl;
                std::cout << "Resolution in Y direction: " << setResy.getValue() << std::endl;
                std::cout << "Set --resx <value>, --resy <value> to modify default resolutions." << std::endl;
                std::cout << std::endl;

                //FOR JSON ...
                paramSurface.type = "QUADMESH";
                paramSurface.resx = setResx.getValue();
                paramSurface.resy = setResy.getValue();
                paramSurface.resz = 0.0;

                MUSE::Quadmesh<> quadmesh (setResx.getValue(), setResy.getValue(), setNewZ.getValue(), boundary);
                quadmesh.save(out_mesh.c_str());

                Surface.setSummary(quadmesh);
            }

            std::cout << std::endl;
            std::cout << "\033[0;32mSaving mesh file: " << out_mesh << "\033[0m" << std::endl;
            std::cout << std::endl;

            Surface.setParameters(paramSurface);

            geometa.setMeshSummary(Surface);
            geometa.setGeospatialData(Geometry);
            geometa.write(out_surf +"/"+ Geometry.name + ".json");

            break;
        }
        case MULTI:
        {
            std::cout << FYEL("MULTI CASE: NOTING TO DO.") << std::endl;
            break;
        }
        case POINT:
        {
            std::vector<Point3D> data, uniq_data;
            load_xyzfile(setPoints.getValue(), data);

            if(data.size() == 0)
            {
                std::cerr << "ERROR on loading points" << std::endl;
                exit(1);
            }

            MUSE::SurfaceMeta::DataSummary dataSummary;
            dataSummary.setDataSummary(data);
            geometa.setDataSummary(dataSummary);

            applyRotation(data);
            geometa.setDataRotation(dataRotation);

            remove_duplicates_test_opt(data, uniq_data, setTolerance.getValue());

            if(proxThreshold.isSet() && setAxis.isSet())
            {
                std::cout << "=== Filter on " << setAxis.getValue() << " values is set. Threshold = " << proxThreshold.getValue() << std::endl;
                std::cout << "=== Intial data dimension: " << uniq_data.size() << std::endl;
                std::vector<Point3D> data_to_filter = uniq_data;
                uniq_data.clear();
                for(size_t i=0; i<data_to_filter.size(); i++)
                {
                    if(setAxis.getValue() == "Z" && data_to_filter.at(i).z <= proxThreshold.getValue())
                        uniq_data.push_back(data_to_filter.at(i));
                }
                std::cout << "=== Filtered data dimension: " << uniq_data.size() << std::endl;
            }

            if(subSet.isSet())
            {
                srand(time(NULL));
                std::vector<size_t> random_id(subSet.getValue());
                for (size_t i = 0; i < subSet.getValue(); i++)
                {
                    random_id[i] = rand() % uniq_data.size();
                    //std::cout << "rand ID: " << random_id[i] << std::endl;
                }

                std::sort(random_id.begin(), random_id.end());
                random_id.erase(std::unique( random_id.begin(), random_id.end() ), random_id.end() );

                std::vector<Point3D> data_rand=uniq_data;
                uniq_data.clear();
                for(int rid:random_id)
                    uniq_data.push_back(data_rand.at(rid));
                std::cout << "### Size of data vector (before random sampling): " << data_rand.size() << std::endl;
                std::cout << "### New size of data vector (after random sampling): " << uniq_data.size() << std::endl;
                std::cout << std::endl;

                std::string filename_rand = "_subset" + ext_txt;
                export3d_xyz(out_surf + "/" + filename_rand, uniq_data);
            }

            /// Questa sezione permette di manipolare un set di dati:
            /// dopo aver caricato una nuvola di punti e sottoposta ad eventuale filtro,
            /// creo un piano passante per quei punti
            /// carico una serie di punti di bordo e li proietto sul piano
            /// salvo il dataset + il bordo proiettato sul piano
            if(setPointManipulation.isSet() && setBoundary.isSet())
            {
                std::cout << "=== Starting point manipulation ... "<< std::endl;
                double plane_min_z = DBL_MAX;
                double plane_max_z = -DBL_MAX;

                std::vector<cinolib::vec3d> points, bpoints, bverts_proj;
                for(auto pd:uniq_data)
                {
                    cinolib::vec3d cp (pd.x, pd.y, pd.z);
                    points.push_back(cp);

                    if(cp.z() < plane_min_z) plane_min_z = cp.z();
                    if(cp.z() > plane_max_z) plane_max_z = cp.z();
                }

                std::cout << "plane_min_z: " << plane_min_z << std::endl;
                std::cout << std::fixed << std::setprecision(setPrecision.getValue()) << "plane_max_z: " << plane_max_z << std::endl;

                std::cout << "Size vector uniq_data: " << uniq_data.size() << std::endl;
                std::cout << "Size vector points: " << points.size() << std::endl;
                std::cout << "=== Point manipulation tools is set." << std::endl;

                cinolib::Plane plane_frompoints (points);
                std::ifstream pboundary;
                pboundary.open(setBoundary.getValue());

                double x,y,z;
                while (pboundary >> x >> y >> z)
                    bpoints.push_back(cinolib::vec3d(x,y,z));

                std::cout << "=== Reading dataset from file: " << setBoundary.getValue() << " COMPLETED." << std::endl;
                std::cout << "Size vector new dataset: " << bpoints.size() << std::endl;

                // Proietta ogni punto del bordo sul piano
                int last_best_idx = 0;

                // for (auto pvid : bpoints)
                // {
                //     cinolib::vec3d intersection;
                //     cinolib::vec3d p2 = pvid;
                //     p2.z() = plane_min_z;

                //     cinolib::Segment segm (0, pvid, p2);
                //     std::cout << "p2: " << p2 << std::endl;
                //     std::cout << "pvid: " << pvid << std::endl;

                //     if (!intersectPlaneSegment(plane_frompoints, segm, intersection)) //.vert(vid));
                //     {
                //         std::cout << FRED("=== ERROR: No intersection is found! ") << intersection << std::endl;
                //         continue;
                //     }
                //     //std::cout << "FOUND INTERSECTION: " << intersection << std::endl;
                //     // Usa il region growing
                //     last_best_idx = find_nearest_in_local_region(intersection, points, last_best_idx, 20);
                //     double z = points[last_best_idx].z();

                //     cinolib::vec3d adjusted(intersection.x(), intersection.y(), z);
                //     bverts_proj.push_back(adjusted);
                // }
                // std::cout << "Size vector points intersection: " << bverts_proj.size() << std::endl;


                for (const auto &pvid : bpoints)
                {
                    cinolib::vec3d projected = projectPointOntoPlane(pvid, plane_frompoints);

                    // Usa il region growing
                    last_best_idx = find_nearest_in_local_region(projected, points, last_best_idx, 20);
                    double z = points[last_best_idx].z();
                    cinolib::vec3d adjusted(projected.x(), projected.y(), z);

                    bverts_proj.push_back(adjusted);
                }

                std::cout << "Size vector points intersection: " << bverts_proj.size() << std::endl;

                std::cout << std::endl;

                for(auto pnew:bverts_proj)
                    points.push_back(pnew);

                std::cout << "Size vector points and intersections: " << points.size() << std::endl;

                std::string filename_new = "_newset" + ext_txt;
                std::ofstream file_out;
                file_out.open(out_surf + "/" + filename_new, std::fstream::out);

                for(size_t i=0; i<points.size(); i++)
                    file_out << std::fixed << std::setprecision(setPrecision.getValue()) << points.at(i).x() << " " << points.at(i).y() << " " << points.at(i).z() << std::endl;

                file_out.close();
            }


            if(triFlag.isSet())
            {
                cinolib::Trimesh<> trimesh;

                //FOR JSON ...
                paramSurface.type = "TRIMESH";

                std::cout << "WARNING: Triangulation is performed on XY plane." << std::endl;

                //Convex hull
                if (convexFlag.isSet())
                {
                    //remove_duplicates_test_opt(data, uniq_data);

                    paramSurface.opt = "c";

                    if(optFlag.isSet())
                        paramSurface.opt = paramSurface.opt + optFlag.getValue();

                    trimesh.clear();
                    trimesh = points_triangulation(uniq_data, paramSurface.opt);

                    if(trimesh.num_verts() > uniq_data.size())
                    {
                        std::cout << std::endl;
                        std::cout << "Restore z for additional points ..." << std::endl;
                        if(setMethodZ.getValue().compare("CONSTANT") == 0)
                        {
                            std::cout << "Constant value " << setNewZ.getValue() << " is set for Z additional points." << std::endl;
                            for(uint vid=uniq_data.size(); vid < trimesh.num_verts(); vid++)
                                trimesh.vert(vid).z() = setNewZ.getValue();
                        }
                        else if (setMethodZ.getValue().compare("MEAN") == 0)
                        {
                            std::cout << "Interpolation method is adopted to set z for additional points ... NOT ACTIVE." << std::endl;
                            // fittedPlane plane = fitPlane(uniq_data);
                            // for(uint vid=uniq_data.size(); vid < trimesh.num_verts(); vid++)
                            //     trimesh.vert(vid).z() = (trimesh.vert(vid).x()-plane.meanX)*plane.meanA0+(trimesh.vert(vid).y()-plane.meanY)*plane.meanA1 + plane.meanZ;
                        }
                        else if (setMethodZ.getValue().compare("NEAR") == 0)
                        {
                            std::cout << "Interpolation method is adopted to set z for additional points" << std::endl;
                            for(uint vid=uniq_data.size(); vid < trimesh.num_verts(); vid++)
                            {
                                std::vector<uint> adj_vert = trimesh.adj_v2v(vid);
                                double mean, sum=0.0;
                                int count=0;
                                for(uint bv:adj_vert)
                                {
                                    if(trimesh.vert(bv).z() != 0.0)
                                    {
                                        sum+=trimesh.vert(bv).z();
                                        count++;
                                    }
                                    //std::cout << sum << std::endl;
                                }
                                mean=sum/count;
                                std::cout << mean << std::endl;
                                trimesh.vert(vid).z() = mean;
                            }
                        }
                        else
                        {
                            std::cout << "No valid interpolation method is set." << std::endl;
                            exit(1);
                        }
                        std::cout << "Restore z for additional points ... COMPLETED." << std::endl;
                        std::cout << std::endl;
                    }

                    remove_isolate_vertices(trimesh);

                    paramSurface.boundary = "CONVEX HULL";

                    std::cout << "\033[0;32mTriangulation with convex hull ... COMPLETED.\033[0m" << std::endl;
                }

                else if (concaveFlag.isSet())
                {
                    //remove_duplicates_test_opt(data, uniq_data);

                    // 1. Calcolo il convex hull (passando per la triangolazione dei punti) e lo trasformo in int da uint
                    trimesh = points_triangulation(uniq_data, "c");
                    std::vector<int> convexhull;
                    std::vector<unsigned int> convex_uint = trimesh.get_ordered_boundary_vertices();
                    // for(int i: convex_uint)
                    //     convexhull.push_back((short) i);
                    for (unsigned int i : convex_uint)
                        convexhull.push_back(static_cast<int>(i));

                    std::vector<int> b_id;
                    std::vector<Point3D> concavehull = computing_concavehull(uniq_data, convexhull, b_id);

                    // 2. Removing points of concavehull (boundary) from datasets
                    std::sort(b_id.begin(), b_id.end());
                    std::vector<Point3D> unique_data;
                    for(size_t i=0; i< uniq_data.size(); i++)
                    {
                        if (!check_index(b_id, i))
                        {
                            Point3D unique_p;
                            unique_p.x = uniq_data.at(i).x;
                            unique_p.y = uniq_data.at(i).y;
                            unique_p.z = uniq_data.at(i).z;

                            unique_data.push_back(unique_p);
                        }
                    }

                    if(optFlag.isSet())
                        paramSurface.opt = paramSurface.opt + optFlag.getValue();

                    trimesh.clear();
                    trimesh = concavehull_triangulation(concavehull, unique_data, paramSurface.opt);
                    remove_isolate_vertices(trimesh);

                    paramSurface.boundary = "CONCAVE HULL";

                    std::cout << "\033[0;32mTriangulation with concave hull ... COMPLETED.\033[0m" << std::endl;
                }

                // External boundary from cmd
                else if (setBoundary.isSet()) //se gli passo da linea di comando un bordo esterno: 1) leggi 2) triangola i punti vincolati al bordo
                {
                    //remove_duplicates_test_opt(data, uniq_data);
                    //uniq_data=data;
                    std::cout << std::endl;

                    std::vector<Point3D> boundary, uniq_boundary;
                    load_xyzfile(setBoundary.getValue(), boundary);
                    remove_duplicates_test_opt(boundary, uniq_boundary, setTolerance.getValue());

                    applyRotation(uniq_boundary);

                    paramSurface.opt = "";

                    if(optFlag.isSet())
                        paramSurface.opt = paramSurface.opt + optFlag.getValue();

                    trimesh.clear();
                    trimesh = constrained_triangulation2(uniq_boundary, uniq_data, paramSurface.opt);

                }
                else
                {
                    std::cerr << "ERROR: Required argument missing: --convex, --concave or --boundary -m <filename>." << std::endl;
                    break;
                }

                Surface.setParameters(paramSurface);
                Surface.setSummary(trimesh);

                trimesh.save(out_mesh.c_str());
            }
            else if (gridFlag.isSet())
            {
                std::cerr << FRED("GRID FLAG IS NOT ACTIVE!!") << std::endl;
                exit(1);
            }

            geometa.setMeshSummary(Surface);
            geometa.setGeospatialData(Geometry);

            geometa.write(out_surf + "/"+ Geometry.name + ".json");

            break;
        }
        case LINESTRING:
            break;
        }
    }


    ///
    /// Creating an offset of a surface
    ///
    if(setOffset.isSet())
    {
        if(!meshFiles.isSet())
        {
            std::cout << FRED("ERROR. Set mesh to extrude by -m or --mesh command.") << std::endl;
            exit(1);
        }

        if(meshFiles.getValue().size() > 1)
        {
            std::cout << FRED("ERROR. Only one mesh is supported for the extrusion.") << std::endl;
            exit(1);
        }

        if (!zOffset.isSet())
        {
            std::cerr << FRED("ERROR: Missing zOffset value. Use -z (or --zoffset) <value>.") << std::endl;
            exit(1);
        }

        std::cout << FMAG("##########################################") << std::endl;
        std::cout << FMAG("WARNING: Only z-direction extrusion is enabled!") << std::endl;
        std::cout << FMAG("##########################################") << std::endl;
        std::cout << std::endl;

        std::vector<std::string> files = meshFiles.getValue();

        MUSE::SurfaceMeta georef;
        //georef.read(get_basename(files.at(0)) + ".json");
        std::string json_path = filesystem::path(files.at(0)).replace_extension(".json").string();
        georef.read(json_path);

        MUSE::SurfaceMeta geometa;
        geometa.setProject(Project);

        std::vector<std::string> excommands;
        excommands.push_back(command);
        geometa.setCommands(excommands);

        std::vector<std::string> deps;
        deps.push_back(filesystem::relative(get_basename(files.at(0)) + ".json", Project.folder));
        geometa.setDependencies(deps);

        geometa.setDataSummary(georef.getDataSummary());
        geometa.setDataRotation(georef.getDataRotation());

        // Loading polygonmesh
        cinolib::Polygonmesh <> mesh;
        mesh.load(files.at(0).c_str());

        std::cout << "\033[0;32m=== Loading mesh file: " << files.at(0) << " ... COMPLETED.\033[0m" << std::endl;

        std::string basename = get_basename(get_filename(files.at(0)));
        std::string out_name = get_basename(files.at(0)); //out_surf +"/" + basename;

        MUSE::SurfaceMeta::Extrusion objinfo;
        if(deltazExtrusion.isSet()) //estrusione con delta costante (z + delta)
        {
            objinfo.type = "delta";
            std::cout << "=== The extrusion is set on delta ... " << std::endl;

            if(zOffset.isSet())
            {
                std::cout << "=== The extrusion value is set on " << zOffset.getValue() << " in z direction ..." << std::endl;
                objinfo.value = zOffset.getValue();
                objinfo.direction = "z";

                for(size_t vid =0 ; vid < mesh.num_verts(); vid++)
                    mesh.vert(vid).z() = mesh.vert(vid).z() + zOffset.getValue();

                //out_name = out_name + "dz" + std::to_string(zOffset.getValue());
                out_name = out_name + "_d" + objinfo.direction;
            }
        }
        else if(abszExtrusion.isSet()) //estrusione fino ad una quota assoluta fissata
        {
            objinfo.type = "absolute elevation";
            std::cout << "=== The extrusion is set on absolute elevation ... " << std::endl;

            if(zOffset.isSet())
            {
                std::cout << "=== The extrusion value is set on " << zOffset.getValue() << " in z direction ..." << std::endl;
                objinfo.value = zOffset.getValue();
                objinfo.direction = "z";

                for(uint pid=0; pid < mesh.num_polys(); pid++)
                {
                    std::vector<cinolib::vec2d> vec2d;
                    for(size_t i=0; i < mesh.poly_verts(pid).size(); i++)
                    {
                        cinolib::vec2d v;
                        v.x() = mesh.poly_verts(pid).at(i).x();
                        v.y() = mesh.poly_verts(pid).at(i).y();
                        vec2d.push_back(v);
                    }
                    //std::cout << polygon_is_CCW(vec2d) << std::endl;
                    //std::cout << polygon_signed_area(vec2d) << std::endl;
                    // se l'area è positiva, quindi ccw = 1 (true) -> normale uscente
                    if(cinolib::polygon_is_CCW(vec2d) == false)
                        mesh.poly_flip_winding_order(pid);
                }
                //extr_trimesh = trimesh;
                for(size_t vid =0 ; vid < mesh.num_verts(); vid++)
                    mesh.vert(vid).z() = zOffset.getValue();

                //out_name = out_name + "absz" + std::to_string(zOffset.getValue());
                out_name = out_name + "_abs" + objinfo.direction;
            }
        }
        else
        {
            std::cerr << "\033[0;31mERROR: Required argument missing: --delta or --abs for surface extrusion in z direction.\033[0m" << std::endl;
            exit(1);
        }

        mesh.save((out_name + ext_surf).c_str());
        std::cout << "\033[0;32mExport mesh file: " << out_name + ext_surf << " ... COMPLETED.\033[0m" << std::endl;


        MUSE::Surface summary;
        summary.setSummary(mesh);

        geometa.setExtrusion(objinfo);
        geometa.setMeshSummary(summary);

        geometa.write(out_name + ".json");
    }


    ///
    /// Closing two (triangular) surface meshes (GEO3D approach)
    ///
    if(createTriObject.isSet() && meshFiles.getValue().size() == 2)
    {
        std::cout << FMAG("##########################################") << std::endl;
        std::cout << FMAG("NOTA BENE: la chiusura delle mesh avviene per ora solo per mesh triangolari e in direzione z") << std::endl;
        std::cout << FMAG("##########################################") << std::endl;
        std::cout << std::endl;

        std::vector<std::string> files = meshFiles.getValue();

        cinolib::Trimesh<> trimesh0;
        trimesh0.load(files.at(0).c_str());
        std::cout << "\033[0;32mLoading mesh file: " << files.at(0) << " ... COMPLETED.\033[0m" << std::endl;
        std::string filename0 = files.at(0).substr(files.at(0).find_last_of("/")+1, files.at(0).length());
        std::string basename0 = get_basename(filename0);

        cinolib::Trimesh<> trimesh1;
        trimesh1.load(files.at(1).c_str());
        std::cout << "\033[0;32mLoading mesh file: " << files.at(1) << " ... COMPLETED.\033[0m" << std::endl;
        std::string filename1 = files.at(1).substr(files.at(1).find_last_of("/")+1, files.at(1).length());
        std::string basename1 = get_basename(filename1);

        std::cout << "bb1 completa: " << trimesh1.bbox()<< std::endl;
        std::cout << "bb0 completa: " << trimesh0.bbox()<< std::endl;


        //Creazione json del triobject
        MUSE::SurfaceMeta geometa;
        geometa.setProject(Project);

        std::vector<std::string> excommands;
        excommands.push_back(command);
        geometa.setCommands(excommands);

        std::vector<std::string> deps;
        deps.push_back(filesystem::relative(get_basename(files.at(0)) + ".json", Project.folder));
        deps.push_back(filesystem::relative(get_basename(files.at(1)) + ".json", Project.folder));

        geometa.setDependencies(deps);


        // 5) Creazione superficie laterale
        // 5.1) Vector indici dei punti sul convex hull nelle rispettive mesh
        std::vector<uint> or_idch0, or_idch1;
        or_idch0 = trimesh0.get_ordered_boundary_vertices();
        or_idch1 = trimesh1.get_ordered_boundary_vertices();

        size_t n = 0; //n vertici uguali
        for(unsigned int vid1 : or_idch1)
        {
            cinolib::vec2d vert2d_1 (trimesh1.vert(vid1).x(), trimesh1.vert(vid1).y());
            for(unsigned int vid0 : or_idch0)
            {
                cinolib::vec2d vert2d_0 (trimesh0.vert(vid0).x(), trimesh0.vert(vid0).y());
                if(vert2d_0.dist(vert2d_1) < 1e-2)
                {
                    n++;
                    break;
                }
            }
        }

        std::cout << std::endl;
        std::cout << or_idch0.size() << " boundary points of: " << basename0 << std::endl;
        std::cout << or_idch1.size() << " boundary points of: " << basename1 << std::endl;
        std::cout << n << " equal boundary points between " << basename0 << " and " << basename1 << std::endl;
        std::cout << FGRN("Check on boundary ... COMPLETED.") << std::endl;
        std::cout << std::endl;


        // INTEGRAZIONE CODICE DI GEO3D -> VALIDA E FUNZIONANTE PER PUNTI TRIANGOLATI CON CONVEX HULL! DA ESTENDERE CON BOUNDARY/CONCAVE
        // TO DO ...
        if(n != or_idch0.size()) //CONDIZIONE SU BORDI UGUALE
        {
            std::cout << "\033[0;31mERROR: Meshes boundaries are different!\033[0m" << std::endl;

            std::cout << FMAG("##########################################") << std::endl;
            std::cout << FMAG("RIFERIMENTO: implementazione di GEO3D") << std::endl;
            std::cout << FMAG("##########################################") << std::endl;
            std::cout << std::endl;


            // Criterio di scelta del bordo: area mesh
            std::vector<double> polygon_area (2);
            polygon_area.at(0) = trimesh0.mesh_area();
            polygon_area.at(1) = trimesh1.mesh_area();

            // Calcolo mesh con area minima
            double min_area = DBL_MAX;
            uint min_index = 0;
            for(uint i=0; i<polygon_area.size(); i++)
            {
                if(polygon_area.at(i) < min_area)
                {
                    min_area = polygon_area.at(i);
                    min_index = i;
                }
            }
            std::cout << "Area minima: " << polygon_area.at(min_index) << std::endl;

            cinolib::Trimesh<> tmp;
            std::string basename_tmp;
            if(min_index != 0)
            {
                tmp = trimesh0;
                trimesh0 = trimesh1;
                trimesh1 = tmp;

                basename_tmp = basename0;
                basename0 = basename1;
                basename1 = basename_tmp;
                //std::cout << "ORDINE MESH INVERTITO" << std::endl;
            }
            tmp.clear();

            //estrazione ch di riferimento
            std::vector<Point2D> ref_ch;
            for(uint i : trimesh0.get_ordered_boundary_vertices())
            {
                Point2D p;
                p.x = trimesh0.vert(i).x();
                p.y = trimesh0.vert(i).y();
                //p.z = trimesh0.vert(i).y();
                ref_ch.push_back(p);
            }

            std::vector<cinolib::vec3d> verts2 = trimesh1.vector_verts();
            std::vector<Point3D> sub_verts2;
            for(unsigned int i=0; i< verts2.size(); i++)
            {
                Point2D p;
                p.x = verts2.at(i).x();
                p.y = verts2.at(i).y();

                bool internal = point_in_polygon(p, ref_ch);
                if(internal)
                {
                    Point3D p3d;
                    p3d.x = verts2.at(i).x();
                    p3d.y = verts2.at(i).y();
                    p3d.z = verts2.at(i).z();

                    sub_verts2.push_back(p3d);
                }
            }
            std::cout << sub_verts2.size() << " vertices into " << basename0 << " mesh convex hull." << std::endl;

            trimesh1.clear();
            trimesh1 = points_triangulation(sub_verts2, "c");


            std::vector<Point3D> ch1_tmp, ch2_tmp, ch;
            for(uint i : trimesh0.get_ordered_boundary_vertices())
            {
                Point3D p;
                p.x = trimesh0.vert(i).x();
                p.y = trimesh0.vert(i).y();
                p.z = trimesh0.vert(i).z();
                ch1_tmp.push_back(p);
            }

            for(uint i : trimesh1.get_ordered_boundary_vertices())
            {
                Point3D p;
                p.x = trimesh1.vert(i).x();
                p.y = trimesh1.vert(i).y();
                p.z = trimesh1.vert(i).z();
                ch2_tmp.push_back(p);
            }

            ch = ch1_tmp;
            ch.insert(ch.end(), ch2_tmp.begin(), ch2_tmp.end());


            // ////////////////////////////////////////////////////////////////////////
            // CH_0
            // ////////////////////////////////////////////////////////////////////////

            std::cout << "Interpolation for added points (related to unique convex hull of first level)" << std::endl;
            std::cout << std::endl;

            std::vector<Point3D> chf = ch;
            std::vector<Point3D> new_points_chf;
            for (unsigned int i=ch1_tmp.size(); i < chf.size(); i++)
                new_points_chf.push_back(chf.at(i));

            std::cout << new_points_chf.size() << " points to estimate z value for first level." << std::endl;
            std::cout << std::endl;


            std::vector<Point3D> verts_1, verts_2;
            for(uint vid=0; vid < trimesh0.num_verts(); vid++)
            {
                Point3D p;
                p.x = trimesh0.vert(vid).x();
                p.y = trimesh0.vert(vid).y();
                p.z = trimesh0.vert(vid).z();
                verts_1.push_back(p);
            }

            fittedPlane planef = fitPlane(verts_1);

            for(uint i=0; i<new_points_chf.size(); i++)
                new_points_chf.at(i).z = (new_points_chf.at(i).x-planef.meanX)*planef.meanA0+(new_points_chf.at(i).y-planef.meanY)*planef.meanA1 + planef.meanZ;

            uint ii=0;
            for (uint i = ii+ch1_tmp.size(); i < chf.size(); i++)
            {
                chf.at(i).z = new_points_chf.at(ii).z;
                ii++;
            }

            std::vector<Point3D> points_exf;
            for(uint vid=0; vid<trimesh0.num_verts(); vid++)
            {
                if(!trimesh0.vert_is_boundary(vid))
                {
                    Point3D p;
                    p.x = trimesh0.vert(vid).x();
                    p.y = trimesh0.vert(vid).y();
                    p.z = trimesh0.vert(vid).z();
                    points_exf.push_back(p);
                }
            }
            points_exf.insert(points_exf.end(), chf.begin(), chf.end());


            // ////////////////////////////////////////////////////////////////////////
            // CH_1
            // ////////////////////////////////////////////////////////////////////////

            std::cout << "Interpolation for added points (related to unique convex hull of second level)" << std::endl;
            std::cout << std::endl;

            std::vector<Point3D> chs = ch;
            std::vector<Point3D> new_points_chs;
            for (uint i=0; i < ch1_tmp.size(); i++)
                new_points_chs.push_back(chs.at(i));

            for (uint i=ch1_tmp.size()+ch2_tmp.size(); i < chs.size(); i++)
                new_points_chs.push_back(chs.at(i));

            std::cout << new_points_chs.size() << " points to estimate z value for second level." << std::endl;
            std::cout << std::endl;

            for(uint vid=0; vid < trimesh1.num_verts(); vid++)
            {
                Point3D p;
                p.x = trimesh1.vert(vid).x();
                p.y = trimesh1.vert(vid).y();
                p.z = trimesh1.vert(vid).z();
                verts_2.push_back(p);
            }

            fittedPlane planes = fitPlane(verts_2);
            for(uint i=0; i<new_points_chs.size(); i++)
                new_points_chs.at(i).z = (new_points_chs.at(i).x-planes.meanX)*planes.meanA0+(new_points_chs.at(i).y-planes.meanY)*planes.meanA1 + planes.meanZ;

            for (uint i=0; i < ch1_tmp.size(); i++)
                chs.at(i).z = new_points_chs.at(i).z;

            uint jj = ch1_tmp.size();
            for (uint i = jj +ch2_tmp.size(); i < chs.size(); i++)
            {
                chs.at(i).z = new_points_chs.at(jj).z;
                jj++;
            }


            std::vector<Point3D> points_exs;
            for(uint vid=0; vid<trimesh1.num_verts(); vid++)
            {
                if(!trimesh1.vert_is_boundary(vid))
                {
                    Point3D p;
                    p.x = trimesh1.vert(vid).x();
                    p.y = trimesh1.vert(vid).y();
                    p.z = trimesh1.vert(vid).z();
                    points_exs.push_back(p);
                }
            }
            points_exs.insert(points_exs.end(), chs.begin(), chs.end());


            trimesh0.clear();
            trimesh1.clear();

            trimesh0 = points_triangulation(points_exf, "c");
            trimesh1 = points_triangulation(points_exs, "c");

            //CHIUSURA MESH
            or_idch0.clear();
            or_idch1.clear();

            //chiusura laterale mesh!
            or_idch0 = trimesh0.get_ordered_boundary_vertices();
            or_idch1 = trimesh1.get_ordered_boundary_vertices();

            std::cerr << or_idch0.size() << " " << or_idch1.size() << std::endl;

            trimesh0.update_bbox();
            trimesh1.update_bbox();
        }

        or_idch0.clear();
        or_idch1.clear();


        // Check on normals
        double offset = trimesh1.bbox().center().z() - trimesh0.bbox().center().z();
        //        double offset = trimesh1.bbox().max.z() - trimesh0.bbox().min.z();
        //        std::cout << "bb1: " << trimesh1.bbox().max.z()<< std::endl;
        //        std::cout << "bb0: " << trimesh0.bbox().min.z()<< std::endl;


        std::cout << "Offset in z direction: " << offset << std::endl;

        if(offset > 0) //check for normals and updated (if necessary)
        {
            for(unsigned int pid=0; pid<trimesh0.num_polys(); pid++)
                trimesh0.poly_flip_winding_order(pid);
        }
        else if(offset < 0)
        {
            for(unsigned int pid=0; pid<trimesh1.num_polys(); pid++)
                trimesh1.poly_flip_winding_order(pid);
        }
        else
        {
            std::cerr << "ERROR: z offset cannot be equal to 0." << std::endl;
            exit(1);
        }
        std::cout << FGRN("Check on normals ... COMPLETED.") << std::endl;
        std::cout << std::endl;


        cinolib::Trimesh<> closed_m;
        double step = trimesh1.edge_avg_length();
        if(step < trimesh0.edge_avg_length())
            step = trimesh0.edge_avg_length();
        std::cout << "Step to discretize lateral gap: " << step << std::endl;

        std::cout << FMAG("##########################################") << std::endl;
        std::cout << FMAG("Lo step di discretizzazione viene definito in base all'edge medio (minimo edge medio tra le due mesh)") << std::endl;
        std::cout << FMAG("##########################################") << std::endl;

        trimesh0.edge_mark_boundaries();
        trimesh1.edge_mark_boundaries();

        if(offset < 0)
            closed_m = closing_2trimeshes(trimesh0, trimesh1, step);
        else
            closed_m = closing_2trimeshes(trimesh1, trimesh0, step);

        if(!check_closing_mesh(closed_m))
        {
            std::cout << "\033[0;31mERROR on surfaces closing!\033[0m" << std::endl;
            exit(1);
        }

        std::string out_closed_mesh = out_surf +"/" + basename0 +"-" + basename1 + ext_surf;
        closed_m.save(out_closed_mesh.c_str());
        std::cout << "\033[0;32mExport mesh file: " << out_closed_mesh << " ... COMPLETED.\033[0m" << std::endl;


        MUSE::Surface summary;
        summary.setSummary(closed_m);
        geometa.setMeshSummary(summary);

        geometa.write(out_surf + "/" + basename0 + "-" + basename1 + ".json");
    }

    if(createTriObject.isSet()  && meshFiles.getValue().size() > 2)
        std::cerr << "ERROR: Unexpected number of input files!" << std::endl;


    ///
    /// Creating volume mesh
    ///
    if(createVolObject.isSet() && meshFiles.getValue().size() == 1)
    {
        // 0) Creazione cartella per il salvataggio delle mesh volumetriche
        if(!filesystem::exists(out_volume))
            filesystem::create_directory(out_volume);

        // 1) Passaggio meshfile
        std::vector<std::string> files = meshFiles.getValue();

        MUSE::VolumeMeta geometa;
        geometa.setProject(Project);

        std::vector<std::string> excommands;
        excommands.push_back(command);
        geometa.setCommands(excommands);

        std::vector<std::string> deps;
        deps.push_back(filesystem::relative(get_basename(files.at(0)) + ".json", Project.folder));
        geometa.setDependencies(deps);


        MUSE::Volume summary;

        std::cout << "\033[0;32m=== Loading mesh file: " << files.at(0) << " ... COMPLETED.\033[0m" << std::endl;
        std::string basename = get_basename(get_filename(files.at(0)));
        std::string out_mesh = out_volume +"/" + basename + ext_vol;

        int flagCount = tetFlag.isSet() + voxFlag.isSet() + hexFlag.isSet();
        if (flagCount > 1)
        {
            std::cerr << "=== ERROR: Only one of tetFlag, voxFlag, or hexFlag can be set at a time." << std::endl;
            exit(1);
        }

        // Se è settato il flag per i tetraedri ...
        if(tetFlag.isSet() && xyzPlane.isSet())
        {
            std::cout << "=== tetFlag is set ... " << std::endl;
            std::cout << "=== xyzPlane is set ... " << std::endl;

            bool shift_plane = false;
            double plane_shift = 0.0;
            if (planeShift.isSet())
            {
                plane_shift = planeShift.getValue();
                shift_plane = true;
            }
            cinolib::Trimesh<> trimesh;
            trimesh.load(files.at(0).c_str());

            double avg_length = trimesh.edge_avg_length();
            std::vector<cinolib::vec3d> ppoints;
            double plane_min_z = DBL_MAX;

            std::ifstream pf;
            pf.open(xyzPlane.getValue());

            double x,y,z;
            while (pf >> x >> y >> z)
            {
                ppoints.push_back(cinolib::vec3d(x,y,z));
                if (z < plane_min_z) plane_min_z = z;
            }
            std::cout << "=== Reading plane from file: " << xyzPlane.getValue() << " COMPLETED." << std::endl;

            cinolib::Plane plane (ppoints);
            const std::vector<uint> bverts = trimesh.get_ordered_boundary_vertices();
            std::vector<cinolib::vec3d> bverts_proj;

            for (uint vid : bverts)
            {
                cinolib::vec3d intersection;
                cinolib::vec3d p2 = trimesh.vert(vid);
                p2.z() = plane_min_z;
                cinolib::Segment segm (0, trimesh.vert(vid), p2);

                if (!intersectPlaneSegment(plane,segm,intersection)) //.vert(vid));
                {
                    std::cerr << "=== ERROR: No intersection is found!" << std::endl;
                }

                bverts_proj.push_back(intersection);
            }

            std::cout << "=== Projected" << std::endl;
            std::vector<uint> segms;

            for (uint i=0; i < bverts_proj.size(); i++)
            {
                segms.push_back(i);
                segms.push_back(i+1);
            }

            segms.push_back(bverts_proj.size()-1);
            segms.push_back(0);

            std::string triangle_flags = "";
            if (optAdditionalFlag.isSet())
                triangle_flags = optAdditionalFlag.getValue();

            if (triangle_flags.length() == 0)
            {
                double avg_area = 0.0;

                for (uint pid=0; pid < trimesh.num_polys(); pid++)
                    avg_area += trimesh.poly_area(pid);

                avg_area /= trimesh.num_polys();

                avg_area *=2.0;

                triangle_flags = "pYa" + std::to_string(avg_area);
            }
            std::cout << "=== Set triangle flags as ... " << triangle_flags << std::endl;

            cinolib::Trimesh<> mbot;
            triangle_wrap(bverts_proj, segms, std::vector<cinolib::vec3d> (), plane_min_z, triangle_flags.c_str() , mbot);

            std::cout << "=======" << std::endl;

            for (uint vid=0 ; vid < mbot.num_verts(); vid++)
            {
                // vec3d pv = plane.project_onto(mbot.vert(vid));

                cinolib::vec3d intersection;
                cinolib::vec3d p2 = mbot.vert(vid);
                p2.z() = trimesh.bbox().max.z();
                cinolib::Segment segm (0, mbot.vert(vid), p2);
                intersectPlaneSegment(plane,segm,intersection); //.vert(vid));

                mbot.vert(vid).z() = intersection.z();
            }

            std::cout << "projected" << std::endl;

            std::cout << trimesh.get_ordered_boundary_vertices().size() << std::endl;
            std::cout << mbot.get_ordered_boundary_vertices().size() << std::endl;

            for (uint pid=0; pid < mbot.num_polys(); pid++)
            {
                const std::vector<uint> vp = mbot.vector_polys().at(pid);
                mbot.vector_polys().at(pid).at(0) = vp.at(2);
                mbot.vector_polys().at(pid).at(2) = vp.at(0);
            }

            mbot.update_p_normals();

            cinolib::Trimesh<> total = trimesh;
            total += mbot;

            for (uint pid=0; pid < total.num_polys(); pid++)
                total.poly_data(pid).label=INT_MAX;

            uint n_bverts = bverts.size();
            uint n_edges = total.num_edges();

            for (uint i=0; i < n_bverts-1; i++)
            {
                total.poly_add(bverts.at(i), trimesh.num_verts() +i, bverts.at(i+1));
                total.poly_add(bverts.at(i+1), trimesh.num_verts() +i, trimesh.num_verts() +i+1);
                // std::cout << total.num_polys() << std::endl;
            }

            total.poly_add(bverts.at(bverts.size()-1), trimesh.num_verts() +bverts.size()-1, bverts.at(0));
            total.poly_add(bverts.at(0), trimesh.num_verts() +bverts.size()-1, trimesh.num_verts() +0);

            bool split = false;
            uint n_splits = 0;

            do
            {
                split = false;
                for (int eid=total.num_edges()-1; eid >=n_edges; eid--)
                    if (total.edge_length(eid) > avg_length)
                    {
                        if (total.poly_data(total.adj_e2p(eid).at(0)).label < INT_MAX && total.poly_data(total.adj_e2p(eid).at(1)).label < INT_MAX)
                        {
                            total.edge_split(eid);
                            split=true;
                        }
                    }
                n_splits++;
            } while (/*split == true*/n_splits < 2);

            cinolib::Tetmesh<> tetmesh, tmp;

            cinolib::vec3d translate_vec = total.bbox().center();
            total.translate(-translate_vec);

            std::string tetgen_flags = "";
            if (optFlag.isSet())
                tetgen_flags = optFlag.getValue();

            if (tetgen_flags.length() == 0)
            {
                tetgen_wrap(total, "Y", tmp);

                double avg_vol=0;

                for (uint pid=0; pid < tmp.num_polys(); pid++)
                {
                    // total.vert_add(tmp.poly_centroid(pid));
                    avg_vol += tmp.poly_volume(pid);
                }

                avg_vol /= tmp.num_polys();

                double vol = avg_vol / 1.1;

                tetgen_flags = "Ya" + std::to_string(vol);
            }
            std::cout << "=== Set tetgen flags as ... " << tetgen_flags << std::endl;
            tetgen_wrap(total, tetgen_flags, tetmesh);
            // tetgen_wrap(total, "Y" , tetmesh);

            tetmesh.translate(translate_vec);

            for (cinolib::vec3d &p : ppoints)
                p.z() += plane_shift;

            cinolib::Plane shifted_plane (ppoints);

            std::vector<uint> tbsplit;
            std::cout << "=== Check on z value: consider only the tetrahedron cutted by the plane by comparing z coordinates." << std::endl;
            for (int pid = 0; pid <= tetmesh.num_polys()-1; pid++)
            {
                double min_z = DBL_MAX;
                double max_z = -DBL_MAX;

                // Calcola la z min/max dei vertici del tet
                for (uint v = 0; v < 4; ++v)
                {
                    double z = tetmesh.poly_vert(pid, v).z();
                    if (z < min_z) min_z = z;
                    if (z > max_z) max_z = z;
                }
                // Se il piano è completamente sopra o sotto il tet, salta
                if (shifted_plane.p.z() < min_z || shifted_plane.p.z() > max_z)
                {
                    //std::cout << "=== The plane is out. Skip tet-pid = " << pid << std::endl;
                    continue;
                }

                // Ora che sappiamo che lo attraversa, procedi con il test dettagliato
                uint v_up = 0, v_down = 0;

                for (uint v=0; v < 4; v++)
                {
                    // Proiezione del segmento lungo z! Generalizzabile ...
                    cinolib::vec3d c = tetmesh.poly_vert(pid, v);
                    cinolib::vec3d s0 = c; s0.z() = tetmesh.bbox().min.z();
                    cinolib::vec3d s1 = c; s1.z() = tetmesh.bbox().max.z();

                    cinolib::Segment segm (0, s0, s1);

                    cinolib::vec3d intersection;
                    intersectPlaneSegment(shifted_plane, segm, intersection);

                    if (intersection.z() > c.z())
                        v_up++;
                    else
                        v_down++;
                }

                if ( v_up > 0 && v_down > 0 )
                    tbsplit.push_back(pid);
            }

            std::cout << "Splitting " << tbsplit.size() << " tets... " << std::endl;

            tetmesh.polys_split(tbsplit);

            std::cout << tetmesh.num_verts() << "V / " << tetmesh.num_polys() << "P" << std::endl;

            tetmesh.save(out_mesh.c_str());

        }

        if(tetFlag.isSet() && !xyzPlane.isSet())
        {
            std::cout << "=== tetFlag is set ... " << std::endl;

            //se voglio i tet ...
            cinolib::Trimesh<> trimesh;
            trimesh.load(files.at(0).c_str());


            // Ratio rappresenta un indice di anisotropia geometrica della mesh:
            // confronta l’estensione massima nel piano orizzontale (XY) rispetto a quella verticale (Z)
            double delta_max = trimesh.bbox().delta_x();
            if(trimesh.bbox().delta_y() >= delta_max)
                delta_max = trimesh.bbox().delta_y();

            std::cout << "max{delta_x,delta_y} = " << delta_max << std::endl;
            std::cout << "delta_z = " << trimesh.bbox().delta_z() << std::endl;

            if(trimesh.bbox().delta_z() == 0.0)
            {
                std::cerr << "=== ERROR: delta_z is zero, cannot compute ratio max{delta_x,delta_y}/delta_z." << std::endl;
                exit(1);
            }
            double ratio = delta_max/trimesh.bbox().delta_z();
            std::cout << "=== Computing anisotropy ratio between max{delta_x,delta_y}/delta_z: " << ratio << std::endl;


            //AGGIUNGERE LA CONDIZIONE PER LA TRASLAZIONE
            cinolib::vec3d center = trimesh.bbox().center();
            std::cout << "=== Translate mesh at BBOX center: " << center << std::endl;
            std::cout << std::endl;
            trimesh.translate(-center);
            if(setSave.isSet())
                trimesh.save((get_basename(files.at(0)) + "_translate"+ext_surf).c_str());

            // Set parameters in opt
            std::string opt = "";
            if(optFlag.isSet())
                opt = opt + optFlag.getValue();

            // Run tetrahedralization by exploting Tetgen Library in Cinolib and create a tetrahedralization mesh (m_tet)
            cinolib::Tetmesh<> volmesh;
            cinolib::tetgen_wrap(trimesh.vector_verts(), trimesh.vector_polys(), trimesh.vector_edges(), opt, volmesh);

            volmesh.translate(center);
            std::cout << "### Restore coordinates mesh from BBOX center: " << center << " COMPLETED." <<std::endl;

            std::cout << std::endl;
            std::cout << "#############################################" << std::endl;
            std::cout << "### Statistics on volume ... " << std::endl;
            std::cout << "### Poly average volume: " << volmesh.mesh_volume()/volmesh.num_polys() << std::endl;
            std::cout << "### Edge average length: " << volmesh.edge_avg_length() << std::endl;
            std::cout << "### Edge max length: " << volmesh.edge_max_length() << std::endl;
            std::cout << "### Edge min length: " << volmesh.edge_min_length() << std::endl;

            if(volmesh.edge_max_length() > volmesh.edge_avg_length() * 1.5)
                std::cout << FYEL("### WARNING: (max) edge length major than 1.5 times edge average length ...") << std::endl;

            std::cout << "#############################################" << std::endl;
            std::cout << std::endl;

            volmesh.save(out_mesh.c_str());
            std::cout << "\033[0;32mExport mesh file: " << out_mesh << " ... COMPLETED.\033[0m" << std::endl;

            summary.setSummary(volmesh);
        }


        if(voxFlag.isSet())
        {
            std::cout << "voxFlag is set ... " << std::endl;

            //MUSE::Quadmesh<> quadmesh;
            cinolib::Polygonmesh<> quadmesh;
            quadmesh.load(files.at(0).c_str());

            uint max_voxels_per_side = setMaxVoxelperSide.getValue();
            cinolib::VoxelGrid grid;
            cinolib::voxelize(quadmesh, max_voxels_per_side, grid);

            std::cout << "Grid dimensions: " << grid.dim[0] << " x " << grid.dim[1] << " x " << grid.dim[2] << std::endl;

            cinolib::Hexmesh<> volmesh;
            voxel_grid_to_hexmesh(grid, volmesh, cinolib::VOXEL_INSIDE);

            volmesh.save(out_mesh.c_str());
            std::cout << "\033[0;32mExport mesh file: " << out_mesh << " ... COMPLETED.\033[0m" << std::endl;

            summary.setSummary(volmesh);
        }


        if(hexFlag.isSet())
        {
            std::cout << "hexFlag is set ... " << std::endl;

            //cinolib::Quadmesh<> quadmesh;
            cinolib::Trimesh<> mesh;
            mesh.load(files.at(0).c_str());

            MUSE::Hexmesh<> hexmesh(setResx.getValue(), setResy.getValue(), setResz.getValue(), mesh);

            hexmesh.save(out_mesh.c_str());
            std::cout << "\033[0;32mExport mesh file: " << out_mesh << " ... COMPLETED.\033[0m" << std::endl;

            summary.setSummary(hexmesh);
        }

        geometa.setMeshSummary(summary);
        geometa.write(out_volume +"/" + basename + ".json");
    }

    if(createVolObject.isSet()  && meshFiles.getValue().size() > 1)
        std::cerr << "ERROR: Unexpected number of input files!" << std::endl;



    ///
    /// Loading and editing surface mesh
    ///
    /// Questo comando permette la lettura di superfici di origine esterna e la creazione del file json corrispondente
    /// Questo comando ha anche la possibilità di effettuare infittimento mediante split su centroide o punto medio edge
    /// Può prendere in input più superfici (-m <filename> -m filename ...)
    ///
    if(loadSurface.isSet() && !setRemeshing.isSet())
    {
        if(!meshFiles.isSet())
        {
            std::cout << FRED("ERROR. Set a surface mesh by -m command") << std::endl;
            exit(1);
        }

        if(meshFiles.getValue().size() >= 2)
        {
            std::cout << FRED("ERROR. Only a surface mesh file is supported.") << std::endl;
            exit(1);
        }

        std::string filename_mesh = meshFiles.getValue().at(0);

        MUSE::SurfaceMeta geometa;
        geometa.setProject(Project);

        std::vector<std::string> excommands;
        excommands.push_back(command);
        geometa.setCommands(excommands);

        if(splitMethod.isSet())
        {
            std::vector<std::string> deps;
            deps.push_back(filesystem::relative(get_basename(filename_mesh) + ".json", Project.folder));
            geometa.setDependencies(deps);
        }

        MUSE::SurfaceMesh<>mesh;
        mesh.load(filename_mesh.c_str());
        std::cout << "\033[0;32m=== Loading mesh file: " << filename_mesh << " ... COMPLETED.\033[0m" << std::endl;

        MeshType type = mesh.set_meshtype();
        std::cout << "=== Check mesh type ... " << std::endl;
        std::cout << "=== Number of verts per poly: " << mesh.verts_per_poly(0) <<  std::endl;

        MUSE::Surface surf;
        MUSE::Surface::Parameters surf_par;

        if(setRotAxis.isSet())
        {
            double rad = (setRotAngle.getValue() * M_PI)/180;
            cinolib::vec3d axis = set_rotation_axis(setRotAxis.getValue());

            cinolib::mat3d R = cinolib::mat3d::ROT_3D(axis, rad);
            cinolib::vec3d rotcenter {setRotCenterX.getValue(), setRotCenterY.getValue(), setRotCenterZ.getValue()};

            for(uint vid=0; vid<mesh.num_verts(); vid++)
            {
                mesh.vert(vid) -= rotcenter;
                mesh.vert(vid) = R*mesh.vert(vid);
                mesh.vert(vid) += rotcenter;
            }
        }

        if(type == MeshType::TRIMESH)
        {
            std::cout << "=== Mesh type: TRIMESH" <<  std::endl;
            std::cout << std::endl;

            surf_par.type = "TRIMESH";

            if(splitMethod.isSet())
            {
                if(splitMethod.getValue().compare("CENTROID") == 0)
                {
                    std::cout << "=== Management of new degree of resolution by poly split at centroid ..." << std::endl;
                    mesh.triangles_split_on_centroid();
                }
                else if(splitMethod.getValue().compare("EDGE") == 0)
                {
                    std::cout << "=== Management of new degree of resolution by poly split at edges middle point ..." << std::endl;
                    mesh.triangles_split_on_edge();
                }
                else if(splitMethod.getValue().compare("BEDGE") == 0)
                {
                    std::cout << "=== Poly split at boundary edges middle point ..." << std::endl;
                    double inedge_avg = 0.0;
                    int n_inedge = 0;
                    for(uint eid=0; eid < mesh.num_edges(); eid++)
                    {
                        if(!mesh.edge_is_boundary(eid))
                        {
                            n_inedge++;
                            inedge_avg += mesh.edge_length(eid);
                        }
                    }
                    inedge_avg = inedge_avg/n_inedge;
                    std::cout << "=== (Internal) edge average lenght: " << inedge_avg << std::endl;
                    std::cout << "=== (Internal) edge average lenght * factor: " << inedge_avg + (inedge_avg * setFactor.getValue()) << std::endl;
                    std::cout << "=== Edge average lenght: " << mesh.edge_avg_length() << std::endl;

                    //DA CONTROLLARE
                    for(uint eid=0; eid < mesh.num_edges(); eid++)
                    {
                        if(!mesh.edge_is_boundary(eid)) continue;

                        if(mesh.edge_length(eid) > inedge_avg + (inedge_avg * setFactor.getValue()))
                        {
                            // Un edge di bordo dovrebbe avere solo un poligono adiacente
                            if(mesh.adj_e2p(eid).size() != 1) continue;

                            std::cout << "=== Edge ID: " << eid << " - Edge lenght: " << mesh.edge_length(eid) << std::endl;

                            //std::vector<uint> pid_eid = mesh.adj_e2p(eid);

                            //std::cout << "=== pid adj edge: " << pid_eid.size() << std::endl;

                            cinolib::vec3d v0 = mesh.edge_vert(eid, 0);
                            cinolib::vec3d v1 = mesh.edge_vert(eid, 1);
                            cinolib::vec3d delta = (v1-v0)/2;

                            cinolib::vec3d v_med (v0.x()+delta.x(), v0.y()+delta.y(), v0.z()+delta.z());
                            mesh.vert_add(v_med);
                        }
                    }
                }
                else
                {
                    std::cout << FRED("ERROR. Split method: ") << splitMethod.getValue()  << FRED(" is not supported.") << std::endl;
                    exit(1);
                }
            }

            if(boundaryExtract.isSet())
            {
                std::vector<Point3D> vec_bv;
                for(uint i:mesh.get_ordered_boundary_vertices())
                {
                    Point3D bv;
                    bv.x = mesh.vert(i).x();
                    bv.y = mesh.vert(i).y();
                    bv.z = mesh.vert(i).z();
                    vec_bv.push_back(bv);
                }
                export3d_xyz(out_surf + "/"+ get_basename(get_filename(filename_mesh)) + "_BP" + ext_txt, vec_bv);
            }

            if(setPerturbation.isSet())
            {
                std::srand(static_cast<unsigned int>(std::time(nullptr)));

                for(uint vid=0; vid < mesh.num_verts(); vid++)
                {
                    float perturbation = static_cast<float>(1000 + (std::rand() % 1001));
                    std::cout << "pertubazione z vertice: " << perturbation << std::endl;
                    mesh.vert(vid).z() += perturbation;
                }
            }
        }
        else if(type == MeshType::QUADMESH)
        {
            std::cout << "Mesh type: QUADMESH" <<  std::endl;
            std::cout << std::endl;

            surf_par.type = "QUADMESH";

            if(splitMethod.isSet())
            {
                std::cout << "### For quads mesh: poly split on edge/centroid corresponds." << std::endl;
                if(splitMethod.getValue().compare("EDGE") == 0 || splitMethod.getValue().compare("CENTROID") == 0)
                {
                    std::cout << "Management of new degree of resolution by poly split ..." << std::endl;
                    mesh.quads_split_on_edge();
                }
                else
                {
                    std::cout << FRED("ERROR. Split method: ") << splitMethod.getValue()  << FRED(" is not supported.") << std::endl;
                    exit(1);
                }
            }
        }
        else
        {
            std::cout << FRED("ERROR. Only triangle/quadrilateral meshes are supported!") << std::endl;
            exit(1);
        }

        surf.setParameters(surf_par);
        surf.setSummary(mesh);
        geometa.setMeshSummary(surf);

        std::string suffix;
        if(splitMethod.isSet()) suffix += "_res";
        if(setRotAxis.isSet()) suffix += "_rot";
        if(setPerturbation) suffix += "_perturb";
        if(setScaleMesh.isSet())
        {
            mesh.scale(setScaleFactorX, setScaleFactorY, setScaleFactorZ);
            suffix += "_scale";
        }

        //std::string out_mesh = get_filename(filename_mesh);
        std::string final_out = get_basename(filename_mesh) + suffix + ext_surf;
        mesh.save(final_out.c_str());
        std::cout << "\033[0;32mSaving mesh file: " << final_out << "\033[0m" << std::endl;


        std::cout << "=================== INFO bbox vertices ===================" << std::endl;
        std::cout << "=== bbox x_min, y_min: " << std::fixed << std::setprecision(setPrecision.getValue()) << mesh.bbox().min.x() << "; " << mesh.bbox().min.y() <<  std::endl;
        std::cout << "=== bbox x_max, y_min: " << std::fixed << std::setprecision(setPrecision.getValue()) << mesh.bbox().max.x() << "; " << mesh.bbox().min.y() <<  std::endl;
        std::cout << "=== bbox x_max, y_max: " << std::fixed << std::setprecision(setPrecision.getValue()) << mesh.bbox().max.x() << "; " << mesh.bbox().max.y() <<  std::endl;
        std::cout << "=== bbox x_min, y_max: " << std::fixed << std::setprecision(setPrecision.getValue()) << mesh.bbox().min.x() << "; " << mesh.bbox().max.y() <<  std::endl;
        std::cout << "=== bbox diag: " << std::fixed << std::setprecision(setPrecision.getValue()) << mesh.bbox().diag() << std::endl;

        geometa.write(get_filename(filename_mesh) + suffix + ".json");
    }



    ///
    /// Merging meshes
    ///
    if(mergeMeshes.isSet() && meshFiles.getValue().size() == 2)
    {
        std::vector<std::string> files = meshFiles.getValue();

        std::string filename_mesh0 = files.at(0);
        std::string ext0 = get_extension(filename_mesh0);

        std::string filename_mesh1 = files.at(1);
        std::string ext1 = get_extension(filename_mesh1);


        if(ext0.compare(".off") == 0 || ext0.compare(".obj") == 0)
        {
            if(ext1.compare(".off") == 0 || ext1.compare(".obj") == 0)
            {
                std::cout << "Meshes are surfaces." << std::endl;
                //Le mesh sono superfici (controllo sull'estensione), quindi le carico come trimesh

                MUSE::SurfaceMesh<> trimesh0;
                trimesh0.load(filename_mesh0.c_str());
                std::cout << "\033[0;32mLoading mesh file: " << filename_mesh0 << " ... COMPLETED.\033[0m" << std::endl;
                std::string basename0 = get_basename(get_filename(filename_mesh0));

                MUSE::SurfaceMesh<> trimesh1;
                trimesh1.load(filename_mesh1.c_str());
                std::cout << "\033[0;32mLoading mesh file: " << filename_mesh1 << " ... COMPLETED.\033[0m" << std::endl;
                std::string basename1 = get_basename(get_filename(filename_mesh1));

                MUSE::SurfaceMesh<> trimesh;
                std::string out_mesh = out_surf +"/" + basename0 + "_" + basename1 + ext_surf;

                if(!trimesh0.check_lateral_closing() && !trimesh1.check_lateral_closing())
                {
                    merge_and_wrap_meshes_old(trimesh0, trimesh1, trimesh);
                    //merge_meshes(trimesh0, trimesh1, trimesh);
                    std::cout << "Meshes merge on boundary ... COMPLETED." << std::endl;
                }
                else if(trimesh0.check_lateral_closing() && trimesh1.check_lateral_closing())
                {
                    cinolib::merge_meshes_at_coincident_vertices(trimesh0, trimesh1, trimesh);
                    std::cout << "Meshes merge at coincident vertices ... COMPLETED." << std::endl;
                }
                else
                {
                    std::cout << "Error on meshes type!" << std::endl;
                    exit(1);
                }
                trimesh.save(out_mesh.c_str());

                std::cout << std::endl;
                std::cout << "P " << trimesh.num_polys() << std::endl;
                std::cout << "E "<< trimesh.num_edges() << std::endl;
                std::cout << "V "<< trimesh.num_verts() << std::endl;
                std::cout << "\033[0;32mSaving mesh file in : " << out_mesh << " ... COMPLETED.\033[0m" << std::endl;
            }
            else
            {
                std::cerr << "ERROR: Meshes format are different!" << std::endl;
                exit(1);
            }
        }
        else if(ext0.compare(".mesh") == 0 || ext0.compare(".vtk") == 0) //caso volumetrico
        {
            if(ext1.compare(".mesh") == 0 || ext0.compare(".vtk") == 0)
            {
                std::cout << "Meshes are volumes." << std::endl;

                MUSE::VolumeMesh<> tetmesh0;
                //cinolib::Hexmesh<> tetmesh0;
                tetmesh0.load(filename_mesh0.c_str());
                MeshType type0 = tetmesh0.set_meshtype();
                std::cout << "\033[0;32mLoading mesh file: " << filename_mesh0 << " ... COMPLETED.\033[0m" << std::endl;
                std::string basename0 = get_basename(get_filename(filename_mesh0));

                MUSE::VolumeMesh<> tetmesh1;
                tetmesh1.load(filename_mesh1.c_str());
                MeshType type1 = tetmesh1.set_meshtype();
                std::cout << "\033[0;32mLoading mesh file: " << filename_mesh1 << " ... COMPLETED.\033[0m" << std::endl;
                std::string basename1 = get_basename(get_filename(filename_mesh1));

                if(type0 != type1)
                {
                    std::cout << FRED("Mesh types are different. Merge not possible!") << std::endl;
                    exit(1);
                }

                std::string out_mesh = out_volume +"/" + basename0 + "_" + basename1 + ext_vol;

                if(type0 == MeshType::HEXMESH)
                {
                    cinolib::Hexmesh<> tetmesh;
                    //cinolib::Tetmesh<> tetmesh;
                    std::cout << "The proximity threshold is set on " << proxThreshold.getValue() << std::endl;
                    merge_meshes_at_coincident_vertices(tetmesh0, tetmesh1, tetmesh, proxThreshold.getValue());
                    std::cout << "Meshes merge at coincident vertices ... COMPLETED." << std::endl;

                    std::cout << std::endl;
                    std::cout << "P " << tetmesh.num_polys() << std::endl;
                    std::cout << "F " << tetmesh.num_faces() << std::endl;
                    std::cout << "E "<< tetmesh.num_edges() << std::endl;
                    std::cout << "V "<< tetmesh.num_verts() << std::endl;

                    tetmesh.save(out_mesh.c_str());
                }
                else
                {
                    //cinolib::Hexmesh<> tetmesh;
                    cinolib::Tetmesh<> tetmesh;
                    std::cout << "The proximity threshold is set on " << proxThreshold.getValue() << std::endl;
                    merge_meshes_at_coincident_vertices(tetmesh0, tetmesh1, tetmesh, proxThreshold.getValue());
                    std::cout << "Meshes merge at coincident vertices ... COMPLETED." << std::endl;

                    std::cout << std::endl;
                    std::cout << "P " << tetmesh.num_polys() << std::endl;
                    std::cout << "F " << tetmesh.num_faces() << std::endl;
                    std::cout << "E "<< tetmesh.num_edges() << std::endl;
                    std::cout << "V "<< tetmesh.num_verts() << std::endl;

                    tetmesh.save(out_mesh.c_str());
                }
                std::cout << "\033[0;32mSaving mesh file in : " << out_mesh << " ... COMPLETED.\033[0m" << std::endl;
            }
            else
            {
                std::cerr << "ERROR: Meshes format are different!" << std::endl;
                exit(1);
            }
        }
        else
        {
            std::cerr << "ERROR: Mesh format is not supported." << std::endl;
            exit(1);
        }
    }

    if(gridData.isSet())
    {
        std::vector<std::string> excommands;
        excommands.push_back(command);

        std::vector<Point3D> boundary;
        if(setBBPoints.isSet())
        {
            std::vector<std::string> bbpoints = setBBPoints.getValue();

            for(uint i=0; i< bbpoints.size(); i++)
                std::cout << "BBPOINTS: " << bbpoints.at(i) << std::endl;
            std::cout << std::endl;

            for(uint i=0; i< bbpoints.size(); i++)
            {
                std::vector<std::string> direc = split_string(bbpoints.at(i), ',');

                Point3D p0;
                p0.x = std::stod(direc.at(0));
                p0.y = std::stod(direc.at(1));
                p0.z = std::stod(direc.at(2));
                boundary.push_back(p0);
            }
        }
        else
        {
            std::cout << FRED("ERROR. Set --bbp boundary points to extract mesh from grid.") << std::endl;
            exit(1);
        }

        //TO DOOOOOOOOOOOOOOOOOO: CONTROLLARE FUNZIONE DI ROTAZIONE PUNTI!!!
        if(setRotAxis.isSet())
        {
            std::vector<Point3D> boundary_tmp = boundary;
            boundary.clear();
            for(size_t i=0; i<boundary_tmp.size(); i++)
            {
                Point3D p_rot = rotPoint(boundary_tmp.at(i), setRotAxis.getValue(), setRotAngle.getValue());
                boundary.push_back(p_rot);
            }
        }

        if(gridFlag.isSet())
        {
            if(!filesystem::exists(out_surf))
                filesystem::create_directory(out_surf);

            MUSE::SurfaceMeta geometa;
            geometa.setProject(Project);
            geometa.setCommands(excommands);

            MUSE::Quadmesh<> quadmesh (setResx.getValue(), setResy.getValue(), setNewZ.getValue(), boundary);

            std::string out_mesh = out_surf + "/grid" + ext_surf;

            quadmesh.save(out_mesh.c_str());

            MUSE::Surface summary;
            MUSE::Surface::Parameters par;
            par.resx = setResx.getValue();
            par.resy = setResy.getValue();
            summary.setParameters(par);
            summary.setSummary(quadmesh);
            geometa.setMeshSummary(summary);

            geometa.write(out_surf + "/grid" + ".json");
        }

        if(hexFlag.isSet())
        {
            if(!filesystem::exists(out_volume))
                filesystem::create_directory(out_volume);

            MUSE::VolumeMeta geometa;
            geometa.setProject(Project);
            geometa.setCommands(excommands);

            MUSE::Hexmesh<> mesh (setResx.getValue(), setResy.getValue(), setResz.getValue(), boundary);

            std::string out_mesh = out_volume + "/grid" + ext_vol;

            mesh.save(out_mesh.c_str());

            MUSE::Volume summary;
            MUSE::Volume::Parameters par;
            par.resx = setResx.getValue();
            par.resy = setResy.getValue();
            par.resz = setResz.getValue();
            summary.setParameters(par);
            summary.setSummary(mesh);
            geometa.setMeshSummary(summary);

            geometa.write(out_volume + "/grid" + ".json");
        }

        std::cout << "Grid creation ... COMPLETED." << std::endl;

    }

    if(extractMeshes.isSet() && meshFiles.getValue().size() == 1)
    {
        std::vector<std::string> excommands;
        excommands.push_back(command);

        std::vector<std::string> deps;

        std::vector<std::string> files = meshFiles.getValue();

        if(gridFlag.isSet())
        {
            MUSE::SurfaceMeta geometa;
            geometa.setProject(Project);
            geometa.setCommands(excommands);

            std::cout << "### Load grid: " << files.at(0) << std::endl;
            MUSE::Quadmesh<> quadmesh;
            quadmesh.load(files.at(0).c_str());
            deps.push_back(filesystem::relative(get_basename(files.at(0)) + ".json", Project.folder));

            std::string name = get_basename(files.at(0));

            std::cout << "### Load mesh for extracting boundary: " << setBoundary.getValue() << std::endl;
            cinolib::Trimesh<> mesh_bound;
            mesh_bound.load(setBoundary.getValue().c_str());
            deps.push_back(filesystem::relative(get_basename(setBoundary.getValue()) + ".json", Project.folder));

            std::string bound_name = setBoundary.getValue().substr(setBoundary.getValue().find_last_of("/")+1, setBoundary.getValue().length());
            bound_name = get_basename(bound_name);

            //split mesh
            std::vector<Point2D> bound_2d;
            for(uint i: mesh_bound.get_ordered_boundary_vertices())
            {
                Point2D p;
                p.x = mesh_bound.vert(i).x();
                p.y = mesh_bound.vert(i).y();
                bound_2d.push_back(p);
            }

            MUSE::Quadmesh<> sub_quadmesh;
            std::map<cinolib::vec3d, uint> verts;

            for(uint pid=0; pid <quadmesh.num_polys(); pid++)
            {
                cinolib::vec3d centr = quadmesh.poly_centroid(pid);

                Point2D c;
                c.x = centr.x();
                c.y = centr.y();

                if(point_in_polygon(c, bound_2d))
                {
                    cinolib::vec3d v0_pos = quadmesh.poly_vert(pid, 0);
                    cinolib::vec3d v1_pos = quadmesh.poly_vert(pid, 1);
                    cinolib::vec3d v2_pos = quadmesh.poly_vert(pid, 2);
                    cinolib::vec3d v3_pos = quadmesh.poly_vert(pid, 3);

                    // Definizione dell'iteratore
                    auto v0_it = verts.find(v0_pos);
                    auto v1_it = verts.find(v1_pos);
                    auto v2_it = verts.find(v2_pos);
                    auto v3_it = verts.find(v3_pos);

                    // Definizione indici vertici
                    uint v0_id = 0;
                    uint v1_id = 0;
                    uint v2_id = 0;
                    uint v3_id = 0;

                    if (v0_it == verts.end()) //se non lo trovo, quindi il vertice non è stato ancora aggiunto
                    {
                        v0_id = sub_quadmesh.vert_add(v0_pos);
                        verts.insert(std::pair<cinolib::vec3d,uint> (v0_pos, v0_id));
                    }
                    else
                        v0_id = v0_it->second;

                    if (v1_it == verts.end())
                    {
                        v1_id = sub_quadmesh.vert_add(v1_pos);
                        verts.insert(std::pair<cinolib::vec3d,uint> (v1_pos, v1_id));
                    }
                    else
                        v1_id = v1_it->second;

                    if (v2_it == verts.end())
                    {
                        v2_id = sub_quadmesh.vert_add(v2_pos);
                        verts.insert(std::pair<cinolib::vec3d,uint> (v2_pos, v2_id));
                    }
                    else
                        v2_id = v2_it->second;

                    if (v3_it == verts.end())
                    {
                        v3_id = sub_quadmesh.vert_add(v3_pos);
                        verts.insert(std::pair<cinolib::vec3d,uint> (v3_pos, v3_id));
                    }
                    else
                        v3_id = v3_it->second;

                    std::vector<uint> vlist;
                    vlist.push_back(v0_id);
                    vlist.push_back(v1_id);
                    vlist.push_back(v2_id);
                    vlist.push_back(v3_id);

                    sub_quadmesh.poly_add(vlist);
                }
            }
            //std::cout << cinolib::connected_components(sub_quadmesh) << std::endl;
            if(cleanPoly.isSet())
                sub_quadmesh.remove_isolate_poly();

            std::string out_mesh = name + "_" + bound_name + ext_surf;
            sub_quadmesh.save(out_mesh.c_str());

            std::cout << std::endl;
            std::cout << "Saving quadmesh: " << out_mesh << std::endl;

            MUSE::Surface summary;
            MUSE::Surface::Parameters par;
            par.resx = setResx.getValue();
            par.resy = setResy.getValue();
            summary.setParameters(par);
            summary.setSummary(quadmesh);
            geometa.setMeshSummary(summary);
            geometa.setDependencies(deps);

            geometa.write(name + "_" + bound_name + ".json");
        }

        if(hexFlag.isSet())
        {
            MUSE::VolumeMeta geometa;
            geometa.setProject(Project);
            geometa.setCommands(excommands);

            std::cout << "### Load grid: " << files.at(0) << std::endl;
            MUSE::Hexmesh<> hexmesh;
            hexmesh.load(files.at(0).c_str());
            deps.push_back(filesystem::relative(get_basename(files.at(0)) + ".json", Project.folder));

            std::string name = files.at(0).substr(files.at(0).find_last_of("/")+1, files.at(0).length());
            name = get_basename(name);

            //std::string ext0 = get_extension(files.at(0));

            std::cout << "### Load mesh for extracting boundary: " << setBoundary.getValue() << std::endl;
            cinolib::Trimesh<> mesh_bound;
            mesh_bound.load(setBoundary.getValue().c_str());
            deps.push_back(filesystem::relative(get_basename(setBoundary.getValue()) + ".json", Project.folder));

            std::string bound_name = setBoundary.getValue().substr(setBoundary.getValue().find_last_of("/")+1, setBoundary.getValue().length());
            bound_name = get_basename(bound_name);

            MUSE::Hexmesh<> sub_hexmesh;
            sub_hexmesh.subHexmesh_from_trimesh(hexmesh, mesh_bound);

            MUSE::Volume summary;
            MUSE::Volume::Parameters par;
            par.resx = setResx.getValue();
            par.resy = setResy.getValue();
            par.resz = setResz.getValue();
            summary.setParameters(par);
            summary.setSummary(sub_hexmesh);
            geometa.setMeshSummary(summary);
            geometa.setDependencies(deps);

            geometa.write(out_volume + "/" + name + "_" + bound_name + ".json");

            std::string out_mesh = out_volume +"/" + name + "_" + bound_name + ext_vol;
            sub_hexmesh.save(out_mesh.c_str());

            std::cout << std::endl;
            std::cout << "Saving hexmesh: " << out_mesh << std::endl;
        }
        std::cout << FGRN("Extracting sub-mesh constrained to boundary ... COMPLETED.") << std::endl;
    }




    if(loadSurface.isSet() && setRemeshing.isSet())
    {
        if(!meshFiles.isSet())
        {
            std::cout << FRED("ERROR. Set a mesh (surface/volume) by -m command") << std::endl;
            exit(1);
        }

        if(meshFiles.getValue().size() >= 2)
        {
            std::cout << FRED("ERROR. Only a mesh (surface/volume) is supported.") << std::endl;
            exit(1);
        }

        std::cout << "############################" << std::endl;
        std::cout << "### REMESHING ALGORITHM" << std::endl;
        std::cout << "### Reference: M.Botsch, L.Kobbelt, A Remeshing Approach to Multiresolution Modeling." << std::endl;
        std::cout << "### Remeshing ONLY accepts triangular meshes as input." << std::endl;
        std::cout << "############################" << std::endl;
        std::cout << std::endl;

        std::string filename_mesh = meshFiles.getValue().at(0);

        MUSE::SurfaceMeta geometa;
        geometa.setProject(Project);

        std::vector<std::string> excommands;
        excommands.push_back(command);
        geometa.setCommands(excommands);

        std::vector<std::string> deps;
        deps.push_back(filesystem::relative(get_basename(filename_mesh) + ".json", Project.folder));
        geometa.setDependencies(deps);

        cinolib::Trimesh<> mesh;
        mesh.load(filename_mesh.c_str());
        std::cout << "\033[0;32mLoading mesh file: " << filename_mesh << " ... COMPLETED.\033[0m" << std::endl;
        std::cout << std::endl;


        MUSE::Surface surf;
        MUSE::Surface::Parameters surf_par;

        if(setRotAxis.isSet())
        {
            double rad = (setRotAngle.getValue() * M_PI)/180;
            cinolib::vec3d axis = set_rotation_axis(setRotAxis.getValue());

            cinolib::mat3d R = cinolib::mat3d::ROT_3D(axis, rad);
            cinolib::vec3d rotcenter {setRotCenterX.getValue(), setRotCenterY.getValue(), setRotCenterZ.getValue()};

            for(uint vid=0; vid<mesh.num_verts(); vid++)
            {
                mesh.vert(vid) -= rotcenter;
                mesh.vert(vid) = R*mesh.vert(vid);
                mesh.vert(vid) += rotcenter;
            }
        }

        if(setMarkedEdge.isSet())
        {
            mesh.edge_mark_boundaries();
            std::cout << "### Remeshing with marked boundary edges." << std::endl;
        }

        std::cout << "### Remeshing fixed on: mean edge." << std::endl;
        remesh_Botsch_Kobbelt_2004(mesh, -1, setMarkedEdge.getValue());

        std::cout << "Remeshing ... " << mesh.num_verts() << "V / " << mesh.num_edges() << "E / " << mesh.num_polys() << "P" << std::endl;

        surf.setParameters(surf_par);
        surf.setSummary(mesh);
        geometa.setMeshSummary(surf);

        std::string out_mesh = out_surf + "/"+ get_basename(get_filename(filename_mesh)) + "_rem";
        mesh.save((out_mesh + ext_surf).c_str());
        std::cout << "\033[0;32mSaving mesh file: " << out_mesh + ext_surf << "\033[0m" << std::endl;

        if(setRotAxis.isSet())
        {
            out_mesh += "_rot";
            mesh.save((out_mesh + ext_surf).c_str());
            std::cout << "\033[0;32mSaving mesh file: " << out_mesh + ext_surf << "\033[0m" << std::endl;
        }

        geometa.write(out_mesh + ".json");
    }

    } catch (ArgException &e)  // catch exceptions
    { std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; }
}
