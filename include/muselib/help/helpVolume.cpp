#include <iostream>
#include "helpVolume.h"

void printHelpVolume()
{
    std::cout << "\n";
    std::cout << "================ Volume Mesh Generation Help ================\n";
    std::cout << "This module allows you to convert a surface mesh into a 3D volume mesh.\n\n";
    std::cout << "Available Flags:\n";
    std::cout << "  --createVolObject       Activate volume mesh generation\n";
    std::cout << "  --tet                   Generate tetrahedral mesh using TetGen\n";
    std::cout << "  --vox                   Generate voxelized hexahedral mesh\n";
    std::cout << "  --hex                   Generate structured hexahedral mesh from a surface\n";
    std::cout << "  --opt <params>          Optional TetGen parameters (e.g., 'pq1.2a0.1')\n";
    std::cout << "  --save                  Save the translated surface mesh (before tetrahedralization)\n";
    std::cout << "  --maxVoxelPerSide <N>   Set voxel resolution per side (for voxel grid)\n";
    std::cout << "  --resX --resY --resZ    Resolution in X/Y/Z (for structured hex mesh)\n";
    std::cout << "  meshFile                Input surface mesh (only one file allowed)\n";
    std::cout << "\n";
    std::cout << "Examples:\n";
    std::cout << "  ./app --createVolObject --tet model.obj --opt 'pq1.2a0.1'\n";
    std::cout << "  ./app --createVolObject --vox model.obj --maxVoxelPerSide 100\n";
    std::cout << "  ./app --createVolObject --hex model.obj --resX 50 --resY 50 --resZ 30\n";
    std::cout << "\n";
    std::cout << "==============================================================\n";
    std::cout << std::endl;
}
