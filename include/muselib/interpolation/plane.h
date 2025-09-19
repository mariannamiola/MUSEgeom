#ifndef PLANE_H
#define PLANE_H

#include <cmath>
#include <string>
#include <vector>
#include <iostream>
#include <algorithm>

#include <cinolib/find_intersections.h>
#include <cinolib/geometry/plane.h>

#include <muselib/data_structures/point.h>

using namespace MUSE;

struct fittedPlane
{
    double meanX;
    double meanY;
    double meanZ;
    double meanA0;
    double meanA1;
    double b;
};


bool intersectPlaneSegment  (const cinolib::Plane& plane, const cinolib::Segment& segment, cinolib::vec3d& intersection);
fittedPlane fitPlane        (const std::vector<Point3D> &points);


#ifndef STATIC_MUSELIB
#include "plane.cpp"
#endif

#endif // PLANE_H
