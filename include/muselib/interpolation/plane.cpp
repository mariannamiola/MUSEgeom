#include "plane.h"

using namespace MUSE;

///
/// \brief intersectPlaneSegment
/// \param plane
/// \param segment
/// \param intersection
/// \return
///
bool intersectPlaneSegment(const cinolib::Plane& plane, const cinolib::Segment& segment, cinolib::vec3d& intersection)
{
    cinolib::vec3d dir = segment.v[1] - segment.v[0]; // Segment direction
    double denom = plane.n.dot(dir);

    //std::cout << "dir: " << dir << std::endl;
    //std::cout << "denom: " << denom << std::endl;


    // Check if segment is parallel to the plane
    if (std::abs(denom) < 1e-6) {
        return false; // No intersection (parallel)
    }

    // Compute intersection parameter t
    double t = plane.n.dot(plane.p - segment.v[0]) / denom;

    // Check if intersection is within segment bounds
    if (t < 0.0 || t > 1.0) {
        return false; // Intersection is outside the segment
    }

    // Compute intersection point
    intersection = segment.v[0] + t * dir;
    return true;
}


///
/// \brief fitPlane
/// \param points
/// \return
///
fittedPlane fitPlane(const std::vector<Point3D> &points)
{
    fittedPlane parameters;

    std::cout << "Fitting points by a Plane in 3 Dimensions ..." << std::endl;

    //Compute the mean of the points
    std::cout << "N. points = " << points.size() << std::endl;

    Point3D mean;
    mean.x = 0.0;
    mean.y = 0.0;
    mean.z = 0.0;

    for(unsigned int i=0; i<points.size(); i++)
    {
        mean.x+=points.at(i).x;
        mean.y+=points.at(i).y;
        mean.z+=points.at(i).z;
    }

    mean.x=mean.x/points.size();
    mean.y=mean.y/points.size();
    mean.z=mean.z/points.size();
    std::cout << "Mean-point = " << mean.x << "; " << mean.y << "; " << mean.z << std::endl;

    //Compute the linear system matrix and vector elements
    double sum_xx=0.0, sum_xy=0.0, sum_xz=0.0, sum_yy=0.0, sum_yz=0.0;

    for(unsigned int i=0; i< points.size(); i++)
    {
        Point3D d;
        d.x = points.at(i).x-mean.x;
        d.y = points.at(i).y-mean.y;
        d.z = points.at(i).z-mean.z;

        sum_xx += d.x*d.x;
        sum_xy += d.x*d.y;
        sum_xz += d.x*d.z;
        sum_yy += d.y*d.y;
        sum_yz += d.y*d.z;
    }

    //Solve the linear system
    double det = sum_xx*sum_yy - sum_xy*sum_xy;
    if(det!=0)
    {
        //Compute the fitted plane h(x,y) = meanZ + meanA0*(x-meanX) + meanA1*(y-meanY)
        parameters.meanX = mean.x;
        parameters.meanY = mean.y;
        parameters.meanZ = mean.z;

        parameters.meanA0 = (sum_yy*sum_xz - sum_xy*sum_yz)/det;
        parameters.meanA1 = (sum_xx*sum_yz - sum_xy*sum_xz)/det;
        parameters.b = parameters.meanZ - parameters.meanA0 * parameters.meanX - parameters.meanA1 * parameters.meanY;
    }
    else
    {
        parameters.meanX = 0;
        parameters.meanY = 0;
        parameters.meanZ = 0;

        parameters.meanA0 = 0;
        parameters.meanA1 = 0;
        parameters.b = 0;
    }

    std::cout << "Plane coefficients: meanA0 = " << parameters.meanA0 << "; meanA1 = " << parameters.meanA1 <<  "; b = " << parameters.b << std::endl;

    std::cout << "Fitting simplificated points by a Plane in 3 Dimensions ... COMPLETED." << std::endl;

    return parameters;
}



