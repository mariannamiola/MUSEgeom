#include "tools.h"

#include <float.h>
#include <math.h>
#include <vector>

#include <deque>
#include <tuple>

#include <cinolib/geometry/aabb.h>
#include <cinolib/geometry/plane.h>

bool check_index (const std::vector<int> &id_dupl, int index)
{
    bool id_found = false;
    for(int i : id_dupl)
    {
        if(i == index)
        {
            //std::cout << "indici uguali" << std::endl;
            id_found = true;
            break;
        }
    }
    return id_found;
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

double point_to_line_distance (const double &x0, const double &y0, const double &angle_rad, const double &x, const double &y)
{
    //https://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line#A_vector_projection_proof
    double dist = abs(cos(angle_rad)*(y0-y) - sin(angle_rad)*(x0-x));
    return dist;
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

bool point_in_polygon (const Point2D p, const std::vector<Point2D> &boundaries)
{
    double min_x =  DBL_MAX;
    double min_y =  DBL_MAX;
    double max_x = -DBL_MAX;
    double max_y = -DBL_MAX;

    //Trova minimo e massimo x,y per i punti di bordo
    for (const Point2D &bp : boundaries)
    {
        if (bp.x < min_x) min_x = bp.x;
        if (bp.x > max_x) max_x = bp.x;
        if (bp.y < min_y) min_y = bp.y;
        if (bp.y > max_y) max_y = bp.y;
    }

    //Check se i punti sono fuori/dentro i min/max
    if (p.x < min_x || p.x > max_x || p.y < min_y || p.y > max_y )
        return false;


    // https://stackoverflow.com/questions/11716268/point-in-polygon-algorithm
    //The first line of the if control if the point's y-coord is within the edge's scope
    //The second line checks whether the test point is to the left of the line
    // -> If that is true the line drawn rightwards from the test point crosses that edge.
    int i, j, c = 0;
    unsigned int nvert = boundaries.size();

    for (i = 0, j = nvert-1; i < nvert; j = i++)
    {
        if ( ((boundaries[i].y>p.y) != (boundaries[j].y>p.y)) &&
                 (p.x < (boundaries[j].x-boundaries[i].x) * (p.y-boundaries[i].y) / (boundaries[j].y-boundaries[i].y) + boundaries[i].x) )
            c = !c;
    }
    return c;
}


void points_in_polygon (const std::vector<Point2D> &points, const std::vector<Point2D> &boundaries, std::vector<unsigned int> &id_in_points)
{
    id_in_points.clear();

    unsigned int id = 0;
    for (const Point2D &p : points)
    {
        if (point_in_polygon(p, boundaries)) //se è vero, salva l'indice del punto nel vettore id_is_points
            id_in_points.push_back(id);
        id++;
    }
}



void string_to_double_conversion_vectors (const std::vector<std::string> &values, const std::vector<std::string> &id, const std::vector<double> &xcoord, const std::vector<double> &ycoord, const std::vector<double> &zcoord,
                                          std::vector<double> &conv_values, std::vector<std::string> &corr_id, std::vector<double> &corr_xcoord, std::vector<double> &corr_ycoord, std::vector<double> &corr_zcoord)
{
    for(size_t i = 0 ; i< values.size(); i++)
    {
        std::string val_tmp = values.at(i);

        double val = 0.0;
        if(!val_tmp.empty() && val_tmp.compare("nd")!=0)
        {
            if(val_tmp.compare("*")!=0)
            {
                if(val_tmp.compare("NA")!=0)
                {
                    val = std::stod(val_tmp);
                    conv_values.push_back(val);

                    if(id.size() > 0)
                        corr_id.push_back(id.at(i));
                    if(xcoord.size() > 0)
                        corr_xcoord.push_back(xcoord.at(i));
                    if(ycoord.size() > 0)
                        corr_ycoord.push_back(ycoord.at(i));
                    if(zcoord.size() > 0)
                        corr_zcoord.push_back(zcoord.at(i));
                }
            }
        }
    }
}

void string_to_double_conversion_vectors (const std::vector<uint> &indices,
                                          const std::vector<std::string> &values,
                                          const std::vector<std::string> &id,
                                          const std::vector<double> &xcoord, const std::vector<double> &ycoord, const std::vector<double> &zcoord,
                                          std::vector<double> &conv_values, std::vector<std::string> &corr_id, std::vector<double> &corr_xcoord, std::vector<double> &corr_ycoord, std::vector<double> &corr_zcoord)
{
    for(size_t i:indices)
    {
        std::string val_tmp = values.at(i);

        double val = 0.0;
        if(!val_tmp.empty() && val_tmp.compare("nd")!=0)
        {
            if(val_tmp.compare("*")!=0)
            {
                if(val_tmp.compare("NA")!=0)
                {
                    val = std::stod(val_tmp);
                    conv_values.push_back(val);

                    if(id.size() > 0)
                        corr_id.push_back(id.at(i));
                    if(xcoord.size() > 0)
                        corr_xcoord.push_back(xcoord.at(i));
                    if(ycoord.size() > 0)
                        corr_ycoord.push_back(ycoord.at(i));
                    if(zcoord.size() > 0)
                        corr_zcoord.push_back(zcoord.at(i));
                }
            }
        }
    }
}


cinolib::vec2d segment_segment_intersection_2d (const cinolib::vec2d &p1, const cinolib::vec2d &p2, const cinolib::vec2d &v0, const cinolib::vec2d &v1)
{
    cinolib::vec2d p_inter;

    double A1 = p2.y() - p1.y();
    double B1 = p1.x() - p2.x();
    double C1 = A1 * p1.x() + B1 * p1.y();

    double A2 = v1.y() - v0.y();
    double B2 = v0.x() - v1.x();
    double C2 = A2 * v0.x() + B2 * v0.y();

    double det = A1 * B2 - A2 * B1;

    //std::cout << "determinant = " <<  det << std::endl;

    if(det != 0)
    {
        p_inter.x() = (B2 * C1 - B1 * C2) / det;
        p_inter.y() = (A1 * C2 - A2 * C1) / det;

        //std::cout << "x_intersection = " << p_inter.x() << std::endl;
        //std::cout << "y_intersection = " << p_inter.y() << std::endl;
    }
    else
    {
        p_inter.x() = DBL_MAX;
        p_inter.y() = DBL_MAX;
    }

    return p_inter;
}

//https://stackoverflow.com/questions/42740765/intersection-between-line-and-triangle-in-3d
cinolib::vec3d segment_triangle_intersection_3d (const cinolib::vec3d &q1, const cinolib::vec3d &q2, const cinolib::vec3d &p1, const cinolib::vec3d &p2, const cinolib::vec3d &p3)
{
    cinolib::vec3d p_inter;

    cinolib::vec3d p21 = p2-p1;
    cinolib::vec3d p31 = p3-p1;

    cinolib::vec3d cross = p21.cross(p31);

    cinolib::vec3d q1p1 = q1-p1;
    cinolib::vec3d q2q1 = q2-q1;

    double dot1 = q1p1.dot(cross);
    double dot2 = q2q1.dot(cross);

    double t = -(dot1/dot2);

    p_inter = q1 + t*(q2-q1);

    return p_inter;
}

cinolib::vec3d set_rotation_axis (const std::string &rot_axis)
{
    cinolib::vec3d axis (0,0,0);

    if(rot_axis.compare("X") == 0)
        axis.x() = 1;
    else if (rot_axis.compare("Y") == 0)
        axis.y() = 1;
    else if (rot_axis.compare("Z") == 0)
        axis.z() = 1;

    return axis;
}


void point_rotation (const double &x, const double &y, const double &z, const cinolib::vec3d &rot_axis, const double &rot_angle, const cinolib::vec3d &rot_center, double &x_rot, double &y_rot, double &z_rot) //angle in degree
{
    double rad = (rot_angle * M_PI)/180;

    cinolib::vec3d point (x,y,z);

    cinolib::mat3d R = cinolib::mat3d::ROT_3D(rot_axis, rad);

    point -= rot_center;
    point = R * point;
    point += rot_center;

    x_rot = point.x();
    y_rot = point.y();
    z_rot = point.z();
}


///
/// \brief point_rotation
/// \param point
/// \param rot_axis
/// \param rot_angle (in degree)
/// \param rot_center
/// \return
///
cinolib::vec3d point_rotation (cinolib::vec3d &point, const cinolib::vec3d &rot_axis, const double &rot_angle, const cinolib::vec3d &rot_center) //angle in degree
{
    double rad = (rot_angle * M_PI)/180;

    cinolib::mat3d R = cinolib::mat3d::ROT_3D(rot_axis, rad);

    point -= rot_center;
    point = R * point;
    point += rot_center;

    return point;
}

std::vector<cinolib::vec3d> points_rotation (std::vector<cinolib::vec3d> &points, const cinolib::vec3d &rot_axis, const double &rot_angle, const cinolib::vec3d &rot_center) //angle in degree
{
    double rad = (rot_angle * M_PI)/180;

    cinolib::mat3d R = cinolib::mat3d::ROT_3D(rot_axis, rad);

    for(uint i=0; i<points.size(); i++)
    {
        points.at(i) -= rot_center;
        points.at(i)  = R * points.at(i);
        points.at(i) += rot_center;
    }

    return points;
}

///
/// \brief align_points_to_xyplane
/// \param points
///
void align_points_to_xyplane (std::vector<Point3D> &points, const double tol = 1e-02)
{
    std::vector<cinolib::vec3d> data_for_plane;
    for(size_t di=0; di < points.size(); di++)
        data_for_plane.push_back(cinolib::vec3d({points.at(di).x, points.at(di).y, points.at(di).z}));

    cinolib::AABB aabb (data_for_plane);
    data_for_plane.clear();

    data_for_plane.push_back(cinolib::vec3d({aabb.min.x(), aabb.min.y(), aabb.min.z()}));
    data_for_plane.push_back(cinolib::vec3d({aabb.max.x(), aabb.min.y(), aabb.min.z()}));
    data_for_plane.push_back(cinolib::vec3d({aabb.max.x(), aabb.max.y(), aabb.min.z()}));
    data_for_plane.push_back(cinolib::vec3d({aabb.min.x(), aabb.max.y(), aabb.min.z()}));
    data_for_plane.push_back(cinolib::vec3d({aabb.min.x(), aabb.min.y(), aabb.max.z()}));
    data_for_plane.push_back(cinolib::vec3d({aabb.max.x(), aabb.min.y(), aabb.max.z()}));
    data_for_plane.push_back(cinolib::vec3d({aabb.max.x(), aabb.max.y(), aabb.max.z()}));
    data_for_plane.push_back(cinolib::vec3d({aabb.min.x(), aabb.max.y(), aabb.max.z()}));
    std::cout << "=== Computing best plane on points bounding box | vector size: " << data_for_plane.size() << std::endl;


    cinolib::Plane plane (data_for_plane);
    std::cout << "=== Plane | normal: " << plane.n << std::endl;
    cinolib::vec3d normal_xy (0,0,1);

    if(plane.n.dist(normal_xy) > tol)
    {
        //std::cout << plane.n.dist(normal_xy) << std::endl;
        std::cout << "=== Rotating points on x-y plane ..." << std::endl;

        /// Computing rotation axis
        cinolib::vec3d rot_axis = plane.n.cross(normal_xy);

        /// Computing angle between two normals
        plane.n /= plane.n.norm();
        normal_xy /= normal_xy.norm();

        double dot = plane.n.dot(normal_xy);
        double angle_rad = std::acos(dot);
        double angle_deg = angle_rad * 180.0 / M_PI;
        std::cout << "=== Rotation axis: " << rot_axis << std::endl;
        //std::cout << "dot = " << dot << std::endl;
        std::cout << "=== Rotation angle = " << angle_rad << " rad | " << angle_deg << " degree" << std::endl;

        cinolib::AABB aabb (data_for_plane);
        cinolib::vec3d center (aabb.delta_x()/2.0, aabb.delta_y()/2.0, aabb.delta_z()/2.0);
        std::cout << "=== Rotation center: " <<center.x() << "; " << center.y() << "; " << center.z() << std::endl;

        rot_axis /= rot_axis.norm();
        for(auto& point : points)
        {
            std::cout << "=== Original point: " <<point.x << "; " << point.y << "; " << point.z << std::endl;
            cinolib::vec3d sample(point.x, point.y, point.z);

            sample = point_rotation(sample, rot_axis, angle_deg, center);
            point.x = sample.x();
            point.y = sample.y();
            point.z = sample.z();
            std::cout << "=== Rotate point: " <<point.x << "; " << point.y << "; " << point.z << std::endl;
        }
    }
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

///
/// \brief resample_elevation_grid
/// \param elevation
/// \param current_res_x
/// \param current_res_y
/// \param res_target_x
/// \param res_target_y
/// \return
///
std::vector<std::vector<float>> resample_elevation_grid(const std::vector<std::vector<float>>& elevation,
                                                        float current_res_x, float current_res_y, float res_target_x, float res_target_y, float XOrigin,
                                                        float YOrigin, float &corrected_XOrigin, float &corrected_YOrigin)
{
    // Calcola i fattori di downsampling (quante celle originali cadono in una nuova:
    //risoluzione originale 1m -> risoluzione target 5m -> fattore 5)
    int factor_x = static_cast<int>(std::floor(res_target_x / std::abs(current_res_x)));
    int factor_y = static_cast<int>(std::floor(res_target_y / std::abs(current_res_y)));

    if (factor_x <= 0 || factor_y <= 0)
    {
        throw std::runtime_error("Invalid resample factor: must be > 0");
    }

    const int original_rows = elevation.size();
    const int original_cols = elevation[0].size();

    const int new_rows = original_rows / factor_y;
    const int new_cols = original_cols / factor_x;

    std::vector<std::vector<float>> downscaled;

    //Downsampling tramite media dei valori delle celle
    for (int r = 0; r < new_rows; ++r)
    {
        std::vector<float> row;
        for (int c = 0; c < new_cols; ++c)
        {
            float sum = 0.0f;
            int count = 0;

            for (int dy = 0; dy < factor_y; ++dy)
            {
                for (int dx = 0; dx < factor_x; ++dx)
                {
                    int orig_r = r * factor_y + dy;
                    int orig_c = c * factor_x + dx;

                    if (orig_r < original_rows && orig_c < original_cols)
                    {
                        sum += elevation[orig_r][orig_c];
                        ++count;
                    }
                }
            }

            float avg = (count > 0) ? (sum / count) : 0.0f;
            row.push_back(avg);
        }
        downscaled.push_back(row);
    }

    // // Calcolo nuova origine centrata sui blocchi
    // corrected_XOrigin = XOrigin + (factor_x * current_res_x) / 2.0f;
    // corrected_YOrigin = YOrigin + (factor_y * current_res_y) / 2.0f;

    // Dimensioni totali delle due griglie
    float original_width = original_cols * current_res_x;
    float original_height = original_rows * current_res_y;
    float new_width = new_cols * res_target_x;
    float new_height = new_rows * res_target_y;

    // Sposta la nuova origine in modo che il centro resti lo stesso
    corrected_XOrigin = XOrigin + (original_width - new_width) / 2.0f;
    //corrected_YOrigin = YOrigin + (original_height - new_height) / 2.0f;
    corrected_YOrigin = YOrigin + ((original_height - new_height) / 2.0f) * (current_res_y < 0 ? -1.0f : 1.0f);

    // Se invece preferisci avere l'origine in basso a sinistra (come nei GIS):
    std::reverse(downscaled.begin(), downscaled.end());
    //corrected_YOrigin = YOrigin + ((original_height - new_height) / 2.0f);  // solo se current_res_y > 0

    std::cout << "=== Resampling factor: X = " << factor_x << ", Y = " << factor_y << std::endl;
    std::cout << "=== Original size: " << elevation.size() << "x" << elevation[0].size() << std::endl;
    std::cout << "=== New size: " << downscaled.size() << "x" << downscaled[0].size() << std::endl;
    std::cout << "=== New Origin: " << corrected_XOrigin << ", " << corrected_YOrigin << std::endl;

    return downscaled;
}


///
///
//Funzione di proiezione ortogonale su un piano
cinolib::vec3d projectPointOntoPlane (const cinolib::vec3d &P, const cinolib::Plane &plane)
{
    cinolib::vec3d diff = P - plane.p;
    double dist = plane.n.dot(diff);

    return P - dist * plane.n;
};

///
/// \brief find_nearest_in_local_region
/// \param target
/// \param points
/// \param last_index
/// \param window
/// \return
///
int find_nearest_in_local_region(const cinolib::vec3d& target,
                                 const std::vector<cinolib::vec3d>& points,
                                 int last_index,
                                 int window)
{
    double min_dist = std::numeric_limits<double>::max();
    int best_idx = last_index;

    int start = std::max(0, last_index - window);
    int end   = std::min(static_cast<int>(points.size()) - 1, last_index + window);

    for (int i = start; i <= end; ++i)
    {
        double dx = target.x() - points[i].x();
        double dy = target.y() - points[i].y();
        double dist = dx * dx + dy * dy;

        if (dist < min_dist)
        {
            min_dist = dist;
            best_idx = i;
        }
    }

    return best_idx;
}

// std::vector<std::vector<std::tuple<double, double>>> buffered_points(const std::vector<std::vector<std::tuple<double, double>>>& poly_list, double distance = 1000)
// {
//     std::vector<std::vector<std::tuple<double, double>>> buff_vertices_l;
//     std::vector<std::vector<std::tuple<double, double>>> buff_vertices_t;

//     for (const auto& poly : poly_list)
//     {
//         size_t size_of_list = poly.size();
//         std::deque<std::tuple<double, double>> poly_coord;

//         for (size_t index = 0; index < size_of_list; ++index)
//         {
//             auto current_vertex = poly[index];
//             std::vector<std::tuple<double, double>> control_vertices =
//             {
//                 poly[(index - 1 + size_of_list) % size_of_list],
//                 current_vertex,
//                 poly[(index + 1) % size_of_list]
//             };
//         }
//     }
//     return buff_vertices_l; // Adjust return as needed
// }


//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::








// Remove duplicates for Point2D, considering a fixed tolerance t = 1e-06
//std::vector<Point3D> remove_duplicates (std::vector<Point3D> &points, const double t)
//{
//    // Fill index
//    for(size_t i=0; i<points.size(); i++)
//        points.at(i).id = i;

//    std::vector<Point3D> sorted_points = points;
//    std::vector<Point3D> unique_points;
//    std::vector<int> id_dupl;

//    std::sort(sorted_points.begin(), sorted_points.end(), comparePoint);

//    for(size_t i=1; i < sorted_points.size(); i++)
//    {
//        Point2D p0, p1;
//        p0.x = sorted_points.at(i-1).x;
//        p0.y = sorted_points.at(i-1).y;
//        p1.x = sorted_points.at(i).x;
//        p1.y = sorted_points.at(i).y;

//        if (dist(p0, p1) <= t)
//            id_dupl.push_back(sorted_points.at(i).id);
//    }

//    if(id_dupl.size() > 0)
//    {
//        std::sort(id_dupl.begin(), id_dupl.end());

//        for(size_t i=0; i< points.size(); i++)
//            if (!check_index(id_dupl, points.at(i).id))
//                unique_points.push_back(points.at(i));
//    }
//    else
//        unique_points = points;

//    std::cout << "Removing duplicated points ... COMPLETED. " << unique_points.size() << " left" << std::endl;
//    return unique_points;
//}

























