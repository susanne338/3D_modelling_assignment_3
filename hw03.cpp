//
// Created by P Riyadi on 03/04/2023.
//

#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <list>
#include <cassert>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polygon_2.h>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef K::Point_2                  Point2;
typedef K::Point_3                  Point3;
typedef CGAL::Polygon_2<K>          Polygon2;
typedef K::Plane_3                  Plane;
//const std::string input_file = Duplex_A_20110907.obj;

std::vector<Point3> points;

//struct Vertex
//{
//    // Position Vector
//    std::vector<Point3> points;
//    // Normal Vector
//    std::vector<Point3> normal;
//    // Texture Coordinate Vector
//    std::vector<Point2> TextureCoordinate;
//};

struct Face {
    std::vector<unsigned long> verticesid; // indices in vector of points
};

struct Shell {
    std::vector<Face> faces;
};

struct Object {
    std::string id;
    std::vector<Shell> shells;
};

struct Material {

};

std::map<std::string, Object> objects;
std::map<std::string, Shell> shells;
std::map<int, Face> faces;

struct VoxelGrid {
    std::vector<unsigned int> voxels;
    unsigned int max_x, max_y, max_z;

    VoxelGrid(unsigned int x, unsigned int y, unsigned int z) {
        max_x = x;
        max_y = y;
        max_z = z;
        unsigned int total_voxels = x*y*z;
        voxels.reserve(total_voxels);
        for (unsigned int i = 0; i < total_voxels; ++i) voxels.push_back(0);
    }

    unsigned int &operator()(const unsigned int &x, const unsigned int &y, const unsigned int &z) {
        assert(x >= 0 && x < max_x);
        assert(y >= 0 && y < max_y);
        assert(z >= 0 && z < max_z);
        return voxels[x + y*max_x + z*max_x*max_y];
    }

    unsigned int operator()(const unsigned int &x, const unsigned int &y, const unsigned int &z) const {
        assert(x >= 0 && x < max_x);
        assert(y >= 0 && y < max_y);
        assert(z >= 0 && z < max_z);
        return voxels[x + y*max_x + z*max_x*max_y];
    }
};

std::vector<double> maxmin_coo(std::vector<Point3> allpoints) {
    std::vector<double> maxmin;
    std::vector<double> points_x, points_y, points_z;
    for (auto &point : allpoints) {
        points_x.push_back(point.x());
        points_y.push_back(point.y());
        points_z.push_back(point.z());
    }
    double max_x = *std::max_element(points_x.begin(), points_x.end());
    double min_x = *std::min_element(points_x.begin(), points_x.end());
    double max_y = *std::max_element(points_y.begin(), points_y.end());
    double min_y = *std::min_element(points_y.begin(), points_y.end());
    double max_z = *std::max_element(points_z.begin(), points_z.end());
    double min_z = *std::min_element(points_z.begin(), points_z.end());
    maxmin.insert(maxmin.end(), {max_x, min_x, max_y, min_y, max_z, min_z});
    return maxmin;
};


void create_voxelgrid(std::vector<Point3> allpoints, int res) {

    // maxminlist = { max_x, min_x, max_y, min_y, max_z, min_z }
    std::vector<double> maxminlist = maxmin_coo(allpoints);

    int rows_x = int(((int(maxminlist[0])+1 - int(maxminlist[1])) / res) + 1) + 2;
    int rows_y = int(((int(maxminlist[2])+1 - int(maxminlist[3])) / res) + 1) + 2;
    int rows_z = int(((int(maxminlist[4])+1 - int(maxminlist[5])) / res) + 1) + 2;
    VoxelGrid voxels(rows_x, rows_y, rows_z);
}

Point3 modelcoo_to_voxcoo(Point3 point, std::vector<Point3> allpoints, int res) {
    std::vector<double> maxminlist = maxmin_coo(allpoints);
    int position_x = int((point.x() - (int(maxminlist[1]) - res)) / res);
    int position_y = int((point.y() - (int(maxminlist[3]) - res)) / res);
    int position_z = int((point.y() - (int(maxminlist[5]) - res)) / res);
    return Point3(position_x, position_y, position_z);
}



int main() {


}


