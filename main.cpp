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

std::vector<Point3> points;
std::vector<K::Point_3> normals;

struct Face {
    std::vector<unsigned long> vertices; // indices in vector of points
};

struct Shell {
    std::vector<Face> faces;
};

struct Object {
    std::string id;
    std::vector<Shell> shells;
};

std::map<std::string, Object> objects;

void loadObjFile(const std::string& filename) {
    std::ifstream infile(filename);
    if (!infile.is_open()) {
        std::cerr << "Could not open file " << filename << " for reading." << std::endl;
        return;
    }
    std::string line;
    Object current_object;
    Shell current_shell;
    int default_object_id = 1;
    bool default_object_created = false;
    while (std::getline(infile, line)) {
        std::stringstream ss(line);
        std::string type;
        ss >> type;

        if (type == "v") {
            // Vertex data
            double x, y, z;
            ss >> x >> y >> z;
            points.emplace_back(x, y, z);
            //std::cout << "Vertex: (" << x << ", " << y << ", " << z << ")" << std::endl;
        }

        if (type == "vn"){
            //Vertex Normals
            double x, y, z;
            ss >> x >> y >> z;
            normals.emplace_back(x, y, z);
            //std::cout << "Normal: (" << x << ", " << y << ", " << z << ")" << std::endl;
        }

        if (type == "g") {
            // Start of a new object, add current shell to current object
            if (!current_shell.faces.empty()) {
                current_object.shells.push_back(current_shell);
                current_shell = Shell();
            }

            // Add current object to objects map
            if (!current_object.id.empty()) {
                objects[current_object.id] = current_object;
            }

            // Start a new object
            current_object = Object();
            ss >> current_object.id;
            //std::cout << "Object: " << current_object.id << std::endl;
        }

        if (type == "f"){
            //Face data
            Face face;
            char ch;
            double x, y, z;
            double a, b, c;
            ss >> x >> ch >> ch >> a >> y >> ch >> ch >> b >> z >> ch >> ch >> c ;
            face.vertices.push_back(x - 1); // .obj file indices are 1-based
            face.vertices.push_back(y - 1);
            face.vertices.push_back(z - 1);
            current_shell.faces.push_back(face);
            //std::cout << "Face: " << x << " " << y << " " << z << std::endl;
        }
    }
    infile.close();

    // Add last shell to current object
    if (!current_shell.faces.empty()) {
        current_object.shells.push_back(current_shell);
    }

    // Add last object to objects map
    if (!current_object.id.empty()) {
        objects[current_object.id] = current_object;
    }

    // If no objects were created, create a default object
    if (objects.empty() && !default_object_created) {
        current_object.id = std::to_string(default_object_id);
        current_object.shells.push_back(current_shell);
        objects[current_object.id] = current_object;
        default_object_created = true;
        std::cout << "Default object created with ID: " << current_object.id << std::endl;
    }


    for (const auto& [id, object] : objects) {
        std::cout << "Object " << id << ":" << std::endl;
        for (const auto& shell : object.shells) {
            std::cout << "  Shell: " << std::endl;
            for (const auto& face : shell.faces) {
                std::cout << "    Face: ";
                for (const auto& vertex : face.vertices) {
                    std::cout << vertex << " ";
                }
                std::cout << std::endl;
            }
        }
    }

}

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


void create_voxelgrid(std::vector<Point3> allpoints, double res) {

    // maxminlist = { max_x, min_x, max_y, min_y, max_z, min_z }
    std::vector<double> maxminlist = maxmin_coo(allpoints);

    int rows_x = int((maxminlist[0] - maxminlist[1] / res) + 1) + 2;
    int rows_y = int((maxminlist[2] - maxminlist[3] / res) + 1) + 2;
    int rows_z = int((maxminlist[4] - maxminlist[5] / res) + 1) + 2;
    VoxelGrid voxels(rows_x, rows_y, rows_z);
}

Point3 modelcoo_to_voxcoo(Point3 point, std::vector<Point3> allpoints, double res) {
    std::vector<double> maxminlist = maxmin_coo(allpoints);
    int position_x = int((point.x() - (int(maxminlist[1]) - res)) / res);
    int position_y = int((point.y() - (int(maxminlist[3]) - res)) / res);
    int position_z = int((point.y() - (int(maxminlist[5]) - res)) / res);
    return Point3(position_x, position_y, position_z);
}

Point3 voxcoo_to_modelcoo(Point3 point, std::vector<Point3> allpoints, int res){
    std::vector<double> maxminlist = maxmin_coo(allpoints);
    double x = (maxminlist[1] - res) + (point[0] * res);
    double y = (maxminlist[3] - res) + (point[0] * res);
    double z = (maxminlist[5] - res)+ (point[0] * res);
    //Gives minimum coordinate of the voxel
    return Point3(x, y, z);
}


int voxelisation(double res) {
    std::cout << "line 181" << std::endl;
    std::cout << "size of object, shell: " << objects.size() << std::endl;
    for (const auto& object : objects) {
        for (const auto &shell: object.second.shells) {
            std::cout << "line 184" << std::endl;
            for (const auto &face: shell.faces) {
                // computing axis-oriented bounding box of a triangle
                // -> the min & max of x, y, z
                std::vector<Point3> tri_coo;
                for (auto faceid: face.vertices) {
                    tri_coo.push_back(points[faceid]);
                }
                std::vector<double> maxmin_xyz = maxmin_coo(tri_coo);
                Point3 min_pt(maxmin_xyz[0], maxmin_xyz[2], maxmin_xyz[4]);
                Point3 min_pt_vox = modelcoo_to_voxcoo(min_pt, points, 0.1);
                Point3 max_pt(maxmin_xyz[1], maxmin_xyz[3], maxmin_xyz[5]);
                Point3 max_pt_vox = modelcoo_to_voxcoo(max_pt, points, 0.1);
                std::vector<Point3> voxel_coordinates;
                std::vector<double> yvoxels;
                std::vector<double> zvoxels;
                std::cout << "line 198" << std::endl;

                for (double i = min_pt_vox.x(); i <= max_pt_vox.x(); ++res) {
                    for (double j = min_pt_vox.y(); j <= max_pt_vox.y(); ++res) {
                        for (double k = min_pt_vox.z(); k <= max_pt_vox.z(); ++res) {
                            Point3 vox(i, j, k);
                            voxel_coordinates.push_back(vox);
                        }
                    }
                }
                std::cout << "size of voxel coord: " << voxel_coordinates.size() << std::endl;

            }
        }
    }
    return 0;
}


int main (int argc, const char * argv[]) {
    const char* filename = (argc > 1) ? argv[1] : "/Users/putiriyadi/Documents/TU/q3/3d/hw03/Duplex_A_20110907.obj";
    loadObjFile(filename);
//    create_voxelgrid(points, 0.1);
    voxelisation(0.1);

    return 0;
}



