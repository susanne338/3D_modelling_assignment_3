#include <iostream>
#include <map>
#include <list>
#include <cassert>

#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polygon_2.h>
#include <utility>
#include <vector>
#include <fstream>
#include <sstream>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/repair_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include "json.hpp"
using json = nlohmann::json;

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef K::Point_2                  Point2;
typedef K::Point_3                  Point3;
typedef CGAL::Polygon_2<K>          Polygon2;
typedef K::Plane_3                  Plane3;
typedef K::Line_3                   Line3;
//const std::string input_file = Duplex_A_20110907.obj;


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

//LOAD OBJ FILE INTO MEMORY-------------------------------------------------------------------------------------------
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

//VOXELGRID------------------------------------------------------------------------------------------------------------
std::vector<double> maxmin_coo(const std::vector<Point3>& allpoints) {
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


VoxelGrid create_voxelgrid(const std::vector<Point3>& allpoints, double res) {

    // maxminlist = { max_x, min_x, max_y, min_y, max_z, min_z }
    std::vector<double> maxminlist = maxmin_coo(allpoints);

    int rows_x = int(((maxminlist[0] - maxminlist[1]) / res) + 1) + 2;
    int rows_y = int(((maxminlist[2] - maxminlist[3]) / res) + 1) + 2;
    int rows_z = int(((maxminlist[4] - maxminlist[5]) / res) + 1) + 2;
    VoxelGrid voxels(rows_x, rows_y, rows_z) ;
    return voxels;
}

std::vector<double> create_voxelgrid2(const std::vector<Point3>& allpoints, double res) {

    // maxminlist = { max_x, min_x, max_y, min_y, max_z, min_z }
    std::vector<double> maxminlist = maxmin_coo(allpoints);
    std::vector<double> rows;

    int rows_x = int(((maxminlist[0]- maxminlist[1]) / res) + 1) + 2;
    int rows_y = int(((maxminlist[2] - maxminlist[3]) / res) + 1) + 2;
    int rows_z = int(((maxminlist[4] - maxminlist[5]) / res) + 1) + 2;
    rows.push_back(rows_x);
    rows.push_back(rows_y);
    rows.push_back(rows_z);
    return rows;
}


Point3 modelcoo_to_voxcoo(Point3 point, const std::vector<Point3>& allpoints, double res) {
    std::vector<double> maxminlist = maxmin_coo(allpoints);
    int position_x = int(abs((point.x() - (maxminlist[1] - res))) / res);
    int position_y = int(abs((point.y() - (maxminlist[3] - res))) / res);
    int position_z = int(abs((point.z() - (maxminlist[5] - res))) / res);
    return {position_x, position_y, position_z};
}

Point3 voxcoo_to_modelcoo(Point3 point, const std::vector<Point3>& allpoints, double res){
    std::vector<double> maxminlist = maxmin_coo(allpoints);
    double x = (maxminlist[1] - res) + (point[0] * res);
    double y = (maxminlist[3] - res) + (point[0] * res);
    double z = (maxminlist[5] - res) + (point[0] * res);
    //Gives minimum coordinate of the voxel
    return {x, y, z};

}

//VOXELISATION---------------------------------------------------------------------------------------------------------
int voxelisation(double res, VoxelGrid grid) {
    std::cout << "size of object, shell: " << objects.size() << std::endl;
    std::map<int, Point3> vertices_dict;
    int vid = 0;
    for (const auto& object : objects) {

        for (const auto &shell: object.second.shells) {
            std::vector<Point3> voxel_intersect;
            for (const auto &face: shell.faces) {
                // computing axis-oriented bounding box of a triangle
                // -> the min & max of x, y, z
                std::vector<Point3> tri_coo;
                for (auto faceid: face.vertices) {
                    tri_coo.push_back(points[faceid]);
                }
//                std::cout <<"triangle coord: " << tri_coo[0] << "; " << tri_coo[1] << "; " << tri_coo[2] << std::endl;
                std::vector<double> maxmin_xyz = maxmin_coo(tri_coo);
                Point3 max_pt(maxmin_xyz[0], maxmin_xyz[2], maxmin_xyz[4]);
                Point3 max_pt_vox = modelcoo_to_voxcoo(max_pt, points, res);
                Point3 min_pt(maxmin_xyz[1], maxmin_xyz[3], maxmin_xyz[5]);
                Point3 min_pt_vox = modelcoo_to_voxcoo(min_pt, points, res);
                std::vector<Point3> voxel_coordinates;
//                std::vector<Point3> voxel_intersect;
//                std::cout << "max pt: " << max_pt << std::endl;
//                std::cout << "min pt: " << min_pt << std::endl;
//                std::cout << "max pt voxel: " << max_pt_vox << std::endl;
//                std::cout << "min pt voxel: " << min_pt_vox << std::endl;
                int step = 1;
                Plane3 triangle(tri_coo[0], tri_coo[1], tri_coo[2]);

                int i = int(min_pt_vox.x());
                while (i <= max_pt_vox.x()) {
                    int j = int(min_pt_vox.y());
                    while (j <= max_pt_vox.y()) {
                        int k = int(min_pt_vox.z());
                        while (k <= max_pt_vox.z()) {
                            Point3 vox(i, j, k);
//                            std::cout << i << " " << j << " " << k << std::endl;
                            voxel_coordinates.push_back(vox);

                            // 26-connectivity
                            // (a) top - down line
                            Point3 a_start(vox.x() + (0.5 * res), vox.y() + (0.5 * res), vox.z());
                            Point3 a_end(vox.x() + (0.5 * res), vox.y() + (0.5 * res), vox.z() + res);
                            Line3 a_line(a_start, a_end);
                            // (b) front - back line
                            Point3 b_start(vox.x() + (0.5 * res), vox.y(), vox.z() + 0.5 * res);
                            Point3 b_end(vox.x() + (0.5 * res), vox.y() + res, vox.z() + 0.5 * res);
                            Line3 b_line(b_start, b_end);
                            // (c) left - right line
                            Point3 c_start(vox.x(), vox.y() + (0.5 * res), vox.z() + 0.5 * res);
                            Point3 c_end(vox.x() + res, vox.y() + (0.5 * res), vox.z() + 0.5 * res);
                            Line3 c_line(c_start, c_end);
                            // check intersection triangle(plane) & line
                            bool intersects_linea = CGAL::do_intersect(a_line, triangle);
                            bool intersects_lineb = CGAL::do_intersect(b_line, triangle);
                            bool intersects_linec = CGAL::do_intersect(c_line, triangle);
                            bool intersects_any = intersects_linea || intersects_lineb || intersects_linec;

                            if (intersects_any) {
                                voxel_intersect.push_back(vox);
                                grid(i, j, k) = 1;
                                Point3 real_coo = voxcoo_to_modelcoo(vox, points, res);
                                vertices_dict[vid] = real_coo;
                            }
                            vid++;
                            k += step;
                        }
                        j += step;
                    }
                    i += step;
                } // end looping the triangle bounding box (max min of x, y, z)
            }
            std::cout << "size of voxel intersect: " << voxel_intersect.size() << std::endl;
            std::vector<Point3> coo_list_intersect;
        }
    }
    return 0;
}


//MARKING----------------------------------------------------------------------------------------------------------------
std::vector<Point3> vector_six_connect(Point3 coordinate, const VoxelGrid& grid){ //neighbouring voxels
    std::vector<Point3> vector;
    vector.push_back(Point3(coordinate.x()-1, coordinate.y(), coordinate.z()));
    vector.push_back(Point3(coordinate.x()+1, coordinate.y(), coordinate.z()));
    vector.push_back(Point3(coordinate.x(), coordinate.y()-1, coordinate.z()));
    vector.push_back(Point3(coordinate.x(), coordinate.y()+1, coordinate.z()));
    vector.push_back(Point3(coordinate.x(), coordinate.y(), coordinate.z()-1));
    vector.push_back(Point3(coordinate.x(), coordinate.y(), coordinate.z()+1));
    std::vector<Point3> new_vector;
    for (auto &p : vector){
        if (  p.x() >grid.max_x-1 or p.x() < 0 or p.y() > grid.max_y-1 or p.y() < 0 or p.z() > grid.max_z-1 or p.z() < 0){
            continue;
        }
        if(grid(p.x(), p.y(), p.z()) != 0){
            continue;
        }
        else{
            new_vector.push_back(p);
        }
    }
    return new_vector;
}

//6CONNECT FOR NON-RECURSIVE FUNCTION
std::vector<Point3> vector_six_connect_2(Point3 coordinate, const VoxelGrid& grid){ //neighbouring voxels
    std::vector<Point3> vector1;
    vector1.push_back(Point3(coordinate.x() - 1, coordinate.y(), coordinate.z()));
    vector1.push_back(Point3(coordinate.x() + 1, coordinate.y(), coordinate.z()));
    vector1.push_back(Point3(coordinate.x(), coordinate.y()-1, coordinate.z()));
    vector1.push_back(Point3(coordinate.x(), coordinate.y()+1, coordinate.z()));
    vector1.push_back(Point3(coordinate.x(), coordinate.y(), coordinate.z()-1));
    vector1.push_back(Point3(coordinate.x(), coordinate.y(), coordinate.z()+1));

    std::vector<Point3> new_vector;
    for (auto &p : vector1){
        if (  p.x() >grid.max_x or p.x() < 0 or p.y() > grid.max_y or p.y() < 0 or p.z() > grid.max_z or p.z() < 0){
            continue;
        }
        else{
            new_vector.push_back(p);
        }
    }
    std::cout << "vector1: " << vector1[0] << ' ' << vector1[1] << ' ' << vector1[2] << std::endl;
    std::cout << "new vector: " << new_vector[0] << ' ' << new_vector[1] << ' ' << new_vector[2] << std::endl;
    return new_vector;
}

//SMALLER RECURSIVE FUNCTION
//void marking2(Point3 coordinate, int id, VoxelGrid& grid){
//    std::vector<Point3> vector = vector_six_connect(coordinate, grid);
//    grid(0,0,0) = 2;
//    for (auto &point : vector) {
//        grid(point.x(), point.y(), point.z()) = id;
//        marking2(point, id, grid);
//
//    }
//}
//MARKING FUNCTION
void marking(Point3 coordinate, int id, VoxelGrid& grid){
    std::vector<Point3> vector = vector_six_connect(coordinate, grid);
    std::cout << "vector size" << vector.size() <<std::endl;
    grid(coordinate.x(), coordinate.y(), coordinate.z()) = id;
    std::cout << " id element" << grid(coordinate.x(), coordinate.y(), coordinate.z()) << std::endl;
    for (auto &point : vector) {
        if (grid(point.x(), point.y(), point.z()) == 0) {
//            grid(point.x(), point.y(), point.z()) = id;
            Point3 new_coord = point;
            std::cout << "coord: " << point.x() << "; " << point.y() << "; " << point.z() << std::endl;
            std::cout << "id of point" << grid(point.x(), point.y(), point.z());
            marking(new_coord, id, grid);
        }
        else {
            continue;
        }

    }
}

//NON-RECURSIVE MARKING EXTERIOR FUNCTION. USES DIFFERENT 6-CONNECT FUNCTION!!
void marking_exterior(VoxelGrid& grid){
    //marking boundary
    for(int i =0; i < grid.max_x; ++i){
        grid(i,0,0) = 2;
    }
    for(int j =0; j < grid.max_y; ++j){
        grid(0,j,0) = 2;
    }
    for(int k =0; k < grid.max_z; ++k){
        grid(0,0,k) = 2;
    }
    for(int n = 1; n < grid.max_x -1; ++n){
        for(int o = 1; o < grid.max_y -1; ++o){
            for(int p = 1; p < grid.max_z-1; ++p){

                if (grid(n,o,p) == 0){
                    std::vector<Point3> conn = vector_six_connect_2(Point3(n,o,p), grid);
                    for(auto&pt:conn){
                        if(grid(pt.x(), pt.y(), pt.z()) == 2){
                            grid(n,o,p) = 2;
                            break;
                        }
                        else {
                            continue;
                        }
                    }
                }
                else {
                    continue;
                }
                std::cout << "id of point" <<n<<o<<p << " is " <<grid(n,o,p) << std::endl;

            }
        }
    }
}


//WRITE CITYJSON--------------------------------------------------------------------------------------------------------
void write_cityjson(std::string outputfile,
                    std::map<int, Point3>& vertices_dict,
                    std::vector<std::vector<int>>& markedsurface,
                    std::vector<std::string> buildingtype){
    /*
     * input:   vertex dictionary
     *              {id -> Point3}
     *          nested list of marked surface
     *              [[1, 2, 3, 4, 5], [3, 4, 5, 6] ...]
     *          list of building type which its position refers to the marked surface list
     *              ["Building", "BuildingRoom", ...]  -> option: "Building"/"BuildingRoom"
     */

    std::vector<Point3> vertices_all;
    nlohmann::json json;
    json["type"] = "CityJSON";
    json["version"] = "1.1";
    json["transform"] = nlohmann::json::object();
    json["transform"]["scale"] = nlohmann::json::array({1.0, 1.0, 1.0});
    json["transform"]["translate"] = nlohmann::json::array({0.0, 0.0, 0.0});
    json["CityObjects"] = nlohmann::json::object();

    for (auto& type : buildingtype) {
        json["vertices"] = nlohmann::json::array({});
        nlohmann::json surface_array = nlohmann::json::array({});
        for (auto& eachsurface : markedsurface) {
            nlohmann::json surface;
            for (auto sid : eachsurface) {
                surface.push_back(sid);
                Point3 vert_pts = vertices_dict[sid];
                nlohmann::json jvertex;
                jvertex.push_back(vert_pts.x());
                jvertex.push_back(vert_pts.y());
                jvertex.push_back(vert_pts.z());
                json["vertices"].push_back(jvertex);
            }
            surface_array.push_back(surface);
        }
        json["CityObjects"][type] = nlohmann::json::object();
        json["CityObjects"][type]["type"] = "Solid";
        json["CityObjects"][type]["geometry"] = nlohmann::json::object();
        json["CityObjects"][type]["geometry"]["type"] = "MultiSurface";
        json["CityObjects"][type]["geometry"]["boundaries"] = nlohmann::json::array();
        json["CityObjects"][type]["geometry"]["boundaries"].push_back(surface_array);
    }

    std::string json_string = json.dump(2);
    std::ofstream out_stream(outputfile);
    out_stream << json_string;
    out_stream.close();
}




//MAIN-----------------------------------------------------------------------------------------------------------------
int main(int argc, const char * argv[]) {
    double res = 1.0;

    //file
    const char* filename = (argc > 1) ? argv[1] : "/Users/putiriyadi/Documents/TU/q3/3d/hw03/Duplex_A_20110907_nofurni.obj";

    //Writing to obj
    loadObjFile(filename);


// -------------------new
    //trying out stuff with POINTERS, this makes the voxelisation run on my computer.
//    double rowx = create_voxelgrid2(points,res)[0];
//    double rowy = create_voxelgrid2(points,res)[1];
//    double rowz = create_voxelgrid2(points,res)[2];
//    VoxelGrid *pointer_vox;
//    VoxelGrid voxels(rowx,rowy,rowz) ;
//    pointer_vox = &voxels;
//
//
//    //Voxelisation
//    voxelisation(res, voxels);
//
//    //Marking
//    //marking exterior
//    marking(Point3(0,0,0), 2, voxels);
// -------------------end of new

    //create voxelgrid
    VoxelGrid grid = create_voxelgrid(points, res);
    std::cout << "max x: " << grid.max_x << std::endl;
    std::cout << "max y: " << grid.max_y << std::endl;
    std::cout << "max z: " << grid.max_z << std::endl;


    //Voxelisation
    voxelisation(res, grid);

//    marking_exterior(grid);
//    std::cout <<"vector six func: " << vector_six_connect(Point3(0, 0, 0), grid)[0][0] << std::endl;

    //Marking
    //marking exterior
    marking(Point3(0, 0, 0), 2, grid); //--> thus this would mark the exterior with id 2, doesnt work
    std::cout << "exterior marked" << std::endl;
    //marking rooms
//    int id_rooms = 3;
//    for (int i = 0; i < grid.max_x; ++i) { //--> marks rooms with each a different id
//        for (int j = 0; j < grid.max_y; ++j){
//            for (int k = 0; k < grid.max_z; ++k){
//                if (grid(i,j,k) == 0) {
//                    marking(Point3(i, j, k), id_rooms, grid);
//                    id_rooms++;
//                }
//                else {
//                    continue;
//                }
//
//            }
//        }
//    }
//    //end of marking rooms
//    std::cout << "does it reach <<" << std::endl;

    return 0;


}



