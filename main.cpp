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
#include <CGAL/Polyhedron_3.h>
#include <CGAL/poisson_surface_reconstruction.h>
#include <CGAL/IO/read_points.h>
#include "json.hpp"
using json = nlohmann::json;

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef K::Point_2                  Point2;
typedef K::Point_3                  Point3;
typedef K::Plane_3                  Plane3;
typedef K::Line_3                   Line3;
typedef K::Vector_3                 vector3;
typedef std::pair<Point3, vector3>  pwn;
typedef CGAL::Polyhedron_3<K> Polyhedron;
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
    double x = ((maxminlist[1] - res) + point.x()) * res;
    double y = ((maxminlist[3] - res) + point.y()) * res;
    double z = ((maxminlist[5] - res) + point.z()) * res;
    //Gives minimum coordinate of the voxel
    return {x, y, z};

}

//VOXELISATION---------------------------------------------------------------------------------------------------------
std::map<int, Point3> voxelisation(double res, VoxelGrid& grid) {
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
                std::vector<double> maxmin_xyz = maxmin_coo(tri_coo);
                Point3 max_pt(maxmin_xyz[0], maxmin_xyz[2], maxmin_xyz[4]);
                Point3 max_pt_vox = modelcoo_to_voxcoo(max_pt, points, res);
                Point3 min_pt(maxmin_xyz[1], maxmin_xyz[3], maxmin_xyz[5]);
                Point3 min_pt_vox = modelcoo_to_voxcoo(min_pt, points, res);
                std::vector<Point3> voxel_coordinates;
                int step = 1;
                Plane3 triangle(tri_coo[0], tri_coo[1], tri_coo[2]);

                int i = int(min_pt_vox.x());
                while (i <= max_pt_vox.x()) {
                    int j = int(min_pt_vox.y());
                    while (j <= max_pt_vox.y()) {
                        int k = int(min_pt_vox.z());
                        while (k <= max_pt_vox.z()) {
                            Point3 vox(i, j, k);
                            voxel_coordinates.push_back(vox);

                            // check intersections
                            // (a) top - down line
                            Point3 a_start(vox.x() + (0.5 ), vox.y() + (0.5 ), vox.z());
                            Point3 a_end(vox.x() + (0.5 ), vox.y() + (0.5 ), vox.z() + res);
                            Line3 a_line(a_start, a_end);
                            // (b) front - back line
                            Point3 b_start(vox.x() + (0.5 ), vox.y(), vox.z() + 0.5 );
                            Point3 b_end(vox.x() + (0.5 ), vox.y() + res, vox.z() + 0.5 );
                            Line3 b_line(b_start, b_end);
                            // (c) left - right line
                            Point3 c_start(vox.x(), vox.y() + (0.5 ), vox.z() + 0.5 );
                            Point3 c_end(vox.x() + res, vox.y() + (0.5 ), vox.z() + 0.5 );
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
    return vertices_dict;
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
    /*
     * input:   vertex dictionary
     *              {id -> Point3}
     *          nested list of marked surface
     *              [ [ [1, 2, 3], [4, 5, 6], ...], [ [10, 11, 12], [13, 14, 15] ...], ... ]
     *          list of building type which its position refers to the marked surface list
     *              ["Building", "BuildingRoom", ...]  -> option: "Building"/"BuildingRoom"
     */
void write_cityjson(std::string outputfile,
                     std::map<int, Point3>& vertices_dict,
                    const std::vector<std::vector<std::vector<int>>>& markedsurface,
                    const std::vector<std::string>& buildingtype){

    std::vector<Point3> vertices_all;
    nlohmann::json json;
    json["type"] = "CityJSON";
    json["version"] = "1.1";
    json["transform"] = nlohmann::json::object();
    json["transform"]["scale"] = nlohmann::json::array({1.0, 1.0, 1.0});
    json["transform"]["translate"] = nlohmann::json::array({0.0, 0.0, 0.0});
    json["CityObjects"] = nlohmann::json::object();
    json["vertices"] = nlohmann::json::array({});
//    nlohmann::json jvertex;
    std::cout << "vertex size" << vertices_dict.size() << std::endl;
    for (int i = 0; i < vertices_dict.size(); i++) {
        nlohmann::json jvertex;
        Point3 vertex_pt = vertices_dict[i];
        jvertex.push_back(vertex_pt.x());
        jvertex.push_back(vertex_pt.y());
        jvertex.push_back(vertex_pt.z());
        json["vertices"].push_back(jvertex);
        std::cout << i << " ";
    }
    std::cout << "vertex done" << std::endl;

    for (int j = 0; j < buildingtype.size(); j++){
        nlohmann::json surface_array = nlohmann::json::array({});
        std::string indexj = std::to_string(j);
        std::string buildingtypej = buildingtype[j];
        json["CityObjects"][indexj] = nlohmann::json::object();
        json["CityObjects"][indexj]["type"] = buildingtypej;
        json["CityObjects"][indexj]["geometry"] = nlohmann::json::object();
        json["CityObjects"][indexj]["geometry"]["type"] = "Solid";
        json["CityObjects"][indexj]["geometry"]["lod"] = "3.0";
        json["CityObjects"][indexj]["geometry"]["boundaries"] = nlohmann::json::array({});
        nlohmann::json surface = json::array({});
        std::vector<std::vector<std::vector<int>>> ssurface;
        for (const auto& tri : markedsurface[j]) {
            nlohmann::json sid;
            std::vector<std::vector<int>> ssid_d;
            std::vector<int> ssid;
            for (const auto& id : tri) {
                sid.push_back(id);
                ssid.push_back(id);
            }
            surface.push_back(sid);
            ssid_d.push_back(ssid);
            ssurface.push_back(ssid_d);
        }
        surface_array.push_back(ssurface);

        json["CityObjects"][indexj]["geometry"]["boundaries"].push_back(surface_array);
    }
    std::cout << "surface done" << std::endl;

    std::string json_string = json.dump(2);
    std::ofstream out_stream(outputfile);
    out_stream << json_string;
    out_stream.close();
}

//ISOSURFACE-----------------------------------------------------------------------------------------------------------
std::vector<std::vector<Point3>> march_cube(Point3 voxel, VoxelGrid& grid, double res, int id);



std::map<int, Point3> surface_extraction(VoxelGrid& grid, double res, int id, const std::vector<Point3>& allpoints){
    std::vector<Point3> surface;
    std::map<int, Point3> vertex_map;
    int vid = 0;
    for(int n = 1; n < grid.max_x -1; ++n){
        for(int o = 1; o < grid.max_y -1; ++o){
            for(int p = 1; p < grid.max_z-1; ++p){
                std::vector<std::vector<Point3>> triangles = march_cube(Point3(n,o,p), grid, res, id);
                for(auto& triangle: triangles){
                    Point3 v1 = voxcoo_to_modelcoo(triangle[0], allpoints, res);
                    Point3 v2 = voxcoo_to_modelcoo(triangle[1], allpoints, res);
                    Point3 v3 = voxcoo_to_modelcoo(triangle[2], allpoints, res);
                    vertex_map[vid] = v1;
                    vertex_map[vid+1] = v2;
                    vertex_map[vid+2] = v3;
                    vid = vid +3;
                }
            }
        }
    }
    std::cout << "surface size: " << vertex_map.size() << std::endl;
    return vertex_map;
}

std::vector<std::vector<int>> compute_surface(VoxelGrid& grid, double res, int id, const std::vector<Point3>& allpoints) {
    std::map<int, Point3> vertexmap = surface_extraction(grid, res, id, allpoints);
    std::vector<std::vector<int>> surface;
    for (int i = 0; i < vertexmap.size(); i=i+3) {
        std::vector<int> triangle;
        triangle.push_back(i);
        triangle.push_back(i+1);
        triangle.push_back(i+2);
        surface.push_back(triangle);
    }
    return surface;
}


std::vector<pwn> poisson_input_point(std::vector<pwn> pwnthing, Point3 voxel, double res){
    std::vector<pwn> points_surface_poisson;
    Point3 center = Point3(voxel.x() +1/2*res, voxel.y() + 1/2*res, voxel.z()+1/2*res);
    double x = center.x();
    double y = center.y();
    double z = center.z();
    Point3 zero(x - 0.5 * res, y, z);
    vector3 n0 = {-1, 0,0 };
    pwnthing.push_back(pwn(zero,n0));
    Point3 one(x, y - 0.5 * res, z);
    vector3 n1 = {0, -1,0 };
    pwnthing.push_back(pwn(one, n1));
    Point3 two(x + 0.5 * res, y, z);
    vector3 n2 = {1, 0,0 };
    pwnthing.push_back(pwn(two, n2));
    Point3 three(x , y +0.5*res, z);
    vector3 n3 = {0, 1,0 };
    pwnthing.push_back(pwn(three, n3));
    Point3 four(x , y, z - 0.5 * res);
    vector3  n4 = {0, 0, -1 };
    pwnthing.push_back(pwn(four, n4));
    Point3 five(x , y , z + 0.5 * res);
    vector3 n5 = {0, 0, 1};
    pwnthing.push_back(pwn(five,n5));
    return pwnthing;
}




void write_toobj(std::string filename, std::vector<std::vector<std::vector<int>>> allsurface, std::map<int, Point3> vertexdict) {
    std::ofstream output_stream;
    output_stream.open(filename, std::fstream::app);
    int i = 1;
    for (auto & surface : allsurface) {
        for (auto & tri : surface) {
            for (auto & vertex : tri) {
                Point3 vertcoo = vertexdict[vertex];
                output_stream << "v " << vertcoo.x() << " " << vertcoo.y() <<" "<< vertcoo.z() << "\n";
            }
            output_stream << "f " << i << " " << i+1 << " " << i+2 << "\n";
            i = i +3;
        }
    }
}


//MAIN-----------------------------------------------------------------------------------------------------------------
int main(int argc, const char * argv[]) {
    double res = 1.0;

    //file
    const std::string filename = "Duplex_A_20110907_nofurni.obj";

    //Writing to obj
    loadObjFile(filename);

    //create voxelgrid
    VoxelGrid grid = create_voxelgrid(points, res);
    std::cout << "max x: " << grid.max_x << std::endl;
    std::cout << "max y: " << grid.max_y << std::endl;
    std::cout << "max z: " << grid.max_z << std::endl;


    //Voxelisation
    voxelisation(res, grid); // this also returns vertex dictionary

    //Marking
    //marking exterior
    marking(Point3(0, 0, 0), 2, grid); //--> thus this would mark the exterior with id 2
    std::cout << "exterior marked" << std::endl;
    //marking rooms
    int id_rooms = 3;
    int voxels_with_1 = 0;
    int voxels_with_2 = 0;
    int voxels_with_3 = 0;
    for (int i = 0; i < grid.max_x; ++i) { //--> marks rooms with each a different id
        for (int j = 0; j < grid.max_y; ++j){
            for (int k = 0; k < grid.max_z; ++k){
                if (grid(i,j,k) == 0) {
                    marking(Point3(i, j, k), id_rooms, grid);
                    id_rooms++;
                }
                else {
                    continue;
                }

            }
        }
    }
    // calculate voxels
    for (int i = 0; i < grid.max_x; ++i) { //--> marks rooms with each a different id
        for (int j = 0; j < grid.max_y; ++j) {
            for (int k = 0; k < grid.max_z; ++k) {
                if (grid(i, j, k) == 1) {
                    voxels_with_1++;
                }
                if (grid(i, j, k) == 2) {
                    voxels_with_2++;
                }
                if (grid(i, j, k) == 3) {
                    voxels_with_3++;
                }
            }
        }
    }

    std::cout << "num of voxels marked with 1: " << voxels_with_1 << std::endl;
    std::cout << "num of voxels marked with 2: " << voxels_with_2 << std::endl;
    std::cout << "num of voxels marked with 3: " << voxels_with_3 << std::endl;
    std::cout << "number of rooms" << id_rooms - 3 << std::endl;

//    //end of marking rooms
    std::cout << "rooms marked" << std::endl;
    //SURFACE
    //outer envelope
    //gives a vector consisting of triangles with each 3 points
    std::vector<std::vector<int>> surface_outer = compute_surface(grid, res, 2, points);
    // convert point3 to id
    std::vector<std::vector<std::vector<int>>> allsurface;
    std::map<int, Point3> vertexdict = surface_extraction(grid, res, 2, points);
    // create list of building type
    std::vector<std::string> buildingtype = {"Building"};
    allsurface.push_back(surface_outer);

    //rooms
    //i is the room number.
    std::vector<std::vector<std::vector<int>>> room_surfaces;
    for (int i = 3; i <= id_rooms; ++i){
        room_surfaces.push_back(compute_surface(grid, res, i, points));
    }
    std::cout << "surface extracted " << std::endl;

    for (auto const& room_sur : room_surfaces) {
        std::vector<std::vector<int>> room_surfaces_ids;
        // result: room surface 1 = [ [1, 2, 3, ...], [4, 5, 6, 7, ...] ]
        for (auto & room_id : room_sur) {
            room_surfaces_ids.push_back(room_id);
        }
        allsurface.push_back(room_surfaces_ids);
        buildingtype.push_back("BuildingRoom");
    }

    // a whole list of surface ids

    std::cout << "building type list size: " << buildingtype.size() << std::endl;
    std::cout << "surface id list size: " << allsurface.size() << std::endl;

    const std::string outputfile = "output.city.json";
    const std::string outputobj = "output4.obj";
    write_cityjson(outputfile, vertexdict, allsurface, buildingtype);
    std::cout << "cityjson done" << std::endl;
    write_toobj(outputobj, allsurface, vertexdict);
    std::cout << "obj done" << std::endl;

    // trying the poisson surface reconstruction --------------------------------------------------
//    std::vector<pwn> surface_pts_poisson;
//    for (int i = 0; i < grid.max_x; ++i) { //--> marks rooms with each a different id
//        for (int j = 0; j < grid.max_y; ++j){
//            for (int k = 0; k < grid.max_z; ++k){
//                if (grid(i,j,k) == 1) {
//                    surface_pts_poisson = poisson_input_point(surface_pts_poisson, Point3(i,j,k), res);
//
//                    }
//                }
//            }
//        }
//
//        Polyhedron  output_mesh;
//        double average_spacing = CGAL::compute_average_spacing<CGAL::Sequential_tag>
//                (surface_pts_poisson, 6, CGAL::parameters::point_map(CGAL::First_of_pair_property_map<pwn>()));
//        if (CGAL::poisson_surface_reconstruction_delaunay
//                (surface_pts_poisson.begin(), surface_pts_poisson.end(),
//                 CGAL::First_of_pair_property_map<pwn>(),
//                 CGAL::Second_of_pair_property_map<pwn>(),
//                 output_mesh, average_spacing))
//        {
//            std::ofstream out("surface_poisson.off");
//            out << output_mesh;
//        }
//        else
//            return EXIT_FAILURE;
//        return EXIT_SUCCESS;
    // -------------------------------------------------------------------------------------------

        return 0;
}


int triTable[256][16] =
        {{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {0, 1, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {1, 8, 3, 9, 8, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {0, 8, 3, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {9, 2, 10, 0, 2, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {2, 8, 3, 2, 10, 8, 10, 9, 8, -1, -1, -1, -1, -1, -1, -1},
         {3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {0, 11, 2, 8, 11, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {1, 9, 0, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {1, 11, 2, 1, 9, 11, 9, 8, 11, -1, -1, -1, -1, -1, -1, -1},
         {3, 10, 1, 11, 10, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {0, 10, 1, 0, 8, 10, 8, 11, 10, -1, -1, -1, -1, -1, -1, -1},
         {3, 9, 0, 3, 11, 9, 11, 10, 9, -1, -1, -1, -1, -1, -1, -1},
         {9, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {4, 3, 0, 7, 3, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {0, 1, 9, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {4, 1, 9, 4, 7, 1, 7, 3, 1, -1, -1, -1, -1, -1, -1, -1},
         {1, 2, 10, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {3, 4, 7, 3, 0, 4, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1},
         {9, 2, 10, 9, 0, 2, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1},
         {2, 10, 9, 2, 9, 7, 2, 7, 3, 7, 9, 4, -1, -1, -1, -1},
         {8, 4, 7, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {11, 4, 7, 11, 2, 4, 2, 0, 4, -1, -1, -1, -1, -1, -1, -1},
         {9, 0, 1, 8, 4, 7, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1},
         {4, 7, 11, 9, 4, 11, 9, 11, 2, 9, 2, 1, -1, -1, -1, -1},
         {3, 10, 1, 3, 11, 10, 7, 8, 4, -1, -1, -1, -1, -1, -1, -1},
         {1, 11, 10, 1, 4, 11, 1, 0, 4, 7, 11, 4, -1, -1, -1, -1},
         {4, 7, 8, 9, 0, 11, 9, 11, 10, 11, 0, 3, -1, -1, -1, -1},
         {4, 7, 11, 4, 11, 9, 9, 11, 10, -1, -1, -1, -1, -1, -1, -1},
         {9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {9, 5, 4, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {0, 5, 4, 1, 5, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {8, 5, 4, 8, 3, 5, 3, 1, 5, -1, -1, -1, -1, -1, -1, -1},
         {1, 2, 10, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {3, 0, 8, 1, 2, 10, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1},
         {5, 2, 10, 5, 4, 2, 4, 0, 2, -1, -1, -1, -1, -1, -1, -1},
         {2, 10, 5, 3, 2, 5, 3, 5, 4, 3, 4, 8, -1, -1, -1, -1},
         {9, 5, 4, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {0, 11, 2, 0, 8, 11, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1},
         {0, 5, 4, 0, 1, 5, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1},
         {2, 1, 5, 2, 5, 8, 2, 8, 11, 4, 8, 5, -1, -1, -1, -1},
         {10, 3, 11, 10, 1, 3, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1},
         {4, 9, 5, 0, 8, 1, 8, 10, 1, 8, 11, 10, -1, -1, -1, -1},
         {5, 4, 0, 5, 0, 11, 5, 11, 10, 11, 0, 3, -1, -1, -1, -1},
         {5, 4, 8, 5, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1},
         {9, 7, 8, 5, 7, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {9, 3, 0, 9, 5, 3, 5, 7, 3, -1, -1, -1, -1, -1, -1, -1},
         {0, 7, 8, 0, 1, 7, 1, 5, 7, -1, -1, -1, -1, -1, -1, -1},
         {1, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {9, 7, 8, 9, 5, 7, 10, 1, 2, -1, -1, -1, -1, -1, -1, -1},
         {10, 1, 2, 9, 5, 0, 5, 3, 0, 5, 7, 3, -1, -1, -1, -1},
         {8, 0, 2, 8, 2, 5, 8, 5, 7, 10, 5, 2, -1, -1, -1, -1},
         {2, 10, 5, 2, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1},
         {7, 9, 5, 7, 8, 9, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1},
         {9, 5, 7, 9, 7, 2, 9, 2, 0, 2, 7, 11, -1, -1, -1, -1},
         {2, 3, 11, 0, 1, 8, 1, 7, 8, 1, 5, 7, -1, -1, -1, -1},
         {11, 2, 1, 11, 1, 7, 7, 1, 5, -1, -1, -1, -1, -1, -1, -1},
         {9, 5, 8, 8, 5, 7, 10, 1, 3, 10, 3, 11, -1, -1, -1, -1},
         {5, 7, 0, 5, 0, 9, 7, 11, 0, 1, 0, 10, 11, 10, 0, -1},
         {11, 10, 0, 11, 0, 3, 10, 5, 0, 8, 0, 7, 5, 7, 0, -1},
         {11, 10, 5, 7, 11, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {0, 8, 3, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {9, 0, 1, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {1, 8, 3, 1, 9, 8, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1},
         {1, 6, 5, 2, 6, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {1, 6, 5, 1, 2, 6, 3, 0, 8, -1, -1, -1, -1, -1, -1, -1},
         {9, 6, 5, 9, 0, 6, 0, 2, 6, -1, -1, -1, -1, -1, -1, -1},
         {5, 9, 8, 5, 8, 2, 5, 2, 6, 3, 2, 8, -1, -1, -1, -1},
         {2, 3, 11, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {11, 0, 8, 11, 2, 0, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1},
         {0, 1, 9, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1},
         {5, 10, 6, 1, 9, 2, 9, 11, 2, 9, 8, 11, -1, -1, -1, -1},
         {6, 3, 11, 6, 5, 3, 5, 1, 3, -1, -1, -1, -1, -1, -1, -1},
         {0, 8, 11, 0, 11, 5, 0, 5, 1, 5, 11, 6, -1, -1, -1, -1},
         {3, 11, 6, 0, 3, 6, 0, 6, 5, 0, 5, 9, -1, -1, -1, -1},
         {6, 5, 9, 6, 9, 11, 11, 9, 8, -1, -1, -1, -1, -1, -1, -1},
         {5, 10, 6, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {4, 3, 0, 4, 7, 3, 6, 5, 10, -1, -1, -1, -1, -1, -1, -1},
         {1, 9, 0, 5, 10, 6, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1},
         {10, 6, 5, 1, 9, 7, 1, 7, 3, 7, 9, 4, -1, -1, -1, -1},
         {6, 1, 2, 6, 5, 1, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1},
         {1, 2, 5, 5, 2, 6, 3, 0, 4, 3, 4, 7, -1, -1, -1, -1},
         {8, 4, 7, 9, 0, 5, 0, 6, 5, 0, 2, 6, -1, -1, -1, -1},
         {7, 3, 9, 7, 9, 4, 3, 2, 9, 5, 9, 6, 2, 6, 9, -1},
         {3, 11, 2, 7, 8, 4, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1},
         {5, 10, 6, 4, 7, 2, 4, 2, 0, 2, 7, 11, -1, -1, -1, -1},
         {0, 1, 9, 4, 7, 8, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1},
         {9, 2, 1, 9, 11, 2, 9, 4, 11, 7, 11, 4, 5, 10, 6, -1},
         {8, 4, 7, 3, 11, 5, 3, 5, 1, 5, 11, 6, -1, -1, -1, -1},
         {5, 1, 11, 5, 11, 6, 1, 0, 11, 7, 11, 4, 0, 4, 11, -1},
         {0, 5, 9, 0, 6, 5, 0, 3, 6, 11, 6, 3, 8, 4, 7, -1},
         {6, 5, 9, 6, 9, 11, 4, 7, 9, 7, 11, 9, -1, -1, -1, -1},
         {10, 4, 9, 6, 4, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {4, 10, 6, 4, 9, 10, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1},
         {10, 0, 1, 10, 6, 0, 6, 4, 0, -1, -1, -1, -1, -1, -1, -1},
         {8, 3, 1, 8, 1, 6, 8, 6, 4, 6, 1, 10, -1, -1, -1, -1},
         {1, 4, 9, 1, 2, 4, 2, 6, 4, -1, -1, -1, -1, -1, -1, -1},
         {3, 0, 8, 1, 2, 9, 2, 4, 9, 2, 6, 4, -1, -1, -1, -1},
         {0, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {8, 3, 2, 8, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1},
         {10, 4, 9, 10, 6, 4, 11, 2, 3, -1, -1, -1, -1, -1, -1, -1},
         {0, 8, 2, 2, 8, 11, 4, 9, 10, 4, 10, 6, -1, -1, -1, -1},
         {3, 11, 2, 0, 1, 6, 0, 6, 4, 6, 1, 10, -1, -1, -1, -1},
         {6, 4, 1, 6, 1, 10, 4, 8, 1, 2, 1, 11, 8, 11, 1, -1},
         {9, 6, 4, 9, 3, 6, 9, 1, 3, 11, 6, 3, -1, -1, -1, -1},
         {8, 11, 1, 8, 1, 0, 11, 6, 1, 9, 1, 4, 6, 4, 1, -1},
         {3, 11, 6, 3, 6, 0, 0, 6, 4, -1, -1, -1, -1, -1, -1, -1},
         {6, 4, 8, 11, 6, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {7, 10, 6, 7, 8, 10, 8, 9, 10, -1, -1, -1, -1, -1, -1, -1},
         {0, 7, 3, 0, 10, 7, 0, 9, 10, 6, 7, 10, -1, -1, -1, -1},
         {10, 6, 7, 1, 10, 7, 1, 7, 8, 1, 8, 0, -1, -1, -1, -1},
         {10, 6, 7, 10, 7, 1, 1, 7, 3, -1, -1, -1, -1, -1, -1, -1},
         {1, 2, 6, 1, 6, 8, 1, 8, 9, 8, 6, 7, -1, -1, -1, -1},
         {2, 6, 9, 2, 9, 1, 6, 7, 9, 0, 9, 3, 7, 3, 9, -1},
         {7, 8, 0, 7, 0, 6, 6, 0, 2, -1, -1, -1, -1, -1, -1, -1},
         {7, 3, 2, 6, 7, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {2, 3, 11, 10, 6, 8, 10, 8, 9, 8, 6, 7, -1, -1, -1, -1},
         {2, 0, 7, 2, 7, 11, 0, 9, 7, 6, 7, 10, 9, 10, 7, -1},
         {1, 8, 0, 1, 7, 8, 1, 10, 7, 6, 7, 10, 2, 3, 11, -1},
         {11, 2, 1, 11, 1, 7, 10, 6, 1, 6, 7, 1, -1, -1, -1, -1},
         {8, 9, 6, 8, 6, 7, 9, 1, 6, 11, 6, 3, 1, 3, 6, -1},
         {0, 9, 1, 11, 6, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {7, 8, 0, 7, 0, 6, 3, 11, 0, 11, 6, 0, -1, -1, -1, -1},
         {7, 11, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {3, 0, 8, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {0, 1, 9, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {8, 1, 9, 8, 3, 1, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1},
         {10, 1, 2, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {1, 2, 10, 3, 0, 8, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1},
         {2, 9, 0, 2, 10, 9, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1},
         {6, 11, 7, 2, 10, 3, 10, 8, 3, 10, 9, 8, -1, -1, -1, -1},
         {7, 2, 3, 6, 2, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {7, 0, 8, 7, 6, 0, 6, 2, 0, -1, -1, -1, -1, -1, -1, -1},
         {2, 7, 6, 2, 3, 7, 0, 1, 9, -1, -1, -1, -1, -1, -1, -1},
         {1, 6, 2, 1, 8, 6, 1, 9, 8, 8, 7, 6, -1, -1, -1, -1},
         {10, 7, 6, 10, 1, 7, 1, 3, 7, -1, -1, -1, -1, -1, -1, -1},
         {10, 7, 6, 1, 7, 10, 1, 8, 7, 1, 0, 8, -1, -1, -1, -1},
         {0, 3, 7, 0, 7, 10, 0, 10, 9, 6, 10, 7, -1, -1, -1, -1},
         {7, 6, 10, 7, 10, 8, 8, 10, 9, -1, -1, -1, -1, -1, -1, -1},
         {6, 8, 4, 11, 8, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {3, 6, 11, 3, 0, 6, 0, 4, 6, -1, -1, -1, -1, -1, -1, -1},
         {8, 6, 11, 8, 4, 6, 9, 0, 1, -1, -1, -1, -1, -1, -1, -1},
         {9, 4, 6, 9, 6, 3, 9, 3, 1, 11, 3, 6, -1, -1, -1, -1},
         {6, 8, 4, 6, 11, 8, 2, 10, 1, -1, -1, -1, -1, -1, -1, -1},
         {1, 2, 10, 3, 0, 11, 0, 6, 11, 0, 4, 6, -1, -1, -1, -1},
         {4, 11, 8, 4, 6, 11, 0, 2, 9, 2, 10, 9, -1, -1, -1, -1},
         {10, 9, 3, 10, 3, 2, 9, 4, 3, 11, 3, 6, 4, 6, 3, -1},
         {8, 2, 3, 8, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1},
         {0, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {1, 9, 0, 2, 3, 4, 2, 4, 6, 4, 3, 8, -1, -1, -1, -1},
         {1, 9, 4, 1, 4, 2, 2, 4, 6, -1, -1, -1, -1, -1, -1, -1},
         {8, 1, 3, 8, 6, 1, 8, 4, 6, 6, 10, 1, -1, -1, -1, -1},
         {10, 1, 0, 10, 0, 6, 6, 0, 4, -1, -1, -1, -1, -1, -1, -1},
         {4, 6, 3, 4, 3, 8, 6, 10, 3, 0, 3, 9, 10, 9, 3, -1},
         {10, 9, 4, 6, 10, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {4, 9, 5, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {0, 8, 3, 4, 9, 5, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1},
         {5, 0, 1, 5, 4, 0, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1},
         {11, 7, 6, 8, 3, 4, 3, 5, 4, 3, 1, 5, -1, -1, -1, -1},
         {9, 5, 4, 10, 1, 2, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1},
         {6, 11, 7, 1, 2, 10, 0, 8, 3, 4, 9, 5, -1, -1, -1, -1},
         {7, 6, 11, 5, 4, 10, 4, 2, 10, 4, 0, 2, -1, -1, -1, -1},
         {3, 4, 8, 3, 5, 4, 3, 2, 5, 10, 5, 2, 11, 7, 6, -1},
         {7, 2, 3, 7, 6, 2, 5, 4, 9, -1, -1, -1, -1, -1, -1, -1},
         {9, 5, 4, 0, 8, 6, 0, 6, 2, 6, 8, 7, -1, -1, -1, -1},
         {3, 6, 2, 3, 7, 6, 1, 5, 0, 5, 4, 0, -1, -1, -1, -1},
         {6, 2, 8, 6, 8, 7, 2, 1, 8, 4, 8, 5, 1, 5, 8, -1},
         {9, 5, 4, 10, 1, 6, 1, 7, 6, 1, 3, 7, -1, -1, -1, -1},
         {1, 6, 10, 1, 7, 6, 1, 0, 7, 8, 7, 0, 9, 5, 4, -1},
         {4, 0, 10, 4, 10, 5, 0, 3, 10, 6, 10, 7, 3, 7, 10, -1},
         {7, 6, 10, 7, 10, 8, 5, 4, 10, 4, 8, 10, -1, -1, -1, -1},
         {6, 9, 5, 6, 11, 9, 11, 8, 9, -1, -1, -1, -1, -1, -1, -1},
         {3, 6, 11, 0, 6, 3, 0, 5, 6, 0, 9, 5, -1, -1, -1, -1},
         {0, 11, 8, 0, 5, 11, 0, 1, 5, 5, 6, 11, -1, -1, -1, -1},
         {6, 11, 3, 6, 3, 5, 5, 3, 1, -1, -1, -1, -1, -1, -1, -1},
         {1, 2, 10, 9, 5, 11, 9, 11, 8, 11, 5, 6, -1, -1, -1, -1},
         {0, 11, 3, 0, 6, 11, 0, 9, 6, 5, 6, 9, 1, 2, 10, -1},
         {11, 8, 5, 11, 5, 6, 8, 0, 5, 10, 5, 2, 0, 2, 5, -1},
         {6, 11, 3, 6, 3, 5, 2, 10, 3, 10, 5, 3, -1, -1, -1, -1},
         {5, 8, 9, 5, 2, 8, 5, 6, 2, 3, 8, 2, -1, -1, -1, -1},
         {9, 5, 6, 9, 6, 0, 0, 6, 2, -1, -1, -1, -1, -1, -1, -1},
         {1, 5, 8, 1, 8, 0, 5, 6, 8, 3, 8, 2, 6, 2, 8, -1},
         {1, 5, 6, 2, 1, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {1, 3, 6, 1, 6, 10, 3, 8, 6, 5, 6, 9, 8, 9, 6, -1},
         {10, 1, 0, 10, 0, 6, 9, 5, 0, 5, 6, 0, -1, -1, -1, -1},
         {0, 3, 8, 5, 6, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {10, 5, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {11, 5, 10, 7, 5, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {11, 5, 10, 11, 7, 5, 8, 3, 0, -1, -1, -1, -1, -1, -1, -1},
         {5, 11, 7, 5, 10, 11, 1, 9, 0, -1, -1, -1, -1, -1, -1, -1},
         {10, 7, 5, 10, 11, 7, 9, 8, 1, 8, 3, 1, -1, -1, -1, -1},
         {11, 1, 2, 11, 7, 1, 7, 5, 1, -1, -1, -1, -1, -1, -1, -1},
         {0, 8, 3, 1, 2, 7, 1, 7, 5, 7, 2, 11, -1, -1, -1, -1},
         {9, 7, 5, 9, 2, 7, 9, 0, 2, 2, 11, 7, -1, -1, -1, -1},
         {7, 5, 2, 7, 2, 11, 5, 9, 2, 3, 2, 8, 9, 8, 2, -1},
         {2, 5, 10, 2, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1},
         {8, 2, 0, 8, 5, 2, 8, 7, 5, 10, 2, 5, -1, -1, -1, -1},
         {9, 0, 1, 5, 10, 3, 5, 3, 7, 3, 10, 2, -1, -1, -1, -1},
         {9, 8, 2, 9, 2, 1, 8, 7, 2, 10, 2, 5, 7, 5, 2, -1},
         {1, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {0, 8, 7, 0, 7, 1, 1, 7, 5, -1, -1, -1, -1, -1, -1, -1},
         {9, 0, 3, 9, 3, 5, 5, 3, 7, -1, -1, -1, -1, -1, -1, -1},
         {9, 8, 7, 5, 9, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {5, 8, 4, 5, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1},
         {5, 0, 4, 5, 11, 0, 5, 10, 11, 11, 3, 0, -1, -1, -1, -1},
         {0, 1, 9, 8, 4, 10, 8, 10, 11, 10, 4, 5, -1, -1, -1, -1},
         {10, 11, 4, 10, 4, 5, 11, 3, 4, 9, 4, 1, 3, 1, 4, -1},
         {2, 5, 1, 2, 8, 5, 2, 11, 8, 4, 5, 8, -1, -1, -1, -1},
         {0, 4, 11, 0, 11, 3, 4, 5, 11, 2, 11, 1, 5, 1, 11, -1},
         {0, 2, 5, 0, 5, 9, 2, 11, 5, 4, 5, 8, 11, 8, 5, -1},
         {9, 4, 5, 2, 11, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {2, 5, 10, 3, 5, 2, 3, 4, 5, 3, 8, 4, -1, -1, -1, -1},
         {5, 10, 2, 5, 2, 4, 4, 2, 0, -1, -1, -1, -1, -1, -1, -1},
         {3, 10, 2, 3, 5, 10, 3, 8, 5, 4, 5, 8, 0, 1, 9, -1},
         {5, 10, 2, 5, 2, 4, 1, 9, 2, 9, 4, 2, -1, -1, -1, -1},
         {8, 4, 5, 8, 5, 3, 3, 5, 1, -1, -1, -1, -1, -1, -1, -1},
         {0, 4, 5, 1, 0, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {8, 4, 5, 8, 5, 3, 9, 0, 5, 0, 3, 5, -1, -1, -1, -1},
         {9, 4, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {4, 11, 7, 4, 9, 11, 9, 10, 11, -1, -1, -1, -1, -1, -1, -1},
         {0, 8, 3, 4, 9, 7, 9, 11, 7, 9, 10, 11, -1, -1, -1, -1},
         {1, 10, 11, 1, 11, 4, 1, 4, 0, 7, 4, 11, -1, -1, -1, -1},
         {3, 1, 4, 3, 4, 8, 1, 10, 4, 7, 4, 11, 10, 11, 4, -1},
         {4, 11, 7, 9, 11, 4, 9, 2, 11, 9, 1, 2, -1, -1, -1, -1},
         {9, 7, 4, 9, 11, 7, 9, 1, 11, 2, 11, 1, 0, 8, 3, -1},
         {11, 7, 4, 11, 4, 2, 2, 4, 0, -1, -1, -1, -1, -1, -1, -1},
         {11, 7, 4, 11, 4, 2, 8, 3, 4, 3, 2, 4, -1, -1, -1, -1},
         {2, 9, 10, 2, 7, 9, 2, 3, 7, 7, 4, 9, -1, -1, -1, -1},
         {9, 10, 7, 9, 7, 4, 10, 2, 7, 8, 7, 0, 2, 0, 7, -1},
         {3, 7, 10, 3, 10, 2, 7, 4, 10, 1, 10, 0, 4, 0, 10, -1},
         {1, 10, 2, 8, 7, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {4, 9, 1, 4, 1, 7, 7, 1, 3, -1, -1, -1, -1, -1, -1, -1},
         {4, 9, 1, 4, 1, 7, 0, 8, 1, 8, 7, 1, -1, -1, -1, -1},
         {4, 0, 3, 7, 4, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {4, 8, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {9, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {3, 0, 9, 3, 9, 11, 11, 9, 10, -1, -1, -1, -1, -1, -1, -1},
         {0, 1, 10, 0, 10, 8, 8, 10, 11, -1, -1, -1, -1, -1, -1, -1},
         {3, 1, 10, 11, 3, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {1, 2, 11, 1, 11, 9, 9, 11, 8, -1, -1, -1, -1, -1, -1, -1},
         {3, 0, 9, 3, 9, 11, 1, 2, 9, 2, 11, 9, -1, -1, -1, -1},
         {0, 2, 11, 8, 0, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {3, 2, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {2, 3, 8, 2, 8, 10, 10, 8, 9, -1, -1, -1, -1, -1, -1, -1},
         {9, 10, 2, 0, 9, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {2, 3, 8, 2, 8, 10, 0, 1, 8, 1, 10, 8, -1, -1, -1, -1},
         {1, 10, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {1, 3, 8, 9, 1, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {0, 9, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {0, 3, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
         {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}};

std::vector<std::vector<Point3>> march_cube(Point3 voxel, VoxelGrid& grid, double res, int id){
    int cubeindex = 0;
    //vertex 0
    if (grid(voxel.x(), voxel.y()+1, voxel.z()) == id) {
        cubeindex += 1;
    }
    //vertex 1
    if (grid(voxel.x() +1, voxel.y() +1, voxel.z()) == id) {
        cubeindex += 2;
    }
    //vertex 2
    if (grid(voxel.x() +1, voxel.y(), voxel.z()) == id) {
        cubeindex += 4;
    }
    //vertex 3
    if (grid(voxel.x() +1, voxel.y(), voxel.z()) == id) {
        cubeindex += 8;
    }
    //vertex 4
    if (grid(voxel.x(), voxel.y() + 1, voxel.z() + 1) == id) {
        cubeindex += 16;
    }
    //vertex 5
    if (grid(voxel.x() +1, voxel.y() + 1, voxel.z() + 1) == id) {
        cubeindex += 32;
    }
    //vertex 6
    if (grid(voxel.x() +1, voxel.y(), voxel.z() + 1) == id) {
        cubeindex += 64;
    }
    //vertex 7
    if (grid(voxel.x(), voxel.y() + 1, voxel.z() + 1) == id) {
        cubeindex += 128;
    }


    Point3 edge_0 = Point3(voxel.x() + 0.5 , voxel.y() + 1, voxel.z());
    Point3 edge_1 = Point3(voxel.x() + 1, voxel.y() + 0.5 , voxel.z());
    Point3 edge_2 = Point3(voxel.x() + 0.5 , voxel.y(), voxel.z());
    Point3 edge_3 = Point3(voxel.x(), voxel.y() + 0.5 , voxel.z());

    Point3 edge_4 = Point3(voxel.x() + 0.5, voxel.y() + 1, voxel.z()+ 1);
    Point3 edge_5 = Point3(voxel.x() + 1, voxel.y() + 0.5 , voxel.z()+ 1);
    Point3 edge_6 = Point3(voxel.x() + 0.5 , voxel.y(), voxel.z()+ 1);
    Point3 edge_7 = Point3(voxel.x() , voxel.y() + 0.5 , voxel.z()+  1);

    Point3 edge_8 = Point3(voxel.x(), voxel.y() + 1, voxel.z()+ 0.5 );
    Point3 edge_9 = Point3(voxel.x() + 1, voxel.y() + 1, voxel.z()+ 0.5);
    Point3 edge_10 = Point3(voxel.x() + 1, voxel.y(), voxel.z()+ 0.5);
    Point3 edge_11 = Point3(voxel.x() , voxel.y(), voxel.z() + 0.5);

    std::vector<Point3> edge_list = {edge_0, edge_1, edge_2, edge_3, edge_4, edge_5, edge_6, edge_7, edge_8, edge_9, edge_10, edge_11};


    std::vector<std::vector<Point3>> triangles;
    for(int i=0; triTable[cubeindex][i] != -1; i+=3){
        std::vector<Point3> triangle;
        triangle.push_back(edge_list[triTable[cubeindex][i]]);
        triangle.push_back(edge_list[triTable[cubeindex][i + 1]]);
        triangle.push_back(edge_list[triTable[cubeindex][i + 2]]);
        triangles.push_back(triangle);
    }
    return triangles;
}

