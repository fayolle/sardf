#include <vector>
#include <iostream>
#include <fstream>
#include <cstdlib>

#include "sardf.h"


const double EPSILON = 1e-7;


struct Point2 {
  Point2() : x(0), y(0) {}
  Point2(double x, double y) : x(x), y(y) {}
  double x;
  double y;
  bool operator==(const Point2& p) const {
    return fabs(x - p.x) <= EPSILON && fabs(y - p.y) <= EPSILON;
  }
};

typedef std::vector<std::vector<Point2>> Grid;

// Make a regular grid of point coordinates.
// In grid[x][y], the first index corresponds to the x coordinate
// being fixed and the second index to the y coordinate being fixed.
// Example:
// grid[x][y0], grid[x][y1], ..., correspond to points on the grid
// with fixed x value and increasing y values
void MakeGrid(const Point2& pmin, const Point2& pmax,
              int number_steps_x, int number_steps_y,
              Grid& grid) {

  grid.clear();

  double delta_x, delta_y;
  delta_x = (pmax.x - pmin.x) / number_steps_x;
  delta_y = (pmax.y - pmin.y) / number_steps_y;

  for (int step_x = 0; step_x < number_steps_x; ++step_x) {
    double x = pmin.x + step_x * delta_x;
    std::vector<Point2> yraw;
    for (int step_y = 0; step_y < number_steps_y; ++step_y) {
      double y = pmin.y + step_y * delta_y;
      yraw.push_back(Point2(x, y));
    }
    grid.push_back(yraw);
  }
}

template<class T>
std::vector<T>
flatten(const std::vector<std::vector<T> >& grid) {
  std::vector<T> flat_grid;
  for (size_t i = 0; i < grid.size(); ++i) {
    for (size_t j = 0; j < grid[i].size(); ++j) {
      flat_grid.push_back(grid[i][j]);
    }
  }
  
  return flat_grid;
}

typedef std::vector<std::vector<double>> ValueGrid;

void create_vtk_file(
                     const std::string& out_filename,
                     const Grid& coordinates,
                     const ValueGrid& values,
                     const std::string& field_title = "Density")
{
  std::ofstream out(out_filename.c_str());

  std::size_t subx = values.size();
  std::size_t suby = values[0].size();
  std::size_t subz = 1;

  // flatten the grids:
  std::vector<Point2> grid = flatten(coordinates);
  std::vector<double> data = flatten(values);

  // header
  out << "# vtk DataFile Version 3.0" << std::endl;
  out << "vtk output" << std::endl;
  out << "ASCII" << std::endl;
  out << "DATASET STRUCTURED_GRID" << std::endl;
  out << "DIMENSIONS " <<
    subx << " " <<
    suby << " " <<
    subz << std::endl;
  out << "POINTS " << subx * suby*subz << " double" << std::endl;

  // structured grid
  std::vector<Point2>::const_iterator it;
  for (it = grid.begin(); it != grid.end(); ++it) {
    Point2 curr = *it;
    out << curr.x << " " << curr.y << " 0.0" << std::endl;
  }
  out << std::endl;

  // data
  // header
  out << std::endl;
  out << "POINT_DATA " << subx * suby*subz << std::endl;
  out << "SCALARS " << field_title.c_str() << " double" << std::endl;
  out << "LOOKUP_TABLE default" << std::endl;

  // data
  std::vector<double>::const_iterator datait;
  for (datait = data.begin(); datait != data.end(); ++datait) {
    out << *datait << std::endl;
  }

  out << std::endl;

  out.close();
}

void sample_sardf(const Grid& p, ValueGrid& f) {
  f.clear();

  for (std::size_t i = 0; i < p.size(); ++i) {
    std::vector<double> temp;
    for (std::size_t j = 0; j < p[i].size(); ++j) {
      Point2 current = p[i][j];
      temp.push_back(uni_sardf_1(current.x, current.y, 0.1, 0.5));
    }
    f.push_back(temp);
  }
}

int main(int argc, char** argv) {
  Grid coordinates;
  Point2 pmin(-5.0,-5.0), pmax(5.0,5.0);

  // Number of steps in each direction for the grid
  int number_steps_x = 128;
  int number_steps_y = 128;
  MakeGrid(pmin, pmax,
           number_steps_x, number_steps_y, coordinates);

  std::string output_base_name = "test.vtk";


  // Save psi1, normalized psi1
  ValueGrid f;
  sample_sardf(coordinates, f);
  create_vtk_file(output_base_name, coordinates, f, "Sardf");

  return EXIT_SUCCESS;
}
