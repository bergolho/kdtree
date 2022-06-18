/*
-----------------------------------------------------------------------------------------------------------------------------
  Program that reads a cloud of points from a surface given in Legacy VTK format and builds a Kdtree using the library
  Author: Lucas Berg @bergolho
-----------------------------------------------------------------------------------------------------------------------------
*/
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include <vector>
#include "kdtree/kdtree.h"

/* returns the distance squared between two dims-dimensional double arrays */
static double dist_sq( double *a1, double *a2, int dims );

static std::vector<std::vector<double>> read_cloud_points_from_vtk (const char filename[]);

static void print_cloud_points (std::vector<std::vector<double>> in);

static void write_points_in_vtk (const char filename[], std::vector<std::vector<double>> in);

int main(int argc, char *argv[]) {
  int i, num_pts;
  void *ptree;
  char *data, *pch;
  struct kdres *presults;
  double pos[3], dist;
  double pt[3] = { 19234.4, 19886, 15900.9 };
  double radius = 50000.0;

  if (argc-1 != 1) {
    printf("----------------------------------------------------------------------------------------------\n");
    printf("Usage:> %s <input_file>\n",argv[0]);
    printf("----------------------------------------------------------------------------------------------\n");
    printf("<input_file> = Input filename with the surface cloud of points in Legacy VTK format\n");
    printf("----------------------------------------------------------------------------------------------\n");
    exit(EXIT_FAILURE);
  }
  char *filename = "inputs/surface_cloud_points.vtk";
  std::vector<std::vector<double>> points = read_cloud_points_from_vtk(filename);
  //print_cloud_points(points);

  /* create a k-d tree for 3-dimensional points */
  ptree = kd_create( 3 );

  /* add some random nodes to the tree (assert nodes are successfully inserted) */
  for(u_int32_t i = 0; i < points.size(); i++ ) {
    assert( 0 == kd_insert3( (kdtree*)ptree, points[i][0], points[i][1], points[i][2], NULL ) );
  }

  /* find points closest to the target and within distance radius */
  //presults = kd_nearest_range( (kdtree*)ptree, pt, radius );

  /* find N points closest to the target */
  presults = kd_nearest_n( (kdtree*)ptree, pt, 40 );

  /* print out all the points found in results */
  printf( "found %d results:\n", kd_res_size(presults) );

  /* Store the nearest points in a vector */
  std::vector<std::vector<double>> nearest_points;

  while( !kd_res_end( presults ) ) {
    /* get the data and position of the current result item */
    pch = (char*)kd_res_item( presults, pos );

    /* compute the distance of the current result from the pt */
    dist = sqrt( dist_sq( pt, pos, 3 ) );

    /* print out the retrieved data */
    printf( "node at (%.3f, %.3f, %.3f) is %.3f away\n", pos[0], pos[1], pos[2], dist );

    /* store the retrieved point in an array */
    nearest_points.push_back({pos[0], pos[1], pos[2]});

    /* go to the next entry */
    kd_res_next( presults );
  }

  /* Write the nearest points in VTK format */
  write_points_in_vtk("outputs/nearest_points.vtk",nearest_points);

  /* free our tree, results set, and other allocated memory */
  kd_res_free( presults );
  kd_free( (kdtree*)ptree );

  return 0;
}

static double dist_sq( double *a1, double *a2, int dims ) {
  double dist_sq = 0, diff;
  while( --dims >= 0 ) {
    diff = (a1[dims] - a2[dims]);
    dist_sq += diff*diff;
  }
  return dist_sq;
}

static std::vector<std::vector<double>> read_cloud_points_from_vtk (const char filename[]) {
  FILE *file = fopen(filename,"r");
  if (!file) {
    fprintf(stderr,"[-] ERROR! Cannot read file '%s'!\n",filename);
    exit(EXIT_FAILURE);
  }

  char str[200];
  while (fscanf(file,"%s",str) != EOF)
    if (strcmp(str,"POINTS") == 0) break;
  
  u_int32_t np;
  fscanf(file,"%u %s",&np,str);
  std::vector<std::vector<double>> out;
  for (u_int32_t i = 0; i < np; i++) {
    double pos[3];
    fscanf(file,"%lf %lf %lf",&pos[0],&pos[1],&pos[2]);
    out.push_back({pos[0],pos[1],pos[2]});
  }
  
  fclose(file);
  return out;
}

static void print_cloud_points (std::vector<std::vector<double>> in) {
  u_int32_t np = in.size();
  for (u_int32_t i = 0; i < np; i++) {
    printf("Point %u = (%g, %g, %g)\n",i,in[i][0],in[i][1],in[i][2]);
  }
}

static void write_points_in_vtk (const char filename[], std::vector<std::vector<double>> in) {
  u_int32_t np = in.size();
  FILE *file = fopen(filename,"w+");
  fprintf(file,"# vtk DataFile Version 4.1\n");
  fprintf(file,"vtk output\n");
  fprintf(file,"ASCII\n");
  fprintf(file,"DATASET POLYDATA\n");
  fprintf(file,"POINTS %u float\n",np);
  for (u_int32_t i = 0; i < np; i++) {
    fprintf(file,"%g %g %g\n",in[i][0],in[i][1],in[i][2]);
  }
  fprintf(file,"VERTICES %u %u\n",np,np*2);
  for (u_int32_t i = 0; i < np; i++) {
    fprintf(file,"1 %u\n",i);
  }

  fclose(file);
}