/*
 * =====================================================================================
 *
 *       Filename:  main.c
 *
 *    Description:  HBFind 
 *
 *        Version:  1.0
 *        Created:  13/09/21 09:43:14 PM IST
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  PARTHAJIT ROY (PR), roy.parthajit@gmail.com
 *   Organization:  The University of Burdwan
 *
 * =====================================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include "geom3d.h"
#include "biodefs.h"
#include "bioio.h"


struct parameter{
      struct{
	    char file_type;
      }cmd;
      struct{
	    double hbdist;
      }bio;
      struct{
	    char os[10];
      }sys;
};
/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  main
 *  Description:  
 * =====================================================================================
 */
int main ( int argc, char *argv[] )
{
      printf("HBFind Starts\n");
      struct parameter param;
      Point3d CA = vec3d_create(35.097, 6.017, 40.600);
      Point3d CB = vec3d_create(36.041, 5.361, 41.614);
      Point3d CG = vec3d_create(35.868, 5.827, 43.041);
      Point3d CD = vec3d_create(34.759, 5.145, 43.815);
      Point3d C  = vec3d_create(35.451,5.552 , 39.197);
      Point3d N  = vec3d_create(33.701,5.765 , 40.912);
      Plane ABG  = plane_create(CA, CB, CG);
      Plane BGD  = plane_create(CB, CG, CD);
      Plane CAB  = plane_create(C, CA, CB);
      Plane NAB  = plane_create(N, CA, CB);
      double dihed = torsion_angle(N, CA, CB, CG);
      printf("The angle is %lf, %lf\n", dihed, 180.0 * dihed / 3.14159);
      return EXIT_SUCCESS;
}				/* ----------  end of function main  ---------- */
