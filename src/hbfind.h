/*
 * =====================================================================================
 *
 *       Filename:  hbfind.h
 *
 *    Description: hydrogen addition. 
 *
 *        Version:  1.0
 *        Created:  13/09/21 09:46:15 PM IST
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  PARTHAJIT ROY (PR), roy.parthajit@gmail.com
 *   Organization:  The University of Burdwan
 *
 * =====================================================================================
 */

 
#ifndef  __hbfind_H__
#define  __hbfind_H__


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
 
point3d find_location(Vector3d pcursor_vector, Vector3d pcursor_normal, Point3d reference_point, double bond_len, double bond_angle, double torsion_angle){
      Vector3d bond_angle_adjusted_vec = vec3d_polar_rotation(precursor_normal, pcursor_vector, PI - bond_angle);
      Vector3d final_position_vec = vec3d_polar_rotation(pcursor_vector, bond_angle_adjusted_vec, torsion_angle);
      Vector3d bondlen_adjusted_vec = vec3d_scal_mult(final_position_vec, bond_len);
      Point3d  scaled_location = (Point3d) vec3d_add(reference_point, bondlen_adjusted_vec);
}

Point3d find_hydro_location(Point3d pcur1, Point3d pcur2, Point3d pcur3, double bond_len, double bond_angle, double torsion_angle){
      Plane pcur_plane = plane_create(pcur1, pcur2, pcur3);
      Vector3d pcur_unit_vector = vec3d_unit(vec3d_sub(pcur3, pcur2));
      find_location(pcur_unit_vector, pcur_plane.normal, pcur3, bond_len, bond_angle, torsion_angle);
}


#endif   /* ----- #ifndef __hbfind_H__  ----- */
