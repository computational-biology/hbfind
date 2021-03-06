/*
 * =====================================================================================
 *
 *       Filename:  hbfind.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  16/09/21 08:05:44 PM IST
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  PARTHAJIT ROY (PR), roy.parthajit@gmail.com
 *                  ANWESHA DUTTA, Masters Student, Dept of Computer Science,
 *                  e-mail: anweshadatta808@gmail.com 
 *   Organization:  The University of Burdwan
 *
 * =====================================================================================
 */

#include "hbfind.h"


static Point3d find_location(Vector3d axis, Vector3d ref_normal, Point3d reference_point, double bond_len, double bond_angle, double torsion_angle)
{
      Vector3d bang_add = vec3d_polar_rotation(ref_normal, axis, PI - bond_angle);
      Vector3d tors_add = vec3d_polar_rotation(axis, bang_add, torsion_angle);
      Vector3d blen_add = vec3d_scal_mult(tors_add, bond_len);
      Point3d point_add = (Point3d)vec3d_add(reference_point, blen_add);
      return point_add;
}

static Point3d fix_atom(Point3d pcur1, Point3d pcur2, Point3d pcur3, double bond_len, double bond_angle, double torsion_angle)
{
      if(torsion_angle < UNDEFINED_VAL + 1.0 || 
		  Point3d_is_uninit(&pcur1) == TRUE ||
		  Point3d_is_uninit(&pcur2) == TRUE ||
		  Point3d_is_uninit(&pcur3) == TRUE
		  ){
	    Point3d point;
	    point3d_uninitialize(&point);
	    return point;
      }
      Plane pcur_plane = plane_create(pcur1, pcur2, pcur3);
      Vector3d pcur_unit_vector = vec3d_unit(vec3d_sub(pcur3, pcur2));
      Point3d hloc = find_location(pcur_unit_vector, pcur_plane.unit_normal, pcur3, bond_len, bond_angle, torsion_angle);
      return hloc;
}

//static void update_atom(struct atom* atom, int atomid, char* loc_name, char* symb, Point3d* point){
//      atom->id = atomid;
//      strcpy(atom->loc, loc_name);
//      strcpy(atom->symbol, symb);
//      atom->center = *point;
//}

void cys_addh(struct residue *res)
{

      // SG incomplete
      Point3d N;
      Point3d CA;
      Point3d CB;
      Point3d SG;
      point3d_uninitialize(&N);
      point3d_uninitialize(&CA);
      point3d_uninitialize(&CB);
      point3d_uninitialize(&SG);
      int Nii =-1;
      int CAii=-1;
      int CBii=-1;
      int SGii=-1;
      int count = 0;
      for (int i = 0; i < res->size; ++i)
      {
	    if (strcmp(res->atoms[i].loc, "N") == 0)
	    {
		  N = res->atoms[i].center;

		  Nii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "CA") == 0)
	    {
		  CA = res->atoms[i].center;

		  CAii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "CB") == 0)
	    {
		  CB = res->atoms[i].center;

		  CBii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "SG") == 0)
	    {
		  SG = res->atoms[i].center;

		  SGii = i;
		  count++;
	    }
      }
      if (count != 4)
      { /* Exception Handling */
	    fprintf(stderr, "Error in function %s()... All required atoms not found %d\n", __func__, count);
	     
      }
      double adjust = torad(120.0);
      double hbangle = torad(106.97);
      double hbdist = 0.97;
      double tors_ang = torsion_angle(N, CA, CB, SG);


      Point3d HB1 = fix_atom(N, CA, CB, hbdist, hbangle, tors_ang + adjust);
      residue_addh(res, CBii, HB1, "HB1");
      Point3d HB2 = fix_atom(N, CA, CB, hbdist, hbangle, tors_ang - adjust);
      residue_addh(res, CBii, HB2, "HB2");
      Point3d HA = fix_atom(SG, CB, CA, 0.97, torad(90.0), tors_ang + adjust);
      residue_addh(res, CAii, HA, "HA");

      Point3d HG = fix_atom(CA, CB, SG, 1.20, torad(90.0), torad(180.0));

      residue_addh(res, SGii, HG, "HG");

}
void ala_addh(struct residue *res)
{
      Point3d N;
      Point3d CA;
      Point3d CB;
      Point3d C;
      Point3d O;
      point3d_uninitialize(& N);
      point3d_uninitialize(&CA);
      point3d_uninitialize(&CB);
      point3d_uninitialize(&C);
      point3d_uninitialize(&O);
      int Nii =-1;
      int CAii=-1;
      int CBii=-1;
      int Cii =-1;
      int Oii =-1;
      int count = 0;
      for (int i = 0; i < res->size; ++i)
      {
	    if (strcmp(res->atoms[i].loc, "N") == 0)
	    {
		  N = res->atoms[i].center;

		  Nii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "CA") == 0)
	    {
		  CA = res->atoms[i].center;

		  CAii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "CB") == 0)
	    {
		  CB = res->atoms[i].center;

		  CBii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "C") == 0)
	    {
		  C = res->atoms[i].center;

		  Cii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "O") == 0)
	    {
		  O= res->atoms[i].center;

		  Oii = i;
		  count++;
	    }
      }
      if (count != 5)
      { /* Exception Handling */
	    fprintf(stderr, "Error in function %s()... All required atoms not found %d\n", __func__, count);
	     
      }
      double adjust = torad(120.0);
      double hbangle = torad(106.97);
      double hbdist = 0.97;
      double tors_ang = torsion_angle(N, CA, C, O);


      Point3d HA = fix_atom(O,C,CA, hbdist, hbangle, tors_ang+torad(240.0) );
      residue_addh(res, CAii, HA, "HA");
      Point3d H = fix_atom(C,CA,N, hbdist, hbangle, torad(118.0)+torad(90.0));
      residue_addh(res, Nii, H, "H");


      Point3d HB1 = fix_atom(C,CA,CB, 0.97, hbangle, torad(120.0));
      residue_addh(res, CBii, HB1, "HB1");
      Point3d HB2 = fix_atom(C,CA,CB, 0.97, hbangle, torad(0.0));
      residue_addh(res, CBii, HB2, "HB2");
      Point3d HB3 = fix_atom(C,CA,CB, 0.97, hbangle, torad(-120.0));
      residue_addh(res, CBii, HB3, "HB3");

}

void asp_addh(struct residue *res)
{
      Point3d N;
      Point3d CA;
      Point3d CB;
      Point3d CG;
      point3d_uninitialize(&N);
      point3d_uninitialize(&CA);
      point3d_uninitialize(&CB);
      point3d_uninitialize(&CG);
      int Nii  = -1;
      int CAii = -1;
      int CBii = -1;
      int CGii = -1;
      int count = 0;
      for (int i = 0; i < res->size; ++i)
      {
	    if (strcmp(res->atoms[i].loc, "N") == 0)
	    {
		  N = res->atoms[i].center;

		  Nii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "CA") == 0)
	    {
		  CA = res->atoms[i].center;

		  CAii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "CB") == 0)
	    {
		  CB = res->atoms[i].center;

		  CBii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "CG") == 0)
	    {
		  CG = res->atoms[i].center;

		  CGii = i;
		  count++;
	    }
      }
      if (count != 4)
      { /* Exception Handling */
	    fprintf(stderr, "Error in function %s()... All required atoms not found %d\n", __func__, count);
	     
      }
      double tors_ang = torsion_angle(N, CA, CB, CG);

      double adjust = torad(120.0);
      double hbangle = torad(106.97);
      double hbdist = 0.97;
      Point3d HB1 = fix_atom(N, CA, CB, hbdist, hbangle, tors_ang + adjust);
      residue_addh(res, CBii, HB1, "HB1");
      Point3d HB2 = fix_atom(N, CA, CB, hbdist, hbangle, tors_ang - adjust);
      residue_addh(res, CBii, HB2, "HB2");

      double tors_ang1 = torsion_angle(N, CA, CB, CG);
      Point3d HA = fix_atom(CG, CB, CA, hbdist, hbangle, tors_ang1 + adjust);
      residue_addh(res, CAii, HA, "HA");
      Point3d HN = fix_atom(CB, CA, N, 0.86, torad(111.37), torad(180.0));

      residue_addh(res, Nii, HN, "H");
}

void glu_addh(struct residue *res)
{
      Point3d CA;
      Point3d CB;
      Point3d CG;
      Point3d CD;
      Point3d N;
      Point3d C;
      point3d_uninitialize(&CA);
      point3d_uninitialize(&CB);
      point3d_uninitialize(&CG);
      point3d_uninitialize(&CD);
      point3d_uninitialize(&N);
      point3d_uninitialize(&C);
      int CAii = -1;
      int CBii = -1;
      int CGii = -1;
      int CDii = -1;
      int Nii =  -1;
      int Cii =  -1;
      // Point3d O;

      int count = 0;
      for (int i = 0; i < res->size; ++i)
      {
	    if (strcmp(res->atoms[i].loc, "CA") == 0)
	    {
		  CA = res->atoms[i].center;

		  CAii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "CB") == 0)
	    {
		  CB = res->atoms[i].center;

		  CBii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "CG") == 0)
	    {
		  CG = res->atoms[i].center;

		  CGii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "CD") == 0)
	    {
		  CD = res->atoms[i].center;

		  CDii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "N") == 0)
	    {
		  N= res->atoms[i].center;

		  Nii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "C") == 0)
	    {
		  C= res->atoms[i].center;

		  Cii = i;
		  count++;
	    }
      }

      if (count != 6)
      { /* Exception Handling */
	    fprintf(stderr, "Error in function %s()... All required atoms not found %d\n", __func__, count);
	     
      }
      double tors_ang = torsion_angle(CA, CB, CG, CD);

      double adjust = torad(120.0);
      double hbangle = torad(106.97);
      double hbdist = 0.97;
      Point3d HB1 = fix_atom(CD, CG, CB, hbdist, hbangle, tors_ang + adjust);
      residue_addh(res, CBii, HB1, "HB1");
      Point3d HB2 = fix_atom(CD, CG, CB, hbdist, hbangle, tors_ang - adjust);
      residue_addh(res, CBii, HB2, "HB2");
      // double tors_ang1 = torsion_angle( CA, CB, CG,CD);
      Point3d HE1 = fix_atom(CA, CB, CG, hbdist, hbangle, tors_ang + adjust);
      residue_addh(res, CGii, HE1, "HE1");
      Point3d HE2 = fix_atom(CA, CB, CG, hbdist, hbangle, tors_ang - adjust);
      residue_addh(res, CGii, HE2, "HE2");
      tors_ang = torsion_angle(C, CA, CB, CG);
      Point3d HA = fix_atom(CG, CB, CA, hbdist, torad(280), tors_ang + adjust);
      residue_addh(res, CAii, HA, "HA");
      Point3d H= fix_atom(CB, CA, N, .86, torad(183.33), torad(90));

      residue_addh(res, Nii, H, "H");
}

void phe_addh(struct residue* res)
{
      Point3d C;
      Point3d N;
      Point3d CA;
      Point3d CB;
      Point3d CG;
      Point3d CD1;
      Point3d CD2;
      Point3d CE1;
      Point3d CE2;
      Point3d CZ;
      point3d_uninitialize(&C);
      point3d_uninitialize(&N);
      point3d_uninitialize(&CA);
      point3d_uninitialize(&CB);
      point3d_uninitialize(&CG);
      point3d_uninitialize(&CD1);
      point3d_uninitialize(&CD2);
      point3d_uninitialize(&CE1);
      point3d_uninitialize(&CE2);
      point3d_uninitialize(&CZ);

      int  Cii  = -1;
      int  Nii  = -1;
      int CAii  = -1;
      int CBii  = -1;
      int CGii  = -1;
      int CD1ii = -1;
      int CD2ii = -1;
      int CE1ii = -1;
      int CE2ii = -1;
      int CZii  = -1;
      int count = 0;
      for(int i=0; i<res->size; ++i){
	    if(strcmp(res->atoms[i].loc, "C") == 0){
		  C = res->atoms[i].center;

		  Cii = i;
		  count++;
	    } 
	    /*  else if(strcmp(res->atoms[i].loc, "O") == 0){
		O= res->atoms[i].center;
		count++;
		} */
	    else if(strcmp(res->atoms[i].loc, "N") == 0){
		  N = res->atoms[i].center;

		  Nii = i;
		  count++;
	    }
	    else if(strcmp(res->atoms[i].loc, "CA") == 0){
		  CA= res->atoms[i].center;

		  CAii = i;
		  count++;
	    }else if(strcmp(res->atoms[i].loc, "CB") == 0){
		  CB = res->atoms[i].center;
		  count++;
	    }else if(strcmp(res->atoms[i].loc, "CG") == 0){
		  CG = res->atoms[i].center;

		  CGii = i;
		  count++;
	    }else if(strcmp(res->atoms[i].loc, "CD1") == 0){
		  CD1= res->atoms[i].center;

		  CD1ii = i;
		  count++;
	    }
	    else if(strcmp(res->atoms[i].loc, "CD2") == 0){
		  CD2= res->atoms[i].center;

		  CD2ii = i;
		  count++;
	    }
	    else if(strcmp(res->atoms[i].loc, "CE1") == 0){
		  CE1= res->atoms[i].center;

		  CE1ii = i;
		  count++;
	    }
	    else if(strcmp(res->atoms[i].loc, "CE2") == 0){
		  CE2= res->atoms[i].center;

		  CE2ii = i;
		  count++;
	    }
	    else if(strcmp(res->atoms[i].loc, "CZ") == 0){
		  CZ= res->atoms[i].center;

		  CZii = i;
		  count++;
	    }
      }
      if( count != 10){    /* Exception Handling */
	    fprintf(stderr, "Error in function %s()... All required atoms not found %d\n", __func__, count);
	     
      }
      double tors_ang = torsion_angle(C ,CA, CB, CG);

      double adjust = torad(120.0);
      double hbangle = torad(106.97);
      double hbdist  = 0.97;
      Point3d HB2 = fix_atom(C,CA,CB,hbdist, hbangle, tors_ang + adjust );
      residue_addh(res, CBii, HB2, "HB2");
      Point3d HB3 = fix_atom(C, CA, CB, hbdist, hbangle, tors_ang - adjust );
      residue_addh(res, CBii, HB3, "HB3");
      Point3d HA = fix_atom(CG,CB,CA,hbdist,hbangle,tors_ang -adjust);
      residue_addh(res, CAii, HA, "HA");

      tors_ang = torsion_angle(CG,CD2,CE2,CZ);
      Point3d HD2 = fix_atom(CZ,CE2,CD2,hbdist, torad(-120.0), tors_ang);
      residue_addh(res, CD2ii, HD2, "HD2");
      Point3d HE2 = fix_atom(CG,CD2,CE2, hbdist, torad(-120), tors_ang);
      residue_addh(res, CE2ii, HE2, "HE2");

      tors_ang= torsion_angle( CE2, CZ, CE1,CD1);
      Point3d HZ = fix_atom(CD1,CE1,CZ ,hbdist, torad(-120.0), tors_ang);
      residue_addh(res, CZii, HZ, "HZ");
      Point3d HE1= fix_atom(CE2,CZ,CE1 ,hbdist, torad(-120.0), tors_ang);
      residue_addh(res, CE1ii, HE1, "HE1");

      tors_ang= torsion_angle(CZ,CE1,CD1,CG);
      Point3d HD1= fix_atom(CZ,CE1 ,CD1,hbdist, torad(-120.0), tors_ang);
      residue_addh(res, CD1ii, HD1, "HD1");

      tors_ang= torsion_angle(CG, CB, CA, C);
      Point3d H= fix_atom(CB, CA, N, .86, hbangle, tors_ang + adjust);
      residue_addh(res, Nii, H, "H");

}

void gly_addh(struct residue *res)
{
      Point3d O;
      Point3d C;
      Point3d CA;
      Point3d N;
      point3d_uninitialize(&O);
      point3d_uninitialize(&C);
      point3d_uninitialize(&CA);
      point3d_uninitialize(&N);

      int Oii  = -1;
      int Cii  = -1;
      int CAii = -1;
      int  Nii  = -1;

      int count = 0;
      for (int i = 0; i < res->size; ++i)
      {
	    if (strcmp(res->atoms[i].loc, "O") == 0)
	    {
		  O = res->atoms[i].center;

		  Oii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "C") == 0)
	    {
		  C = res->atoms[i].center;

		  Cii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "CA") == 0)
	    {
		  CA = res->atoms[i].center;

		  CAii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "N") == 0)
	    {
		  N = res->atoms[i].center;

		  Nii = i;
		  count++;
	    }
      }
      if (count != 4)
      { /* Exception Handling */
	    fprintf(stderr, "Error in function %s()... All required atoms not found %d\n", __func__, count);
	     
      }
      double tors_ang = torsion_angle(O, C, CA, N);

      double adjust = torad(120.0);
      double hbangle = torad(106.97);
      double hbdist = 0.97;
      Point3d HB1 = fix_atom(O, C, CA, hbdist, hbangle, tors_ang + adjust);
      residue_addh(res, CAii, HB1, "HB1");
      Point3d HB2 = fix_atom(O, C, CA, hbdist, hbangle, tors_ang - adjust);
      residue_addh(res, CAii, HB2, "HB2");
      Point3d H = fix_atom(C, CA, N, 0.86, hbangle, -adjust);

      residue_addh(res, Nii, H, "H");
}
void his_addh(struct residue* res)
{
      Point3d O;
      Point3d C;
      Point3d CA;
      Point3d N;
      Point3d CB;
      Point3d CG;
      Point3d ND1;
      Point3d CE1;
      Point3d NE2;
      Point3d CD2;

      point3d_uninitialize(&O);
      point3d_uninitialize(&C);
      point3d_uninitialize(&CA);
      point3d_uninitialize(&N);
      point3d_uninitialize(&CB);
      point3d_uninitialize(&CG);
      point3d_uninitialize(&ND1);
      point3d_uninitialize(&CE1);
      point3d_uninitialize(&NE2);
      point3d_uninitialize(&CD2);

      int  Oii  = -1;
      int  Cii  = -1;
      int CAii  = -1;
      int  Nii  = -1;
      int CBii  = -1;
      int CGii  = -1;
      int ND1ii = -1;
      int CE1ii = -1;
      int NE2ii = -1;
      int CD2ii = -1;


      int count = 0;
      for(int i=0; i<res->size; ++i){
	    if(strcmp(res->atoms[i].loc, "O") == 0){
		  O = res->atoms[i].center;

		  Oii = i;

		  count++;
	    }else if(strcmp(res->atoms[i].loc, "C") == 0){
		  C = res->atoms[i].center;

		  Cii = i;
		  count++;
	    }else if(strcmp(res->atoms[i].loc, "CA") == 0){
		  CA = res->atoms[i].center;

		  CAii = i;
		  count++;
	    }else if(strcmp(res->atoms[i].loc, "N") == 0){
		  N = res->atoms[i].center;

		  Nii = i;
		  count++;
	    }
	    else if(strcmp(res->atoms[i].loc, "CB") == 0){
		  CB = res->atoms[i].center;

		  CBii = i;
		  count++;
	    }
	    else if(strcmp(res->atoms[i].loc, "CG") == 0){
		  CG= res->atoms[i].center;

		  CGii = i;
		  count++;
	    }
	    else if(strcmp(res->atoms[i].loc, "ND1") == 0){
		  ND1 = res->atoms[i].center;

		  ND1ii = i;
		  count++;
	    }
	    else if(strcmp(res->atoms[i].loc, "CE1") == 0){
		  CE1= res->atoms[i].center;

		  CE1ii = i;
		  count++;
	    }
	    else if(strcmp(res->atoms[i].loc, "NE2") == 0){
		  NE2 = res->atoms[i].center;

		  NE2ii = i;
		  count++;
	    }
	    else if(strcmp(res->atoms[i].loc, "CD2") == 0){
		  CD2 = res->atoms[i].center;

		  CD2ii = i;
		  count++;
	    }

      }
      if( count != 10 ){ /* Exception Handling */
	    fprintf(stderr, "Error in function %s()... All required atoms not found %d\n", __func__, count);
	     
      }
      double tors_ang = torsion_angle(CG, CB, CA, C);

      double adjust = torad(120.0);
      double hbangle = torad(106.97);
      double hbdist  = 0.97;
      Point3d HA= fix_atom(O,C ,CA, hbdist, hbangle, tors_ang-adjust);
      residue_addh(res, CAii, HA, "HA");

      Point3d H  = fix_atom(C,CA,N ,0.86, hbangle, adjust);

      residue_addh(res, Nii, H, "H");

      double tors_ang1 = torsion_angle(C,CA,CB,CG);
      Point3d HB1 = fix_atom(C ,CA,CB, hbdist, hbangle, tors_ang1 + adjust );
      residue_addh(res, CBii, HB1, "HB1");
      Point3d HB2 = fix_atom(C, CA,CB, hbdist, hbangle, tors_ang1 - adjust );
      residue_addh(res, CBii, HB2, "HB2");

      tors_ang = torsion_angle(CG,ND1,CE1,NE2);
      Point3d HE1= fix_atom(CG,ND1,CE1, hbdist, torad(-126.0), tors_ang);
      residue_addh(res, CE1ii, HE1, "HE1");

      tors_ang = torsion_angle(CB,CG,CD2,NE2);
      Point3d HD2= fix_atom(CB,CG,CD2, hbdist, torad(-126.0), tors_ang);
      residue_addh(res, CD2ii, HD2, "HD2");
}
void ile_addh(struct residue *res)
{
      Point3d N;
      Point3d CA;
      Point3d CB;
      Point3d CG1;
      Point3d CG2;
      Point3d CD1;

      point3d_uninitialize(&N);
      point3d_uninitialize(&CA);
      point3d_uninitialize(&CB);
      point3d_uninitialize(&CG1);
      point3d_uninitialize(&CG2);
      point3d_uninitialize(&CD1);

      int  Nii  = -1;
      int CAii  = -1;
      int CBii  = -1;
      int CG1ii = -1;
      int CG2ii = -1;
      int CD1ii = -1;

      int count = 0;
      for (int i = 0; i < res->size; ++i)
      {
	    if (strcmp(res->atoms[i].loc, "N") == 0)
	    {
		  N = res->atoms[i].center;

		  Nii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "CA") == 0)
	    {
		  CA = res->atoms[i].center;

		  CAii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "CB") == 0)
	    {
		  CB = res->atoms[i].center;

		  CBii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "CG1") == 0)
	    {
		  CG1 = res->atoms[i].center;

		  CG1ii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "CG2") == 0)
	    {
		  CG2 = res->atoms[i].center;

		  CG2ii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "CD1") == 0)
	    {
		  CD1 = res->atoms[i].center;

		  CD1ii = i;
		  count++;
	    }
      }
      if (count != 6)
      { /* Exception Handling */
	    fprintf(stderr, "Error in function %s()... All required atoms not found %d\n", __func__, count);
	     
      }
      double tors_ang = torsion_angle(CG2,CB,CA,N);

      double adjust = torad(120.0);
      double hbangle = torad(106.97);
      double hbdist = 0.97;
      Point3d HA = fix_atom(CG2, CB, CA, 0.97, hbangle, tors_ang + adjust);
      residue_addh(res, CAii, HA, "HA");
      Point3d HB = fix_atom(N, CA, CB, hbdist, hbangle, tors_ang - adjust);
      residue_addh(res, CBii, HB, "HB");
      Point3d H22 = fix_atom(CA, CB, CG2, hbdist, hbangle, -adjust);
      residue_addh(res, CG2ii, H22, "H22");
      Point3d H21 = fix_atom(CA, CB, CG2, hbdist, hbangle, torad(0.0));
      residue_addh(res, CG2ii, H21, "H21");
      Point3d H23= fix_atom(CA, CB, CG2, hbdist, hbangle, adjust);
      residue_addh(res, CG2ii, H23, "H23");
      Point3d H = fix_atom(CB, CA, N, 0.86, hbangle, tors_ang + adjust);
      residue_addh(res, Nii, H, "H");
      tors_ang = torsion_angle(CA,CB,CG1,CD1);
      Point3d H11 = fix_atom(CA,CB,CG1, hbdist, hbangle, tors_ang + adjust);
      residue_addh(res, CG1ii, H11, "H11");
      Point3d H13 = fix_atom(CA,CB,CG1, hbdist, hbangle, tors_ang - adjust);
      residue_addh(res, CG1ii, H13, "H13");
      Point3d HD1 = fix_atom(CB, CG1, CD1,hbdist, hbangle, adjust);
      residue_addh(res, CD1ii, HD1, "HD1");
      Point3d HD2 = fix_atom( CB, CG1,CD1,hbdist, hbangle, torad(0.0));
      residue_addh(res, CD1ii, HD2, "HD2");
      Point3d HD3= fix_atom( CB, CG1,CD1, hbdist, hbangle, -adjust);
      residue_addh(res, CD1ii, HD3, "HD3");


}
void asn_addh(struct residue *res)
{
      Point3d N;
      Point3d CA;
      Point3d CB;
      Point3d CG;
      Point3d C;
      Point3d ND2;

      point3d_uninitialize(&N);
      point3d_uninitialize(&CA);
      point3d_uninitialize(&CB);
      point3d_uninitialize(&CG);
      point3d_uninitialize(&C);
      point3d_uninitialize(&ND2);

      int  Nii  = -1;
      int CAii  = -1;
      int CBii  = -1;
      int CGii  = -1;
      int Cii   = -1;
      int ND2ii = -1;

      int count = 0;
      for (int i = 0; i < res->size; ++i)
      {
	    if (strcmp(res->atoms[i].loc, "N") == 0)
	    {
		  N = res->atoms[i].center;

		  Nii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "CA") == 0)
	    {
		  CA = res->atoms[i].center;

		  CAii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "CB") == 0)
	    {
		  CB = res->atoms[i].center;

		  CBii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "CG") == 0)
	    {
		  CG = res->atoms[i].center;

		  CGii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "C") == 0)
	    {
		  C = res->atoms[i].center;

		  Cii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "ND2") == 0)
	    {
		  ND2= res->atoms[i].center;

		  ND2ii = i;
		  count++;
	    }
      }
      if (count != 6)
      { /* Exception Handling */
	    fprintf(stderr, "Error in function %s()... All required atoms not found %d\n", __func__, count);
	     
      }
      double tors_ang = torsion_angle(C, CA, CB, CG);

      double adjust = torad(120.0);
      double hbangle = torad(106.97);
      double hbdist = 0.97;
      Point3d HB1 = fix_atom(C, CA, CB, hbdist, hbangle, tors_ang + adjust);
      residue_addh(res, CBii, HB1, "HB1");
      Point3d HB2 = fix_atom(C, CA, CB, hbdist, hbangle, tors_ang - adjust);
      residue_addh(res, CBii, HB2, "HB2");
      Point3d HA = fix_atom(CG, CB, CA, hbdist, hbangle, tors_ang - adjust);
      residue_addh(res, CAii, HA, "HA");
      Point3d H = fix_atom(CB, CA, N, 0.86, hbangle, adjust);
      residue_addh(res, Nii, H, "H");
      Point3d H21 = fix_atom(CB, CG, ND2, 0.86, torad(120.0), torad(0.0));

      residue_addh(res, ND2ii, H21, "H21");
      Point3d H22 = fix_atom(CB, CG, ND2, 0.86, torad(-120.0), torad(0.0));

      residue_addh(res, ND2ii, H22, "H22");

}
void lys_addh(struct residue *res)
{
      Point3d N;
      Point3d C;
      Point3d CA;
      Point3d CB;
      Point3d CG;
      Point3d CD;
      Point3d CE;
      Point3d NZ;

      point3d_uninitialize(&N);
      point3d_uninitialize(&C);
      point3d_uninitialize(&CA);
      point3d_uninitialize(&CB);
      point3d_uninitialize(&CG);
      point3d_uninitialize(&CD);
      point3d_uninitialize(&CE);
      point3d_uninitialize(&NZ);

      int  Nii  = -1;
      int  Cii  = -1;
      int CAii  = -1;
      int CBii  = -1;
      int CGii  = -1;
      int CDii  = -1;
      int CEii  = -1;
      int NZii  = -1;

      int count = 0;
      for (int i = 0; i < res->size; ++i)
      {
	    if (strcmp(res->atoms[i].loc, "N") == 0)
	    {
		  N= res->atoms[i].center;

		  Nii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "C") == 0)
	    {
		  C = res->atoms[i].center;

		  Cii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "CA") == 0)
	    {
		  CA = res->atoms[i].center;

		  CAii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "CB") == 0)
	    {
		  CB = res->atoms[i].center;

		  CBii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "CG") == 0)
	    {
		  CG = res->atoms[i].center;

		  CGii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "CD") == 0)
	    {
		  CD = res->atoms[i].center;

		  CDii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "CE") == 0)
	    {
		  CE = res->atoms[i].center;

		  CEii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "NZ") == 0)
	    {
		  NZ = res->atoms[i].center;

		  NZii = i;
		  count++;
	    }
      }
      if (count != 8)
      { /* Exception Handling */
	    fprintf(stderr, "Error in function %s()... All required atoms not found %d\n", __func__, count);
	    return ;
	     
      }
      double tors_ang = torsion_angle(C,CA, CB, CG);

      double adjust = torad(120.0);
      double hbangle = torad(106.97);
      double hbdist = 0.97;
      Point3d HA = fix_atom(CG, CB,CA, hbdist, hbangle, tors_ang - adjust);
      residue_addh(res, CAii, HA, "HA");
      Point3d HB1 = fix_atom(C, CA, CB, hbdist, hbangle, tors_ang + adjust);
      residue_addh(res, CBii, HB1, "HB1");
      Point3d HB2 = fix_atom(C, CA, CB, hbdist, hbangle, tors_ang - adjust);
      residue_addh(res, CBii, HB2, "HB2");

      tors_ang = torsion_angle(CB, CG, CD, CE);
      Point3d HG2 = fix_atom(CE,CD,CG, hbdist, hbangle, tors_ang + adjust);
      residue_addh(res, CGii, HG2, "HG2");
      Point3d HG3 = fix_atom(CE, CD, CG, hbdist, hbangle, tors_ang - adjust);
      residue_addh(res, CGii, HG3, "HG3");
      Point3d HD2 = fix_atom(CB, CG, CD, hbdist, hbangle, tors_ang + adjust);
      residue_addh(res, CDii, HD2, "HD2");
      Point3d HD3 = fix_atom(CB, CG, CD, hbdist, hbangle, tors_ang - adjust);
      residue_addh(res, CDii, HD3, "HD3");

      tors_ang = torsion_angle(CG, CD, CE, NZ);
      Point3d HE1 = fix_atom(CG,CD,CE, hbdist, hbangle, tors_ang + adjust);
      residue_addh(res, CEii, HE1, "HE1");
      Point3d HE2 = fix_atom(CG,CD,CE, hbdist, hbangle, tors_ang - adjust);
      residue_addh(res, CEii, HE2, "HE2");

      Point3d HZ1 = fix_atom(CD, CE, NZ,0.89, hbangle, torad(0.0)+torad(120.0));
      residue_addh(res, NZii, HZ1, "HZ1");
      Point3d HZ2 = fix_atom(CD, CE, NZ, 0.89, hbangle, torad(0.0)+torad(0.0));
      residue_addh(res, NZii, HZ2, "HZ2");
      Point3d HZ3 = fix_atom(CD, CE, NZ, 0.89, hbangle, torad(0.0)-torad(120.0));
      residue_addh(res, NZii, HZ3, "HZ3");
      Point3d H = fix_atom(CB, CA, N, 0.86, hbangle, torad(180.0)+adjust);
      residue_addh(res, Nii, H, "H");


}
void leu_addh(struct residue *res)
{
      Point3d C;
      Point3d N;
      Point3d CD2;
      Point3d CD1;
      Point3d CA;
      Point3d CB;
      Point3d CG;

      point3d_uninitialize(&C);
      point3d_uninitialize(&N);
      point3d_uninitialize(&CD2);
      point3d_uninitialize(&CD1);
      point3d_uninitialize(&CA);
      point3d_uninitialize(&CB);
      point3d_uninitialize(&CG);

      int  Cii  = -1;
      int  Nii  = -1;
      int CD2ii = -1;
      int CD1ii = -1;
      int CAii  = -1;
      int CBii  = -1;
      int CGii  = -1;

      int count = 0;
      for (int i = 0; i < res->size; ++i)
      {
	    if (strcmp(res->atoms[i].loc, "N") == 0)
	    {
		  N = res->atoms[i].center;

		  Nii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "C") == 0)
	    {
		  C = res->atoms[i].center;

		  Cii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "CA") == 0)
	    {
		  CA = res->atoms[i].center;

		  CAii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "CB") == 0)
	    {
		  CB = res->atoms[i].center;

		  CBii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "CG") == 0)
	    {
		  CG = res->atoms[i].center;

		  CGii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "CD1") == 0)
	    {
		  CD1 = res->atoms[i].center;

		  CD1ii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "CD2") == 0)
	    {
		  CD2 = res->atoms[i].center;

		  CD2ii = i;
		  count++;
	    }
      }      
      if (count != 7)
      { /* Exception Handling */
	    fprintf(stderr, "Error in function %s()... All required atoms not found %d\n", __func__, count);
	     
      }
      double tors_ang = torsion_angle(C, CA, CB, CG);

      double adjust = torad(120.0);
      double hbangle = torad(106.97);
      double hbdist = 0.97;
      Point3d HA = fix_atom(CG, CB, CA, hbdist, hbangle, tors_ang - adjust);
      residue_addh(res, CAii, HA, "HA");
      Point3d HB1 = fix_atom(C, CA, CB, hbdist, hbangle, tors_ang + adjust);
      residue_addh(res, CBii, HB1, "HB1");
      Point3d HB2 = fix_atom(C, CA, CB, hbdist, hbangle, tors_ang - adjust);
      residue_addh(res, CBii, HB2, "HB2");

      tors_ang = torsion_angle(CA, CB, CG,CD1);
      Point3d HG = fix_atom(CA, CB, CG,hbdist, hbangle, tors_ang - adjust);
      residue_addh(res, CGii, HG, "HG");
      Point3d H11 = fix_atom(CB, CG,CD1, 0.97, hbangle, adjust);
      residue_addh(res, CD1ii, H11, "H11");
      Point3d H12 = fix_atom(CB, CG,CD1, 0.97, hbangle, torad(0.0));
      residue_addh(res, CD1ii, H12, "H12");
      Point3d H13 = fix_atom(CB, CG,CD1, 0.97, hbangle, -adjust);
      residue_addh(res, CD1ii, H13, "H13");
      Point3d H21 = fix_atom(CB,CG,CD2, 0.97, hbangle, adjust);
      residue_addh(res, CD2ii, H21, "H21");
      Point3d H22 = fix_atom(CB, CG,CD2, 0.97, hbangle, torad(0.0));
      residue_addh(res, CD2ii, H22, "H22");
      Point3d H23 = fix_atom(CB, CG,CD2, 0.97, hbangle, - adjust);
      residue_addh(res, CD2ii, H23, "H23");
      Point3d H = fix_atom(C, CA, N, 0.86, hbangle, torad(0.0));
      residue_addh(res, Nii, H, "H");

}
void met_addh(struct residue *res)
{
      Point3d C;
      Point3d N;
      Point3d CA;
      Point3d CB;
      Point3d CG;
      Point3d SD;
      Point3d CE;

      point3d_uninitialize(&C);
      point3d_uninitialize(&N);
      point3d_uninitialize(&CA);
      point3d_uninitialize(&CB);
      point3d_uninitialize(&CG);
      point3d_uninitialize(&SD);
      point3d_uninitialize(&CE);

      int  Cii  = -1;
      int  Nii  = -1;
      int CAii  = -1;
      int CBii  = -1;
      int CGii  = -1;
      int SDii  = -1;
      int CEii  = -1;
      int count = 0;
      for (int i = 0; i < res->size; ++i)
      {
	    if (strcmp(res->atoms[i].loc, "N") == 0)
	    {
		  N = res->atoms[i].center;

		  Nii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "C") == 0)
	    {
		  C = res->atoms[i].center;

		  Cii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "CA") == 0)
	    {
		  CA = res->atoms[i].center;

		  CAii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "CB") == 0)
	    {
		  CB = res->atoms[i].center;

		  CBii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "CG") == 0)
	    {
		  CG = res->atoms[i].center;

		  CGii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "SD") == 0)
	    {
		  SD = res->atoms[i].center;

		  SDii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "CE") == 0)
	    {
		  CE = res->atoms[i].center;

		  CEii = i;
		  count++;
	    }

      }
      if (count != 7)
      { /* Exception Handling */
	    fprintf(stderr, "Error in function %s()... All required atoms not found %d\n", __func__, count);
	     
      }
      double tors_ang = torsion_angle(C, CA, CB, CG);

      double adjust = torad(120.0);
      double hbangle = torad(106.97);
      double hbdist = 0.97;
      Point3d HA =fix_atom(CG, CB, CA, hbdist, hbangle, tors_ang -adjust);
      residue_addh(res, CAii, HA, "HA");
      Point3d HB1 = fix_atom(C, CA, CB, hbdist, hbangle, tors_ang + adjust);
      residue_addh(res, CBii, HB1, "HB1");
      Point3d HB2 = fix_atom(C, CA, CB, hbdist, hbangle, tors_ang - adjust);
      residue_addh(res, CBii, HB2, "HB2");

      tors_ang = torsion_angle(CB, CG, SD, CE);
      Point3d HG1 = fix_atom(CE, SD, CG, hbdist, hbangle, tors_ang + adjust);
      residue_addh(res, CGii, HG1, "HG1");
      Point3d HG2 = fix_atom(CE,SD, CG, hbdist, hbangle, tors_ang - adjust);
      residue_addh(res, CGii, HG2, "HG2");

      Point3d H = fix_atom(CB, CA, N, 0.86, torad(105.90), torad(180.0)+torad(150.0));
      residue_addh(res, Nii, H, "H");
      Point3d HE1 = fix_atom(CG,SD,CE, 0.97, hbangle, adjust);
      residue_addh(res, CEii, HE1, "HE1");
      Point3d HE2 = fix_atom(CG,SD,CE ,0.97, hbangle, torad(0.0));
      residue_addh(res, CEii, HE2, "HE2");
      Point3d HE3 = fix_atom(CG, SD,CE, 0.97, hbangle, -adjust);
      residue_addh(res, CEii, HE3, "HE3");
}  
void pro_addh(struct residue *res)
{
      Point3d N;
      Point3d CA;
      Point3d CB;
      Point3d CG;
      Point3d CD;

      point3d_uninitialize(&N);
      point3d_uninitialize(&CA);
      point3d_uninitialize(&CB);
      point3d_uninitialize(&CG);
      point3d_uninitialize(&CD);

      int  Nii  = -1;
      int CAii  = -1;
      int CBii  = -1;
      int CGii  = -1;
      int CDii  = -1;

      int count = 0;
      for (int i = 0; i < res->size; ++i)
      {
	    if (strcmp(res->atoms[i].loc, "N") == 0)
	    {
		  N = res->atoms[i].center;

		  Nii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "CA") == 0)
	    {
		  CA = res->atoms[i].center;

		  CAii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "CB") == 0)
	    {
		  CB = res->atoms[i].center;

		  CBii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "CG") == 0)
	    {
		  CG = res->atoms[i].center;

		  CGii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "CD") == 0)
	    {
		  CD = res->atoms[i].center;

		  CDii = i;
		  count++;
	    }
      }
      if (count != 5)
      { /* Exception Handling */
	    fprintf(stderr, "Error in function %s()... All required atoms not found %d\n", __func__, count);
	     
      }
      double tors_ang = torsion_angle(CB,CG,CD,N);

      double adjust = torad(120.0);
      double hbangle = torad(106.97);
      double hbdist = 0.97;
      Point3d HD2 = fix_atom(CB,CG,CD, hbdist, hbangle, tors_ang + adjust);
      residue_addh(res, CDii, HD2, "HD2");
      Point3d HD3 = fix_atom(CB,CG,CD, hbdist, hbangle, tors_ang - adjust);
      residue_addh(res, CDii, HD3, "HD3");
      Point3d HG2 = fix_atom(N,CD,CG, hbdist, hbangle, tors_ang + adjust);
      residue_addh(res, CGii, HG2, "HG2");
      Point3d HG3 = fix_atom(N,CD,CG ,hbdist, hbangle, tors_ang - adjust);
      residue_addh(res, CGii, HG3, "HG3");
      tors_ang = torsion_angle(CG,CB,CA,N);
      Point3d HB2 = fix_atom(N,CA,CB, hbdist, hbangle, tors_ang + adjust);
      residue_addh(res, CBii, HB2, "HB2");
      Point3d HB3 = fix_atom(N,CA,CB, hbdist, hbangle, tors_ang - adjust);
      residue_addh(res, CBii, HB3, "HB3");
      Point3d HA = fix_atom(CG,CB,CA, hbdist, hbangle, tors_ang + adjust);
      residue_addh(res, CAii, HA, "HA");
}
void gln_addh(struct residue *res)
{
      Point3d O;
      Point3d C;
      Point3d N;
      Point3d CA;
      Point3d CB;
      Point3d CG;
      Point3d CD;
      Point3d NE2;

      point3d_uninitialize(&O);
      point3d_uninitialize(&C);
      point3d_uninitialize(&N);
      point3d_uninitialize(&CA);
      point3d_uninitialize(&CB);
      point3d_uninitialize(&CG);
      point3d_uninitialize(&CD);
      point3d_uninitialize(&NE2);

      int  Oii  = -1;
      int  Cii  = -1;
      int  Nii  = -1;
      int CAii  = -1;
      int CBii  = -1;
      int CGii  = -1;
      int CDii  = -1;
      int NE2ii = -1;

      int count = 0;
      for (int i = 0; i < res->size; ++i)
      {
	    if (strcmp(res->atoms[i].loc, "O") == 0)
	    {
		  O = res->atoms[i].center;

		  Oii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "C") == 0)
	    {
		  C = res->atoms[i].center;

		  Cii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "N") == 0)
	    {
		  N= res->atoms[i].center;

		  Nii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "CA") == 0)
	    {
		  CA = res->atoms[i].center;

		  CAii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "CB") == 0)
	    {
		  CB = res->atoms[i].center;

		  CBii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "CG") == 0)
	    {
		  CG = res->atoms[i].center;

		  CGii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "CD") == 0)
	    {
		  CD = res->atoms[i].center;

		  CDii = i;
		  count++;
	    }

	    else if (strcmp(res->atoms[i].loc, "NE2") == 0)
	    {
		  NE2 = res->atoms[i].center;

		  NE2ii = i;
		  count++;
	    }
      }
      if (count != 8)
      { /* Exception Handling */
	    fprintf(stderr, "Error in function %s()... All required atoms not found %d\n", __func__, count);
	    return ;
	     
      }
      double tors_ang = torsion_angle( CA, CB, CG,CD);

      double adjust = torad(120.0);
      double hbangle = torad(106.97);
      double hbdist = 0.97;
      Point3d HB2 =fix_atom(CD,CG, CB, hbdist, hbangle, tors_ang +adjust);
      residue_addh(res, CBii, HB2, "HB2");
      Point3d HB3 =fix_atom(CD,CG, CB, hbdist, hbangle, tors_ang -adjust);
      residue_addh(res, CBii, HB3, "HB3");
      Point3d HG2 = fix_atom(CA, CB,CG, hbdist, hbangle, tors_ang + adjust);
      residue_addh(res, CGii, HG2, "HG2");
      Point3d HG3 = fix_atom(CA, CB, CG, hbdist, hbangle, tors_ang - adjust);
      residue_addh(res, CGii, HG3, "HG3");

      tors_ang = torsion_angle(O,C,CA,N);
      Point3d HA = fix_atom(O,C,CA, hbdist, hbangle, tors_ang - adjust);
      residue_addh(res, CAii, HA, "HA");
      Point3d H = fix_atom(C, CA, N, 0.86, torad(115.90), torad(170.0)+torad(150.0));
      residue_addh(res, Nii, H, "H");
      Point3d H21 = fix_atom(CG,CD,NE2, 0.86, hbangle, torad(0.0));
      residue_addh(res, NE2ii, H21, "H21");
      Point3d H22 = fix_atom(CG,CD,NE2 ,0.86, -hbangle ,torad(0.0));
      residue_addh(res, NE2ii, H22, "H22");

} 
void arg_addh(struct residue *res)
{
      Point3d C;
      Point3d CA;
      Point3d N;
      Point3d CB;
      Point3d CG;
      Point3d CD;
      Point3d NE;
      Point3d CZ;
      Point3d NH1;
      Point3d NH2;

      point3d_uninitialize(&C);
      point3d_uninitialize (&CA);
      point3d_uninitialize (&N);
      point3d_uninitialize (&CB);
      point3d_uninitialize (&CG);
      point3d_uninitialize (&CD);
      point3d_uninitialize (&NE);
      point3d_uninitialize (&CZ);
      point3d_uninitialize (&NH1);
      point3d_uninitialize (&NH2);

      int  Cii  = -1;
      int CAii  = -1;
      int  Nii  = -1;
      int CBii  = -1;
      int CGii  = -1;
      int CDii  = -1;
      int NEii  = -1;
      int CZii  = -1;
      int NH1ii = -1;
      int NH2ii = -1;



      int count = 0;
      for (int i = 0; i < res->size; ++i)
      {
	    if (strcmp(res->atoms[i].loc, "C") == 0)
	    {
		  C = res->atoms[i].center;

		  Cii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "CA") == 0)
	    {
		  CA = res->atoms[i].center;

		  CAii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "N") == 0)
	    {
		  N= res->atoms[i].center;

		  Nii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "CB") == 0)
	    {
		  CB = res->atoms[i].center;

		  CAii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "CG") == 0)
	    {
		  CG = res->atoms[i].center;

		  CGii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "CD") == 0)
	    {
		  CD= res->atoms[i].center;

		  CDii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "NE") == 0)
	    {
		  NE= res->atoms[i].center;

		  NEii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "CZ") == 0)
	    {
		  CZ = res->atoms[i].center;

		  CZii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "NH1") == 0)
	    {
		  NH1 = res->atoms[i].center;

		  NH1ii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "NH2") == 0)
	    {
		  NH2 = res->atoms[i].center;

		  NH2ii = i;
		  count++;
	    }
      }
      if (count != 10)
      { /* Exception Handling */
	    fprintf(stderr, "Error in function %s()... All required atoms not found %d\n", __func__, count);
	    return ;
	     
      }
      double tors_ang = torsion_angle(C, CA, CB, CG);

      double adjust = torad(120.0);
      double hbangle = torad(106.97);
      double hbdist = 0.97;
      Point3d HB2 = fix_atom(C, CA, CB, hbdist, hbangle, tors_ang + adjust);
      residue_addh(res, CBii, HB2, "HB2");
      Point3d HB3 = fix_atom(C, CA, CB, hbdist, hbangle, tors_ang - adjust);
      residue_addh(res, CBii, HB3, "HB3");
      Point3d HA = fix_atom(CG, CB, CA, hbdist, hbangle, tors_ang -adjust);
      residue_addh(res, CAii, HA, "HA");

      tors_ang = torsion_angle(CB,CG,CD,NE);
      Point3d HG2 = fix_atom(NE,CD,CG, hbdist, hbangle, tors_ang + adjust);
      residue_addh(res, CGii, HG2, "HG2");
      Point3d HG3 = fix_atom(NE,CD,CG, hbdist, hbangle, tors_ang - adjust);
      residue_addh(res, CGii, HG3, "HG3");
      Point3d HD2 = fix_atom(CB,CG,CD, hbdist, hbangle, tors_ang + adjust);
      residue_addh(res, CDii, HD2, "HD2");
      Point3d HD3 = fix_atom(CB,CG,CD, hbdist, hbangle, tors_ang - adjust);
      residue_addh(res, CDii, HD3, "HD3");

      tors_ang = torsion_angle(CD,NE,CZ,NH2);
      Point3d HE = fix_atom(NH2,CZ,NE, hbdist, -hbangle, tors_ang);
      residue_addh(res, NEii, HE,"HE");

      Point3d H21 = fix_atom(NE,CZ,NH2, 0.86, hbangle, tors_ang);
      residue_addh(res, NH2ii, H21, "H21");
      Point3d H22 = fix_atom(NE,CZ,NH2, 0.86, -hbangle, tors_ang);
      residue_addh(res, NH2ii, H22, "H22");
      tors_ang = torsion_angle(CD,NE,CZ,NH1);
      Point3d H11 = fix_atom(NE,CZ,NH1, 0.86, hbangle, tors_ang);
      residue_addh(res, NH1ii, H11, "H11");
      Point3d H12 = fix_atom(NE,CZ,NH1, 0.86, -hbangle, torad(0.0));
      residue_addh(res, NH1ii, H12, "H12");

      Point3d H = fix_atom(C, CA,N, 0.86, torad(112.76), torad(115.92) + adjust);
      residue_addh(res, Nii, H, "H");
}   
void ser_addh(struct residue *res)
{
      Point3d C;
      Point3d N;
      Point3d CA;
      Point3d CB;
      Point3d OG;

      point3d_uninitialize (&C);
      point3d_uninitialize (&N);
      point3d_uninitialize (&CA);
      point3d_uninitialize (&CB);
      point3d_uninitialize (&OG);

      int  Cii  = -1;
      int  Nii  = -1;
      int CAii  = -1;
      int CBii  = -1;
      int OGii  = -1;
      int count = 0;
      for (int i = 0; i < res->size; ++i)
      {
	    if (strcmp(res->atoms[i].loc, "N") == 0)
	    {
		  N = res->atoms[i].center;

		  Nii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "C") == 0)
	    {
		  C = res->atoms[i].center;

		  Cii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "CA") == 0)
	    {
		  CA = res->atoms[i].center;

		  CAii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "CB") == 0)
	    {
		  CB = res->atoms[i].center;

		  CBii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "OG") == 0)
	    {
		  OG = res->atoms[i].center;

		  OGii = i;
		  count++;
	    }
      }
      if (count != 5)
      { /* Exception Handling */
	    fprintf(stderr, "Error in function %s()... All required atoms not found %d\n", __func__, count);
	     
      }
      double tors_ang = torsion_angle(OG,CB,CA,C);

      double adjust = torad(120.0);
      double hbangle = torad(106.97);
      double hbdist = 0.97;
      Point3d HB1 = fix_atom(C, CA, CB, hbdist, hbangle, tors_ang + adjust);
      residue_addh(res, CBii, HB1, "HB1");
      Point3d HB2 = fix_atom(C, CA, CB, hbdist, hbangle, tors_ang - adjust);
      residue_addh(res, CBii, HB2, "HB2");
      Point3d HA = fix_atom(OG,CB,CA, hbdist, hbangle, tors_ang - adjust);
      residue_addh(res, CAii, HA, "HA");
      Point3d HG = fix_atom(CA, CB, OG, 0.84, hbangle, torad(0.0));
      residue_addh(res, OGii, HG, "HG");
      Point3d H = fix_atom(C, CA,N, 0.86, hbangle, torad(0.0));
      residue_addh(res, Nii, H, "H");
} 
void thr_addh(struct residue *res)
{
      Point3d CG2;
      Point3d N;
      Point3d CA;
      Point3d CB;
      Point3d OG1;

      point3d_uninitialize (&CG2);
      point3d_uninitialize (&N);
      point3d_uninitialize (&CA);
      point3d_uninitialize (&CB);
      point3d_uninitialize (&OG1);

      int CG2ii = -1;
      int  Nii  = -1;
      int CAii  = -1;
      int CBii  = -1;
      int OG1ii = -1;

      int count = 0;
      for (int i = 0; i < res->size; ++i)
      {
	    if (strcmp(res->atoms[i].loc, "N") == 0)
	    {
		  N = res->atoms[i].center;

		  Nii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "CG2") == 0)
	    {
		  CG2 = res->atoms[i].center;

		  CG2ii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "CA") == 0)
	    {
		  CA = res->atoms[i].center;

		  CAii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "CB") == 0)
	    {
		  CB = res->atoms[i].center;

		  CBii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "OG1") == 0)
	    {
		  OG1 = res->atoms[i].center;

		  OG1ii = i;
		  count++;
	    }
      }
      if (count != 5)
      { /* Exception Handling */
	    fprintf(stderr, "Error in function %s()... All required atoms not found %d\n", __func__, count);
	     
      }
      double tors_ang = torsion_angle(CG2,CB,CA,N);

      double adjust = torad(120.0);
      double hbangle = torad(106.97);
      double hbdist = 0.97;
      Point3d HB = fix_atom(N,CA,CB, hbdist, hbangle, tors_ang - adjust);
      residue_addh(res, CBii, HB, "HB");
      Point3d HA = fix_atom(CG2,CB,CA, hbdist, hbangle, tors_ang +adjust);
      residue_addh(res, CAii, HA, "HA");

      Point3d HG21 = fix_atom(CA,CB,CG2, hbdist, hbangle, torad(120.0));
      residue_addh(res, CG2ii, HG21, "H21");
      Point3d HG22 = fix_atom(CA, CB,CG2, hbdist, hbangle, torad(0.0));
      residue_addh(res, CG2ii, HG22, "H22");
      Point3d HG23 = fix_atom(CA, CB,CG2, hbdist, hbangle, torad(-120.0));
      residue_addh(res, CG2ii, HG23, "H23");

      Point3d HG = fix_atom(CG2, CB, OG1, 0.84, torad(109.0), torad(109.0)+adjust);
      residue_addh(res, OG1ii, HG, "HG");
      Point3d H = fix_atom(CB, CA,N, 0.86, torad(115.0), torad(171.0)+ adjust);
      residue_addh(res, Nii, H, "H");
} 
void val_addh(struct residue *res)
{
      Point3d CG2;
      Point3d N;
      Point3d CA;
      Point3d CB;
      Point3d CG1;

      point3d_uninitialize (&CG2);
      point3d_uninitialize (&N);
      point3d_uninitialize (&CA);
      point3d_uninitialize (&CB);
      point3d_uninitialize (&CG1);

      int CG2ii = -1;
      int  Nii  = -1;
      int CAii  = -1;
      int CBii  = -1;
      int CG1ii = -1;

      int count = 0;
      for (int i = 0; i < res->size; ++i)
      {
	    if (strcmp(res->atoms[i].loc, "N") == 0)
	    {
		  N = res->atoms[i].center;

		  Nii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "CG2") == 0)
	    {
		  CG2 = res->atoms[i].center;

		  CG2ii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "CA") == 0)
	    {
		  CA = res->atoms[i].center;

		  CAii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "CB") == 0)
	    {
		  CB = res->atoms[i].center;

		  CBii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "CG1") == 0)
	    {
		  CG1 = res->atoms[i].center;

		  CG1ii = i;
		  count++;
	    }
      }
      if (count != 5)
      { /* Exception Handling */
	    fprintf(stderr, "Error in function %s()... All required atoms not found %d\n", __func__, count);
	     
      }
      double tors_ang = torsion_angle(CG2,CB,CA,N);

      double adjust = torad(120.0);
      double hbangle = torad(106.97);
      double hbdist = 0.97;
      Point3d HB = fix_atom(N,CA,CB, hbdist, hbangle, tors_ang + adjust);
      residue_addh(res, CBii, HB, "HB");
      Point3d HA = fix_atom(CG2,CB,CA, hbdist, hbangle, tors_ang + adjust);
      residue_addh(res, CAii, HA, "HA");

      Point3d HG11 = fix_atom(CA,CB,CG1, hbdist, hbangle, torad(120.0));
      residue_addh(res, CG1ii, HG11, "H11");
      Point3d HG12 = fix_atom(CA, CB,CG1, hbdist, hbangle, torad(0.0));
      residue_addh(res, CG1ii, HG12, "H12");
      Point3d HG13 = fix_atom(CA, CB,CG1, hbdist, hbangle, torad(-120.0));
      residue_addh(res, CG1ii, HG13, "H13");

      Point3d HG21 = fix_atom(CA,CB,CG2, hbdist, hbangle, torad(120.0));
      residue_addh(res, CG2ii, HG21, "H21");
      Point3d HG22 = fix_atom(CA, CB,CG2, hbdist, hbangle, torad(0.0));
      residue_addh(res, CG2ii, HG22, "H22");
      Point3d HG23 = fix_atom(CA, CB,CG2, hbdist, hbangle, torad(-120.0));
      residue_addh(res, CG2ii, HG23, "H23");
      Point3d H = fix_atom(CB, CA,N, 0.86, hbangle, tors_ang + adjust);
      residue_addh(res, Nii, H, "H");

} 
void trp_addh(struct residue *res)
{
      Point3d C;
      Point3d CA;
      Point3d CB;
      Point3d CG;
      Point3d CD1;
      Point3d NE1;
      Point3d CD2;
      Point3d CE3;
      Point3d CZ3;
      Point3d CH2;
      Point3d CE2;
      Point3d CZ2;

      point3d_uninitialize (&C);
      point3d_uninitialize (&CA);
      point3d_uninitialize (&CB);
      point3d_uninitialize (&CG);
      point3d_uninitialize (&CD1);
      point3d_uninitialize (&NE1);
      point3d_uninitialize (&CD2);
      point3d_uninitialize (&CE3);
      point3d_uninitialize (&CZ3);
      point3d_uninitialize (&CH2);
      point3d_uninitialize (&CE2);
      point3d_uninitialize (&CZ2);

      int  Cii  = -1;
      int CAii  = -1;
      int CBii  = -1;
      int CGii  = -1;
      int CD1ii = -1;
      int NE1ii = -1;
      int CD2ii = -1;
      int CE3ii = -1;
      int CZ3ii = -1;
      int CH2ii = -1;
      int CE2ii = -1;
      int CZ2ii = -1;


      int count = 0;
      for (int i = 0; i < res->size; ++i)
      {
	    if (strcmp(res->atoms[i].loc, "CD1") == 0)
	    {
		  CD1 = res->atoms[i].center;

		  CD1ii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "C") == 0)
	    {
		  C = res->atoms[i].center;

		  Cii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "CA") == 0)
	    {
		  CA = res->atoms[i].center;

		  CAii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "CB") == 0)
	    {
		  CB = res->atoms[i].center;

		  CBii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "CG") == 0)
	    {
		  CG = res->atoms[i].center;

		  CGii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "NE1") == 0)
	    {
		  NE1 = res->atoms[i].center;

		  NE1ii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "CD2") == 0)
	    {
		  CD2 = res->atoms[i].center;

		  CD2ii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "CE3") == 0)
	    {
		  CE3 = res->atoms[i].center;

		  CE3ii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "CZ3") == 0)
	    {
		  CZ3 = res->atoms[i].center;

		  CZ3ii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "CH2") == 0)
	    {
		  CH2 = res->atoms[i].center;

		  CH2ii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "CE2") == 0)
	    {
		  CE2 = res->atoms[i].center;

		  CE2ii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "CZ2") == 0)
	    {
		  CZ2 = res->atoms[i].center;

		  CZ2ii = i;
		  count++;
	    }
      }
      if (count != 12)
      { /* Exception Handling */
	    fprintf(stderr, "Error in function %s()... All required atoms not found %d\n", __func__, count);
	     
      }
      double tors_ang = torsion_angle(C,CA,CB,CG);

      double adjust = torad(120.0);
      double hbangle = torad(106.97);
      double hbdist = 0.97;
      Point3d HB1 = fix_atom(C, CA, CB, hbdist, hbangle, tors_ang + adjust);
      residue_addh(res, CBii, HB1, "HB1");
      Point3d HB2 = fix_atom(C, CA, CB, hbdist, hbangle, tors_ang - adjust);
      residue_addh(res, CBii, HB2, "HB2");
      Point3d HA = fix_atom(CG,CB,CA, hbdist, hbangle, tors_ang - adjust );
      residue_addh(res, CAii, HA, "HA");

      tors_ang = torsion_angle(CG,CD1,NE1,CE2);
      Point3d HD = fix_atom(CE2,NE1,CD1, hbdist, torad(-126.0), tors_ang);
      residue_addh(res, CD1ii, HD, "HD");
      Point3d HN = fix_atom(CG,CD1,NE1, 0.86, torad(-126.0), tors_ang);
      residue_addh(res, NE1ii, HN, "HN");

      tors_ang = torsion_angle(CD2,CE3,CZ3,CH2);
      Point3d HZ = fix_atom(CD2,CE3,CZ3, hbdist, torad(-120.0), tors_ang);
      residue_addh(res, CZ3ii, HZ, "HZ");
      Point3d HE = fix_atom(CH2,CZ3,CE3,hbdist, torad(-120.0), tors_ang);
      residue_addh(res, CE3ii, HE, "HE");

      tors_ang = torsion_angle(CE2,CZ2,CH2,CZ3);
      Point3d HZ2 = fix_atom(CZ3,CH2,CZ2, hbdist, torad(-120.0), tors_ang);
      residue_addh(res, CZ2ii, HZ2, "HZ2");
      Point3d H2 = fix_atom(CE2,CZ2,CH2,hbdist, torad(-120.0), tors_ang);
      residue_addh(res, CH2ii, H2, "H2");

} 
void tyr_addh(struct residue *res)
{
      Point3d CZ;
      Point3d CE2;
      Point3d CD2;
      Point3d CG;
      Point3d CE1;
      Point3d CD1;
      Point3d OH;
      Point3d CB;
      Point3d CA;
      Point3d N;

      point3d_uninitialize (&CZ);
      point3d_uninitialize (&CE2);
      point3d_uninitialize (&CD2);
      point3d_uninitialize (&CG);
      point3d_uninitialize (&CE1);
      point3d_uninitialize (&CD1);
      point3d_uninitialize (&OH);
      point3d_uninitialize (&CB);
      point3d_uninitialize (&CA);
      point3d_uninitialize (&N);

      int  CZii  = -1;
      int CE2ii = -1;
      int CD2ii = -1;
      int CGii  = -1;
      int CE1ii = -1;
      int CD1ii = -1;
      int  OHii  = -1;
      int CBii  = -1;
      int CAii  = -1;     
      int Nii   = -1;

      int count = 0;
      for (int i = 0; i < res->size; ++i)
      {
	    if (strcmp(res->atoms[i].loc, "N") == 0)
	    {
		  N = res->atoms[i].center;

		  Nii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "CG") == 0)
	    {
		  CG = res->atoms[i].center;

		  CGii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "CA") == 0)
	    {
		  CA = res->atoms[i].center;

		  CAii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "CB") == 0)
	    {
		  CB = res->atoms[i].center;

		  CBii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "CE1") == 0)
	    {
		  CE1 = res->atoms[i].center;

		  CE1ii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "CD1") == 0)
	    {
		  CD1 = res->atoms[i].center;

		  CD1ii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "CE2") == 0)
	    {
		  CE2 = res->atoms[i].center;

		  CE2ii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "CD2") == 0)
	    {
		  CD2 = res->atoms[i].center;

		  CD2ii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "OH") == 0)
	    {
		  OH = res->atoms[i].center;

		  OHii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "CZ") == 0)
	    {
		  CZ = res->atoms[i].center;
		  count++;
	    }
      }

      if (count != 10)
      { /* Exception Handling */
	    fprintf(stderr, "Error in function %s()... All required atoms not found %d\n", __func__, count);
	     
      }
      double tors_ang = torsion_angle(CZ,CE2,CD2,CG);

      double adjust = torad(120.0);
      double hbangle = torad(106.97);
      double hbdist = 0.97;
      Point3d HD2 = fix_atom(CZ,CE2,CD2, hbdist, torad(-120.0), tors_ang);
      residue_addh(res, CD2ii, HD2, "HD2");
      Point3d HE2 = fix_atom(CG,CD2,CE2, hbdist, torad(-120.0), tors_ang);
      residue_addh(res, CE2ii, HE2, "HE2");

      tors_ang = torsion_angle(CZ,CE1,CD1,CG);
      Point3d HD1 = fix_atom(CZ,CE1,CD1, hbdist, torad(-120.0), tors_ang);
      residue_addh(res, CD1ii, HD1, "HD1");
      Point3d HE1 = fix_atom(CG,CD1,CE1, hbdist, torad(-120.0), tors_ang);
      residue_addh(res, CE1ii, HE1, "HE1");

      tors_ang = torsion_angle(CG,CB,CA,N);
      Point3d HB1 = fix_atom(N, CA, CB, hbdist, hbangle, tors_ang + adjust);
      residue_addh(res, CBii, HB1, "HB1");
      Point3d HB2 = fix_atom(N, CA, CB, hbdist, hbangle, tors_ang - adjust);
      residue_addh(res, CBii, HB2, "HB2");
      Point3d HA = fix_atom(CG, CB, CA, hbdist, hbangle, tors_ang + adjust);
      residue_addh(res, CAii, HA, "HA");

      Point3d HH = fix_atom(CE2, CZ, OH, 0.84, hbangle, torad(0.0));
      residue_addh(res, OHii, HH, "HH");
      Point3d H = fix_atom(CB,CA,N, 0.86, hbangle, tors_ang + adjust);
      residue_addh(res, Nii, H, "H"); 

} 

void pentose_sugar_addh(struct residue* res){
      Point3d C1_P;
      Point3d C2_P;
      Point3d C3_P;
      Point3d O3_P;
      Point3d C4_P;
      Point3d C5_P;
      Point3d O2_P;
      Point3d O5_P;
      Point3d O4_P;
      Point3d P;

      point3d_uninitialize (&C1_P);
      point3d_uninitialize (&C2_P);
      point3d_uninitialize (&C3_P);
      point3d_uninitialize (&O3_P);
      point3d_uninitialize (&C4_P);
      point3d_uninitialize (&C5_P);
      point3d_uninitialize(&O2_P);
      point3d_uninitialize (&O5_P);
      point3d_uninitialize (&O4_P);
      point3d_uninitialize (&P);
      int C1_Pii = -1;
      int C2_Pii = -1;
      int C3_Pii = -1;
      int O3_Pii = -1;
      int C4_Pii = -1;
      int C5_Pii = -1;
      int O2_Pii = -1;
      int O5_Pii = -1;
      int O4_Pii = -1;
      int Pii    = -1;
      int count = 0;
      for (int i = 0; i < res->size; ++i)
      {

	    if (strcmp(res->atoms[i].loc, "C1'") == 0)
	    {
		  C1_P = res->atoms[i].center;

		  C1_Pii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "C2'") == 0)
	    {
		  C2_P = res->atoms[i].center;

		  C2_Pii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "C3'") == 0)
	    {
		  C3_P = res->atoms[i].center;

		  C3_Pii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "O3'") == 0)
	    {
		  O3_P = res->atoms[i].center;

		  O3_Pii = i;
		  count++;
	    }

	    else if (strcmp(res->atoms[i].loc, "C4'") == 0)
	    {
		  C4_P = res->atoms[i].center;
		  C4_Pii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "C5'") == 0)
	    {
		  C5_P = res->atoms[i].center;
		  C5_Pii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "O5'") == 0)
	    {
		  O5_P = res->atoms[i].center;
		  O5_Pii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "O2'") == 0)
	    {
		  O2_P = res->atoms[i].center;
		  O2_Pii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "O4'") == 0)
	    {
		  O4_P = res->atoms[i].center;
		  O4_Pii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "P") == 0)
	    {
		  P = res->atoms[i].center;
		  Pii = i;
		  count++;
	    }
      }
      if (count != 10)
      { /* Exception Handling */
	    fprintf(stderr, "Error in function %s()... All required atoms not found in %s (%d-%s%c)\n", __func__, res->name, res->atoms[0].resid, res->atoms[0].chain, res->atoms[0].ins[0]);
	    //	    
      }

      double adjust = torad(120.0);
      double hbangle = torad(106.97);
      double hbdist = 0.97;
      double tors_ang = torsion_angle(C2_P, C3_P, C4_P, C5_P); //added by parthajit roy to check.
      Point3d H4_P = fix_atom(C2_P,C3_P,C4_P, hbdist, hbangle, tors_ang - adjust);
      residue_addh(res, C4_Pii, H4_P, "H4'");
      ///      tors_ang = torsion_angle(C2_P, C3_P, C4_P, C5_P); //added by parthajit roy to check.
      Point3d H3_P = fix_atom(C5_P,C4_P,C3_P, hbdist, hbangle, tors_ang+adjust );
      residue_addh(res, C3_Pii, H3_P, "H3'");

      tors_ang = torsion_angle(C4_P,C5_P,O5_P,P);
      Point3d H5_P1 = fix_atom(P, O5_P, C5_P, hbdist, hbangle, tors_ang+adjust);
      residue_addh(res, C5_Pii, H5_P1, "H5'");
      Point3d H5_P2 = fix_atom(P, O5_P, C5_P, hbdist, hbangle, tors_ang - adjust);
      residue_addh(res, C5_Pii, H5_P2, "H5''");

      Point3d HO_P = fix_atom(C2_P,C3_P,O2_P, 0.86, hbangle, torad(180.0)+adjust);
      residue_addh(res, O2_Pii, HO_P, "HO'");

      tors_ang = torsion_angle(O4_P,C1_P,C2_P,C3_P);
      Point3d H2_P = fix_atom(O4_P,C1_P,C2_P, hbdist, hbangle, tors_ang -adjust);
      residue_addh(res, C2_Pii, H2_P, "H2'");
      Point3d H1_P = fix_atom(C3_P,C2_P,C1_P, hbdist, hbangle, tors_ang - adjust);
      residue_addh(res, C1_Pii, H1_P, "H1'");

}
void ade_addh(struct residue *res)
{
      Point3d N1;
      Point3d C2;
      Point3d N3;
      Point3d C4;
      Point3d C5;
      Point3d N7;
      Point3d C8;
      Point3d N9;
      Point3d C6;
      Point3d N6;

      point3d_uninitialize (&N1);
      point3d_uninitialize (&C2);
      point3d_uninitialize (&N3);
      point3d_uninitialize (&C4);
      point3d_uninitialize (&C5);
      point3d_uninitialize (&N7);
      point3d_uninitialize (&C8);
      point3d_uninitialize (&N9);
      point3d_uninitialize (&C6);
      point3d_uninitialize (&N6);

      int  N1ii = -1;
      int  C2ii = -1;
      int  N3ii = -1;
      int  C4ii = -1;
      int  C5ii = -1;
      int  N7ii = -1;
      int  C8ii = -1;
      int  N9ii = -1;
      int  C6ii = -1;
      int  N6ii = -1;

      int count = 0;
      for (int i = 0; i < res->size; ++i)
      {
	    if (strcmp(res->atoms[i].loc, "N1") == 0)
	    {
		  N1 = res->atoms[i].center;

		  N1ii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "C2") == 0)
	    {
		  C2 = res->atoms[i].center;

		  C2ii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "N3") == 0)
	    {
		  N3 = res->atoms[i].center;

		  N3ii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "C4") == 0)
	    {
		  C4 = res->atoms[i].center;

		  C4ii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "C5") == 0)
	    {
		  C5 = res->atoms[i].center;

		  C5ii = i;
		  count++;

	    }
	    else if (strcmp(res->atoms[i].loc, "N7") == 0)
	    {
		  N7 = res->atoms[i].center;

		  N7ii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "C8") == 0)
	    {
		  C8 = res->atoms[i].center;

		  C8ii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "N9") == 0)
	    {
		  N9 = res->atoms[i].center;

		  N9ii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "C6") == 0)
	    {
		  C6 = res->atoms[i].center;

		  C6ii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "N6") == 0)
	    {
		  N6 = res->atoms[i].center;

		  N6ii = i;
		  count++;
	    }
      }
      if (count != 10)
      { /* Exception Handling */
	    fprintf(stderr, "Error in function %s()... All required atoms not found %d\n", __func__, count);
	     
      }
      double tors_ang = torsion_angle(N1,C2,N3,C4);

      double adjust = torad(120.0);
      double hbangle = torad(106.97);
      double hbdist = 0.93;
      Point3d H2 = fix_atom(C4, N3, C2, hbdist, torad(-120.0), tors_ang);
      residue_addh(res, C2ii, H2, "H2");
      tors_ang = torsion_angle(C5,N7,C8,N9);
      Point3d H8 = fix_atom(C5, N7, C8, hbdist, torad(-120.0), tors_ang);
      residue_addh(res, C8ii, H8, "H8");
      Point3d H61 = fix_atom(N1,C6,N6,  0.86, torad(-120.0), torad(180.0));
      residue_addh(res, N6ii, H61, "H61");
      Point3d H62 = fix_atom(N1,C6,N6, 0.86,  torad(120.0), torad(180.0));
      residue_addh(res, N6ii, H62, "H62");

      pentose_sugar_addh(res);


} 
void gua_addh(struct residue *res)
{
      Point3d N1;
      Point3d N2;
      Point3d N7;
      Point3d N9;
      Point3d C2;
      Point3d C5;
      Point3d C6;
      Point3d C8;

      point3d_uninitialize (&N1);
      point3d_uninitialize (&N2);
      point3d_uninitialize (&N7);
      point3d_uninitialize (&N9);
      point3d_uninitialize (&C2);
      point3d_uninitialize (&C5);
      point3d_uninitialize (&C6);
      point3d_uninitialize (&C8);


      int  N1ii = -1;
      int  N2ii = -1;
      int  N7ii = -1;
      int  N9ii = -1;
      int  C2ii = -1;
      int  C5ii = -1;
      int  C6ii = -1;
      int  C8ii = -1;


      int count = 0;
      for (int i = 0; i < res->size; ++i)
      {
	    if (strcmp(res->atoms[i].loc, "N9") == 0)
	    {
		  N9= res->atoms[i].center;

		  N9ii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "C8") == 0)
	    {
		  C8= res->atoms[i].center;

		  C8ii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "N7") == 0)
	    {
		  N7= res->atoms[i].center;

		  N7ii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "C5") == 0)
	    {
		  C5= res->atoms[i].center;

		  C5ii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "C6") == 0)
	    {
		  C6= res->atoms[i].center;

		  C6ii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "N1") == 0)
	    {
		  N1= res->atoms[i].center;

		  N1ii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "C2") == 0)
	    {
		  C2= res->atoms[i].center;

		  C2ii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "N2") == 0)
	    {
		  N2= res->atoms[i].center;

		  N2ii = i;
		  count++;
	    }

      }
      if (count != 8)
      { /* Exception Handling */
	    fprintf(stderr, "Error in function %s()... All required atoms not found %d\n", __func__, count);
	     
      }

      double adjust = torad(120.0);
      double hbangle = torad(106.97);
      double hbdist = 0.93;
      double tors_ang = torsion_angle(N9,C8,N7,C5);
      Point3d H8 = fix_atom(C5,N7,C8, hbdist, torad(-120.0), tors_ang);
      residue_addh(res, C8ii, H8, "H8");

      tors_ang = torsion_angle(C5, C6, N1, C2);
      Point3d H1 = fix_atom(C5, C6, N1, 0.86, torad(-120.0), tors_ang);
      residue_addh(res, N1ii, H1, "H1");

      Point3d H2_1 = fix_atom(N1, C2, N2, 0.86, torad(120.0), 0.0);
      residue_addh(res, N2ii, H2_1, "H21");
      Point3d H2_2 = fix_atom(N1, C2, N2, 0.86, torad(-120), 0.0);
      residue_addh(res, N2ii, H2_2, "H22");

      pentose_sugar_addh(res);

}
void cyt_addh(struct residue *res)
{

      Point3d N1;
      Point3d N3;
      Point3d N4;
      Point3d C5;
      Point3d C6;
      Point3d C4;

      point3d_uninitialize (&N1);
      point3d_uninitialize (&N3);
      point3d_uninitialize (&N4);
      point3d_uninitialize (&C5);
      point3d_uninitialize (&C6);
      point3d_uninitialize (&C4);

      int N1ii = -1;
      int  N3ii = -1;
      int  N4ii = -1;
      int  C5ii = -1;
      int  C6ii = -1;
      int  C4ii = -1;



      int count = 0;
      for (int i = 0; i < res->size; ++i)
      {
	    if (strcmp(res->atoms[i].loc, "N1") == 0)
	    {
		  N1= res->atoms[i].center;

		  N1ii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "N3") == 0)
	    {
		  N3= res->atoms[i].center;

		  N3ii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "N4") == 0)
	    {
		  N4= res->atoms[i].center;

		  N4ii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "C4") == 0)
	    {
		  C4= res->atoms[i].center;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "C5") == 0)
	    {
		  C5= res->atoms[i].center;

		  C5ii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "C6") == 0)
	    {
		  C6= res->atoms[i].center;

		  C6ii = i;
		  count++;
	    }


      }
      if (count != 6)
      { /* Exception Handling */
	    fprintf(stderr, "Error in function %s()... All required atoms not found %d\n", __func__, count);
	    return;
	     
      }

      double adjust = torad(120.0);
      double hbangle = torad(106.97);
      double hbdist = 0.97;

      double tors_ang = torsion_angle(C4,C5,C6,N1);
      Point3d H5 = fix_atom(N1,C6,C5, hbdist, torad(-120.0), tors_ang);
      residue_addh(res, C5ii, H5, "H5");
      Point3d H6 = fix_atom(C4,C5,C6, hbdist, torad(-120.0), tors_ang);
      residue_addh(res, C6ii, H6, "H6");

      Point3d H4_1 = fix_atom(N3,C4,N4, 0.86, torad(-120.0), torad(180.0));
      residue_addh(res, N4ii, H4_1, "H41");
      Point3d H4_2 = fix_atom(N3,C4,N4, 0.86, torad(120.0), torad(180.0));
      residue_addh(res, N4ii, H4_2, "H42");

      pentose_sugar_addh(res);
}

void thy_addh(struct residue *res)
{
      Point3d N1;
      Point3d N3;
      Point3d C5;
      Point3d C2;
      Point3d C6;
      Point3d C4;
      Point3d C7;

      point3d_uninitialize (&N1);
      point3d_uninitialize (&N3);
      point3d_uninitialize (&C5);
      point3d_uninitialize (&C2);
      point3d_uninitialize (&C6);
      point3d_uninitialize (&C4);
      point3d_uninitialize (&C7);

      int  N1ii = -1;
      int  N3ii = -1;
      int  C5ii = -1;
      int  C2ii = -1;
      int  C6ii = -1;
      int  C4ii = -1;
      int  C7ii = -1;


      int count = 0;
      for (int i = 0; i < res->size; ++i)
      {
	    if (strcmp(res->atoms[i].loc, "N1") == 0)
	    {
		  N1= res->atoms[i].center;

		  N1ii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "N3") == 0)
	    {
		  N3= res->atoms[i].center;

		  N3ii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "C2") == 0)
	    {
		  C2= res->atoms[i].center;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "C4") == 0)
	    {
		  C4 = res->atoms[i].center;

		  C4ii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "C5") == 0)
	    {
		  C5 = res->atoms[i].center;

		  C4ii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "C6") == 0)
	    {
		  C6 = res->atoms[i].center;

		  C6ii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "C7") == 0)
	    {
		  C7= res->atoms[i].center;

		  C7ii = i;
		  count++;
	    }


      }
      if (count != 7)
      { /* Exception Handling */
	    fprintf(stderr, "Error in function %s()... All required atoms not found %d\n", __func__, count);
	    //  return;
	     
      }
      double adjust = torad(120.0);
      double hbangle = torad(106.97);
      double hbdist = 0.97;
      double tors_ang = torsion_angle(C4,C5,C6,N1);
      Point3d H6 = fix_atom(C4,C5,C6, hbdist, torad(-120.0), tors_ang);
      residue_addh(res, C6ii, H6, "H6");

      tors_ang = torsion_angle(C4,N3,C2,N1);
      Point3d H3 = fix_atom(N1,C2,N3, 0.86, torad(-120.0), tors_ang);
      residue_addh(res, N3ii, H3, "H3");

      Point3d H71 = fix_atom(C4,C5,C7, hbdist, hbangle, torad(30.0));
      residue_addh(res, C7ii, H71, "H71");
      Point3d H72 = fix_atom(C4,C5,C7, hbdist, hbangle, torad(150.0));
      residue_addh(res, C7ii, H72, "H72");
      Point3d H73 = fix_atom(C4,C5,C7, hbdist, hbangle, torad(-90.0));
      residue_addh(res, C7ii, H73, "H73");

      pentose_sugar_addh(res);
}
void ura_addh(struct residue *res)
{

      Point3d N1;
      Point3d N3;
      Point3d C5;
      Point3d C2;
      Point3d C6;
      Point3d C4;

      point3d_uninitialize (&N1);
      point3d_uninitialize (&N3);
      point3d_uninitialize (&C5);
      point3d_uninitialize (&C2);
      point3d_uninitialize (&C6);
      point3d_uninitialize (&C4);

      int  N1ii = -1;
      int  N3ii = -1;
      int  C5ii = -1;
      int  C2ii = -1;
      int  C6ii = -1;
      int  C4ii = -1;



      int count = 0;
      for (int i = 0; i < res->size; ++i)
      {
	    if (strcmp(res->atoms[i].loc, "N1") == 0)
	    {
		  N1= res->atoms[i].center;

		  N1ii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "N3") == 0)
	    {
		  N3= res->atoms[i].center;

		  N3ii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "C2") == 0)
	    {
		  C2= res->atoms[i].center;

		  C2ii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "C4") == 0)
	    {
		  C4= res->atoms[i].center;

		  C4ii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "C5") == 0)
	    {
		  C5= res->atoms[i].center;

		  C5ii = i;
		  count++;
	    }
	    else if (strcmp(res->atoms[i].loc, "C6") == 0)
	    {
		  C6= res->atoms[i].center;

		  C6ii = i;
		  count++;
	    }


      }
      if (count != 6)
      { /* Exception Handling */
	    fprintf(stderr, "Error in function %s()... All required atoms not found %d\n", __func__, count);
	    //  return;
	     
      }

      double adjust = torad(120.0);
      double hbangle = torad(106.97);
      double hbdist = 0.97;
      double tors_ang = torsion_angle(C4,C5,C6,N1);
      Point3d H5 = fix_atom(N1,C6,C5, hbdist, torad(-120.0), tors_ang);
      residue_addh(res, C5ii, H5, "H5");
      Point3d H6 = fix_atom(C4,C5,C6, hbdist, torad(120.0), tors_ang+torad(180.0));
      residue_addh(res, C6ii, H6, "H6");

      tors_ang = torsion_angle(C6,N3,C2,N1);
      Point3d H3 = fix_atom(N1,C2,N3, 0.86, torad(-120.0), tors_ang);
      residue_addh(res, N3ii, H3, "H3");

      pentose_sugar_addh(res);
}
