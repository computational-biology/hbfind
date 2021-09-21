/*
 * =====================================================================================
 *
 *       Filename:  controller.c
 *
 *    Description:   
 *
 *        Version:  1.0
 *        Created:  20/09/21 04:32:08 PM IST
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  PARTHAJIT ROY (PR), roy.parthajit@gmail.com
 *   Organization:  The University of Burdwan
 *
 * =====================================================================================
 */

#include "controller.h"


void site_setsrc(struct site* site, struct residue* src)
{
      site->src = src;
      site->numnn = 0;
}

void site_fill_neighbor(struct site* site, struct kdtree* tree, double incldist)
{
      kdtree_neighbors((void*)&site->nnres, &site->numnn, tree, site->src, incldist *incldist);
}

void exec_hbfind(struct atom* atoms, int numatom)
{
      struct polymer poly;
      polymer_create(&poly, atoms, numatom);
      struct kdtree kdtree;
      kdtree_init(&kdtree, poly.residues, sizeof(struct residue*), poly.numres, res_comp, res_value_at, res_distsqr);
      kdtree_build(&kdtree);

      struct site site;
      for(int i=0; i<poly.numres; ++i){
	    site_setsrc(&site, poly.residues+i);
	    site_fill_neighbor(&site, &kdtree, 20.0);
	    for(int j=0; j<site.numnn; ++j){
		  print_pdb_line(stdout, site.nnres[i]->atoms + 0);
	    }
	    exit(1);
      }

}

