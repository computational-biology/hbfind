/*
 * =====================================================================================
 *
 *       Filename:  polymer.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  16/09/21 04:00:40 PM IST
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  PARTHAJIT ROY (PR), roy.parthajit@gmail.com
 *   Organization:  The University of Burdwan
 *
 * =====================================================================================
 */

#include "polymer.h"
/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  partition_into_residues
 *  Description:  This finction returns an array containing the partition 
 *  		  information.
 *  		  Ex:    |0|15|20|30|35|
 *  		  here   num_residues is 4. first residue starts at 0 and goes up to 14
 *  		  i.e. it contains 15 elements. Second starts at 15 and goes up to 19.
 *  		  The last element starts from 30 goes up to 34. 
 * =====================================================================================
 */
static void partition_into_residues(const struct atom* atom_array, const int num_atoms, const int init_numres_guess, int** partition_array, int* num_residues)
{ 
     int max_size = init_numres_guess+1;

      int* partition = (int*) malloc ( max_size * sizeof(int) );
      if ( partition == NULL ) {
	    fprintf ( stderr, "\ndynamic memory allocation failed\n" );
	    exit (EXIT_FAILURE);
      }

      int index = 1;
      partition[0] = 0;

      for(int i=1; i<num_atoms; ++i){
	    if(atom_array[i].resid != atom_array[i-1].resid ||
			strcmp(atom_array[i].chain, atom_array[i-1].chain) != 0||
			atom_array[i].model != atom_array[i-1].model){

		  if(index == max_size){
			max_size *= 2;
			partition = (int*) realloc ( partition, max_size * sizeof(int) );
			if ( partition == NULL ) {
			      fprintf ( stderr, "\ndynamic memory reallocation failed\n" );
			      exit (EXIT_FAILURE);
			}
		  }
		  partition[index] = i;
		  index ++;
	    }
      }
      if(index == max_size){
	    max_size += 1;
	    partition = (int*) realloc ( partition, max_size * sizeof(int) );
	    if ( partition == NULL ) {
		  fprintf ( stderr, "\ndynamic memory reallocation failed\n" );
		  exit (EXIT_FAILURE);
	    }
      }
      partition[index] = num_atoms;
      index++;
      *partition_array = partition;
      *num_residues = index-1;
}
void polymer_create(struct polymer* polymer, struct atom* atoms, int numatom){

      polymer->atoms = atoms;
      polymer->numatom = numatom;
      partition_into_residues(atoms, numatom, (numatom/5), &polymer->residue, &polymer->numres);
}

struct atom* polymer_resbeg(struct polymer* poly, int resindex)
{
      return &poly->atoms[poly->residue[resindex]];
} 

struct atom* polymer_resend(struct polymer* poly, int resindex)
{
      return &poly->atoms[poly->residue[resindex + 1] - 1];
} 

int polymer_ressize(struct polymer* poly, int resindex)
{
      return (poly->residue[resindex + 1] - poly->residue[resindex]);
} 

struct residue residue_at(struct polymer* poly, int resindex){
      struct residue res;
      res.atoms = polymer_resbeg(poly, resindex);
      res.size  = polymer_ressize(poly, resindex);
      strcpy(res.name, res.atoms->resname);
      return res;
}