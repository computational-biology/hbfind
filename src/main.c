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
#include "polymer.h"
#include "hbfind.h"
#include "controller.h"
#include "parameter.h"




/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  main
 *  Description:  
 * =====================================================================================
 */

int main ( int argc, char *argv[] )
{
      printf("HBFind Starts\n");
      struct parameter args;
      int* file_index = (int*) malloc(argc * sizeof(int));

      if ( file_index==NULL ) {
	    fprintf ( stderr, "\ndynamic memory allocation failed in function %s()\n" , __func__);
	    exit (EXIT_FAILURE);
      }

      int file_count = 0;

      param_init(&args);

      
      process_argv(argc, argv, &args, file_index, &file_count);
      
      
      struct atom* atoms;
      int numatoms;
      for(int i=0; i<file_count; ++i){
	    strcpy(args.file.full_name, argv[file_index[i]]);

	    
	    fname_split(args.file.path, args.file.basename, args.file.ext, args.file.full_name);
	    //strcpy(args.file.ext, ".pdb");
	    if(strcmp(args.file.ext, ".cif") == 0){
		  scancif(args.file.full_name, all_residues, NULL, NULL, &atoms, &numatoms, ALL_TYPE, "label", args.bio.occu);
	    }else if(strcmp(args.file.ext, ".pdb") ==0 ){
		  scanpdb(args.file.full_name, all_residues, NULL, NULL, &atoms, &numatoms, ALL_TYPE, args.bio.occu);
	    }else{
		  fprintf(stderr, "Error in function %s()... Unrecognized file type supplied.\n", __func__);
		  exit(EXIT_FAILURE);
	    }
	    
	    struct polymer polymer;
	    polymer_create(&polymer, atoms, numatoms);
	    exec_hbfind(&polymer);

	    struct site* sites;
	    sites_init(&sites, polymer.numres);
	    sites_build(sites, &polymer);
	    find_hydro_clash(sites, polymer.numres);
	    sites_free(&sites);
	    FILE *outfp;										/* output-file pointer */
	    char outfp_file_name[512];		/* output-file name    */
	    fname_join(outfp_file_name, args.file.path, args.file.basename, "_h.pdb"); 
	    

	    

	    outfp = fopen( outfp_file_name, "w" );
	    if ( outfp == NULL ) {
		  fprintf ( stderr, "couldn't open file '%s'; %s\n",
			      outfp_file_name, strerror(errno) );
		  exit (EXIT_FAILURE);
	    }
	    polymer_printpdb(outfp, &polymer);
	    
	    if( fclose(outfp) == EOF ) {			/* close output file   */
		  fprintf ( stderr, "couldn't close file '%s'; %s\n",
			      outfp_file_name, strerror(errno) );
		  exit (EXIT_FAILURE);
	    }

	    polymer_free(&polymer);
	    free(atoms);
	    atoms = NULL;
	    
      }
      free(file_index);
      file_index = NULL;
      return EXIT_SUCCESS;
}				/* ----------  end of function main  ---------- */
