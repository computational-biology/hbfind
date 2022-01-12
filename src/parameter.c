/*
 * =====================================================================================
 *
 *       Filename:  parameter.c
 *
 *    Description:  Parameters descriptions.
 *
 *        Version:  1.0
 *        Created:  Friday 31 December 2021 04:36:14  IST
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  PARTHAJIT ROY (PR), roy.parthajit@gmail.com
 *   Organization:  The University of Burdwan
 *
 * =====================================================================================
 */


#include "parameter.h"

void param_init(struct parameter* args)
{
      args->bio.occu = 'B';
      args->bio.HH_distsqr = 6.00;
      strcpy(args->sys.os, "linux");
}

void process_argv(int argc, char* argv[], struct parameter* args, int file_index[], int* file_count)
{
      for(int i=1; i<argc; ++i){
	    if(strncmp(argv[i], "-occu=", 6) == 0){
	        if(strcmp(argv[i] + 6, "std") == 0){
	               args->bio.occu = 'S';
	        }else if(strcmp(argv[i] + 6, "best") == 0){
	               args->bio.occu = 'B';
	        }else if(strcmp(argv[i] + 6, "first") == 0){
	               args->bio.occu = 'S';
	        }else{
	               /* Exception Handling */ 
	    			fprintf(stderr, "Error in function %s()  .... wrong argument supplied.\n", __func__);
	    			exit(EXIT_FAILURE);
	        }
			
	    }else if(strncmp(argv[i], "-hhdist=", 8) == 0){
		  double tmp1 = atof(argv[i]+ 8);
		  args->bio.HH_distsqr = tmp1 * tmp1;
	    }else{
		  file_index[*file_count] = i;
		  *file_count += 1;
	    }
      }
}

