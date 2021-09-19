/*
 * =====================================================================================
 *
 *       Filename:  polymer.h
 *
 *    Description:  Biopolymer 
 *
 *        Version:  1.0
 *        Created:  05/09/21 08:54:44 PM IST
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  PARTHAJIT ROY (PR), roy.parthajit@gmail.com
 *   Organization:  The University of Burdwan
 *
 * =====================================================================================
 */

#ifndef  __polymer_H__
#define  __polymer_H__

#include "biodefs.h"

struct polymer{
      struct atom* atoms;
      int numatom;
      int* residue;
      int numres;
};

struct residue{
      struct atom* atoms;
      int size;
      char name[4];
};

struct residue residue_at(struct polymer*, int resindex);
void polymer_create(struct polymer* polymer, struct atom* atoms, int numatom);
struct atom* polymer_resbeg(struct polymer* poly, int resindex);

struct atom* polymer_resend(struct polymer* poly, int resindex);
int polymer_ressize(struct polymer* poly, int resindex);
 

#endif   /* ----- #ifndef __polymer_H__  ----- */
