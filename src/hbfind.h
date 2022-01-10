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

#include "biodefs.h"
#include "polymer.h"


 
void cys_addh(struct residue* res);
void asp_addh(struct residue* res);
void ala_addh(struct residue* res);
void glu_addh(struct residue* res);
void phe_addh(struct residue* res);
void gly_addh(struct residue* res);
void his_addh(struct residue* res);
void ile_addh(struct residue* res);
void asn_addh(struct residue* res);
void lys_addh(struct residue* res);
void leu_addh(struct residue* res);
void met_addh(struct residue* res);
void pro_addh(struct residue* res);
void gln_addh(struct residue* res);
void arg_addh(struct residue* res);
void ser_addh(struct residue* res);
void thr_addh(struct residue* res);
void val_addh(struct residue* res);
void trp_addh(struct residue* res);
void tyr_addh(struct residue* res);
void ade_addh(struct residue* res);
void gua_addh(struct residue* res);
void cyt_addh(struct residue* res);
void thy_addh(struct residue* res);
void ura_addh(struct residue* res);
#endif   /* ----- #ifndef __hbfind_H__  ----- */
