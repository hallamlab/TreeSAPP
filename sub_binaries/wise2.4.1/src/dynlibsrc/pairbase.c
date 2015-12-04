#ifdef _cplusplus
extern "C" {
#endif
#include "pairbase.h"

/* Function:  diagonal_tweak_PairBaseCodonModel(m,ratio_on,ratio_off_positive,ratio_off_negative)
 *
 * Descrip:    Tweaks a PairBaseCodonModel with on and off diagonal ratios
 *
 *
 * Arg:                         m [UNKN ] Undocumented argument [PairBaseCodonModel *]
 * Arg:                  ratio_on [UNKN ] Undocumented argument [double]
 * Arg:        ratio_off_positive [UNKN ] Undocumented argument [double]
 * Arg:        ratio_off_negative [UNKN ] Undocumented argument [double]
 *
 */
# line 43 "pairbase.dy"
void diagonal_tweak_PairBaseCodonModel(PairBaseCodonModel * m,double ratio_on,double ratio_off_positive,double ratio_off_negative)
{
  int a,b,c,x,y,z;
  int p;
  int codon_a;
  int codon_b;
  pairbase_type seq[3];

  for(a=0;a<5;a++) {
    for(b=0;b<5;b++) {
      for(c=0;c<5;c++) {
	for(x=0;x<5;x++) {
	  for(y=0;y<5;y++) {
	    for(z=0;z<5;z++) {

	      /* build the sequence */
	      seq[0] = MAKE_PAIRBASE(a,x);
	      seq[1] = MAKE_PAIRBASE(b,y);
	      seq[2] = MAKE_PAIRBASE(c,z);
 
	      p = pairbase_codon_from_seq(seq);



	      if( a == x && b == y && c == z ) {
		if( a == 4 && b == 4 && c== 4 ) {
		  m->codon[p] *= 0.5;
		} else {
		  m->codon[p] *= ratio_on;
		}
	      } else if ( m->codon[p] > 1.0 ) {
		/*	fprintf(stderr,"For %d %d,%d,%d vs %d,%d,%d\n",p,a,b,c,x,y,z); */
		m->codon[p] *= ratio_off_positive;
	      } else {
		m->codon[p] *= ratio_off_negative;
	      }
	    }
	  }
	}
      }
    }
  }
}


/* Function:  flatten_diagonal_PairBaseCodonModel(m,ct)
 *
 * Descrip:    flattens out the diagonal signal
 *
 *
 * Arg:         m [UNKN ] Undocumented argument [PairBaseCodonModel *]
 * Arg:        ct [UNKN ] Undocumented argument [CodonTable *]
 *
 */
# line 91 "pairbase.dy"
void flatten_diagonal_PairBaseCodonModel(PairBaseCodonModel * m,CodonTable * ct)
{
  int a,b,c;
  pairbase_type seq[3];
  int codon_a;
  int p;


  for(a=0;a<5;a++) {
    for(b=0;b<5;b++) {
      for(c=0;c<5;c++) {

	if( a < 4 && b < 4 && c < 4 ) {
	  codon_a = (25 * a) + (5 * b) + c;
	  if( is_stop_codon(codon_a,ct) ) {
	    m->codon[p] = 0.0;
	    continue;
	  }
	}

	/* build the sequence */
	seq[0] = MAKE_PAIRBASE(a,a);
	seq[1] = MAKE_PAIRBASE(b,b);
	seq[2] = MAKE_PAIRBASE(c,c);
	
	
	p = pairbase_codon_from_seq(seq);

	if( a == 4 && b == 4 && c== 4 ) {
	  m->codon[p] = 0.5;
	} else {
	  m->codon[p] = 1.0;
	}
      }
    }
  }

}

/* Function:  flatten_diagonal_PairBaseCodonModelScore(m,ct)
 *
 * Descrip:    flattens out the diagonal signal - for scores!
 *
 *
 * Arg:         m [UNKN ] Undocumented argument [PairBaseCodonModelScore *]
 * Arg:        ct [UNKN ] Undocumented argument [CodonTable *]
 *
 */
# line 133 "pairbase.dy"
void flatten_diagonal_PairBaseCodonModelScore(PairBaseCodonModelScore * m,CodonTable * ct)
{
  int a,b,c;
  pairbase_type seq[3];
  int codon_a;
  int p;


  for(a=0;a<5;a++) {
    for(b=0;b<5;b++) {
      for(c=0;c<5;c++) {

	if( a < 4 && b < 4 && c < 4 ) {
	  codon_a = (25 * a) + (5 * b) + c;
	  if( is_stop_codon(codon_a,ct) ) {
	    m->codon[p] = NEGI;
	    continue;
	  }
	}

	/* build the sequence */
	seq[0] = MAKE_PAIRBASE(a,a);
	seq[1] = MAKE_PAIRBASE(b,b);
	seq[2] = MAKE_PAIRBASE(c,c);
	
	
	p = pairbase_codon_from_seq(seq);

	if( a == 4 && b == 4 && c== 4 ) {
	  m->codon[p] = -2;
	} else {
	  m->codon[p] = 0;
	}
      }
    }
  }

}




/* Function:  zero_PairBaseModelScore(void)
 *
 * Descrip:    a 0 pairbasemodel score 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [PairBaseModelScore *]
 *
 */
# line 178 "pairbase.dy"
PairBaseModelScore * zero_PairBaseModelScore(void)
{
  PairBaseModelScore * out;
  int i;

  out = PairBaseModelScore_alloc();
  
  for(i=0;i<PAIRBASE_LENGTH;i++) {
    out->base[i] = 0;
  }

  return out;
}
      


/* Function:  very_simple_PairBaseCodonModel(id,rnd,nonm,gap,ct)
 *
 * Descrip:    Makes a PairBaseCodonModel from just a one parameter! Wow!
 *
 *
 * Arg:          id [UNKN ] Undocumented argument [Probability]
 * Arg:         rnd [UNKN ] Undocumented argument [Probability]
 * Arg:        nonm [UNKN ] Undocumented argument [Probability]
 * Arg:         gap [UNKN ] Undocumented argument [Probability]
 * Arg:          ct [UNKN ] Undocumented argument [CodonTable *]
 *
 * Return [UNKN ]  Undocumented return value [PairBaseCodonModel *]
 *
 */
# line 197 "pairbase.dy"
PairBaseCodonModel * very_simple_PairBaseCodonModel(Probability id,Probability rnd,Probability nonm,Probability gap,CodonTable * ct)
{
  CompProb * p;
  PairBaseCodonModel * out;

  p = simple_aa_CompProb(id,id,rnd);

  out = make_flat_PairBaseCodonModel(p,nonm,gap,ct);
  
  free_CompProb(p);

  return out;
}


/* Function:  make_flat_PairBaseCodonModel(cp,nonm,gap,ct)
 *
 * Descrip:    Makes a PairBaseCodonModel from a protein matrix - assumming a flat 
 *             mapping to CodonMatrix
 *
 *
 * Arg:          cp [UNKN ] Undocumented argument [CompProb *]
 * Arg:        nonm [UNKN ] Undocumented argument [Probability]
 * Arg:         gap [UNKN ] Undocumented argument [Probability]
 * Arg:          ct [UNKN ] Undocumented argument [CodonTable *]
 *
 * Return [UNKN ]  Undocumented return value [PairBaseCodonModel *]
 *
 */
# line 216 "pairbase.dy"
PairBaseCodonModel * make_flat_PairBaseCodonModel(CompProb * cp,Probability nonm,Probability gap,CodonTable * ct)
{
  CodonMatrix * cm;
  PairBaseCodonModel * out;

  cm = naive_CodonMatrix(ct,cp);

  out = make_PairBaseCodonModel(cm,nonm,gap,ct);

  free_CodonMatrix(cm);

  return out;
}

/* Function:  make_start_PairBaseCodonModelScore(ct)
 *
 * Descrip:    Makes a PairBaseCodonModel score for start codon
 *
 *
 * Arg:        ct [UNKN ] Undocumented argument [CodonTable *]
 *
 * Return [UNKN ]  Undocumented return value [PairBaseCodonModelScore *]
 *
 */
# line 233 "pairbase.dy"
PairBaseCodonModelScore * make_start_PairBaseCodonModelScore(CodonTable * ct)
{
  PairBaseCodonModel * model;
  PairBaseCodonModelScore * out;

  model = make_conserved_PairBaseCodonModel(100.0,0.000000000001,'M',ct);
  
  out = new_PairBaseCodonModelScore(model);

  free_PairBaseCodonModel(model);
  
  return out;

}

/* Function:  make_stop_PairBaseCodonModelScore(ct)
 *
 * Descrip:    Makes a PairBaseCodonModel score for start codon
 *
 *
 * Arg:        ct [UNKN ] Undocumented argument [CodonTable *]
 *
 * Return [UNKN ]  Undocumented return value [PairBaseCodonModelScore *]
 *
 */
# line 251 "pairbase.dy"
PairBaseCodonModelScore * make_stop_PairBaseCodonModelScore(CodonTable * ct)
{
  PairBaseCodonModel * model;
  PairBaseCodonModelScore * out;

  model = make_conserved_PairBaseCodonModel(100.0,0.000000000001,'X',ct);
  
  out = new_PairBaseCodonModelScore(model);

  free_PairBaseCodonModel(model);
  
  return out;

}


/* Function:  make_conserved_PairBaseCodonModel(cons,non_cons,aa_var,ct)
 *
 * Descrip:    Makes a PairBaseCodonModel for a particular character in the CodonTable
 *
 *
 * Arg:            cons [UNKN ] Undocumented argument [Probability]
 * Arg:        non_cons [UNKN ] Undocumented argument [Probability]
 * Arg:          aa_var [UNKN ] Undocumented argument [char]
 * Arg:              ct [UNKN ] Undocumented argument [CodonTable *]
 *
 * Return [UNKN ]  Undocumented return value [PairBaseCodonModel *]
 *
 */
# line 270 "pairbase.dy"
PairBaseCodonModel * make_conserved_PairBaseCodonModel(Probability cons,Probability non_cons,char aa_var,CodonTable * ct)
{
  PairBaseCodonModel * out;
  int i;
  base base_a,base_b,base_c;
  pairbase_type seq[3];
  int p;

  assert(ct);

  out = PairBaseCodonModel_alloc();

  for(i=0;i<PAIRBASE_CODON_LENGTH;i++) {
    out->codon[i] = non_cons;
  }

  for(i=0;i<125;i++) {
    if( ct->codon_str[i] == aa_var ) {
      fprintf(stderr,"Assinging %d with %c\n",i,aa_var);
      all_bases_from_codon(i,&base_a,&base_b,&base_c);
      seq[0] = MAKE_PAIRBASE(base_a,base_a);
      seq[1] = MAKE_PAIRBASE(base_b,base_b);
      seq[2] = MAKE_PAIRBASE(base_c,base_c);
      
      p = pairbase_codon_from_seq(seq);
      out->codon[p] = cons;
    }
  }

  return out;
}


/* Function:  make_PairBaseCodonModel(codon_matrix,nonm,gap,ct)
 *
 * Descrip:    Makes a PairBaseCodonModel from a CodonMatrix and parameters
 *
 *
 * Arg:        codon_matrix [UNKN ] Undocumented argument [CodonMatrix *]
 * Arg:                nonm [UNKN ] Undocumented argument [Probability]
 * Arg:                 gap [UNKN ] Undocumented argument [Probability]
 * Arg:                  ct [UNKN ] Undocumented argument [CodonTable *]
 *
 * Return [UNKN ]  Undocumented return value [PairBaseCodonModel *]
 *
 */
# line 306 "pairbase.dy"
PairBaseCodonModel * make_PairBaseCodonModel(CodonMatrix * codon_matrix,Probability nonm,Probability gap,CodonTable * ct)
{
  PairBaseCodonModel * out;
  int a,b,c,x,y,z;
  int i;
  int codon_a;
  int codon_b;
  int p;

  pairbase_type seq[3];


  assert(codon_matrix);
  assert(ct);

  out = PairBaseCodonModel_alloc();

  for(i=0;i<PAIRBASE_CODON_LENGTH;i++) {
    out->codon[i] = 0.0;
  }

  for(a=0;a<5;a++) {
    for(b=0;b<5;b++) {
      for(c=0;c<5;c++) {
	for(x=0;x<5;x++) {
	  for(y=0;y<5;y++) {
	    for(z=0;z<5;z++) {

	      /* build the sequence */
	      seq[0] = MAKE_PAIRBASE(a,x);
	      seq[1] = MAKE_PAIRBASE(b,y);
	      seq[2] = MAKE_PAIRBASE(c,z);
 
	      p = pairbase_codon_from_seq(seq);

	      codon_a = (a * 25) + (b * 5) + c;
	      codon_b = (x * 25) + (y * 5) + z;
 
	      if( is_stop_codon(codon_a,ct) || is_stop_codon(codon_b,ct) ) {
		out->codon[p] = 0.0;
		continue;
	      }

	      /* else */
	      out->codon[p] = codon_matrix->prob[codon_a][codon_b];
	      
	    }
	  }
	}
      }
    }
  }

  /* now to do blank and gap scores */

  for(a=0;a<5;a++) {
    for(b=0;b<5;b++) {
      for(c=0;c<5;c++) {
  
	codon_a = (a * 25) + (b * 5) + c;
	      
	seq[0] = MAKE_PAIRBASE(a,BASE_GAP);
	seq[1] = MAKE_PAIRBASE(b,BASE_GAP);
	seq[2] = MAKE_PAIRBASE(c,BASE_GAP);
 
	p = pairbase_codon_from_seq(seq);
       
	if( is_stop_codon(codon_a,ct) ) {
	  out->codon[p] = 0.0;
	} else { 
	  out->codon[p] = gap;
	}
	      
	seq[0] = MAKE_PAIRBASE(BASE_GAP,a);
	seq[1] = MAKE_PAIRBASE(BASE_GAP,b);
	seq[2] = MAKE_PAIRBASE(BASE_GAP,c);
 
	p = pairbase_codon_from_seq(seq);
       
	if( is_stop_codon(codon_a,ct) ) {
	  out->codon[p] = 0.0;
	} else { 
	  out->codon[p] = gap;
	}

	      
	seq[0] = MAKE_PAIRBASE(a,BASE_OPEN);
	seq[1] = MAKE_PAIRBASE(b,BASE_OPEN);
	seq[2] = MAKE_PAIRBASE(c,BASE_OPEN);
 
	p = pairbase_codon_from_seq(seq);
       
	if( is_stop_codon(codon_a,ct) ) {
	  out->codon[p] = 0.0;
	} else { 
	  out->codon[p] = nonm;
	}
	      
	seq[0] = MAKE_PAIRBASE(BASE_OPEN,a);
	seq[1] = MAKE_PAIRBASE(BASE_OPEN,b);
	seq[2] = MAKE_PAIRBASE(BASE_OPEN,c);
 
	p = pairbase_codon_from_seq(seq);
       
	if( is_stop_codon(codon_a,ct) ) {
	  out->codon[p] = 0.0;
	} else { 
	  out->codon[p] = nonm;
	}
      }
    }
  }

  return out;
}



/* Function:  simple_PairBaseModel(iden,other,gap)
 *
 * Descrip:    Makes a pair base model from simple leading diagonal
 *
 *
 * Arg:         iden [UNKN ] Undocumented argument [Probability]
 * Arg:        other [UNKN ] Undocumented argument [Probability]
 * Arg:          gap [UNKN ] Undocumented argument [Probability]
 *
 * Return [UNKN ]  Undocumented return value [PairBaseModel *]
 *
 */
# line 427 "pairbase.dy"
PairBaseModel * simple_PairBaseModel(Probability iden,Probability other,Probability gap)
{
  PairBaseModel * out;
  int i;
  int j;
  int base;


  out = PairBaseModel_alloc();

  for(i=0;i<7;i++) {
    for(j=0;j<7;j++) {
      base = (i*7)+j;
      if( i == 5 || j == 5 ) {
	out->base[base]= gap/0.25;
      } else if( i == 6 || j == 6 ) {
	out->base[base]= 0.0;
      } else if( i == 4 || j == 4 ) {
	out->base[base]= 1.0;
      } else if( i == j ) {
	out->base[base]= iden / 0.25;
      } else {
	out->base[base]= other / 0.25;
      }
    }
  }

  return out;
}



/* Function:  new_PairBaseCodonModelScore(pbcm)
 *
 * Descrip:    Makes a codon score from a codon model
 *
 *
 * Arg:        pbcm [UNKN ] Undocumented argument [PairBaseCodonModel *]
 *
 * Return [UNKN ]  Undocumented return value [PairBaseCodonModelScore *]
 *
 */
# line 462 "pairbase.dy"
PairBaseCodonModelScore * new_PairBaseCodonModelScore(PairBaseCodonModel * pbcm)
{
  PairBaseCodonModelScore * out;

  out = PairBaseCodonModelScore_alloc();

  Probability2Score_move(pbcm->codon,out->codon,PAIRBASE_CODON_LENGTH);

  return out;
}

/* Function:  new_PairBaseModelScore(pbm)
 *
 * Descrip:    Makes a base score from a base model
 *
 *
 * Arg:        pbm [UNKN ] Undocumented argument [PairBaseModel *]
 *
 * Return [UNKN ]  Undocumented return value [PairBaseModelScore *]
 *
 */
# line 476 "pairbase.dy"
PairBaseModelScore * new_PairBaseModelScore(PairBaseModel * pbm)
{
  PairBaseModelScore * out;

  out = PairBaseModelScore_alloc();

  Probability2Score_move(pbm->base,out->base,PAIRBASE_LENGTH);

  return out;
}

/* Function:  show_PairBaseModelScore(sc,ofp)
 *
 * Descrip:    Debugging
 *
 *
 * Arg:         sc [UNKN ] Undocumented argument [PairBaseModelScore *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
# line 490 "pairbase.dy"
void show_PairBaseModelScore(PairBaseModelScore * sc,FILE * ofp)
{
  int i;
  int anchor;
  int informant;
    

  for(i=0;i<PAIRBASE_LENGTH;i++) {
    anchor =anchor_base_from_pairbase(i);
    informant = informant_base_from_pairbase(i);

    fprintf(ofp," %2d %c %c %d\n",i,char_for_base(anchor),char_for_base(informant),sc->base[i]);
  }

}

/* Function:  show_PairBaseCodonModelScore(sc,ct,ofp)
 *
 * Descrip:    Debugging
 *
 *
 * Arg:         sc [UNKN ] Undocumented argument [PairBaseCodonModelScore *]
 * Arg:         ct [UNKN ] Undocumented argument [CodonTable *]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
# line 509 "pairbase.dy"
void show_PairBaseCodonModelScore(PairBaseCodonModelScore * sc,CodonTable * ct,FILE * ofp)
{
  int i;
  pairbase_type a;
  pairbase_type b;
  pairbase_type c;

  int anchor_a;
  int anchor_b;
  int anchor_c;

  int informant_a;
  int informant_b;
  int informant_c;

  char seq1[4];
  char seq2[4];

  seq1[3] = seq2[3] = '\0';

  for(i=0;i<PAIRBASE_CODON_LENGTH;i++) {
    decompose_pairbase_codon(i,&a,&b,&c);

    anchor_a = anchor_base_from_pairbase(a);
    anchor_b = anchor_base_from_pairbase(b);
    anchor_c = anchor_base_from_pairbase(c);

    informant_a = informant_base_from_pairbase(a);
    informant_b = informant_base_from_pairbase(b);
    informant_c = informant_base_from_pairbase(c);

    seq1[0] = char_for_base(anchor_a);
    seq1[1] = char_for_base(anchor_b);
    seq1[2] = char_for_base(anchor_c);

    seq2[0] = char_for_base(informant_a);
    seq2[1] = char_for_base(informant_b);
    seq2[2] = char_for_base(informant_c);


    fprintf(ofp,"%9d %s[%c] %s[%c] : %d\n",i,seq1,aminoacid_from_seq(ct,seq1),seq2,aminoacid_from_seq(ct,seq2),sc->codon[i]);

  }
  
}

/* Function:  reverse_pairbase_codon(codon)
 *
 * Descrip:    Inverts pairbase codon
 *
 *
 * Arg:        codon [UNKN ] Undocumented argument [pairbase_codon_type]
 *
 * Return [UNKN ]  Undocumented return value [pairbase_codon_type]
 *
 */
# line 558 "pairbase.dy"
pairbase_codon_type reverse_pairbase_codon(pairbase_codon_type codon)
{
  pairbase_type a;
  pairbase_type b;
  pairbase_type c;
  
  decompose_pairbase_codon(codon,&a,&b,&c);

  a = complement_pairbase(a);
  b = complement_pairbase(b);
  c = complement_pairbase(c);

  return (c*(PAIRBASE_LENGTH*PAIRBASE_LENGTH))+(b*PAIRBASE_LENGTH)+a;
}


/* Function:  complement_pairbase(b)
 *
 * Descrip:    complements a pairbase 
 *
 *
 * Arg:        b [UNKN ] Undocumented argument [pairbase_type]
 *
 * Return [UNKN ]  Undocumented return value [pairbase_type]
 *
 */
# line 577 "pairbase.dy"
pairbase_type complement_pairbase(pairbase_type b)
{
  pairbase_type anchor;
  pairbase_type informant;

  anchor = anchor_base_from_pairbase(b);
  informant = informant_base_from_pairbase(b);

  /* we reverse completement anchor */
  
  anchor = complement_base(anchor);

  if( informant != BASE_GAP && informant != BASE_OPEN ) {
    informant = complement_base(informant);
  }

  return MAKE_PAIRBASE(anchor,informant);
}


/* Function:  pairbase_codon_from_seq(seq)
 *
 * Descrip:    Makes a pairbase_codon from a pairbase_sequence
 *
 *
 * Arg:        seq [UNKN ] Undocumented argument [pairbase_type *]
 *
 * Return [UNKN ]  Undocumented return value [pairbase_codon_type]
 *
 */
# line 600 "pairbase.dy"
pairbase_codon_type pairbase_codon_from_seq(pairbase_type * seq)
{
  pairbase_type one;
  pairbase_type two;
  pairbase_type three;

  one   = (*seq);
  two   = (*(seq+1));
  three = (*(seq+2));

  return (one*(PAIRBASE_LENGTH*PAIRBASE_LENGTH))+(two*PAIRBASE_LENGTH)+three;
}



/* Function:  decompose_pairbase_codon(t,a,b,c)
 *
 * Descrip:    Decomposes a pairbase codon
 *
 *
 * Arg:        t [UNKN ] Undocumented argument [pairbase_codon_type]
 * Arg:        a [UNKN ] Undocumented argument [pairbase_type *]
 * Arg:        b [UNKN ] Undocumented argument [pairbase_type *]
 * Arg:        c [UNKN ] Undocumented argument [pairbase_type *]
 *
 */
# line 618 "pairbase.dy"
void decompose_pairbase_codon(pairbase_codon_type t,pairbase_type * a,pairbase_type * b,pairbase_type * c)
{
  assert(a);
  assert(b);
  assert(c);
    
  *a = t/(PAIRBASE_LENGTH*PAIRBASE_LENGTH);
  t -= (*a) * (PAIRBASE_LENGTH*PAIRBASE_LENGTH);

  *b = t/PAIRBASE_LENGTH;
  t -= (*b) * PAIRBASE_LENGTH;

  *c = t;

}


/* Function:  anchor_base_from_pairbase(pairbase)
 *
 * Descrip:    Finds the anchor base from a pair base
 *
 *
 * Arg:        pairbase [UNKN ] Undocumented argument [pairbase_type]
 *
 * Return [UNKN ]  Undocumented return value [base]
 *
 */
# line 638 "pairbase.dy"
base anchor_base_from_pairbase(pairbase_type pairbase)
{
  int top;

  top = (int)pairbase / 7;

  return top;
}

/* Function:  informant_base_from_pairbase(pairbase)
 *
 * Descrip:    Finds the informant base from a pair base
 *
 *
 * Arg:        pairbase [UNKN ] Undocumented argument [pairbase_type]
 *
 * Return [UNKN ]  Undocumented return value [base]
 *
 */
# line 650 "pairbase.dy"
base informant_base_from_pairbase(pairbase_type pairbase)
{
  int top;

  top = (int) pairbase /7;
  pairbase -= top*7;

  return pairbase;
}


/* Function:  char_for_base(base)
 *
 * Descrip:    gives back the character for the base
 *
 *
 * Arg:        base [UNKN ] Undocumented argument [int]
 *
 * Return [UNKN ]  Undocumented return value [char]
 *
 */
# line 664 "pairbase.dy"
char char_for_base(int base)
{
  if( base < 5 ) {
    return char_from_base(base);
  } 
  if( base == BASE_GAP ) {
    return '-';
  } else {
    return ' ';
  }
}
  









# line 795 "pairbase.c"
/* Function:  hard_link_PairBaseModel(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [PairBaseModel *]
 *
 * Return [UNKN ]  Undocumented return value [PairBaseModel *]
 *
 */
PairBaseModel * hard_link_PairBaseModel(PairBaseModel * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a PairBaseModel object: passed a NULL object");   
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  PairBaseModel_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [PairBaseModel *]
 *
 */
PairBaseModel * PairBaseModel_alloc(void) 
{
    PairBaseModel * out;/* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(PairBaseModel *) ckalloc (sizeof(PairBaseModel))) == NULL)  {  
      warn("PairBaseModel_alloc failed ");   
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    /* base[PAIRBASE_LENGTH] is an array: no default possible */ 


    return out;  
}    


/* Function:  free_PairBaseModel(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [PairBaseModel *]
 *
 * Return [UNKN ]  Undocumented return value [PairBaseModel *]
 *
 */
PairBaseModel * free_PairBaseModel(PairBaseModel * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a PairBaseModel obj. Should be trappable"); 
      return NULL;   
      }  


#ifdef PTHREAD   
    assert(pthread_mutex_lock(&(obj->dynamite_mutex)) == 0); 
#endif   
    if( obj->dynamite_hard_link > 1)     {  
      return_early = 1;  
      obj->dynamite_hard_link--; 
      }  
#ifdef PTHREAD   
    assert(pthread_mutex_unlock(&(obj->dynamite_mutex)) == 0);   
#endif   
    if( return_early == 1)   
      return NULL;   


    ckfree(obj); 
    return NULL; 
}    


/* Function:  hard_link_PairBaseCodonModel(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [PairBaseCodonModel *]
 *
 * Return [UNKN ]  Undocumented return value [PairBaseCodonModel *]
 *
 */
PairBaseCodonModel * hard_link_PairBaseCodonModel(PairBaseCodonModel * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a PairBaseCodonModel object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  PairBaseCodonModel_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [PairBaseCodonModel *]
 *
 */
PairBaseCodonModel * PairBaseCodonModel_alloc(void) 
{
    PairBaseCodonModel * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(PairBaseCodonModel *) ckalloc (sizeof(PairBaseCodonModel))) == NULL)    {  
      warn("PairBaseCodonModel_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    /* codon[PAIRBASE_CODON_LENGTH] is an array: no default possible */ 


    return out;  
}    


/* Function:  free_PairBaseCodonModel(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [PairBaseCodonModel *]
 *
 * Return [UNKN ]  Undocumented return value [PairBaseCodonModel *]
 *
 */
PairBaseCodonModel * free_PairBaseCodonModel(PairBaseCodonModel * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a PairBaseCodonModel obj. Should be trappable");    
      return NULL;   
      }  


#ifdef PTHREAD   
    assert(pthread_mutex_lock(&(obj->dynamite_mutex)) == 0); 
#endif   
    if( obj->dynamite_hard_link > 1)     {  
      return_early = 1;  
      obj->dynamite_hard_link--; 
      }  
#ifdef PTHREAD   
    assert(pthread_mutex_unlock(&(obj->dynamite_mutex)) == 0);   
#endif   
    if( return_early == 1)   
      return NULL;   


    ckfree(obj); 
    return NULL; 
}    


/* Function:  hard_link_PairBaseModelScore(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [PairBaseModelScore *]
 *
 * Return [UNKN ]  Undocumented return value [PairBaseModelScore *]
 *
 */
PairBaseModelScore * hard_link_PairBaseModelScore(PairBaseModelScore * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a PairBaseModelScore object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  PairBaseModelScore_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [PairBaseModelScore *]
 *
 */
PairBaseModelScore * PairBaseModelScore_alloc(void) 
{
    PairBaseModelScore * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(PairBaseModelScore *) ckalloc (sizeof(PairBaseModelScore))) == NULL)    {  
      warn("PairBaseModelScore_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    /* base[PAIRBASE_LENGTH] is an array: no default possible */ 


    return out;  
}    


/* Function:  free_PairBaseModelScore(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [PairBaseModelScore *]
 *
 * Return [UNKN ]  Undocumented return value [PairBaseModelScore *]
 *
 */
PairBaseModelScore * free_PairBaseModelScore(PairBaseModelScore * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a PairBaseModelScore obj. Should be trappable");    
      return NULL;   
      }  


#ifdef PTHREAD   
    assert(pthread_mutex_lock(&(obj->dynamite_mutex)) == 0); 
#endif   
    if( obj->dynamite_hard_link > 1)     {  
      return_early = 1;  
      obj->dynamite_hard_link--; 
      }  
#ifdef PTHREAD   
    assert(pthread_mutex_unlock(&(obj->dynamite_mutex)) == 0);   
#endif   
    if( return_early == 1)   
      return NULL;   


    ckfree(obj); 
    return NULL; 
}    


/* Function:  hard_link_PairBaseCodonModelScore(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [PairBaseCodonModelScore *]
 *
 * Return [UNKN ]  Undocumented return value [PairBaseCodonModelScore *]
 *
 */
PairBaseCodonModelScore * hard_link_PairBaseCodonModelScore(PairBaseCodonModelScore * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a PairBaseCodonModelScore object: passed a NULL object"); 
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  PairBaseCodonModelScore_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [PairBaseCodonModelScore *]
 *
 */
PairBaseCodonModelScore * PairBaseCodonModelScore_alloc(void) 
{
    PairBaseCodonModelScore * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(PairBaseCodonModelScore *) ckalloc (sizeof(PairBaseCodonModelScore))) == NULL)  {  
      warn("PairBaseCodonModelScore_alloc failed "); 
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    /* codon[PAIRBASE_CODON_LENGTH] is an array: no default possible */ 


    return out;  
}    


/* Function:  free_PairBaseCodonModelScore(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [PairBaseCodonModelScore *]
 *
 * Return [UNKN ]  Undocumented return value [PairBaseCodonModelScore *]
 *
 */
PairBaseCodonModelScore * free_PairBaseCodonModelScore(PairBaseCodonModelScore * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a PairBaseCodonModelScore obj. Should be trappable");   
      return NULL;   
      }  


#ifdef PTHREAD   
    assert(pthread_mutex_lock(&(obj->dynamite_mutex)) == 0); 
#endif   
    if( obj->dynamite_hard_link > 1)     {  
      return_early = 1;  
      obj->dynamite_hard_link--; 
      }  
#ifdef PTHREAD   
    assert(pthread_mutex_unlock(&(obj->dynamite_mutex)) == 0);   
#endif   
    if( return_early == 1)   
      return NULL;   


    ckfree(obj); 
    return NULL; 
}    



#ifdef _cplusplus
}
#endif
