#ifdef _cplusplus
extern "C" {
#endif
#include "transfactor.h"


# line 78 "transfactor.dy"
void show_context_match_TransFactorMatchSetCompara(TransFactorMatchSetCompara * tfmsc,int context,FILE *ofp)
{
  int i;
  int j;
  int k;
  int p;
  char buffer[MAXLINE];

  assert(tfmsc != NULL);
  assert(tfmsc->overall != NULL);
  
  for(i=0;i<tfmsc->overall->len;i++) {
    fprintf(ofp,"Match %d Factor %s %s %d %d %d %.2f %.2f\n",i,tfmsc->overall->match[i]->factor->name,tfmsc->overall->target->name,
	    tfmsc->overall->match[i]->start,
	    tfmsc->overall->match[i]->end,
	    tfmsc->overall->match[i]->strand,
	    tfmsc->overall->match[i]->bit_score,
	    Probability2Bits(tfmsc->overall->match[i]->factor->max_prob)
	);
    for(j=0;j<tfmsc->sa->len;j++) {
      k = tfmsc->overall->match[i]->start - context;
      if( k <= 0 ) {
	k = 0;
      }
      for(p=0;k<tfmsc->overall->match[i]->start;k++,p++) {
	buffer[p] = tolower(tfmsc->sa->seq[j]->seq[k]);
      }
      for(;k<tfmsc->overall->match[i]->end;k++,p++) {
	buffer[p] = toupper(tfmsc->sa->seq[j]->seq[k]);
      }
      for(;k<tfmsc->overall->match[i]->end+context;k++,p++) {
	buffer[p] = tolower(tfmsc->sa->seq[j]->seq[k]);
      }
      buffer[p] = '\0';
      fprintf(ofp,"%40s %s\n",tfmsc->sa->seq[j]->name,buffer);
    }
    fprintf(ofp,"EndMatch\n");
  }
  

}



# line 122 "transfactor.dy"
TransFactorSet * circular_permuted_TransFactorSet(TransFactorSet * in,int rotate_number)
{
  TransFactorSet * out;
  int i;
  TransFactor * f;

  assert(in != NULL);

  out = TransFactorSet_alloc_len(in->len);

  for(i=0;i<in->len;i++) {
    f = TransFactor_alloc();
    f->name = stringalloc(in->factor[i]->name);
    f->pwm  = circular_permuted_pwmDNA(in->factor[i]->pwm,rotate_number);
    f->seed = hard_link_SeqAlign(in->factor[i]->seed);
    add_TransFactorSet(out,f);
  }

  return out;
}
  


# line 145 "transfactor.dy"
void show_TransFactorMatchSet(TransFactorMatchSet * tfms,FILE * ofp)
{
  int i;
  

  assert(tfms != NULL);
  assert(ofp != NULL);
  assert(tfms->target != NULL);

  for(i=0;i<tfms->len;i++) {
    
    fprintf(ofp,"Motif\t%s\t%d\t%d\t%d\t%s\t%.2f\t%.*s\n",tfms->target->name,tfms->match[i]->start+1,tfms->match[i]->start+tfms->match[i]->factor->seed->seq[0]->len,tfms->match[i]->strand,tfms->match[i]->factor->name,tfms->match[i]->bit_score,tfms->match[i]->factor->seed->seq[0]->len,tfms->target->seq+tfms->match[i]->start);
  }

}

# line 161 "transfactor.dy"
void show_help_TransFactorMatchPara(FILE * ofp)
{
  fprintf(ofp,"TransFactor Match Parameters\n");
  fprintf(ofp,"  -tfm_type [abs/rel/relmix] type of cutoff: absolute, relative, relative mixed\n");
  fprintf(ofp,"  -tfm_cutoff  (abs) bits cutoff for absolute matches, default 11.0\n");
  fprintf(ofp,"  -tfm_rel     [0.95] (rel/relmix) Relative to best possible score, accept if above irregardless of score\nn");
  fprintf(ofp,"  -tfm_relsoft [0.9] (relmix) Relative to best possible score, accept if above this relative and bit score\n");
  fprintf(ofp,"  -tfm_relbits [11.0] (relmix) If above relsoft and above this bits score, accept\n");
}

# line 171 "transfactor.dy"
void show_help_TransFactorComparaPara(FILE * ofp)
{
  fprintf(ofp,"TransFactor comparative filter parameters\n");
  fprintf(ofp,"  -tfc_type [precise/overlap/region] overlap criteria\n");
  fprintf(ofp,"  -tfc_region  [20] region length\n");
  fprintf(ofp,"  -tfc_missing [0]  number of sequences motif can be missing in\n");
}

# line 179 "transfactor.dy"
TransFactorComparaPara * new_TransFactorComparaPara_from_argv(int * argc,char ** argv)
{
  char * temp;
  TransFactorComparaPara * out;

  out = TransFactorComparaPara_alloc();
  
  temp = strip_out_assigned_argument(argc,argv,"tfc_type");
  if( temp != NULL ) {
    if( strcmp(temp,"precise") == 0 ) {
      out->overlap_type = TFCOMPARA_OVERLAP_PRECISE;
    } else if ( strcmp(temp,"overlap") == 0 ) {
      out->overlap_type = TFCOMPARA_OVERLAP_OVERLAP;
    } else if( strcmp(temp,"region") == 0 ) {
      out->overlap_type = TFCOMPARA_OVERLAP_REGION;
    } else {
      warn("string %s is not a valid tf_compara overlap type",temp);
      return NULL;
    }
  }

  strip_out_integer_argument(argc,argv,"tfc_region",&out->overlap_len);
  strip_out_integer_argument(argc,argv,"tfc_missing",&out->missing_seq);

  return out;
}

# line 206 "transfactor.dy"
TransFactorMatchSetCompara * calculate_TransFactorMatchSetCompara(SeqAlign * sa,TransFactorSet * tfs,TransFactorMatchPara * match_para,TransFactorComparaPara * para)
{
  int i;
  int j;
  int k;
  int count;
  int found;
  TransFactorMatchSetCompara * out;

  out = TransFactorMatchSetCompara_alloc_len(sa->len);
  out->sa = hard_link_SeqAlign(sa);
  
  for(i=0;i<sa->len;i++) {
    add_TransFactorMatchSetCompara(out,calculate_TransFactorMatchSet(sa->seq[i],tfs,match_para));
  }

  out->overall = TransFactorMatchSet_alloc_std();
  out->overall->target = hard_link_Sequence(out->tfms[0]->target);

  for(i=0;i<out->tfms[0]->len;i++) {
    count = 0;
    
    for(k=1;k<out->len;k++) {
      found = 0;
      for(j=0;j<out->tfms[k]->len;j++) {
	if( out->tfms[0]->match[i]->factor == out->tfms[k]->match[j]->factor ) {
	  if( para->overlap_type == TFCOMPARA_OVERLAP_PRECISE ) {
	    if( out->tfms[0]->match[i]->start == out->tfms[k]->match[j]->start ) {
	      found = 1;
	      break;
	    }
	  } else if ( para->overlap_type == TFCOMPARA_OVERLAP_OVERLAP ) {
	    if( !(out->tfms[0]->match[i]->end < out->tfms[k]->match[j]->start ||
		  out->tfms[0]->match[i]->start > out->tfms[k]->match[j]->end ) ) {
	      found = 1;
	      break;
	    }
	  } else if ( para->overlap_type == TFCOMPARA_OVERLAP_REGION ) {
	    if( !(out->tfms[0]->match[i]->end - para->overlap_len < out->tfms[k]->match[j]->start ||
		  out->tfms[0]->match[i]->start + para->overlap_len > out->tfms[k]->match[j]->end ) ) {
	      found = 1;
	      break;
	    }
	  } else {
	    fatal("Bad overlap type in compara tf matching");
	  }
	}
      }
      if( found == 0 ) {
	count++;
      } 
    }


    if( count <= para->missing_seq ) {
      add_TransFactorMatchSet(out->overall,hard_link_TransFactorMatch(out->tfms[0]->match[i]));
    }
  }

  return out;
								      
}

# line 269 "transfactor.dy"
void show_help_TransFactorBuildPara(FILE * ofp)
{
  fprintf(ofp,"TransFactor Build Parameters\n");
  fprintf(ofp,"  -tfb_pseudo   simple pseudo count, default 0.3\n");
  fprintf(ofp,"  -[no]tfb_warn warn on small sequence number [default yes]\n");
}

# line 276 "transfactor.dy"
TransFactorBuildPara * new_TransFactorBuildPara_from_argv(int * argc,char ** argv)
{
  TransFactorBuildPara * out;

  out = TransFactorBuildPara_alloc();

  out->pseudo_count = 0.3;
  out->rnd_dna = RandomModelDNA_std();

  strip_out_float_argument(argc,argv,"tfb_pseudo",&out->pseudo_count);

  return out;
}

# line 290 "transfactor.dy"
TransFactorMatchPara * new_TransFactorMatchPara_from_argv(int * argc,char ** argv)
{
  char * temp;
  TransFactorMatchPara * out;

  out = TransFactorMatchPara_alloc();
  out->type = TFM_RELATIVE_MIXED;
  out->relative_prob = 0.95;
  out->relative_prob_bits = 0.9;
  out->min_bits = 11.0;
  out->min_relative = 11.0;
  out->allow_N = 0;

  temp = strip_out_assigned_argument(argc,argv,"tfm_type");
  if( temp != NULL ) {
    if( strcmp(temp,"abs") == 0 ) {
      out->type = TFM_ABSOLUTE;
    } else if ( strcmp(temp,"rel") == 0 ) {
      out->type = TFM_RELATIVE;
    } else if ( strcmp(temp,"relmix") == 0 ) {
      out->type = TFM_RELATIVE_MIXED;
    } else {
      warn("Could not understand %s as a match para type",temp);
      return NULL;
    }
  }

  strip_out_float_argument(argc,argv,"tfm_cutoff",&out->min_bits);
  strip_out_float_argument(argc,argv,"tfm_rel",&out->relative_prob);
  strip_out_float_argument(argc,argv,"tfm_relsoft",&out->relative_prob_bits);
  strip_out_float_argument(argc,argv,"tfm_relbits",&out->min_relative);

  assert(out->relative_prob <= 1.0 );
  assert(out->relative_prob_bits <= 1.0 );

  return out;
}

# line 328 "transfactor.dy"
boolean build_TransFactorSet(TransFactorSet * tfs,TransFactorBuildPara * p)
{
  int i;

  assert(tfs);
  assert(p);

  for(i=0;i<tfs->len;i++) {
    build_pwm_TransFactor(tfs->factor[i],p);
  }

}


# line 342 "transfactor.dy"
boolean build_pwm_TransFactor(TransFactor * tf,TransFactorBuildPara * p)
{
  assert(tf);
  assert(tf->seed);
  assert(p);

  tf->pwm = pwmDNA_from_SeqAlign(tf->seed,p->pseudo_count);

  fold_randommodel_pwmDNA(tf->pwm,p->rnd_dna);

  tf->max_prob = max_prob_TransFactor(tf);
  tf->min_prob = min_prob_TransFactor(tf); 

  return TRUE;
}

# line 358 "transfactor.dy"
TransFactorMatchSet * calculate_TransFactorMatchSet(Sequence * seq,TransFactorSet * tfs,TransFactorMatchPara * p)
{
  int i;
  int j;
  int k;
  int len;
  double prob;
  double rev_prob;
  double cutoff;
  double t;
  TransFactorMatchSet * out;
  TransFactorMatch * m;

  Sequence * comp;

  assert(seq);
  assert(tfs);
  assert(p);


  comp = reverse_complement_Sequence(seq);

  assert(comp);


  out = TransFactorMatchSet_alloc_std();
  out->target = hard_link_Sequence(seq);

  for(i=0;i<tfs->len;i++) {
    len = seq->len - tfs->factor[i]->pwm->len;
    
    /* adjust cutoff on a per motif basis */
    switch(p->type ) {
    case TFM_ABSOLUTE :
      cutoff = Bits2Probability(p->min_bits);
      break;
    case TFM_RELATIVE :
      cutoff = Bits2Probability(Probability2Bits(tfs->factor[i]->max_prob)*p->relative_prob);
      break;
    case TFM_RELATIVE_MIXED :
      cutoff = Bits2Probability(Probability2Bits(tfs->factor[i]->max_prob)*p->relative_prob);
      t = Bits2Probability(Probability2Bits(tfs->factor[i]->max_prob)*p->relative_prob_bits);
      if( t > p->min_relative && t < cutoff) {
	cutoff = t;
      }
      break;
    default :
      fatal("Impossible - bad match parameter passed in");
    }

    /*   fprintf(stderr,"Cutoff is %.2f for factor %d (%.2f) vs %.2f (bits) and %.2f cutoff %.2f\n",cutoff,i,tfs->factor[i]->max_prob,Probability2Bits(tfs->factor[i]->max_prob),Probability2Bits(tfs->factor[i]->max_prob)*p->relative_prob,p->relative_prob); */

    for(j=0;j<len;j++) {
      prob = 1.0;
      rev_prob = 1.0;
      for(k=0;k<tfs->factor[i]->pwm->len;k++) {
	/*	fprintf(stderr,"position is %d, chr %c base %d\n",j+k,seq->seq[j+k],base_from_char(seq->seq[j+k])); */

	if( p->allow_N == 0 && base_from_char(seq->seq[j+k]) == BASE_N ) {
	  prob *= 0.0;
	}
	if( p->allow_N == 0 && base_from_char(comp->seq[j+k]) == BASE_N ) {
	  rev_prob *= 0.0;
	}

	prob *= tfs->factor[i]->pwm->pos[k]->emit[base_from_char(seq->seq[j+k])];
	rev_prob *= tfs->factor[i]->pwm->pos[k]->emit[base_from_char(comp->seq[j+k])];
      }

      if( prob > cutoff ) {
	m = TransFactorMatch_alloc();
	m->start = j;
	m->end   = j+tfs->factor[i]->seed->seq[0]->len;
	m->strand = 1;
	m->bit_score = Score2Bits(Probability2Score(prob));
	m->factor = hard_link_TransFactor(tfs->factor[i]);
	add_TransFactorMatchSet(out,m);
      }

      if( rev_prob > cutoff ) {
	m = TransFactorMatch_alloc();
	m->start = seq->len - j - tfs->factor[i]->seed->seq[0]->len;;
	m->end   = seq->len - j;
	m->strand = -1;
	m->bit_score = Score2Bits(Probability2Score(rev_prob));
	m->factor = hard_link_TransFactor(tfs->factor[i]);
	add_TransFactorMatchSet(out,m);
      }


    }
  }

  free_Sequence(comp);


  return out;

}


# line 459 "transfactor.dy"
double min_prob_TransFactor(TransFactor * tf)
{
  double sc;
  int i,j;
  double min;
  assert(tf != NULL);
  assert(tf->pwm != NULL);

  sc = 1.0;

  for(i=0;i<tf->pwm->len;i++) {
    min = tf->pwm->pos[i]->emit[0];
    for(j=1;j<4;j++) {
      if( min > tf->pwm->pos[i]->emit[j] ) {
	min = tf->pwm->pos[i]->emit[j];
      }
    }
    sc *= min;
  }

  return sc;
}


# line 483 "transfactor.dy"
double max_prob_TransFactor(TransFactor * tf)
{
  double sc;
  int i,j;
  int maxj;
  double max;
  assert(tf != NULL);
  assert(tf->pwm != NULL);

  sc = 1.0;

  for(i=0;i<tf->pwm->len;i++) {
    max = tf->pwm->pos[i]->emit[0];
    maxj = 0;
    for(j=1;j<4;j++) {
      /*      fprintf(stderr,"%d, comparing %.2f with %.2f\n",i,max,tf->pwm->pos[i]->emit[j]); */
      if( max < tf->pwm->pos[i]->emit[j] ) {
	max = tf->pwm->pos[i]->emit[j];
	maxj = j;
      }
    }
    /*fprintf(stderr,"At position %d, so far %.2f with %.2f (%d)\n",i,sc,max,maxj);*/
    sc *= max;
  }

  /*fprintf(stderr,"Returning %.2f as max\n",sc);*/
  return sc;
}



# line 514 "transfactor.dy"
void write_TransFactorSet(TransFactorSet * tfs,FILE * ofp)
{
  int i;

  for(i=0;i<tfs->len;i++) {
    write_TransFactor(tfs->factor[i],ofp);
  }

}

# line 524 "transfactor.dy"
void write_TransFactor(TransFactor * tf,FILE * ofp)
{
  double sc = 0.0;
  int i;
  int j;

  assert(tf);
  assert(tf->seed);
  
  if( tf->pwm != NULL ) {
    for(i=0;i<tf->seed->seq[0]->len;i++) {
      for(j=0;j<4;j++) {
	sc += 0.25 * tf->pwm->pos[i]->emit[j];
      }
    }
  } else {
    sc = 0.0;
  }

  fprintf(ofp,"factor %s\n",tf->name);
  fprintf(ofp,"seed\n");
  fprintf(ofp,"expected %.2f\n",sc);
  write_selex_SeqAlign(tf->seed,15,100,ofp);
  fprintf(ofp,"//\n");
  fprintf(ofp,"end factor\n");
}

# line 551 "transfactor.dy"
TransFactorSet * read_TransFactorSet_file(char * filename)
{
  TransFactorSet * tfs;
  FILE * ifp;

  assert(filename);
  
  ifp = openfile(filename,"r");
  if( ifp == NULL ) {
    warn("Could not open %s for TransFactorSet reading",filename);
    return NULL;
  }

  tfs = read_TransFactorSet(ifp);

  fclose(ifp);

  return tfs;

}


# line 573 "transfactor.dy"
TransFactorSet * read_ben_IUPAC_TransFactorSet_file(char * filename)
{
  TransFactorSet * tfs;
  FILE * ifp;

  assert(filename);
  
  ifp = openfile(filename,"r");
  if( ifp == NULL ) {
    warn("Could not open %s for TransFactorSet reading",filename);
    return NULL;
  }

  tfs = read_ben_IUPAC_TransFactorSet(ifp);

  fclose(ifp);

  return tfs;

}


 static char iupac_map[17][5] = { {'A','A','.','.','.'},
				 {'C','C','.','.','.'}, 
				 {'G','G','.','.','.'}, 
				 {'T','T','.','.','.'}, 
				 {'M','A','C','.','.'}, 
				 {'R','A','G','.','.'}, 
				 {'W','A','T','.','.'}, 
				 {'S','C','G','.','.'}, 
				 {'Y','C','T','.','.'}, 
				 {'K','G','T','.','.'}, 
				 {'V','A','C','G','.'}, 
				 {'H','A','C','T','.'}, 
				 {'D','A','G','T','.'}, 
				 {'B','C','G','T','.'}, 
				 {'X','A','G','T','C'}, 
				 {'N','A','G','T','C'}, 
    };

# line 613 "transfactor.dy"
TransFactorSet * read_ben_IUPAC_TransFactorSet(FILE * ifp)
{
  char buffer[MAXLINE];
  TransFactorSet * out;
  TransFactor * temp;
  SeqAlign * align;
  Sequence * line;

  char lines [12][40];

  char sbuffer[MAXLINE];
  int motif_no = 1;
  int seq_no = 1;
  int i;
  int j;
  int k;
  int l;

  out = TransFactorSet_alloc_std();

  while( fgets(buffer,MAXLINE,ifp) != NULL ) {

    if( buffer[0] == ' ' || buffer[0] == '#' ) {
      continue;
    }

    for(i=0;!isspace(buffer[i]);i++) {
      for(l=0;l<17;l++) {
	if( iupac_map[l][0] == buffer[i] ) {
	  break;
	}
      }

      if( l == 17 ) {
	warn("Weird non IUPAC code [%c] in %s",buffer[i],buffer);
	break;
      }

      /* we make 12 fake sequences, using the IUPAC map, 
	 moving the k, the index in the possible nuc in the map
	 when we hit . - this means A becomes 12 A's, Y becomes
	 6 C's and 6 T's etc */
      for(j=0,k=1;j<12;j++,k++) {
	for(;iupac_map[l][k] == '.' && k < 5;k++) 
	  ;
	if( k == 5 ) {
	  k = 1;
	}
	if( buffer[i] == 'N' ) {
	  lines[j][i] = 'N';
	} else {
	  lines[j][i] = iupac_map[l][k];
	}
      }
    }


    if( !isspace(buffer[i]) ) {
      continue; /* error scenario */
    } else {
      buffer[i] = '\0';
    }

    for(j=0;j<12;j++) {
      lines[j][i] = '\0';
    }

    sprintf(sbuffer,"motif_%d_%s",motif_no,buffer);
    temp = TransFactor_alloc();
    temp->name = stringalloc(sbuffer);
    

    align = SeqAlign_alloc_std();      
    for(j=0;j<12;j++) {
      sprintf(sbuffer,"fake_%d_%d",motif_no,j);
      line = Sequence_from_static_memory(sbuffer,lines[j]);
      add_SeqAlign(align,line);
    }


    temp->seed = trim_from_N_SeqAlign(align);
    free_SeqAlign(align);


    add_TransFactorSet(out,temp);
    motif_no++;
  }
  

  return out;

}

# line 706 "transfactor.dy"
TransFactorSet * read_laurence_TransFactorSet_file(char * filename)
{
  TransFactorSet * tfs;
  FILE * ifp;

  assert(filename);
  
  ifp = openfile(filename,"r");
  if( ifp == NULL ) {
    warn("Could not open %s for TransFactorSet reading",filename);
    return NULL;
  }

  tfs = read_laurence_TransFactorSet(ifp);

  fclose(ifp);

  return tfs;

}

# line 727 "transfactor.dy"
TransFactorSet * read_laurence_TransFactorSet(FILE * ifp)
{
  char buffer[MAXLINE];
  TransFactorSet * out;
  TransFactor * temp;
  SeqAlign * align;
  Sequence * line;

  char sbuffer[MAXLINE];
  int motif_no = 1;
  int seq_no = 1;
  int i;

  out = TransFactorSet_alloc_std();
  
  sprintf(sbuffer,"motif_%d",motif_no);

  temp = TransFactor_alloc();
  temp->name = stringalloc(sbuffer);
  align = SeqAlign_alloc_std();


  while( fgets(buffer,MAXLINE,ifp) ) {
    if( buffer[0] == '#' ) {
      continue;
    }
    if( strstr(buffer,"degenerate") != NULL ) {

      temp->seed = trim_from_N_SeqAlign(align);
      free_SeqAlign(align);

      add_TransFactorSet(out,temp);
      seq_no = 1;
      motif_no++;

      sprintf(sbuffer,"motif_%d",motif_no);
      
      temp = TransFactor_alloc();
      temp->name = stringalloc(sbuffer);
      align = SeqAlign_alloc_std();      
      continue;
    }
    if( buffer[0] == '=' ) {
      continue;
    }

    for(i=0;buffer[i] != '\0';i++) {
      if( buffer[i] != 'A' && buffer[i] != 'T' && buffer[i] != 'G' && buffer[i] != 'C' ) {
	buffer[i] = 'N';
      }
    }

    line = Sequence_alloc();
    sprintf(sbuffer,"motif_%d_seq_%d",motif_no,seq_no);
    seq_no++;
    line->name = stringalloc(sbuffer);
    line->seq  = stringalloc(buffer);
    line->type = SEQUENCE_DNA;
    line->len  = strlen(buffer);
    

    add_SeqAlign(align,line);
  }
      
  

  return out;

}


# line 798 "transfactor.dy"
TransFactorSet * read_TransFactorSet(FILE * ifp)
{
  char buffer[MAXLINE];
  TransFactorSet * out;
  TransFactor * temp;

  out = TransFactorSet_alloc_std();
  
  while( fgets(buffer,MAXLINE,ifp) ) {
    if( buffer[0] == '#' ) {
      continue;
    }
    if( strstartcmp(buffer,"factor") == 0 ) {
      temp = read_TransFactor(buffer,ifp);
      if( temp == NULL ) {
	warn("No transfactor; skipping; parsing may have fallen down");
	continue;
      } else {
	add_TransFactorSet(out,temp);
      }
      continue;
    }
    if( isalpha(buffer[0]) ) {
      warn("could not interpret in reading transfactorset %s",buffer);
    }
     
  }

  return out;


}


# line 832 "transfactor.dy"
TransFactor * read_TransFactor(char * line,FILE * ifp)
{
  char buffer[MAXLINE];
  char * name;
  char * t;
  TransFactor * out;
  SeqAlign * temp;

  if( strstartcmp(line,"factor") != 0 ) {
    warn("Passed in a line [%s] which has no factor to the factor reader. Not good!",line);
    return NULL;
  }


  name = line + strlen("factor")+1;
  for(;isspace(*name);name++)
    ;

  if( *name == '\0' ) {
    warn("No name in factor - which is bad! %s\n",line);
    return NULL;
  }
  for(t=name;*t && !isspace(*t);t++)
    ;
  *t = '\0';

  while( fgets(buffer,MAXLINE,ifp) != NULL ) {
    if( buffer[0] == '#' ) {
      continue;
    }

    if( strstartcmp(buffer,"seed") == 0 ) {
      temp = read_selex_SeqAlign(ifp);

      if( temp == NULL ) {
	warn("No sequence alignment for this factor. Skipping");
	return NULL;
      }
      continue;
    }
    if( strstartcmp(buffer,"end") == 0 ) {
      break;
    }
    if( isalpha(buffer[0]) ) {
      warn("Uninterpretable line in transfactor %s",buffer);
    }

    
  }

  if( temp == NULL ) {
    warn("No seed alignment read. Impossible for %s",name);
    return NULL;
  }

  out = TransFactor_alloc();
  out->seed = temp;
  out->name = stringalloc(name);
  
  return out;

}


# line 896 "transfactor.dy"
int compare_start_TransFactorMatch(TransFactorMatch * a,TransFactorMatch * b)
{
  return a->start - b->start;
}

# line 901 "transfactor.dy"
void sort_by_start_TransFactorMatchSet(TransFactorMatchSet * t)
{

  sort_TransFactorMatchSet(t,compare_start_TransFactorMatch);


}
# line 862 "transfactor.c"
/* Function:  hard_link_TransFactor(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [TransFactor *]
 *
 * Return [UNKN ]  Undocumented return value [TransFactor *]
 *
 */
TransFactor * hard_link_TransFactor(TransFactor * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a TransFactor object: passed a NULL object"); 
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  TransFactor_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [TransFactor *]
 *
 */
TransFactor * TransFactor_alloc(void) 
{
    TransFactor * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(TransFactor *) ckalloc (sizeof(TransFactor))) == NULL)  {  
      warn("TransFactor_alloc failed "); 
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->name = NULL;    
    out->seed = NULL;    
    out->pwm = NULL; 
    out->full = NULL;    
    out->max_prob = 0;   
    out->min_prob = 0;   


    return out;  
}    


/* Function:  free_TransFactor(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [TransFactor *]
 *
 * Return [UNKN ]  Undocumented return value [TransFactor *]
 *
 */
TransFactor * free_TransFactor(TransFactor * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a TransFactor obj. Should be trappable");   
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
    if( obj->name != NULL)   
      ckfree(obj->name);     
    if( obj->seed != NULL)   
      free_SeqAlign(obj->seed);  
    if( obj->pwm != NULL)    
      free_pwmDNA(obj->pwm);     
    if( obj->full != NULL)   
      free_SeqAlign(obj->full);  


    ckfree(obj); 
    return NULL; 
}    


/* Function:  swap_TransFactorSet(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_TransFactorSet
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [TransFactor **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_TransFactorSet(TransFactor ** list,int i,int j)  
{
    TransFactor * temp;  
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_TransFactorSet(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_TransFactorSet which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [TransFactor **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_TransFactorSet(TransFactor ** list,int left,int right,int (*comp)(TransFactor * ,TransFactor * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_TransFactorSet(list,left,(left+right)/2);   
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_TransFactorSet (list,++last,i); 
      }  
    swap_TransFactorSet (list,left,last);    
    qsort_TransFactorSet(list,left,last-1,comp); 
    qsort_TransFactorSet(list,last+1,right,comp);    
}    


/* Function:  sort_TransFactorSet(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_TransFactorSet
 *
 *
 * Arg:         obj [UNKN ] Object containing list [TransFactorSet *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_TransFactorSet(TransFactorSet * obj,int (*comp)(TransFactor *, TransFactor *)) 
{
    qsort_TransFactorSet(obj->factor,0,obj->len-1,comp); 
    return;  
}    


/* Function:  expand_TransFactorSet(obj,len)
 *
 * Descrip:    Really an internal function for add_TransFactorSet
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [TransFactorSet *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_TransFactorSet(TransFactorSet * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_TransFactorSet called with no need"); 
      return TRUE;   
      }  


    if( (obj->factor = (TransFactor ** ) ckrealloc (obj->factor,sizeof(TransFactor *)*len)) == NULL)     {  
      warn("ckrealloc failed for expand_TransFactorSet, returning FALSE");   
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_TransFactorSet(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [TransFactorSet *]
 * Arg:        add [OWNER] Object to add to the list [TransFactor *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_TransFactorSet(TransFactorSet * obj,TransFactor * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_TransFactorSet(obj,obj->len + TransFactorSetLISTLENGTH) == FALSE)   
        return FALSE;    
      }  


    obj->factor[obj->len++]=add; 
    return TRUE; 
}    


/* Function:  flush_TransFactorSet(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [TransFactorSet *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_TransFactorSet(TransFactorSet * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->factor[i] != NULL)    {  
        free_TransFactor(obj->factor[i]);    
        obj->factor[i] = NULL;   
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  TransFactorSet_alloc_std(void)
 *
 * Descrip:    Equivalent to TransFactorSet_alloc_len(TransFactorSetLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [TransFactorSet *]
 *
 */
TransFactorSet * TransFactorSet_alloc_std(void) 
{
    return TransFactorSet_alloc_len(TransFactorSetLISTLENGTH);   
}    


/* Function:  TransFactorSet_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [TransFactorSet *]
 *
 */
TransFactorSet * TransFactorSet_alloc_len(int len) 
{
    TransFactorSet * out;   /* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = TransFactorSet_alloc()) == NULL)   
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->factor = (TransFactor ** ) ckcalloc (len,sizeof(TransFactor *))) == NULL)   {  
      warn("Warning, ckcalloc failed in TransFactorSet_alloc_len");  
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_TransFactorSet(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [TransFactorSet *]
 *
 * Return [UNKN ]  Undocumented return value [TransFactorSet *]
 *
 */
TransFactorSet * hard_link_TransFactorSet(TransFactorSet * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a TransFactorSet object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  TransFactorSet_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [TransFactorSet *]
 *
 */
TransFactorSet * TransFactorSet_alloc(void) 
{
    TransFactorSet * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(TransFactorSet *) ckalloc (sizeof(TransFactorSet))) == NULL)    {  
      warn("TransFactorSet_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->factor = NULL;  
    out->len = out->maxlen = 0;  


    return out;  
}    


/* Function:  free_TransFactorSet(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [TransFactorSet *]
 *
 * Return [UNKN ]  Undocumented return value [TransFactorSet *]
 *
 */
TransFactorSet * free_TransFactorSet(TransFactorSet * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a TransFactorSet obj. Should be trappable");    
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
    if( obj->factor != NULL) {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->factor[i] != NULL)  
          free_TransFactor(obj->factor[i]);  
        }  
      ckfree(obj->factor);   
      }  


    ckfree(obj); 
    return NULL; 
}    


/* Function:  hard_link_TransFactorMatch(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [TransFactorMatch *]
 *
 * Return [UNKN ]  Undocumented return value [TransFactorMatch *]
 *
 */
TransFactorMatch * hard_link_TransFactorMatch(TransFactorMatch * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a TransFactorMatch object: passed a NULL object");    
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  TransFactorMatch_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [TransFactorMatch *]
 *
 */
TransFactorMatch * TransFactorMatch_alloc(void) 
{
    TransFactorMatch * out; /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(TransFactorMatch *) ckalloc (sizeof(TransFactorMatch))) == NULL)    {  
      warn("TransFactorMatch_alloc failed ");    
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->start = 0;  
    out->end = 0;    
    out->strand = 1; 
    out->bit_score = 0;  
    out->factor = NULL;  


    return out;  
}    


/* Function:  free_TransFactorMatch(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [TransFactorMatch *]
 *
 * Return [UNKN ]  Undocumented return value [TransFactorMatch *]
 *
 */
TransFactorMatch * free_TransFactorMatch(TransFactorMatch * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a TransFactorMatch obj. Should be trappable");  
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
    if( obj->factor != NULL) 
      free_TransFactor(obj->factor);     


    ckfree(obj); 
    return NULL; 
}    


/* Function:  swap_TransFactorMatchSet(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_TransFactorMatchSet
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [TransFactorMatch **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_TransFactorMatchSet(TransFactorMatch ** list,int i,int j)  
{
    TransFactorMatch * temp; 
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_TransFactorMatchSet(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_TransFactorMatchSet which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [TransFactorMatch **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_TransFactorMatchSet(TransFactorMatch ** list,int left,int right,int (*comp)(TransFactorMatch * ,TransFactorMatch * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_TransFactorMatchSet(list,left,(left+right)/2);  
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_TransFactorMatchSet (list,++last,i);    
      }  
    swap_TransFactorMatchSet (list,left,last);   
    qsort_TransFactorMatchSet(list,left,last-1,comp);    
    qsort_TransFactorMatchSet(list,last+1,right,comp);   
}    


/* Function:  sort_TransFactorMatchSet(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_TransFactorMatchSet
 *
 *
 * Arg:         obj [UNKN ] Object containing list [TransFactorMatchSet *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_TransFactorMatchSet(TransFactorMatchSet * obj,int (*comp)(TransFactorMatch *, TransFactorMatch *)) 
{
    qsort_TransFactorMatchSet(obj->match,0,obj->len-1,comp); 
    return;  
}    


/* Function:  expand_TransFactorMatchSet(obj,len)
 *
 * Descrip:    Really an internal function for add_TransFactorMatchSet
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [TransFactorMatchSet *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_TransFactorMatchSet(TransFactorMatchSet * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_TransFactorMatchSet called with no need");    
      return TRUE;   
      }  


    if( (obj->match = (TransFactorMatch ** ) ckrealloc (obj->match,sizeof(TransFactorMatch *)*len)) == NULL)     {  
      warn("ckrealloc failed for expand_TransFactorMatchSet, returning FALSE");  
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_TransFactorMatchSet(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [TransFactorMatchSet *]
 * Arg:        add [OWNER] Object to add to the list [TransFactorMatch *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_TransFactorMatchSet(TransFactorMatchSet * obj,TransFactorMatch * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_TransFactorMatchSet(obj,obj->len + TransFactorMatchSetLISTLENGTH) == FALSE) 
        return FALSE;    
      }  


    obj->match[obj->len++]=add;  
    return TRUE; 
}    


/* Function:  flush_TransFactorMatchSet(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [TransFactorMatchSet *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_TransFactorMatchSet(TransFactorMatchSet * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->match[i] != NULL) {  
        free_TransFactorMatch(obj->match[i]);    
        obj->match[i] = NULL;    
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  TransFactorMatchSet_alloc_std(void)
 *
 * Descrip:    Equivalent to TransFactorMatchSet_alloc_len(TransFactorMatchSetLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [TransFactorMatchSet *]
 *
 */
TransFactorMatchSet * TransFactorMatchSet_alloc_std(void) 
{
    return TransFactorMatchSet_alloc_len(TransFactorMatchSetLISTLENGTH); 
}    


/* Function:  TransFactorMatchSet_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [TransFactorMatchSet *]
 *
 */
TransFactorMatchSet * TransFactorMatchSet_alloc_len(int len) 
{
    TransFactorMatchSet * out;  /* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = TransFactorMatchSet_alloc()) == NULL)  
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->match = (TransFactorMatch ** ) ckcalloc (len,sizeof(TransFactorMatch *))) == NULL)  {  
      warn("Warning, ckcalloc failed in TransFactorMatchSet_alloc_len"); 
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_TransFactorMatchSet(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [TransFactorMatchSet *]
 *
 * Return [UNKN ]  Undocumented return value [TransFactorMatchSet *]
 *
 */
TransFactorMatchSet * hard_link_TransFactorMatchSet(TransFactorMatchSet * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a TransFactorMatchSet object: passed a NULL object"); 
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  TransFactorMatchSet_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [TransFactorMatchSet *]
 *
 */
TransFactorMatchSet * TransFactorMatchSet_alloc(void) 
{
    TransFactorMatchSet * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(TransFactorMatchSet *) ckalloc (sizeof(TransFactorMatchSet))) == NULL)  {  
      warn("TransFactorMatchSet_alloc failed "); 
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->match = NULL;   
    out->len = out->maxlen = 0;  
    out->target = NULL;  


    return out;  
}    


/* Function:  free_TransFactorMatchSet(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [TransFactorMatchSet *]
 *
 * Return [UNKN ]  Undocumented return value [TransFactorMatchSet *]
 *
 */
TransFactorMatchSet * free_TransFactorMatchSet(TransFactorMatchSet * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a TransFactorMatchSet obj. Should be trappable");   
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
    if( obj->match != NULL)  {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->match[i] != NULL)   
          free_TransFactorMatch(obj->match[i]);  
        }  
      ckfree(obj->match);    
      }  
    if( obj->target != NULL) 
      free_Sequence(obj->target);    


    ckfree(obj); 
    return NULL; 
}    


/* Function:  swap_TransFactorMatchSetCompara(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_TransFactorMatchSetCompara
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [TransFactorMatchSet **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_TransFactorMatchSetCompara(TransFactorMatchSet ** list,int i,int j)  
{
    TransFactorMatchSet * temp;  
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_TransFactorMatchSetCompara(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_TransFactorMatchSetCompara which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [TransFactorMatchSet **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_TransFactorMatchSetCompara(TransFactorMatchSet ** list,int left,int right,int (*comp)(TransFactorMatchSet * ,TransFactorMatchSet * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_TransFactorMatchSetCompara(list,left,(left+right)/2);   
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_TransFactorMatchSetCompara (list,++last,i); 
      }  
    swap_TransFactorMatchSetCompara (list,left,last);    
    qsort_TransFactorMatchSetCompara(list,left,last-1,comp); 
    qsort_TransFactorMatchSetCompara(list,last+1,right,comp);    
}    


/* Function:  sort_TransFactorMatchSetCompara(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_TransFactorMatchSetCompara
 *
 *
 * Arg:         obj [UNKN ] Object containing list [TransFactorMatchSetCompara *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_TransFactorMatchSetCompara(TransFactorMatchSetCompara * obj,int (*comp)(TransFactorMatchSet *, TransFactorMatchSet *)) 
{
    qsort_TransFactorMatchSetCompara(obj->tfms,0,obj->len-1,comp);   
    return;  
}    


/* Function:  expand_TransFactorMatchSetCompara(obj,len)
 *
 * Descrip:    Really an internal function for add_TransFactorMatchSetCompara
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [TransFactorMatchSetCompara *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_TransFactorMatchSetCompara(TransFactorMatchSetCompara * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_TransFactorMatchSetCompara called with no need"); 
      return TRUE;   
      }  


    if( (obj->tfms = (TransFactorMatchSet ** ) ckrealloc (obj->tfms,sizeof(TransFactorMatchSet *)*len)) == NULL)     {  
      warn("ckrealloc failed for expand_TransFactorMatchSetCompara, returning FALSE");   
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_TransFactorMatchSetCompara(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [TransFactorMatchSetCompara *]
 * Arg:        add [OWNER] Object to add to the list [TransFactorMatchSet *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_TransFactorMatchSetCompara(TransFactorMatchSetCompara * obj,TransFactorMatchSet * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_TransFactorMatchSetCompara(obj,obj->len + TransFactorMatchSetComparaLISTLENGTH) == FALSE)   
        return FALSE;    
      }  


    obj->tfms[obj->len++]=add;   
    return TRUE; 
}    


/* Function:  flush_TransFactorMatchSetCompara(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [TransFactorMatchSetCompara *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_TransFactorMatchSetCompara(TransFactorMatchSetCompara * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->tfms[i] != NULL)  {  
        free_TransFactorMatchSet(obj->tfms[i]);  
        obj->tfms[i] = NULL; 
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  TransFactorMatchSetCompara_alloc_std(void)
 *
 * Descrip:    Equivalent to TransFactorMatchSetCompara_alloc_len(TransFactorMatchSetComparaLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [TransFactorMatchSetCompara *]
 *
 */
TransFactorMatchSetCompara * TransFactorMatchSetCompara_alloc_std(void) 
{
    return TransFactorMatchSetCompara_alloc_len(TransFactorMatchSetComparaLISTLENGTH);   
}    


/* Function:  TransFactorMatchSetCompara_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [TransFactorMatchSetCompara *]
 *
 */
TransFactorMatchSetCompara * TransFactorMatchSetCompara_alloc_len(int len) 
{
    TransFactorMatchSetCompara * out;   /* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = TransFactorMatchSetCompara_alloc()) == NULL)   
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->tfms = (TransFactorMatchSet ** ) ckcalloc (len,sizeof(TransFactorMatchSet *))) == NULL) {  
      warn("Warning, ckcalloc failed in TransFactorMatchSetCompara_alloc_len");  
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_TransFactorMatchSetCompara(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [TransFactorMatchSetCompara *]
 *
 * Return [UNKN ]  Undocumented return value [TransFactorMatchSetCompara *]
 *
 */
TransFactorMatchSetCompara * hard_link_TransFactorMatchSetCompara(TransFactorMatchSetCompara * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a TransFactorMatchSetCompara object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  TransFactorMatchSetCompara_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [TransFactorMatchSetCompara *]
 *
 */
TransFactorMatchSetCompara * TransFactorMatchSetCompara_alloc(void) 
{
    TransFactorMatchSetCompara * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(TransFactorMatchSetCompara *) ckalloc (sizeof(TransFactorMatchSetCompara))) == NULL)    {  
      warn("TransFactorMatchSetCompara_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->sa = NULL;  
    out->tfms = NULL;    
    out->len = out->maxlen = 0;  
    out->overall = NULL; 


    return out;  
}    


/* Function:  free_TransFactorMatchSetCompara(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [TransFactorMatchSetCompara *]
 *
 * Return [UNKN ]  Undocumented return value [TransFactorMatchSetCompara *]
 *
 */
TransFactorMatchSetCompara * free_TransFactorMatchSetCompara(TransFactorMatchSetCompara * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a TransFactorMatchSetCompara obj. Should be trappable");    
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
    if( obj->sa != NULL) 
      free_SeqAlign(obj->sa);    
    if( obj->tfms != NULL)   {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->tfms[i] != NULL)    
          free_TransFactorMatchSet(obj->tfms[i]);    
        }  
      ckfree(obj->tfms); 
      }  
    if( obj->overall != NULL)    
      free_TransFactorMatchSet(obj->overall);    


    ckfree(obj); 
    return NULL; 
}    


/* Function:  hard_link_TransFactorComparaPara(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [TransFactorComparaPara *]
 *
 * Return [UNKN ]  Undocumented return value [TransFactorComparaPara *]
 *
 */
TransFactorComparaPara * hard_link_TransFactorComparaPara(TransFactorComparaPara * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a TransFactorComparaPara object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  TransFactorComparaPara_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [TransFactorComparaPara *]
 *
 */
TransFactorComparaPara * TransFactorComparaPara_alloc(void) 
{
    TransFactorComparaPara * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(TransFactorComparaPara *) ckalloc (sizeof(TransFactorComparaPara))) == NULL)    {  
      warn("TransFactorComparaPara_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->overlap_type = TFCOMPARA_OVERLAP_PRECISE;   
    out->overlap_len = 0;    
    out->missing_seq = 0;    


    return out;  
}    


/* Function:  free_TransFactorComparaPara(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [TransFactorComparaPara *]
 *
 * Return [UNKN ]  Undocumented return value [TransFactorComparaPara *]
 *
 */
TransFactorComparaPara * free_TransFactorComparaPara(TransFactorComparaPara * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a TransFactorComparaPara obj. Should be trappable");    
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


/* Function:  hard_link_TransFactorBuildPara(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [TransFactorBuildPara *]
 *
 * Return [UNKN ]  Undocumented return value [TransFactorBuildPara *]
 *
 */
TransFactorBuildPara * hard_link_TransFactorBuildPara(TransFactorBuildPara * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a TransFactorBuildPara object: passed a NULL object");    
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  TransFactorBuildPara_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [TransFactorBuildPara *]
 *
 */
TransFactorBuildPara * TransFactorBuildPara_alloc(void) 
{
    TransFactorBuildPara * out; /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(TransFactorBuildPara *) ckalloc (sizeof(TransFactorBuildPara))) == NULL)    {  
      warn("TransFactorBuildPara_alloc failed ");    
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->rnd_dna = NULL; 
    out->pseudo_count = 0.3; 
    out->warn_on_small_seq = TRUE;   


    return out;  
}    


/* Function:  free_TransFactorBuildPara(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [TransFactorBuildPara *]
 *
 * Return [UNKN ]  Undocumented return value [TransFactorBuildPara *]
 *
 */
TransFactorBuildPara * free_TransFactorBuildPara(TransFactorBuildPara * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a TransFactorBuildPara obj. Should be trappable");  
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
    if( obj->rnd_dna != NULL)    
      free_RandomModelDNA(obj->rnd_dna);     


    ckfree(obj); 
    return NULL; 
}    


/* Function:  hard_link_TransFactorMatchPara(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [TransFactorMatchPara *]
 *
 * Return [UNKN ]  Undocumented return value [TransFactorMatchPara *]
 *
 */
TransFactorMatchPara * hard_link_TransFactorMatchPara(TransFactorMatchPara * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a TransFactorMatchPara object: passed a NULL object");    
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  TransFactorMatchPara_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [TransFactorMatchPara *]
 *
 */
TransFactorMatchPara * TransFactorMatchPara_alloc(void) 
{
    TransFactorMatchPara * out; /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(TransFactorMatchPara *) ckalloc (sizeof(TransFactorMatchPara))) == NULL)    {  
      warn("TransFactorMatchPara_alloc failed ");    
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->min_bits = 11.0;    
    out->min_relative = 11.0;    
    out->relative_prob = 0;  
    out->relative_prob_bits = 0; 
    out->type = TFM_RELATIVE_MIXED;  
    out->allow_N = 0;    


    return out;  
}    


/* Function:  free_TransFactorMatchPara(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [TransFactorMatchPara *]
 *
 * Return [UNKN ]  Undocumented return value [TransFactorMatchPara *]
 *
 */
TransFactorMatchPara * free_TransFactorMatchPara(TransFactorMatchPara * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a TransFactorMatchPara obj. Should be trappable");  
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
