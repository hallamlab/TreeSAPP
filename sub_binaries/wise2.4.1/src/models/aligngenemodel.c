#ifdef _cplusplus
extern "C" {
#endif
#include "aligngenemodel.h"




# line 77 "aligngenemodel.dy"
Probability weighted_coding_AlignGeneColumn(AlignGeneColumnStore * store,AlignGeneModelParam * agmp)
{
  int aa;
  int seq_no;
  Probability out;
  Probability aa_prob;

  int codon;


  out = 0.0;

  for(aa=0;aa<26;aa++) {
    if( agmp->rm->aminoacid[aa] < 0.0000000001 ) {
      continue;
    }

    aa_prob = agmp->rm->aminoacid[aa];

    for(seq_no=0;seq_no<store->len;seq_no++) {

      if( !isalpha(store->codon[seq_no].seq[0]) ||
	  !isalpha(store->codon[seq_no].seq[1]) ||
	  !isalpha(store->codon[seq_no].seq[2]) ) {
	return 0.5;
      }
      
      codon = codon_from_seq(store->codon[seq_no].seq);
      if( is_stop_codon(codon,agmp->ct) ) {
	if( seq_no > 0 && agmp->tolerate_nonanchor_stops ) {
	  return agmp->nonanchor_stop;
	} else {
	  return 0.0;
	}
      }
      
      aa_prob *= agmp->protein->comp[aa][aminoacid_from_codon(agmp->ct,codon)-'A'] * store->codon[seq_no].weight;
    }

    out += aa_prob;
  }

  return out;
}

# line 122 "aligngenemodel.dy"
void window_AlignGeneModel(SeqAlign * sal,AlignGeneModel * agm,AlignGeneModelParam * agmp)
{
  int i;
  int j;
  int k;
  int count;
  int win_count;
  double average_change;
  int change;

  double * change_array;

  change_array = calloc(agm->len,sizeof(double));


  fprintf(stderr,"Got Window at %f and non window at %f\n",agmp->coding_window_bonus,agmp->noncoding_window_pen);

  for(i=1;i<sal->len;i++) {

    average_change = 0.0;
    count = 0;

    for(j=0;j<sal->seq[0]->len-32;j++) {
      change = 0;
      win_count = 0;
      for(k = 0;k<31;k++) {
	if( sal->seq[0]->seq[j+k] != sal->seq[i]->seq[j+k] && sal->seq[i]->seq[j+k] != '~') {
	  change++;
	}
	if(  sal->seq[i]->seq[j+k] != '~' ) {
	  win_count++;
	}
      }
      if( win_count == 0 ) {
	change_array[j+15] = 0.9;
	continue;
      }

      change_array[j+15] = ((double) change) / win_count;
      average_change += change_array[j+15];
      count++;
      /*      fprintf(stdout,"Got %f, running average %f\n",change_array[j+15],average_change/(double)count);*/
    }

    average_change = average_change / (double) count;

    fprintf(stdout,"For %s, Average change is %f\n",sal->seq[i]->name,average_change);

    for(j=15;j<sal->seq[0]->len-16;j++) {

      if( change_array[j] < (average_change/agmp->coding_window_thres) ) {
	/*	fprintf(stdout,"At position %d in %s giving coding bonus\n",j,sal->seq[i]->name);*/

	agm->forward_coding[j] *= agmp->coding_window_bonus;
	agm->reverse_coding[j] *= agmp->coding_window_bonus;
      }
      if( change_array[j] > (average_change * agmp->noncoding_window_thres) ) {
	/*	fprintf(stdout,"At position %d in %s giving penalty\n",j,sal->seq[i]->name);*/
	agm->forward_coding[j] *= agmp->noncoding_window_pen;
	agm->reverse_coding[j] *= agmp->noncoding_window_pen;
      }
    }
  }

  free(change_array);

  return;
}


# line 192 "aligngenemodel.dy"
NonCodingSimpleModel * create_NonCodingSimpleModel(Probability change)
{
  NonCodingSimpleModel * out;

  out = NonCodingSimpleModel_alloc();

  out->one_off = change * (1.0-change) * (1.0-change) * 3;
  out->two_off = change * change * (1.0-change) * 3;
  out->identical = 1.0 - (out->one_off + out->two_off);

  return out;
}

# line 205 "aligngenemodel.dy"
Probability simple_non_coding_AlignGeneColumn(AlignGeneColumnStore * store,AlignGeneModelParam * agmp)
{
  int i;
  int j;
  int k;
  Probability out = 1.0;
  Probability col = 0.0;
  Probability c   = 0.0;


  for(k=0;k<3;k++) {
      col = 0.0;
      for(j=0;j<4;j++) {
	c= 0.25;
	for(i=0;i<store->len;i++) {
	  if( store->codon[i].seq[2-k] == '-' || store->codon[i].seq[2-k] == '~') {
	    return 0.5;
	  }
	  c *= agmp->dm->prob[j][base_from_char(store->codon[i].seq[2-k])];
	}
	col += c;
      }
      
      out *= col;
  }
  
  return out;
}

# line 234 "aligngenemodel.dy"
Probability weighted_non_coding_AlignGeneColumn(AlignGeneColumnStore * store,AlignGeneModelParam * agmp)
{
  int i;
  int j;
  int k;


  Probability out = 1.0;
  Probability col = 0.0;
  Probability c   = 0.0;


  for(k=0;k<3;k++) {
      col = 0.0;
      for(j=0;j<4;j++) {
	c= 0.25;
	for(i=0;i<store->len;i++) {
	  if( store->codon[i].seq[2-k] == '-'  || store->codon[i].seq[2-k] == '~') {
	    return 0.5;
	  }
	  c *= agmp->dm->prob[j][base_from_char(store->codon[i].seq[2-k])];
	}
	col += c;
      }
      
      out *= col;
  }
  
  return out;


}


# line 268 "aligngenemodel.dy"
Probability weighted_simple_non_coding_AlignGeneColumn(AlignGeneColumnStore * store,AlignGeneModelParam * agmp)
{
  int i;
  int j;
  int k;


  Probability out = 1.0;
  Probability col = 0.0;
  Probability c   = 0.0;


  for(k=0;k<3;k++) {
      col = 0.0;
      for(j=0;j<4;j++) {
	c= 0.25;
	for(i=0;i<store->len;i++) {
	  if( store->codon[i].seq[2-k] == '-' ) {
	    return 0.5;
	  }
	  c *= agmp->dm->prob[j][base_from_char(store->codon[i].seq[2-k])];
	}
	col += c;
      }
      
      out *= col;
  }
  
  return out;


}

# line 301 "aligngenemodel.dy"
Probability simple_coding_AlignGeneColumn(AlignGeneColumnStore * store,AlignGeneModelParam * agmp)
{
  int i;
  int j;
  int codon;
  Probability out = 0.0;
  Probability c = 0.0;
  Probability f;

  assert(store);
  assert(agmp);

  for(j=0;j<26;j++) {
    if( agmp->rm->aminoacid[j] < 0.0000000001 ) {
      continue;
    }
    c = agmp->rm->aminoacid[j];
    for(i=0;i<store->len;i++) {
      if( !isalpha(store->codon[i].seq[0]) ||
	  !isalpha(store->codon[i].seq[1]) ||
	  !isalpha(store->codon[i].seq[2]) ) {
	return 0.0;
      }

      codon = codon_from_seq(store->codon[i].seq);
      if( is_stop_codon(codon,agmp->ct) ) {
	return 0.0;
      }
      f = agmp->protein->comp[j][aminoacid_from_codon(agmp->ct,codon)-'A'];
    
      /*
      fprintf(stdout,"For amino acid %c, match %.4f cum %.4f\n",j+'A',
	      agmp->protein->comp[j][aminoacid_from_codon(agmp->ct,codon)-'A'],
	      c);
      */

      /*     c = c*f;*/

      c = c*f;
      /*      fprintf(stdout,"Refactoring by %f\n",store->codon[i].weight); */
      /* c = c*100; */
    }
    out += c;
  }
    
  /*  fprintf(stdout,"final %.4f\n",out);*/
  return out;
}

# line 350 "aligngenemodel.dy"
boolean fill_reverse_AlignGeneColumn(AlignGeneColumnStore * store,SeqAlign * sa,int forward_codon_end)
{
  int i;

  assert(store);
  assert(sa);
  assert(store->len == sa->len);
  assert(forward_codon_end >= 0 && forward_codon_end < sa->seq[0]->len);

  for(i=0;i<sa->len;i++) {
    store->codon[i].seq[0] = char_complement_base(sa->seq[i]->seq[forward_codon_end]);
    store->codon[i].seq[1] = char_complement_base(sa->seq[i]->seq[forward_codon_end-1]);
    store->codon[i].seq[2] = char_complement_base(sa->seq[i]->seq[forward_codon_end-2]);
    store->codon[i].weight = sa->seq[i]->weight;
  }

  return TRUE;

}

# line 370 "aligngenemodel.dy"
boolean fill_forward_AlignGeneColumn(AlignGeneColumnStore * store,SeqAlign * sa,int codon_end)
{
  int i;

  assert(store);
  assert(sa);
  assert(store->len == sa->len);
  assert(codon_end >= 0 && codon_end < sa->seq[0]->len);

  for(i=0;i<sa->len;i++) {
    store->codon[i].seq[0] = sa->seq[i]->seq[codon_end-2];
    store->codon[i].seq[1] = sa->seq[i]->seq[codon_end-1];
    store->codon[i].seq[2] = sa->seq[i]->seq[codon_end];
    store->codon[i].weight = sa->seq[i]->weight;
  }

  return TRUE;
}

# line 389 "aligngenemodel.dy"
AlignGeneColumnStore * new_empty_AlignGeneColumnStore(int len)
{
  AlignGeneColumnStore * out;

  out = AlignGeneColumnStore_alloc();

  out->codon = calloc(len,sizeof(AlignGeneCodon));
  out->len = len;

  return out;
}

# line 401 "aligngenemodel.dy"
AlignGeneCodon * free_AlignGeneCodon(AlignGeneCodon * c)
{
  ckfree(c);
  return NULL;
}


# line 408 "aligngenemodel.dy"
AlignGeneModel * new_AlignGeneModel(int len)
{
  AlignGeneModel * out;

  out = AlignGeneModel_alloc();

  out->len = len;
  
  out->forward_coding = calloc(len,sizeof(Probability));
  out->reverse_coding = calloc(len,sizeof(Probability));
  out->splice5_forward = calloc(len,sizeof(Probability));
  out->splice3_forward = calloc(len,sizeof(Probability));
  out->splice5_reverse = calloc(len,sizeof(Probability));
  out->splice3_reverse = calloc(len,sizeof(Probability));
  out->change_rate = calloc(len,sizeof(double));

  out->align  = NULL;
  out->anchor = NULL;

  return out;
}

# line 430 "aligngenemodel.dy"
void show_AlignGeneModel(AlignGeneModel * agm,SeqAlign * sal,CodonTable * ct,GenomicRegion * gr,FILE * ofp,AlignGeneModelParam * agmp)
{
  int i;
  int j;
  int * should_show;
  int ex_len;
  int k;
  AlignGeneColumnStore * store;
  int tr_len;
  int phase_score[3];
  int total_phase_0 = 0;
  int phase;
  int * phase_array;
  int pos;
  int count = 0;

  fprintf(stderr,"Gene model at %d with %d\n",gr,gr->len);

  should_show = calloc(agm->len,sizeof(int));
  phase_array = calloc(agm->len,sizeof(int));

  fprintf(stderr,"Gene model at %d with %d\n",gr,gr->len);

  if( gr != NULL ) {
    if( gr->len == 0 ) {
      warn("genomic region given, but no genes on region!");
    }

    for(i=0;i<agm->len;i++) 
      should_show[i] = 0;
  } else {
    for(i=0;i<agm->len;i++) 
      should_show[i] = 1;
  }

 

  if( gr != NULL ) {
    for( i=0;i<gr->len;i++) {
      tr_len = 3;
      for( j=0;j<gr->gene[i]->transcript[0]->ex_len; j++) {
	ex_len = gr->gene[i]->transcript[0]->exon[j]->end - 
	  gr->gene[i]->transcript[0]->exon[j]->start;

	phase_score[0] = 0;
	phase_score[1] = 0;
	phase_score[2] = 0;


	/* deal with splice sites etc */
	pos = gr->gene[i]->start+gr->gene[i]->transcript[0]->exon[j]->start-4;
	should_show[pos] = 1;
	pos = gr->gene[i]->start+gr->gene[i]->transcript[0]->exon[j]->start-3;
	should_show[pos] = 1;
	pos = gr->gene[i]->start+gr->gene[i]->transcript[0]->exon[j]->start-2;
	should_show[pos] = 1;
	pos = gr->gene[i]->start+gr->gene[i]->transcript[0]->exon[j]->start-1;
	should_show[pos] = 1;

	pos = gr->gene[i]->start+gr->gene[i]->transcript[0]->exon[j]->start+ex_len+1;
	should_show[pos] = 1;
	pos = gr->gene[i]->start+gr->gene[i]->transcript[0]->exon[j]->start+ex_len;
	should_show[pos] = 1;
	pos = gr->gene[i]->start+gr->gene[i]->transcript[0]->exon[j]->start+ex_len+2;
	should_show[pos] = 1;
	pos = gr->gene[i]->start+gr->gene[i]->transcript[0]->exon[j]->start+ex_len+3;
	should_show[pos] = 1;
	
	for(k=0;k<ex_len;k++) {
	  phase = ((tr_len -2) %3);
	  pos = gr->gene[i]->start+gr->gene[i]->transcript[0]->exon[j]->start+k;

	  phase_score[phase] += Probability2Score(agm->forward_coding[pos]);
	  phase_array[pos] = phase_score[phase];
	  /*  
	  printf("%d Adding %d in to %d with %f\n",phase,Probability2Score(agm->forward_coding[pos]),phase_score[phase],Score2Bits(phase_score[phase])); 
	  */
	  
	  if( phase == 0 ) {
	    total_phase_0 += Probability2Score(agm->forward_coding[i]);
	    should_show[pos] = 2;
	  } else {
	    should_show[pos] = 1;
	  }
	  tr_len++;
	}
      }
    }
  }

  store = new_empty_AlignGeneColumnStore(sal->len);
  

  for(i=0;i<agm->len;i++) {
    if( should_show[i] == 0 ) {
      continue;
    }
  
    if( i > 0 && should_show[i-1] == 0 ) {
      /* just had a discontinuity */
      fprintf(ofp,"\n ---\n");
    }
    
    fill_forward_AlignGeneColumn(store,sal,i);
    
    fprintf(ofp,"%c %4d %2.4f Bits: %2.4f [%1.6f vs %1.6f cum % 8.1f] %2.4f %2.4f | %2.4f %2.4f %2.4f ",
	    should_show[i] == 2 ? '*' : '-',
	    i,
	    agm->forward_coding[i],
	    Score2Bits(Probability2Score(agm->forward_coding[i])),
	    weighted_coding_AlignGeneColumn(store,agmp),
	    weighted_non_coding_AlignGeneColumn(store,agmp),
	    gr == NULL ? 0.0 : Score2Bits(phase_array[i]),
	    agm->splice5_forward[i],
	    agm->splice3_forward[i],
	    agm->reverse_coding[i],agm->splice5_reverse[i],
	    agm->splice3_reverse[i]
	    );
    for(j=0;j<sal->len;j++) {
      fprintf(ofp,"%c",sal->seq[j]->seq[i]);
    }
    fprintf(ofp,"|");

    for(j=0;j<sal->len;j++) {
      fprintf(ofp,"%c",char_complement_base(sal->seq[j]->seq[i]));
    }

    fprintf(ofp," ");
    if( i > 3 ) {
      for(j=0;j<sal->len;j++) {
	if( sal->seq[j]->seq[i-2] == '~' || sal->seq[j]->seq[i-2] == '-' ||
	    sal->seq[j]->seq[i-1] == '~' || sal->seq[j]->seq[i-1] == '-' ||
	    sal->seq[j]->seq[i] == '~' || sal->seq[j]->seq[i] == '-' ) {
	  fprintf(ofp,"?");
	} else {
	  if( sal->seq[0]->seq[i-2] != sal->seq[1]->seq[i-2] ||
	      sal->seq[0]->seq[i-1] != sal->seq[1]->seq[i-1] ||
	      sal->seq[0]->seq[i] != sal->seq[1]->seq[i]) {
	    fprintf(ofp,"%c",aminoacid_from_seq(ct,sal->seq[j]->seq+i-2));
	  } else {
	    fprintf(ofp,"%c",tolower(aminoacid_from_seq(ct,sal->seq[j]->seq+i-2)));
	  }
	}
      }
      if( agm->forward_coding[i] > 15.0 ) {
	fprintf(ofp,"*");
      }
    }
    fprintf(ofp,"\n");
  }
}

# line 582 "aligngenemodel.dy"
AlignGeneModel * create_AlignGeneModel(SeqAlign * sal,AlignGeneModelParam * agmp)
{
  AlignGeneModel * out;
  int i;
  int j;
  Probability ss5p;
  Probability ss3p;
  AlignGeneColumnStore * store;
  Sequence * rev;

  assert(agmp);
  assert(agmp->ct);
  assert(agmp->protein);
  assert(agmp->dm);


  out = new_AlignGeneModel(sal->seq[0]->len);

  reweight_SeqAlign(sal);

  for(i=0;i<sal->len;i++) {
    fprintf(stdout,"%s weight %f\n",sal->seq[i]->name,sal->seq[i]->weight);
  }

  store = new_empty_AlignGeneColumnStore(sal->len);

  for(i=0;i<sal->seq[0]->len;i++) {
    if( i > 2 ) {

      fill_forward_AlignGeneColumn(store,sal,i);

      /*
      fprintf(stdout,"%d %f vs %f [%f] | %f vs %f %c%c%c %c%c%c\n",i,
	      simple_coding_AlignGeneColumn(store,agmp),
	      weighted_coding_AlignGeneColumn(store,agmp),

	      simple_coding_AlignGeneColumn(store,agmp)/weighted_coding_AlignGeneColumn(store,agmp),

	      simple_non_coding_AlignGeneColumn(store,agmp),
	      non_coding_probability_AlignGeneModel(sal,i,agmp),
	      store->codon[0].seq[0],
	      store->codon[0].seq[1],
	      store->codon[0].seq[2],
	      store->codon[1].seq[0],
	      store->codon[1].seq[1],
	      store->codon[1].seq[2]);
      
      
      */

      out->forward_coding[i] = 
	weighted_coding_AlignGeneColumn(store,agmp)*agmp->total_weight/
	(weighted_non_coding_AlignGeneColumn(store,agmp));
      
    }
    
    ss5p = 1.0;
    ss3p = 1.0;
    for(j=0;j<1;j++) {
      if( sal->seq[j]->seq[i] == 'G' && sal->seq[j]->seq[i+1] == 'T' && i < sal->seq[j]->len-2) 
	ss5p *= prob_SpliceSiteProb(agmp->ss5,sal->seq[j],i);
      else 
	ss5p = 0.00000001;

      if( sal->seq[j]->seq[i] == 'G' && sal->seq[j]->seq[i-1] == 'A' && i > 2 ) 
	ss3p *= prob_SpliceSiteProb(agmp->ss3,sal->seq[j],i);
      else
	ss3p = 0.000000001;
    }
    if( ss5p > 4 ) {
      ss5p = 4.0;
    }
    if( ss3p > 4 ) {
      ss3p = 4.0;
    }

    out->splice5_forward[i] = ss5p;
    out->splice3_forward[i] = ss3p;
  }

  rev = reverse_complement_Sequence(sal->seq[0]);

  for(i=0;i<sal->seq[0]->len;i++) {
    if( i > 2 ) {
      fill_reverse_AlignGeneColumn(store,sal,i);

      out->reverse_coding[i] = 
	weighted_coding_AlignGeneColumn(store,agmp)*agmp->total_weight/
	(weighted_non_coding_AlignGeneColumn(store,agmp));
    }
    
    if( rev->seq[rev->len-i] == 'G' && rev->seq[rev->len-i-1] == 'T' && i > 2) 
      out->splice5_reverse[i] = prob_SpliceSiteProb(agmp->ss5,rev,rev->len-i);
    else
      out->splice5_reverse[i] = 0.00000001;

    if( rev->seq[rev->len-i] == 'G' && rev->seq[rev->len-i+1] == 'A' && i < rev->len -2) 
      out->splice3_reverse[i] = prob_SpliceSiteProb(agmp->ss3,rev,rev->len-i);
    else 
      out->splice3_reverse[i] = 0.00000001;


    if( out->splice5_reverse[i] > 4 ) {
      out->splice5_reverse[i] = 4.0;
    } 

    if( out->splice3_reverse[i] > 4 ) {
      out->splice3_reverse[i] = 4.0;
    } 

    /*
    out->splice5_reverse[i] = 1.0;
    out->splice3_reverse[i] = 1.0;
	*/
  }

  free_AlignGeneColumnStore(store);
  free_Sequence(rev);
  
  window_AlignGeneModel(sal,out,agmp);
  
  out->align = hard_link_SeqAlign(sal);
  out->anchor = hard_link_Sequence(sal->seq[0]);

  return out;
}

# line 709 "aligngenemodel.dy"
Probability prob_SpliceSiteProb(SpliceSiteProb * ssp,Sequence * seq,int pos)
{
  assert(ssp);
  assert(seq);
 
  if( pos<= ssp->offset ) {
    return 0.0;
  }

  return prob_pwmDNA_string(ssp->pwm,seq->seq-ssp->offset+1+pos);
}




# line 724 "aligngenemodel.dy"
Probability non_coding_probability_AlignGeneModel(SeqAlign * sal,int codon_end_pos,AlignGeneModelParam * agmp)
{
  int i;
  int j;
  int k;
  Probability out = 1.0;
  Probability col = 0.0;
  Probability c   = 0.0;


  for(k=0;k<3;k++) {
    for(j=0;j<4;j++) {
      c= 0.25;
      for(i=0;i<sal->len;i++) {
	if( sal->seq[i]->seq[codon_end_pos-k] == '-' ) {
	  return 0.5;
	}
	c *= agmp->dm->prob[j][base_from_char(sal->seq[i]->seq[codon_end_pos-k])];
      }
      col += c;
    }
    
    out *= col;
  }

  return out;

}

# line 753 "aligngenemodel.dy"
Probability coding_probability_AlignGeneModel(SeqAlign * sal,int codon_end_pos,AlignGeneModelParam * agmp)
{
  int i;
  int j;
  int codon;
  Probability out = 0.0;
  Probability c = 0.0;
  Probability f;
  
  assert(sal);
  assert(agmp);
  assert(agmp->rm);
  assert(agmp->protein);

  /*
  fprintf(stdout,"Testing %c%c%c [%c] vs %c%c%c [%c]\n",
	  sal->seq[0]->seq[codon_end_pos-2],
	  sal->seq[0]->seq[codon_end_pos-1],
	  sal->seq[0]->seq[codon_end_pos],
	  aminoacid_from_codon(agmp->ct,codon_from_seq(sal->seq[0]->seq+codon_end_pos-2)),
	  sal->seq[1]->seq[codon_end_pos-2],
	  sal->seq[1]->seq[codon_end_pos-1],
	  sal->seq[1]->seq[codon_end_pos],
	  aminoacid_from_codon(agmp->ct,codon_from_seq(sal->seq[1]->seq+codon_end_pos-2)) );
  */

  for(j=0;j<26;j++) {
    if( agmp->rm->aminoacid[j] < 0.0000000001 ) {
      continue;
    }
    c = agmp->rm->aminoacid[j];
    for(i=0;i<sal->len;i++) {
      if( !isalpha(sal->seq[i]->seq[codon_end_pos-2]) ||
	  !isalpha(sal->seq[i]->seq[codon_end_pos-1]) ||
	  !isalpha(sal->seq[i]->seq[codon_end_pos]) ) {
	return 0.0;
      }

      codon = codon_from_seq(sal->seq[i]->seq+codon_end_pos-2);
      if( is_stop_codon(codon,agmp->ct) ) {
	return 0.0;
      }
      f = agmp->protein->comp[j][aminoacid_from_codon(agmp->ct,codon)-'A'];
    
      /*
      fprintf(stdout,"For amino acid %c, match %.4f cum %.4f\n",j+'A',
	      agmp->protein->comp[j][aminoacid_from_codon(agmp->ct,codon)-'A'],
	      c);
      */

      c = c*f;
    }
    out += c;
  }
    
  /*  fprintf(stdout,"final %.4f\n",out);*/
  return out;

}


# line 814 "aligngenemodel.dy"
AlignGeneModelParam * new_AlignGeneModelParam_from_argv(int * argc,char ** argv)
{
  AlignGeneModelParam * out;
  CompProb * cp;
  DnaProbMatrix * dm;

  GeneModelParam * gmp;
  GeneStats * gs;

  char * temp;
  double temp_prob;

  out = AlignGeneModelParam_alloc();

  if( (temp = strip_out_assigned_argument(argc,argv,"am_codon")) == NULL ) {
    temp = "codon.table";
  }
  out->ct = read_CodonTable_file(temp);
  assert(out->ct);

  if( (temp = strip_out_assigned_argument(argc,argv,"am_protein")) == NULL ) {
    temp = "wag85";
  }
  out->protein = read_Blast_file_CompProb(temp);
  assert(out->protein);

  temp_prob = 0.85;
  strip_out_float_argument(argc,argv,"am_dna",&temp_prob);

  out->dm = DnaProbMatrix_from_match(0.85,NMaskType_VARIABLE);

  gmp = new_GeneModelParam_from_argv(argc,argv);

  if((gs=GeneStats_from_GeneModelParam(gmp)) == NULL ) {
    fatal("Could not build gene stats");
  }

  out->gs = gs;
  out->rm      = default_RandomModel();
  
  out->ss5 = SpliceSiteProb_alloc();
  out->ss5->pwm = pwmDNA_from_SeqAlign(gs->splice5,1);
 
  out->ss5->offset = gs->splice5_offset;
  fold_randommodel_pwmDNA(out->ss5->pwm,gs->rnd);

  out->ss3 = SpliceSiteProb_alloc();
  out->ss3->pwm = pwmDNA_from_SeqAlign(gs->splice3,1);

  out->ss3->offset = gs->splice3_offset;
  fold_randommodel_pwmDNA(out->ss3->pwm,gs->rnd);


  out->nonanchor_stop = 0.001;
  strip_out_float_argument(argc,argv,"am_nonanchor_stop_pen",&out->nonanchor_stop);

  out->tolerate_nonanchor_stops = 1;
  strip_out_boolean_def_argument(argc,argv,"am_nonanchor_stop",&out->tolerate_nonanchor_stops);

  out->coding_window_thres = 3.0;
  strip_out_float_argument(argc,argv,"am_coding_bonus_thres",&out->coding_window_thres);
  

  out->noncoding_window_thres = 2.0;
  strip_out_float_argument(argc,argv,"am_noncoding_bonus_thres",&out->noncoding_window_thres);


  out->coding_window_bonus = 4.0;
  strip_out_float_argument(argc,argv,"am_coding_bonus_factor",&out->coding_window_bonus);

  out->noncoding_window_pen = 0.25;
  strip_out_float_argument(argc,argv,"am_noncoding_bonus_factor",&out->noncoding_window_pen);


  out->total_weight = 0.5;
  strip_out_float_argument(argc,argv,"am_overall_coding",&out->total_weight);


  return out;

}

# line 896 "aligngenemodel.dy"
void show_help_AlignGeneModelParam(FILE * ofp)
{
  fprintf(ofp,"Align Gene Model Parameters\n");
  fprintf(ofp,"  -am_codon [codon.table]    codon table\n");
  fprintf(ofp,"  -am_protein [wag85]        protein comparison probabilities\n");
  fprintf(ofp,"  -am_dna [0.85]             DNA non coding model match probability\n");
  fprintf(ofp,"  -[no]am_nonanchor_stop     allow stops in non anchor sequences [true]\n");
  fprintf(ofp,"  -am_nonanchor_stop_pen     penalty for non anchor stops when allowed [0.001]\n");
  fprintf(ofp,"  -am_coding_bonus_thres     coding threshold to use conservation window [3.0]\n");
  fprintf(ofp,"  -am_noncoding_bonus_thres  nonc   threshold to use conservation window [2.0]\n");
  fprintf(ofp,"  -am_coding_bonus_factor    coding bonus factor when windowed [4.0]\n");
  fprintf(ofp,"  -am_noncoding_bonus_factor nonc   bonus factor when windowed [0.25]\n");
  fprintf(ofp,"  -am_overall_coding         overall per column \"odds-prior\" for coding [0.5]\n");

  fprintf(ofp,"Gene Model statistics used by AlignModel\n");
  show_help_GeneModelParam(ofp);

  return;
  
}


# line 918 "aligngenemodel.dy"
AlignGeneModelParam * std_AlignGeneModelParam(CompProb * cp,DnaProbMatrix * dm,CodonTable * ct,GeneStats * gs,Probability change)
{
  AlignGeneModelParam * out;

  out = AlignGeneModelParam_alloc();

  assert(cp);
  assert(dm);
  assert(ct);
  assert(gs);

  assert(gs->splice5);
  assert(gs->splice5->len > 1);
  assert(gs->splice5->seq[0]->len > 0);

  assert(gs->splice3);
  assert(gs->splice3->len > 1);
  assert(gs->splice3->seq[0]->len > 0);

  out->ncsm    = create_NonCodingSimpleModel(change);
  out->protein = hard_link_CompProb(cp);
  out->rm      = default_RandomModel();
  out->ct      = hard_link_CodonTable(ct);
  out->dm      = hard_link_DnaProbMatrix(dm);
  
  out->ss5 = SpliceSiteProb_alloc();
  out->ss5->pwm = pwmDNA_from_SeqAlign(gs->splice5,1);
 
  out->ss5->offset = gs->splice5_offset;
  fold_randommodel_pwmDNA(out->ss5->pwm,gs->rnd);

  /*  show_pwmDNA_col(out->ss5->pwm,stdout);
      fprintf(stdout,"\n");
  */
  out->ss3 = SpliceSiteProb_alloc();
  out->ss3->pwm = pwmDNA_from_SeqAlign(gs->splice3,1);



  out->ss3->offset = gs->splice3_offset;
  fold_randommodel_pwmDNA(out->ss3->pwm,gs->rnd);
  /*
  show_pwmDNA_col(out->ss3->pwm,stdout);
  fprintf(stdout,"\n");
  */

  out->nonanchor_stop = 0.001;
  out->coding_window_thres = 3.0;
  out->noncoding_window_thres = 2.0;
  out->coding_window_bonus = 4.0;
  out->noncoding_window_pen = 0.25;
  out->tolerate_nonanchor_stops = 1;
  out->total_weight = 0.5;

  return out;
}

# line 975 "aligngenemodel.dy"
AlignGeneModelScore * AlignGeneModelScore_from_AlignGeneModel(AlignGeneModel * agm)
{
  AlignGeneModelScore * out;

  out = AlignGeneModelScore_alloc();

  out->len = agm->len;
  
  out->forward_coding = calloc(agm->len,sizeof(Score));
  out->reverse_coding = calloc(agm->len,sizeof(Score));
  out->splice5_forward = calloc(agm->len,sizeof(Score));
  out->splice3_forward = calloc(agm->len,sizeof(Score));
  out->splice5_reverse = calloc(agm->len,sizeof(Score));
  out->splice3_reverse = calloc(agm->len,sizeof(Score));

  Probability2Score_move(agm->forward_coding,out->forward_coding,agm->len);
  Probability2Score_move(agm->reverse_coding,out->reverse_coding,agm->len);
  Probability2Score_move(agm->splice5_forward,out->splice5_forward,agm->len);
  Probability2Score_move(agm->splice3_forward,out->splice3_forward,agm->len);
  Probability2Score_move(agm->splice5_reverse,out->splice5_reverse,agm->len);
  Probability2Score_move(agm->splice3_reverse,out->splice3_reverse,agm->len);
  

  out->align = hard_link_SeqAlign(agm->align);
  out->anchor = hard_link_Sequence(agm->anchor);

  return out;
}

# line 1004 "aligngenemodel.dy"
Probability * free_Probability(Probability * p)
{
  ckfree(p);
  return NULL;
}

# line 1010 "aligngenemodel.dy"
Score * free_Score(Score * p)
{
  ckfree(p);
  return NULL;
}
# line 969 "aligngenemodel.c"
/* Function:  hard_link_NonCodingSimpleModel(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [NonCodingSimpleModel *]
 *
 * Return [UNKN ]  Undocumented return value [NonCodingSimpleModel *]
 *
 */
NonCodingSimpleModel * hard_link_NonCodingSimpleModel(NonCodingSimpleModel * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a NonCodingSimpleModel object: passed a NULL object");    
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  NonCodingSimpleModel_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [NonCodingSimpleModel *]
 *
 */
NonCodingSimpleModel * NonCodingSimpleModel_alloc(void) 
{
    NonCodingSimpleModel * out; /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(NonCodingSimpleModel *) ckalloc (sizeof(NonCodingSimpleModel))) == NULL)    {  
      warn("NonCodingSimpleModel_alloc failed ");    
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->identical = 0.0;    
    out->one_off = 0.0;  
    out->two_off = 0.0;  


    return out;  
}    


/* Function:  free_NonCodingSimpleModel(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [NonCodingSimpleModel *]
 *
 * Return [UNKN ]  Undocumented return value [NonCodingSimpleModel *]
 *
 */
NonCodingSimpleModel * free_NonCodingSimpleModel(NonCodingSimpleModel * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a NonCodingSimpleModel obj. Should be trappable");  
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


/* Function:  hard_link_SpliceSiteProb(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [SpliceSiteProb *]
 *
 * Return [UNKN ]  Undocumented return value [SpliceSiteProb *]
 *
 */
SpliceSiteProb * hard_link_SpliceSiteProb(SpliceSiteProb * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a SpliceSiteProb object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  SpliceSiteProb_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [SpliceSiteProb *]
 *
 */
SpliceSiteProb * SpliceSiteProb_alloc(void) 
{
    SpliceSiteProb * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(SpliceSiteProb *) ckalloc (sizeof(SpliceSiteProb))) == NULL)    {  
      warn("SpliceSiteProb_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->pwm = NULL; 
    out->offset = 0; 


    return out;  
}    


/* Function:  free_SpliceSiteProb(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [SpliceSiteProb *]
 *
 * Return [UNKN ]  Undocumented return value [SpliceSiteProb *]
 *
 */
SpliceSiteProb * free_SpliceSiteProb(SpliceSiteProb * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a SpliceSiteProb obj. Should be trappable");    
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
    if( obj->pwm != NULL)    
      free_pwmDNA(obj->pwm);     


    ckfree(obj); 
    return NULL; 
}    


/* Function:  hard_link_AlignGeneColumnStore(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [AlignGeneColumnStore *]
 *
 * Return [UNKN ]  Undocumented return value [AlignGeneColumnStore *]
 *
 */
AlignGeneColumnStore * hard_link_AlignGeneColumnStore(AlignGeneColumnStore * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a AlignGeneColumnStore object: passed a NULL object");    
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  AlignGeneColumnStore_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [AlignGeneColumnStore *]
 *
 */
AlignGeneColumnStore * AlignGeneColumnStore_alloc(void) 
{
    AlignGeneColumnStore * out; /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(AlignGeneColumnStore *) ckalloc (sizeof(AlignGeneColumnStore))) == NULL)    {  
      warn("AlignGeneColumnStore_alloc failed ");    
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->codon = NULL;   
    out->len = 0;    


    return out;  
}    


/* Function:  free_AlignGeneColumnStore(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [AlignGeneColumnStore *]
 *
 * Return [UNKN ]  Undocumented return value [AlignGeneColumnStore *]
 *
 */
AlignGeneColumnStore * free_AlignGeneColumnStore(AlignGeneColumnStore * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a AlignGeneColumnStore obj. Should be trappable");  
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
    if( obj->codon != NULL)  
      free_AlignGeneCodon(obj->codon);   


    ckfree(obj); 
    return NULL; 
}    


/* Function:  hard_link_AlignGeneModel(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [AlignGeneModel *]
 *
 * Return [UNKN ]  Undocumented return value [AlignGeneModel *]
 *
 */
AlignGeneModel * hard_link_AlignGeneModel(AlignGeneModel * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a AlignGeneModel object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  AlignGeneModel_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [AlignGeneModel *]
 *
 */
AlignGeneModel * AlignGeneModel_alloc(void) 
{
    AlignGeneModel * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(AlignGeneModel *) ckalloc (sizeof(AlignGeneModel))) == NULL)    {  
      warn("AlignGeneModel_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->len = 0;    
    out->forward_coding = NULL;  
    out->reverse_coding = NULL;  
    out->splice5_forward = NULL; 
    out->splice3_forward = NULL; 
    out->splice5_reverse = NULL; 
    out->splice3_reverse = NULL; 
    out->align = NULL;   
    out->anchor = NULL;  
    out->change_rate = NULL; 


    return out;  
}    


/* Function:  free_AlignGeneModel(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [AlignGeneModel *]
 *
 * Return [UNKN ]  Undocumented return value [AlignGeneModel *]
 *
 */
AlignGeneModel * free_AlignGeneModel(AlignGeneModel * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a AlignGeneModel obj. Should be trappable");    
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
    if( obj->forward_coding != NULL) 
      free_Probability(obj->forward_coding);     
    if( obj->reverse_coding != NULL) 
      free_Probability(obj->reverse_coding);     
    if( obj->splice5_forward != NULL)    
      free_Probability(obj->splice5_forward);    
    if( obj->splice3_forward != NULL)    
      free_Probability(obj->splice3_forward);    
    if( obj->splice5_reverse != NULL)    
      free_Probability(obj->splice5_reverse);    
    if( obj->splice3_reverse != NULL)    
      free_Probability(obj->splice3_reverse);    
    if( obj->align != NULL)  
      free_SeqAlign(obj->align);     
    if( obj->anchor != NULL) 
      free_Sequence(obj->anchor);    
    if( obj->change_rate != NULL)    
      ckfree(obj->change_rate);  


    ckfree(obj); 
    return NULL; 
}    


/* Function:  hard_link_AlignGeneModelScore(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [AlignGeneModelScore *]
 *
 * Return [UNKN ]  Undocumented return value [AlignGeneModelScore *]
 *
 */
AlignGeneModelScore * hard_link_AlignGeneModelScore(AlignGeneModelScore * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a AlignGeneModelScore object: passed a NULL object"); 
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  AlignGeneModelScore_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [AlignGeneModelScore *]
 *
 */
AlignGeneModelScore * AlignGeneModelScore_alloc(void) 
{
    AlignGeneModelScore * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(AlignGeneModelScore *) ckalloc (sizeof(AlignGeneModelScore))) == NULL)  {  
      warn("AlignGeneModelScore_alloc failed "); 
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->len = 0;    
    out->forward_coding = NULL;  
    out->reverse_coding = NULL;  
    out->splice5_forward = NULL; 
    out->splice3_forward = NULL; 
    out->splice5_reverse = NULL; 
    out->splice3_reverse = NULL; 
    out->align = NULL;   
    out->anchor = NULL;  


    return out;  
}    


/* Function:  free_AlignGeneModelScore(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [AlignGeneModelScore *]
 *
 * Return [UNKN ]  Undocumented return value [AlignGeneModelScore *]
 *
 */
AlignGeneModelScore * free_AlignGeneModelScore(AlignGeneModelScore * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a AlignGeneModelScore obj. Should be trappable");   
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
    if( obj->forward_coding != NULL) 
      free_Score(obj->forward_coding);   
    if( obj->reverse_coding != NULL) 
      free_Score(obj->reverse_coding);   
    if( obj->splice5_forward != NULL)    
      free_Score(obj->splice5_forward);  
    if( obj->splice3_forward != NULL)    
      free_Score(obj->splice3_forward);  
    if( obj->splice5_reverse != NULL)    
      free_Score(obj->splice5_reverse);  
    if( obj->splice3_reverse != NULL)    
      free_Score(obj->splice3_reverse);  
    if( obj->align != NULL)  
      free_SeqAlign(obj->align);     
    if( obj->anchor != NULL) 
      free_Sequence(obj->anchor);    


    ckfree(obj); 
    return NULL; 
}    


/* Function:  hard_link_AlignGeneModelParam(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [AlignGeneModelParam *]
 *
 * Return [UNKN ]  Undocumented return value [AlignGeneModelParam *]
 *
 */
AlignGeneModelParam * hard_link_AlignGeneModelParam(AlignGeneModelParam * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a AlignGeneModelParam object: passed a NULL object"); 
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  AlignGeneModelParam_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [AlignGeneModelParam *]
 *
 */
AlignGeneModelParam * AlignGeneModelParam_alloc(void) 
{
    AlignGeneModelParam * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(AlignGeneModelParam *) ckalloc (sizeof(AlignGeneModelParam))) == NULL)  {  
      warn("AlignGeneModelParam_alloc failed "); 
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->protein = NULL; 
    out->rm = NULL;  
    out->cm = NULL;  
    out->ct = NULL;  
    out->dm = NULL;  
    out->ss5 = NULL; 
    out->ss3 = NULL; 
    out->gs = NULL;  
    out->total_weight = 2;   
    out->tolerate_nonanchor_stops = TRUE;    
    out->nonanchor_stop = 0.0;   
    out->ncsm = NULL;    
    out->coding_window_thres = 0;    
    out->noncoding_window_thres = 0; 
    out->coding_window_bonus = 0;    
    out->noncoding_window_pen = 0;   


    return out;  
}    


/* Function:  free_AlignGeneModelParam(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [AlignGeneModelParam *]
 *
 * Return [UNKN ]  Undocumented return value [AlignGeneModelParam *]
 *
 */
AlignGeneModelParam * free_AlignGeneModelParam(AlignGeneModelParam * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a AlignGeneModelParam obj. Should be trappable");   
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
    if( obj->protein != NULL)    
      free_CompProb(obj->protein);   
    if( obj->rm != NULL) 
      free_RandomModel(obj->rm);     
    if( obj->cm != NULL) 
      free_CodonMapper(obj->cm);     
    if( obj->ct != NULL) 
      free_CodonTable(obj->ct);  
    if( obj->dm != NULL) 
      free_DnaProbMatrix(obj->dm);   
    if( obj->ss5 != NULL)    
      free_SpliceSiteProb(obj->ss5);     
    if( obj->ss3 != NULL)    
      free_SpliceSiteProb(obj->ss3);     
    if( obj->gs != NULL) 
      free_GeneStats(obj->gs);   
    if( obj->ncsm != NULL)   
      free_NonCodingSimpleModel(obj->ncsm);  


    ckfree(obj); 
    return NULL; 
}    



#ifdef _cplusplus
}
#endif
