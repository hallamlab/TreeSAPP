#ifdef _cplusplus
extern "C" {
#endif
#include "dnaprofileengine.h"


# line 66 "dnaprofileengine.dy"
TransFactorSet * TransFactorSet_from_DnaProfileSet(DnaProfileSet * in)
{
  int i;
  char name_buf[512];
  TransFactorSet * out;

  out = TransFactorSet_alloc_len(in->len);

  


}

# line 79 "dnaprofileengine.dy"
DnaProfileNode * balanced_4_Sequence_fasta_stream(FILE * ifp)
{
  Sequence * one;
  Sequence * two;
  Sequence * three;
  Sequence * four;

  DnaProfileNode * leafone;
  DnaProfileNode * leaftwo;
  DnaProfileNode * leafthree;
  DnaProfileNode * leaffour;

  DnaProfileNode * midonetwo;
  DnaProfileNode * midthreefour;

  DnaProfileNode * root;

  one = read_fasta_Sequence(ifp);
  fprintf(stderr,"Got 1 %s\n",one->name);
  two = read_fasta_Sequence(ifp);
  fprintf(stderr,"Got 2 %s\n",two->name);
  three = read_fasta_Sequence(ifp);
  fprintf(stderr,"Got 3 %s\n",three->name);
  four = read_fasta_Sequence(ifp);
  fprintf(stderr,"Got 4 %s\n",four->name);

  assert(four != NULL);

  leafone   = new_leaf_DnaProfileNode(one);
  leaftwo   = new_leaf_DnaProfileNode(two);
  leafthree = new_leaf_DnaProfileNode(three);
  leaffour  = new_leaf_DnaProfileNode(four);

  midonetwo = DnaProfileNode_alloc();
  midonetwo->type  = DnaProfileNode_SET;
  midonetwo->left  = leafone;
  midonetwo->right = leaftwo;

  midthreefour = DnaProfileNode_alloc();
  midthreefour->type  = DnaProfileNode_SET;
  midthreefour->left  = leafthree;
  midthreefour->right = leaffour;

  root = DnaProfileNode_alloc();
  root->type  = DnaProfileNode_SET;
  root->left  = midonetwo;
  root->right = midthreefour;

  return root;
}


# line 131 "dnaprofileengine.dy"
DnaProfileNode * simple_cascade_Sequence_fasta_stream(FILE * ifp)
{
  DnaProfileNode * head;
  DnaProfileNode * leaf;
  DnaProfileNode * first;
  DnaProfileNode * temp;

  Sequence * read;

  read = read_fasta_Sequence(ifp);
  uppercase_Sequence(read);

  if( read == NULL ) {
    fatal("Attempting to build a cascade with only no sequences! impossible!");
  }
  first = new_leaf_DnaProfileNode(read);

  read = read_fasta_Sequence(ifp);
  uppercase_Sequence(read);
  if( read == NULL ) {
    fatal("Attempting to build a cascade with only one sequence! impossible!");
  }

  leaf = new_leaf_DnaProfileNode(read);

  head = DnaProfileNode_alloc();
  head->type = DnaProfileNode_SET;
  head->left = first;
  head->right = leaf;


  /* now loop over all remaining sequences */

  while( (read = read_fasta_Sequence(ifp)) != NULL ) {
    uppercase_Sequence(read);
    leaf = new_leaf_DnaProfileNode(read);
    temp = DnaProfileNode_alloc();
    temp->type  = DnaProfileNode_SET;
    temp->left  = head;
    temp->right = leaf;
    
    head = temp;
  }

  return head;
}


# line 179 "dnaprofileengine.dy"
DnaProfileNode * new_leaf_DnaProfileNode(Sequence * seq)
{
  DnaProfileNode * out;

  assert(seq != NULL);

  out = DnaProfileNode_alloc();
  out->type = DnaProfileNode_LEAF;
  out->leaf = seq;

  return out;
}





# line 196 "dnaprofileengine.dy"
void populate_DnaProfileNode_from_root(DnaProfileNode * root,DnaProfileEnginePara * dpep)
{
  assert(root != NULL);
  if( root->left->type == DnaProfileNode_SET ) {
    populate_DnaProfileNode_from_root(root->left,dpep);
  }
  if( root->right->type == DnaProfileNode_SET ) {
    populate_DnaProfileNode_from_root(root->right,dpep);
  }

  /* left and right now populated */

  root->set = join_two_DnaProfileNode(root->left,root->right,dpep);
}

# line 211 "dnaprofileengine.dy"
DnaProfileSet * join_two_DnaProfileNode(DnaProfileNode * left,DnaProfileNode * right,DnaProfileEnginePara * dpep)
{
  /* big switch around the types of left and right */

  fprintf(stderr,"Entering join with %d vs %d\n",left->type,right->type);

  if( left->type == DnaProfileNode_LEAF && right->type == DnaProfileNode_LEAF ) {
    assert(left->leaf);
    assert(right->leaf);
    return DnaProfileSet_from_leaf_leaf(left->leaf,right->leaf,dpep);
  }

  if( left->type == DnaProfileNode_SET && right->type == DnaProfileNode_LEAF ) {
    assert(left->set);
    assert(right->leaf);
    return DnaProfileSet_from_leaf_node(right->leaf,left->set,dpep);
  }

  if( left->type == DnaProfileNode_LEAF && right->type == DnaProfileNode_SET ) {
    assert(left->leaf);
    assert(right->set);
    return DnaProfileSet_from_leaf_node(left->leaf,right->set,dpep);
  }

  if( left->type == DnaProfileNode_SET && right->type == DnaProfileNode_SET ) {
    assert(left->set);
    assert(right->set);
    return DnaProfileSet_from_node_node(left->set,right->set,dpep);
  }


  fatal("Should not get here. Weird no leaf/node case");

  return NULL;

}




# line 251 "dnaprofileengine.dy"
DnaProfileEnginePara * new_DnaProfileEnginePara_from_argv(int * argc,char ** argv)
{
  DnaProfileEnginePara * out;

  out = DnaProfileEnginePara_alloc();

  out->dpri = new_DPRunImpl_from_argv(argc,argv);

  out->setpara = new_LocalCisHitSetPara_from_argv(argc,argv);
  out->lchs    = standard_LocalCisHitScore(NMaskType_VARIABLE);

  out->rm = RandomModelDNA_std();

  out->pseudo = 0.5;
  out->open_unmatched = 0.001;
  out->ext_unmatched = 0.8;
  out->gap_unmatched = 0.5;
  out->seq_id = 0.8;
  out->m2i = 0.1;
  out->m2d = 0.1;
  out->i2i = 0.8;
  out->d2d = 0.8;
  out->min_seq_prof = 400;

  strip_out_float_argument(argc,argv,"dnap_pseudo",&out->pseudo);
  strip_out_float_argument(argc,argv,"dnap_open_un",&out->open_unmatched);
  strip_out_float_argument(argc,argv,"dnap_ext_un",&out->ext_unmatched);
  strip_out_float_argument(argc,argv,"dnap_gap_un",&out->ext_unmatched);
  strip_out_float_argument(argc,argv,"dnap_seq_self",&out->seq_id);
  strip_out_float_argument(argc,argv,"dnap_m2i",&out->m2i);
  strip_out_float_argument(argc,argv,"dnap_m2d",&out->m2d);
  strip_out_float_argument(argc,argv,"dnap_i2i",&out->i2i);
  strip_out_float_argument(argc,argv,"dnap_d2d",&out->d2d);
  strip_out_integer_argument(argc,argv,"dnap_min_seq_prof",&out->min_seq_prof);
  
  return out;

}

# line 290 "dnaprofileengine.dy"
void show_help_DnaProfileEnginePara(FILE * ofp)
{
  fprintf(ofp,"DnaProfile build/matching parameters\n");
  fprintf(ofp,"  -dnap_pseudo   [0.1]  pseudo count used in profile construction\n");
  fprintf(ofp,"  -dnap_open_un  [0.4]  unmatched probability open\n");
  fprintf(ofp,"  -dnap_ext_un   [0.95] unmatched extend probability\n");
  fprintf(ofp,"  -dnap_ext_un   [0.5]  unmatched gap    probability\n");
  fprintf(ofp,"  -dnap_seq_self [0.8]  %% identity for pure sequence matching\n");
  fprintf(ofp,"  -dnap_m2i      [0.1]  Match 2 insert transition in dnaprofiles\n");
  fprintf(ofp,"  -dnap_m2d      [0.1]  Match 2 delete transitions in dnaprofiles\n");
  fprintf(ofp,"  -dnap_i2i      [0.8]  Insert 2 Insert transition in dnaprofiles\n");
  fprintf(ofp,"  -dnap_d2d      [0.8]  Delete 2 Delete transition in dnaprofiles\n");
  fprintf(ofp,"  -dnap_min_seq_prof [400] minimum score for sequence profile matching\n");

  fprintf(ofp,"Local CisHit para for sequence to sequnence matching\n");
  show_help_LocalCisHitSetPara(ofp);

  show_help_DPRunImpl(ofp);

}

# line 311 "dnaprofileengine.dy"
DnaProfileSet * filter_DnaProfileSet(DnaProfileSet * in,int min_length,int min_score)
{
  int i;
  DnaProfileSet * out;

  out = DnaProfileSet_alloc_std();

  for(i=0;i<in->len;i++) {
    if( in->dnap[i]->sa->seq[0]->len <= min_length ) {
      continue;
    }

    add_DnaProfileSet(out,hard_link_DnaProfile(in->dnap[i]));
  }

  return out;
		      
}

# line 330 "dnaprofileengine.dy"
DnaProfileSet * DnaProfileSet_from_leaf_leaf(Sequence * one,Sequence * two,DnaProfileEnginePara * dpep)
{
  DnaProfileSet * out;
  DnaMatrix * dm;
  DnaProbMatrix * dmp;
  PairwiseShortDna * psd;
  LocalCisHitSet * set;
  Sequence * two_rev;
  DnaProfile * dp;
  SeqAlign * sa;

  Sequence * temp1;
  Sequence * temp2;

  char * temp_seq1;
  char * temp_seq2;
  
  int unmatched;
  int seq1_i,seq2_i;

  AlnColumn * alc;
  int i;

  two_rev = reverse_complement_Sequence(two);

  
  dmp = DnaProbMatrix_from_match(0.65,NMaskType_BANNED);  
  assert(dmp);
  flat_null_DnaProbMatrix(dmp);  

  dm = DnaMatrix_from_DnaProbMatrix(dmp);

  show_DnaMatrix(dm,stderr);

  psd = query_to_reverse_target(one,two,dm,0,one->len,0,two->len);

  
  set = make_LocalCisHitSet(one,two,two_rev,psd->forward,psd->reverse,dpep->setpara,dpep->lchs,NULL,NULL,NULL,NULL,0,dpep->dpri);

  temp_seq1 = calloc(one->len > two->len ? one->len : two->len,sizeof(char));
  temp_seq2 = calloc(one->len > two->len ? one->len : two->len,sizeof(char));
  
  out = DnaProfileSet_alloc_std();

  for(i=0;i<set->len;i++) {
    unmatched = 1;
    sa = NULL;

    /*
     * Main loop over DBA style alignment. We need to make one
     * DnaProfile per matching block, which are separated by unmatched
     * blocks. Could potentially be no blocks.
     *
     * Extra annoyance provided by the "wrong" convention being used in
     * DBA alignments, meaning that "inserts" label the "sequence" containing
     * strand, not the non-sequence containing strand. Stupid, but dbadisplay
     * uses this convention, so if we changed, would have to fix lots of exisiting
     * code. Not ideal.
     *
     */


    for(alc=set->lch[i]->alb->start;alc != NULL;alc=alc->next) {
      
      /* hitting an unmatched block */
      if( unmatched == 0 && (strcmp(alc->alu[0]->text_label,"UM") == 0 ||
	  strcmp(alc->alu[0]->text_label,"UI") == 0 || strcmp(alc->alu[0]->text_label,"END") == 0) ) {
	/* if we have an alignment, put it away now */
	if( sa != NULL ) {
	  temp_seq1[seq1_i] = '\0';
	  temp_seq2[seq2_i] = '\0';

	  temp1 = Sequence_from_static_memory(one->name,temp_seq1);
	  temp2 = Sequence_from_static_memory(two->name,temp_seq2);
	  
	  add_SeqAlign(sa,temp1);
	  add_SeqAlign(sa,temp2);
	  
	  dp = naive_DnaProfile_from_SeqAlign(sa,0.15,0.1,0.1,0.8,0.8);
	  fold_RandomModel_DnaProfile(dp,dpep->rm);

	  add_DnaProfileSet(out,dp);
	  free_SeqAlign(sa); /* hard linked inside DP */
	  sa = NULL;
	}

	continue;
      } else if( unmatched == 1 && (strstartcmp(alc->alu[0]->text_label,"MM") == 0 ||
				    strstartcmp(alc->alu[0]->text_label,"MI") == 0 ) ) {
	unmatched = 0;
	  
	sa = SeqAlign_alloc_len(2);
	seq1_i = 0;
	seq2_i = 0;
      }

      /* only if we are in a matched block */
      if( unmatched == 0 ) {
	/* Bloody twisted DBA convention - Niclas has alot to answer for.
	   Evil stuff -- MI is on the wrong strand! */
	if( strstartcmp(alc->alu[0]->text_label,"MI") == 0 ) {
	  /* means 0 has sequence, other has gap */
	  temp_seq1[seq1_i++] = one->seq[alc->alu[0]->end];
	  temp_seq2[seq2_i++] = '-';
	} else if ( strstartcmp(alc->alu[1]->text_label,"MI") == 0 ) {
	  temp_seq1[seq1_i++] = '-';
	  temp_seq2[seq2_i++] = two->seq[alc->alu[1]->end];
	} else if ( strstartcmp(alc->alu[0]->text_label,"MM") == 0 &&
		    strstartcmp(alc->alu[1]->text_label,"MM") == 0 ) {
	  temp_seq1[seq1_i++] = one->seq[alc->alu[0]->end];
	  temp_seq2[seq2_i++] = two->seq[alc->alu[1]->end];
	} else {
	  warn("Impossible label pair reached in matched block local cis hit stuff, %s,%s",alc->alu[0]->text_label,alc->alu[1]->text_label);
	}
      }
	  
    }
  }


  free(temp_seq1);
  free(temp_seq2);
  free_PairwiseShortDna(psd);
  free_LocalCisHitSet(set);
  free_DnaMatrix(dm);
  free_DnaProbMatrix(dmp);

  return out;
}



# line 462 "dnaprofileengine.dy"
DnaProfileSet * DnaProfileSet_from_node_node(DnaProfileSet * one,DnaProfileSet * two,DnaProfileEnginePara * dpep)
{
  DnaProfile * new_dnap;
  DnaProfileSet * out;
  DnaProfileMatchPairSet * dpmps;
  SeqAlign * sa;
  int i;
  int j;

  dpmps = DnaProfileMatchPairSet_alloc_std();

  for(i=0;i<one->len;i++) {
    for(j=0;j<two->len;j++) {
	add_DnaProfileMatchPairSet(dpmps,DnaProfileMatchPair_from_DnaProfile(one->dnap[i],two->dnap[j],dpep));
    }
  }

  sort_DnaProfileMatchPairSet_by_score(dpmps);


  out = DnaProfileSet_alloc_std();

  for(i=0;i<dpmps->len;i++) {
    /* check this profile has not already been used */
    /* not done yet */

    if( dpmps->pair[i]->score < dpep->min_seq_prof ) {
      fprintf(stderr,"Warning... rejecting match due to score %d vs %d\n",dpmps->pair[i]->score,dpep->min_seq_prof);
      break;
    }
    

    sa = merged_SeqAlign(dpmps->pair[i]->query,dpmps->pair[i]->target,dpmps->pair[i]->alb);

    fprintf(stderr,"Node/Node Accepting score at %d length %d\n",dpmps->pair[i]->score,sa->seq[0]->len);

    new_dnap = naive_DnaProfile_from_SeqAlign(sa,dpep->pseudo,dpep->m2i,dpep->m2d,dpep->i2i,dpep->d2d);

    assert(new_dnap != NULL);
    /* need to log-odds dnap here */

    fold_RandomModel_DnaProfile(new_dnap,dpep->rm);


    add_DnaProfileSet(out,new_dnap);
  }

    
  fprintf(stderr,"Returing %d profiles\n",out->len);
  return out;
}

# line 514 "dnaprofileengine.dy"
DnaProfileSet * DnaProfileSet_from_leaf_node(Sequence * one,DnaProfileSet * two,DnaProfileEnginePara * dpep)
{
  DnaProfileSet * out;
  DnaProfile * dnap;
  DnaProfile * dnapr;
  DnaProfileMatchPairSet * dpmps;
  Sequence * rev;
  SeqAlign * sa;
  DnaProfile * new_dnap;
  int i;
  int j;

  
  dpmps = DnaProfileMatchPairSet_alloc_std();

  out = DnaProfileSet_alloc_std();

  rev = reverse_complement_Sequence(one);

  dnap = naive_DnaProfile_from_Sequence(one,dpep->seq_id,dpep->m2i,dpep->m2d,dpep->i2i,dpep->d2d);
  dnapr = naive_DnaProfile_from_Sequence(rev,dpep->seq_id,dpep->m2i,dpep->m2d,dpep->i2i,dpep->d2d);

  fold_RandomModel_DnaProfile(dnap,dpep->rm);
  fold_RandomModel_DnaProfile(dnapr,dpep->rm);


  for(i=0;i<two->len;i++) {
    fprintf(stderr,"Processing %d\n",i);
    add_DnaProfileMatchPairSet(dpmps,DnaProfileMatchPair_from_DnaProfile(dnap,two->dnap[i],dpep));    
    add_DnaProfileMatchPairSet(dpmps,DnaProfileMatchPair_from_DnaProfile(dnapr,two->dnap[i],dpep));    
  }

  fprintf(stderr,"Sorting....\n");

  sort_DnaProfileMatchPairSet_by_score(dpmps);

  for(i=0;i<dpmps->len;i++) {
    /* check this profile has not already been used */
    /* not done yet */
    
    if( dpmps->pair[i]->score < dpep->min_seq_prof ) {
      fprintf(stderr,"Warning... rejecting match due to score %d vs %d\n",dpmps->pair[i]->score,dpep->min_seq_prof);
      break;    }
    
    sa = merged_SeqAlign(dpmps->pair[i]->query,dpmps->pair[i]->target,dpmps->pair[i]->alb);
    

    new_dnap = naive_DnaProfile_from_SeqAlign(sa,dpep->pseudo,dpep->m2i,dpep->m2d,dpep->i2i,dpep->d2d);

    /* need to log-odds dnap here */

    fold_RandomModel_DnaProfile(new_dnap,dpep->rm);


    add_DnaProfileSet(out,new_dnap);
  }

    
  fprintf(stderr,"Freeing DNA profiles...\n");

  free_DnaProfile(dnap);
  free_DnaProfile(dnapr);

  fprintf(stderr,"Freeing sequences\n");

  free_Sequence(rev);
  


  return out;
}


# line 587 "dnaprofileengine.dy"
DnaProfileMatchPair * DnaProfileMatchPair_from_DnaProfile(DnaProfile * query,DnaProfile * target,DnaProfileEnginePara * dpep)
{
  DnaProfileMatchPair * out;
  DnaProfileScore * query_s;
  DnaProfileScore * target_s;
  DnaProfileMatchScore * match;

  PackAln * pal;


  assert(query != NULL);
  assert(target != NULL);
  /*
  assert(query->len > 4 );
  assert(target->len > 4);
  */

  out = DnaProfileMatchPair_alloc();
  out->query = hard_link_DnaProfile(query);
  out->target = hard_link_DnaProfile(target);


  query_s  = DnaProfileScore_from_DnaProfile(query);
  target_s = DnaProfileScore_from_DnaProfile(target);

  fprintf(stderr,"Matching %d to %d\n",query->len,target->len);

  match= new_ALLR_DnaProfileMatchScore(query,target);
  


  pal = PackAln_bestmemory_DnaProfileMat(query_s,target_s,match,Probability2Score(dpep->open_unmatched),Probability2Score(dpep->ext_unmatched),Probability2Score(dpep->gap_unmatched),NULL,dpep->dpri);
  
  fprintf(stderr,"...Made pal %d\n",pal);
  out->alb = convert_PackAln_to_AlnBlock_DnaProfileMat(pal);
  out->score = pal->score;

  fprintf(stderr,"...freeing pal\n");
  free_PackAln(pal);


  fprintf(stderr,"...freeing match\n");
  free_DnaProfileMatchScore(match);
  fprintf(stderr,"...freeing query\n");
  free_DnaProfileScore(query_s);
  fprintf(stderr,"...freeing target\n");
  free_DnaProfileScore(target_s);


  return out;
}


# line 640 "dnaprofileengine.dy"
void sort_DnaProfileMatchPairSet_by_score(DnaProfileMatchPairSet * set)
{
  sort_DnaProfileMatchPairSet(set,compare_DnaProfileMatchPair);
}

# line 645 "dnaprofileengine.dy"
int compare_DnaProfileMatchPair(DnaProfileMatchPair * one,DnaProfileMatchPair * two)
{
  return two->score - one->score;
}


# line 651 "dnaprofileengine.dy"
void show_DnaProfileSet(DnaProfileSet * dnaps,RandomModelDNA * rm,FILE * ofp)
{
  int i;
  
  for(i=0;i<dnaps->len;i++) {
    show_DnaProfile(dnaps->dnap[i],rm,ofp);
  }
}


# line 617 "dnaprofileengine.c"
/* Function:  hard_link_DnaProfileEnginePara(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [DnaProfileEnginePara *]
 *
 * Return [UNKN ]  Undocumented return value [DnaProfileEnginePara *]
 *
 */
DnaProfileEnginePara * hard_link_DnaProfileEnginePara(DnaProfileEnginePara * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a DnaProfileEnginePara object: passed a NULL object");    
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  DnaProfileEnginePara_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [DnaProfileEnginePara *]
 *
 */
DnaProfileEnginePara * DnaProfileEnginePara_alloc(void) 
{
    DnaProfileEnginePara * out; /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(DnaProfileEnginePara *) ckalloc (sizeof(DnaProfileEnginePara))) == NULL)    {  
      warn("DnaProfileEnginePara_alloc failed ");    
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->dpri = NULL;    
    out->rm = NULL;  
    out->setpara = NULL; 
    out->lchs = NULL;    
    out->pseudo = 0.0;   
    out->open_unmatched = 0.0;   
    out->ext_unmatched = 0.0;    
    out->gap_unmatched = 0.0;    
    out->seq_id = 0.0;   
    out->m2i = 0.0;  
    out->m2d = 0.0;  
    out->i2i = 0.0;  
    out->d2d = 0.0;  
    out->min_seq_prof = 0;   


    return out;  
}    


/* Function:  free_DnaProfileEnginePara(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [DnaProfileEnginePara *]
 *
 * Return [UNKN ]  Undocumented return value [DnaProfileEnginePara *]
 *
 */
DnaProfileEnginePara * free_DnaProfileEnginePara(DnaProfileEnginePara * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a DnaProfileEnginePara obj. Should be trappable");  
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
    if( obj->dpri != NULL)   
      free_DPRunImpl(obj->dpri);     
    if( obj->rm != NULL) 
      free_RandomModelDNA(obj->rm);  
    if( obj->setpara != NULL)    
      free_LocalCisHitSetPara(obj->setpara);     
    if( obj->lchs != NULL)   
      free_LocalCisHitScore(obj->lchs);  


    ckfree(obj); 
    return NULL; 
}    


/* Function:  swap_DnaProfileSet(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_DnaProfileSet
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [DnaProfile **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_DnaProfileSet(DnaProfile ** list,int i,int j)  
{
    DnaProfile * temp;   
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_DnaProfileSet(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_DnaProfileSet which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [DnaProfile **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_DnaProfileSet(DnaProfile ** list,int left,int right,int (*comp)(DnaProfile * ,DnaProfile * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_DnaProfileSet(list,left,(left+right)/2);    
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_DnaProfileSet (list,++last,i);  
      }  
    swap_DnaProfileSet (list,left,last); 
    qsort_DnaProfileSet(list,left,last-1,comp);  
    qsort_DnaProfileSet(list,last+1,right,comp); 
}    


/* Function:  sort_DnaProfileSet(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_DnaProfileSet
 *
 *
 * Arg:         obj [UNKN ] Object containing list [DnaProfileSet *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_DnaProfileSet(DnaProfileSet * obj,int (*comp)(DnaProfile *, DnaProfile *)) 
{
    qsort_DnaProfileSet(obj->dnap,0,obj->len-1,comp);    
    return;  
}    


/* Function:  expand_DnaProfileSet(obj,len)
 *
 * Descrip:    Really an internal function for add_DnaProfileSet
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [DnaProfileSet *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_DnaProfileSet(DnaProfileSet * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_DnaProfileSet called with no need");  
      return TRUE;   
      }  


    if( (obj->dnap = (DnaProfile ** ) ckrealloc (obj->dnap,sizeof(DnaProfile *)*len)) == NULL)   {  
      warn("ckrealloc failed for expand_DnaProfileSet, returning FALSE");    
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_DnaProfileSet(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [DnaProfileSet *]
 * Arg:        add [OWNER] Object to add to the list [DnaProfile *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_DnaProfileSet(DnaProfileSet * obj,DnaProfile * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_DnaProfileSet(obj,obj->len + DnaProfileSetLISTLENGTH) == FALSE) 
        return FALSE;    
      }  


    obj->dnap[obj->len++]=add;   
    return TRUE; 
}    


/* Function:  flush_DnaProfileSet(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [DnaProfileSet *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_DnaProfileSet(DnaProfileSet * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->dnap[i] != NULL)  {  
        free_DnaProfile(obj->dnap[i]);   
        obj->dnap[i] = NULL; 
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  DnaProfileSet_alloc_std(void)
 *
 * Descrip:    Equivalent to DnaProfileSet_alloc_len(DnaProfileSetLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [DnaProfileSet *]
 *
 */
DnaProfileSet * DnaProfileSet_alloc_std(void) 
{
    return DnaProfileSet_alloc_len(DnaProfileSetLISTLENGTH); 
}    


/* Function:  DnaProfileSet_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [DnaProfileSet *]
 *
 */
DnaProfileSet * DnaProfileSet_alloc_len(int len) 
{
    DnaProfileSet * out;/* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = DnaProfileSet_alloc()) == NULL)    
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->dnap = (DnaProfile ** ) ckcalloc (len,sizeof(DnaProfile *))) == NULL)   {  
      warn("Warning, ckcalloc failed in DnaProfileSet_alloc_len");   
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_DnaProfileSet(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [DnaProfileSet *]
 *
 * Return [UNKN ]  Undocumented return value [DnaProfileSet *]
 *
 */
DnaProfileSet * hard_link_DnaProfileSet(DnaProfileSet * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a DnaProfileSet object: passed a NULL object");   
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  DnaProfileSet_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [DnaProfileSet *]
 *
 */
DnaProfileSet * DnaProfileSet_alloc(void) 
{
    DnaProfileSet * out;/* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(DnaProfileSet *) ckalloc (sizeof(DnaProfileSet))) == NULL)  {  
      warn("DnaProfileSet_alloc failed ");   
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->dnap = NULL;    
    out->len = out->maxlen = 0;  


    return out;  
}    


/* Function:  free_DnaProfileSet(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [DnaProfileSet *]
 *
 * Return [UNKN ]  Undocumented return value [DnaProfileSet *]
 *
 */
DnaProfileSet * free_DnaProfileSet(DnaProfileSet * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a DnaProfileSet obj. Should be trappable"); 
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
    if( obj->dnap != NULL)   {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->dnap[i] != NULL)    
          free_DnaProfile(obj->dnap[i]); 
        }  
      ckfree(obj->dnap); 
      }  


    ckfree(obj); 
    return NULL; 
}    


/* Function:  hard_link_DnaProfileNode(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [DnaProfileNode *]
 *
 * Return [UNKN ]  Undocumented return value [DnaProfileNode *]
 *
 */
DnaProfileNode * hard_link_DnaProfileNode(DnaProfileNode * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a DnaProfileNode object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  DnaProfileNode_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [DnaProfileNode *]
 *
 */
DnaProfileNode * DnaProfileNode_alloc(void) 
{
    DnaProfileNode * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(DnaProfileNode *) ckalloc (sizeof(DnaProfileNode))) == NULL)    {  
      warn("DnaProfileNode_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->type = 0;   
    out->leaf = NULL;    
    out->set = NULL; 
    out->left = NULL;    
    out->right = NULL;   


    return out;  
}    


/* Function:  free_DnaProfileNode(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [DnaProfileNode *]
 *
 * Return [UNKN ]  Undocumented return value [DnaProfileNode *]
 *
 */
DnaProfileNode * free_DnaProfileNode(DnaProfileNode * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a DnaProfileNode obj. Should be trappable");    
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
    if( obj->leaf != NULL)   
      free_Sequence(obj->leaf);  
    if( obj->set != NULL)    
      free_DnaProfileSet(obj->set);  
    if( obj->left != NULL)   
      free_DnaProfileNode(obj->left);    
    if( obj->right != NULL)  
      free_DnaProfileNode(obj->right);   
    /* obj->parent is linked in */ 


    ckfree(obj); 
    return NULL; 
}    


/* Function:  hard_link_DnaProfileMatchPair(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [DnaProfileMatchPair *]
 *
 * Return [UNKN ]  Undocumented return value [DnaProfileMatchPair *]
 *
 */
DnaProfileMatchPair * hard_link_DnaProfileMatchPair(DnaProfileMatchPair * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a DnaProfileMatchPair object: passed a NULL object"); 
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  DnaProfileMatchPair_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [DnaProfileMatchPair *]
 *
 */
DnaProfileMatchPair * DnaProfileMatchPair_alloc(void) 
{
    DnaProfileMatchPair * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(DnaProfileMatchPair *) ckalloc (sizeof(DnaProfileMatchPair))) == NULL)  {  
      warn("DnaProfileMatchPair_alloc failed "); 
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->query = NULL;   
    out->target = NULL;  
    out->alb = NULL; 
    out->score = 0;  
    out->accepted = 0;   


    return out;  
}    


/* Function:  free_DnaProfileMatchPair(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [DnaProfileMatchPair *]
 *
 * Return [UNKN ]  Undocumented return value [DnaProfileMatchPair *]
 *
 */
DnaProfileMatchPair * free_DnaProfileMatchPair(DnaProfileMatchPair * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a DnaProfileMatchPair obj. Should be trappable");   
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
    if( obj->query != NULL)  
      free_DnaProfile(obj->query);   
    if( obj->target != NULL) 
      free_DnaProfile(obj->target);  
    if( obj->alb != NULL)    
      free_AlnBlock(obj->alb);   


    ckfree(obj); 
    return NULL; 
}    


/* Function:  swap_DnaProfileMatchPairSet(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_DnaProfileMatchPairSet
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [DnaProfileMatchPair **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_DnaProfileMatchPairSet(DnaProfileMatchPair ** list,int i,int j)  
{
    DnaProfileMatchPair * temp;  
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_DnaProfileMatchPairSet(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_DnaProfileMatchPairSet which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [DnaProfileMatchPair **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_DnaProfileMatchPairSet(DnaProfileMatchPair ** list,int left,int right,int (*comp)(DnaProfileMatchPair * ,DnaProfileMatchPair * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_DnaProfileMatchPairSet(list,left,(left+right)/2);   
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_DnaProfileMatchPairSet (list,++last,i); 
      }  
    swap_DnaProfileMatchPairSet (list,left,last);    
    qsort_DnaProfileMatchPairSet(list,left,last-1,comp); 
    qsort_DnaProfileMatchPairSet(list,last+1,right,comp);    
}    


/* Function:  sort_DnaProfileMatchPairSet(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_DnaProfileMatchPairSet
 *
 *
 * Arg:         obj [UNKN ] Object containing list [DnaProfileMatchPairSet *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_DnaProfileMatchPairSet(DnaProfileMatchPairSet * obj,int (*comp)(DnaProfileMatchPair *, DnaProfileMatchPair *)) 
{
    qsort_DnaProfileMatchPairSet(obj->pair,0,obj->len-1,comp);   
    return;  
}    


/* Function:  expand_DnaProfileMatchPairSet(obj,len)
 *
 * Descrip:    Really an internal function for add_DnaProfileMatchPairSet
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [DnaProfileMatchPairSet *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_DnaProfileMatchPairSet(DnaProfileMatchPairSet * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_DnaProfileMatchPairSet called with no need"); 
      return TRUE;   
      }  


    if( (obj->pair = (DnaProfileMatchPair ** ) ckrealloc (obj->pair,sizeof(DnaProfileMatchPair *)*len)) == NULL)     {  
      warn("ckrealloc failed for expand_DnaProfileMatchPairSet, returning FALSE");   
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_DnaProfileMatchPairSet(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [DnaProfileMatchPairSet *]
 * Arg:        add [OWNER] Object to add to the list [DnaProfileMatchPair *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_DnaProfileMatchPairSet(DnaProfileMatchPairSet * obj,DnaProfileMatchPair * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_DnaProfileMatchPairSet(obj,obj->len + DnaProfileMatchPairSetLISTLENGTH) == FALSE)   
        return FALSE;    
      }  


    obj->pair[obj->len++]=add;   
    return TRUE; 
}    


/* Function:  flush_DnaProfileMatchPairSet(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [DnaProfileMatchPairSet *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_DnaProfileMatchPairSet(DnaProfileMatchPairSet * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->pair[i] != NULL)  {  
        free_DnaProfileMatchPair(obj->pair[i]);  
        obj->pair[i] = NULL; 
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  DnaProfileMatchPairSet_alloc_std(void)
 *
 * Descrip:    Equivalent to DnaProfileMatchPairSet_alloc_len(DnaProfileMatchPairSetLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [DnaProfileMatchPairSet *]
 *
 */
DnaProfileMatchPairSet * DnaProfileMatchPairSet_alloc_std(void) 
{
    return DnaProfileMatchPairSet_alloc_len(DnaProfileMatchPairSetLISTLENGTH);   
}    


/* Function:  DnaProfileMatchPairSet_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [DnaProfileMatchPairSet *]
 *
 */
DnaProfileMatchPairSet * DnaProfileMatchPairSet_alloc_len(int len) 
{
    DnaProfileMatchPairSet * out;   /* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = DnaProfileMatchPairSet_alloc()) == NULL)   
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->pair = (DnaProfileMatchPair ** ) ckcalloc (len,sizeof(DnaProfileMatchPair *))) == NULL) {  
      warn("Warning, ckcalloc failed in DnaProfileMatchPairSet_alloc_len");  
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_DnaProfileMatchPairSet(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [DnaProfileMatchPairSet *]
 *
 * Return [UNKN ]  Undocumented return value [DnaProfileMatchPairSet *]
 *
 */
DnaProfileMatchPairSet * hard_link_DnaProfileMatchPairSet(DnaProfileMatchPairSet * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a DnaProfileMatchPairSet object: passed a NULL object");  
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  DnaProfileMatchPairSet_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [DnaProfileMatchPairSet *]
 *
 */
DnaProfileMatchPairSet * DnaProfileMatchPairSet_alloc(void) 
{
    DnaProfileMatchPairSet * out;   /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(DnaProfileMatchPairSet *) ckalloc (sizeof(DnaProfileMatchPairSet))) == NULL)    {  
      warn("DnaProfileMatchPairSet_alloc failed ");  
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->pair = NULL;    
    out->len = out->maxlen = 0;  


    return out;  
}    


/* Function:  free_DnaProfileMatchPairSet(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [DnaProfileMatchPairSet *]
 *
 * Return [UNKN ]  Undocumented return value [DnaProfileMatchPairSet *]
 *
 */
DnaProfileMatchPairSet * free_DnaProfileMatchPairSet(DnaProfileMatchPairSet * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a DnaProfileMatchPairSet obj. Should be trappable");    
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
    if( obj->pair != NULL)   {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->pair[i] != NULL)    
          free_DnaProfileMatchPair(obj->pair[i]);    
        }  
      ckfree(obj->pair); 
      }  


    ckfree(obj); 
    return NULL; 
}    



#ifdef _cplusplus
}
#endif
