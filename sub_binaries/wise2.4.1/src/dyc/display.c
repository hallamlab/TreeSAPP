#ifdef _cplusplus
extern "C" {
#endif
#include "display.h"





# line 97 "display.dy"
void write_Aln2Display_header(DYNFILE * dfp,Aln2Display * a2d)
{
  return;
}


# line 103 "display.dy"
void write_Aln2Display_function(DYNFILE * dfp,Aln2Display * a2d)
{

  write_Aln2Display_convert_func(dfp,a2d);

}


# line 111 "display.dy"
void write_Aln2Display_convert_func(DYNFILE * dfp,Aln2Display * a2d)
{
  register int i;
  register int j;
  register int k;
  char buffer[MAXLINE];


  /*** assume file is already loaded into writec ***/


  sprintf(buffer,"AlnDisplay * convert_AlnBlock_to_AlnDisplay_%s(AlnBlock * alb ",a2d->name);
  for(i=0;i<a2d->res_len;i++) {
    strcat(buffer,",");

    strcat(buffer,a2d->a2dr[i]->type);
    strcat(buffer,a2d->a2dr[i]->name);
  }
  strcat(buffer,")");

  start_function(dfp,buffer);
  expr(dfp,"AlnDisplay * out");
  expr(dfp,"AlnDisplayColumn * adc;");
  expr(dfp,"AlnDisplayUnit * adu;");
  expr(dfp,"AlnDisplayField * adf;");
  expr(dfp,"AlnDisplayNameBlock * adnb;");
  expr(dfp,"AlnDisplayNameUnit * adnu;");
  expr(dfp,"AlnIndexField * aif;");
  expr(dfp,"AlnDisplayIndexUnit * adiu;");
  expr(dfp,"AlnDisplayIndex * adi;");
  expr(dfp,"AlnColumn * alc;");
  expr(dfp,"AlnUnit   * alu;");
  expr(dfp,"char tempbuf[512]");

  add_break(dfp);

  expr(dfp,"out = new_AlnDisplay(alb->len)");
  add_break(dfp);
  
  for(i=0;i<a2d->len;i++) {
    expr(dfp,"adnb = AlnDisplayNameBlock_no_get(out,%d);",i);
    expr(dfp,"adnb->direction = ALN_DISPLAY_%s",a2d->s2d[i]->direction);
    expr(dfp,"adnu = get_next_AlnDisplayNameUnit(adnb);");
    if( a2d->s2d[i]->name_str != NULL) {
      expr(dfp,"adnu->name = %s;",a2d->s2d[i]->name_str);
      expr(dfp,"adnu->should_free = %s;",a2d->s2d[i]->is_static == TRUE ? "FALSE" : "TRUE");
    }
    else {
      expr(dfp,"adnu->name = \"NoName\";");
      expr(dfp,"adnu->should_free = FALSE;");
    }

    if( a2d->s2d[i]->ind_len > 0 ) {
      /* get the index  */
      add_break(dfp);
      add_block_comment(dfp,"This sequence has got an index, get the alndisplayindex");
      expr(dfp,"adi = AlnDisplayIndex_no_get(out,%d);",i);
      expr(dfp,"adi->direction = ALN_DISPLAY_%s",a2d->s2d[i]->direction);
      for(j=0;j< a2d->s2d[i]->ind_len;j++) {
	auto Index2Display  * i2d;
	i2d = a2d->s2d[i]->i2d[j];
	
	expr(dfp,"adiu = get_next_AlnDisplayIndexUnit(adi);");
	expr(dfp,"adiu->start = %s",i2d->eval == NULL ? "0" : i2d->eval);
	add_break(dfp);
      }
    } /* end of if index*/
  } /* end of for */

  


  expr(dfp,"for(alc = alb->start;alc != NULL;alc = alc->next)");
  startbrace_tag(dfp,"all columns in AlnBlock");
  expr(dfp,"adc = new_terminal_AlnDisplayColumn(out);");
  
  /*** go over each sequence, and move the labels across ***/

  for(i=0;i<a2d->len;i++) {
    j = 0;
    expr(dfp,"adu = AlnDisplayUnit_no_get(adc,%d)",i);
    expr(dfp,"alu = alc->alu[%d]",i);

    

    /*** start of the labels ***/

    for(j=0;j<a2d->s2d[i]->len;j++) {
      expr(dfp,"%s if(alu->text_label != NULL && strcmp(alu->text_label,\"%s\") == 0) ",j == 0 ? "" : "else",a2d->s2d[i]->l2d[j]->label);
      startbrace_tag(dfp,"if this label");


      /*** add the index if it is there ***/

      for(k=0;k<a2d->s2d[i]->l2d[j]->ind_len;k++) {
	expr(dfp,"aif = get_next_AlnIndexField(adu)");
	expr(dfp,"aif->index_no = %d;",a2d->s2d[i]->l2d[j]->i2d[k]->number);
	expr(dfp,"aif->length = %s;",a2d->s2d[i]->l2d[j]->i2d[k]->eval);
      }


      if(a2d->s2d[i]->l2d[j]->is_lone == TRUE ) {
	expr(dfp,"adc->is_lone = TRUE;");
      }

      expr(dfp,"adu->direction = ALN_DISPLAY_%s;",a2d->s2d[i]->l2d[j]->direction == NULL ? "DOWN" : a2d->s2d[i]->l2d[j]->direction);
      for(k=0;k<a2d->s2d[i]->l2d[j]->len;k++) {
	auto Aln2DisplayField * a2df;
	a2df = a2d->s2d[i]->l2d[j]->a2df[k];
	expr(dfp,"adf = get_next_AlnDisplayField(adu);");


	if( a2df->is_sub == TRUE ) {
	  expr(dfp,"sprintf(tempbuf,\"%%d\",alu->start);");
	  expr(dfp,"push_scan_and_replace_pair(\"%%START\",tempbuf);");
	  expr(dfp,"sprintf(tempbuf,\"%%d\",alu->end);");
	  expr(dfp,"push_scan_and_replace_pair(\"%%END\",tempbuf);");
	}

	if( a2df->is_string == TRUE )
	  expr(dfp,"adf->characters = \"%s\";",a2df->string);
	else if( a2df->is_single == TRUE ) 
	  expr(dfp,"adf->single = %s;",a2df->string);
	else if( a2df->is_sub == TRUE) { 
	  expr(dfp,"strcpy(tempbuf,\"%s\");",a2df->string);
	  expr(dfp,"adf->characters = stringalloc(scan_and_replace_line(tempbuf))");
	  expr(dfp,"pop_scan_and_replace_pair();");
	  expr(dfp,"pop_scan_and_replace_pair();");
	}
	else expr(dfp,"adf->characters = %s",a2df->string);


	if( a2df->is_sub == TRUE ) 
	  expr(dfp,"adf->length = strlen(adf->characters)");
	else expr(dfp,"adf->length = %s",a2df->length == NULL ? "1" : a2df->length);
	expr(dfp,"adf->direction = ALN_DISPLAY_%s",a2df->direction == NULL ? "RIGHT" : a2df->direction);
	expr(dfp,"adf->should_free = %s",a2df->is_static == TRUE ? "FALSE" : "TRUE");
	if( a2df->convert != NULL ) 
	  expr(dfp,"adf->convert_func = %s;",a2df->convert);
      }
      closebrace(dfp);
    }
    expr(dfp,"else ");
    startbrace(dfp);
    expr(dfp,"warn(\"In AlnBlock to AlnDisplay %s, got unintretable label [%%s] for sequence %d\",alc->alu[%d]->text_label);",a2d->name,i,i);
    expr(dfp,"adf = get_next_AlnDisplayField(adu);");
    expr(dfp,"adf->characters = \"?\";");
    expr(dfp,"adf->length = 1");
    expr(dfp,"adf->direction = ALN_DISPLAY_RIGHT");
    expr(dfp,"adf->should_free = FALSE");
    closebrace(dfp);
  }

  closebrace(dfp);

  expr(dfp,"return out;");

  close_function(dfp);
  add_break(dfp);

}
      

 /*************************/
 /* Access/checking func  */
 /*************************/

# line 278 "display.dy"
boolean prepare_Aln2Display(Aln2Display * a2d)
{
  boolean ret = TRUE;

  number_up_Sequence_Index2Display(a2d);

  if( crosslink_Index(a2d) == FALSE ){
    warn("Unable to crosslink Display %s: probably you have not got one of the indexes in the sequence block",a2d->name);
    ret = FALSE;
  }


  return ret;
}



# line 295 "display.dy"
void number_up_Sequence_Index2Display(Aln2Display * a2d)
{
  register int i;
  register int j;

  for(i=0;i<a2d->len;i++) 
    for(j=0;j<a2d->s2d[i]->ind_len;j++) 
      a2d->s2d[i]->i2d[j]->number = j;
}

# line 305 "display.dy"
boolean crosslink_Index(Aln2Display * a2d)
{
  register int i;
  register int j;
  boolean ret = TRUE;

  for(i=0;i<a2d->len;i++)
    for(j=0;j<a2d->s2d[i]->len;j++) {
      if( crosslink_Label2Display_Index(a2d->s2d[i]->l2d[j],a2d->s2d[i]) == FALSE ) 
	ret = FALSE;
    }


  return ret;
}



# line 323 "display.dy"
boolean crosslink_Label2Display_Index(Label2Display * l2d,Sequence2Display * s2d)
{
  register int i;
  Index2Display * temp;
  boolean ret = TRUE;

  for(i=0;i<l2d->ind_len;i++) {
    temp = Index2Display_from_name(s2d,l2d->i2d[i]->name);
    if( temp == NULL ) {
      warn("In Label %s of Sequence number %d, could not crosslink index name %s",l2d->label,s2d->number,l2d->i2d[i]->name);
      ret = FALSE;
    }

    else l2d->i2d[i]->number = temp->number;
  }
  
  return ret;
}

# line 342 "display.dy"
Index2Display * Index2Display_from_name(Sequence2Display * s2d,char * name)
{
  register int i;

  for(i=0;i<s2d->ind_len;i++) {
    if( strcmp(s2d->i2d[i]->name,name) == 0 )
      return s2d->i2d[i];
  }

  return NULL;
}



 /**************************/
 /* I/O functions          */
 /**************************/





# line 364 "display.dy"
Aln2Display         * read_Aln2Display_line(char * line,FILE * ifp)
{
  Sequence2Display * temp;
  Aln2DisplayResource * tempres;
  Aln2Display * out;
  char buffer[MAXLINE];
  char * runner;


  runner = strtok(line,spacestr);
  if( runner == NULL ) {
    warn("In read_Aln2Display_line, got completely blank line!");
    return NULL;
  }

  if( strcmp(runner,"display") != 0 ) {
    warn("In read_Aln2Display_line, got non display line... [%s]",line);
    return NULL;
  }

  runner = strtok(NULL,spacestr);
  if( runner == NULL ) {
    warn("In read_Sequence2Display_line, got no display name... cannot process!");
    return NULL;
  }
  

  out = Aln2Display_alloc_std();
  if( out == NULL)
    return NULL;

  out->name = stringalloc(runner);


  while( fgets(buffer,MAXLINE,ifp) != NULL ) {
    if( strwhitestartcmp(buffer,"enddisplay",spacestr) == 0 ) 
      break;
    if( strwhitestartcmp(buffer,"sequence",spacestr) == 0 ) {
      temp = read_Sequence2Display_line(buffer,ifp);
      if( temp == NULL ) 
	warn("In display number %d.. got bad sequence",out->name);
      else add_Aln2Display(out,temp);
    }
    else if ( strwhitestartcmp(buffer,"resource",spacestr) == 0 ) {
      tempres = read_Aln2DisplayResource_line(buffer);
     if( tempres == NULL ) 
	warn("In display %s.. got bad resource",out->name);
      else add_res_Aln2Display(out,tempres);
    }
    else {
      warn("In aln2display read.. got uninterpretable line [%s]",buffer);
    }

  }
  
  return out;
}


# line 423 "display.dy"
Sequence2Display    * read_Sequence2Display_line(char * line,FILE * ifp)
{
  Sequence2Display * out;
  Label2Display * temp;
  Index2Display * ind;
  char buffer[MAXLINE];
  char * runner;


  runner = strtok(line,spacestr);
  if( runner == NULL ) {
    warn("In read_Label2Display_line, got completely blank line!");
    return NULL;
  }

  if( strcmp(runner,"sequence") != 0 ) {
    warn("In read_Sequecne2Display_line, got non sequence line... [%s]",line);
    return NULL;
  }

  runner = strtok(NULL,spacestr);
  if( runner == NULL ) {
    warn("In read_Sequence2Display_line, got no  sequence number... cannot process!");
    return NULL;
  }
  

  out = Sequence2Display_alloc_std();
  if( out == NULL)
    return NULL;

  if( strcmp(runner,"ALLOTHERS") == 0)
    out->number = ALL_OTHER_NUMBER;
  else out->number = atoi(runner);


  while( fgets(buffer,MAXLINE,ifp) != NULL ) {
    if( strwhitestartcmp(buffer,"endsequence",spacestr) == 0 ) 
      break;
    else if( strwhitestartcmp(buffer,"name",spacestr) == 0 )
      read_name_line(out,buffer);
    else if( strwhitestartcmp(buffer,"label",spacestr) == 0 ) {
      temp = read_Label2Display_line(buffer,ifp);
      if( temp == NULL ) 
	warn("In sequence number %d.. got bad label",out->number);
      else add_Sequence2Display(out,temp);
    }
    else if ( strwhitestartcmp(buffer,"index",spacestr) == 0 ) {
      ind = read_Index2Display_line(buffer);
      if( ind == NULL ) 
	warn("In sequence number %d, unable to read Index line ",out->number);
      else add_ind_Sequence2Display(out,ind);
    }
      
    else {
      warn("In sequence2display.. got uninterpretable line [%s]",buffer);
    }

  }
  
  return out;
}

# line 486 "display.dy"
Label2Display       * read_Label2Display_line(char * line,FILE * ifp)
{
  Label2Display * out;
  Aln2DisplayField * temp;
  Index2Display * ind;
  char * runner;
  char buffer[MAXLINE];



  runner = strtok(line,spacestr);
  if( runner == NULL ) {
    warn("In read_Label2Display_line, got completely blank line!");
    return NULL;
  }

  if( strcmp(runner,"label") != 0 ) {
    warn("In read_Label2Display_line, got non label line... [%s]",line);
    return NULL;
  }

  runner = strtok(NULL,spacestr);
  if( runner == NULL ) {
    warn("In read_Label2Display_line, got no label name... cannot process!");
    return NULL;
  }

  


  out = Label2Display_alloc_std();

  if( out == NULL)
    return NULL;

  out->label = stringalloc(runner);


  runner = strtok(NULL,spacestr);
  while( runner != NULL ) {
    if( strcmp(runner,"!lone") == 0 ) {
      out->is_lone = TRUE;
    }

    else {
      warn("Could not understand label modifier %s",runner);
    }
    runner = strtok(NULL,spacestr);
  }


  while( fgets(buffer,MAXLINE,ifp) != NULL ) {
    if( strwhitestartcmp(buffer,"endlabel",spacestr) == 0)
      break;
    if( strwhitestartcmp(buffer,"field",spacestr) == 0 ) {
      temp = read_Aln2DisplayField_line(buffer,ifp);
      if( temp == NULL ) 
	warn("In label %s, unable to read field lines...",out->label);
      else add_Label2Display(out,temp);
    }
    else if( strwhitestartcmp(buffer,"index",spacestr) == 0 ) {
      ind = read_Index2Display_line(buffer);
      if( ind == NULL ) 
	warn("In label %s, unable to read Index line ",out->label);
      else add_ind_Label2Display(out,ind);
    }

    else if( strwhitestartcmp(buffer,"direction",spacestr) == 0 ) {
      runner = strtok(buffer,spacestr);
      out->direction = string_from_quoted_equality(runner);

    }
    else {
      warn("In reading label [%s], unable to interpret [%s]",out->label,buffer);
    }

  }

  return out;
}
    


# line 569 "display.dy"
Aln2DisplayField    * read_Aln2DisplayField_line(char * line,FILE * ifp)
{
  Aln2DisplayField * out;
  char ** base;
  char ** bkstr;
  char buffer[MAXLINE];

  out = Aln2DisplayField_alloc();
  if( out == NULL)
    return NULL;



  base = bkstr = breakstring(line,spacestr);

  for(bkstr++;*bkstr != NULL;bkstr++) {
    put_away_Aln2DisplayField_strpair(out,*bkstr);
  }
  ckfree(base);

  while( fgets(buffer,MAXLINE,ifp) != NULL ) {
    if( buffer[0] == '#')
      continue;
    if( strwhitestartcmp(buffer,"endfield",spacestr) == 0)
      break;
    base = bkstr = breakstring(buffer,spacestr);
    for(;*bkstr != NULL;bkstr++) {
      put_away_Aln2DisplayField_strpair(out,*bkstr);
    }


    ckfree(base);
  }


  return out;

}

# line 608 "display.dy"
boolean read_name_line(Sequence2Display * s2d,char * name_line)
{
  char ** base;
  char ** bkstr;
 
  if( strwhitestartcmp(name_line,"name",spacestr) != 0 ) {
    warn("Passed a non name line to read name line... [%s]",name_line);
    return FALSE;
  }
  

  base = bkstr = breakstring(name_line,spacestr);

  for(bkstr++;*bkstr != NULL;bkstr++) {
    if( strstartcmp(*bkstr,"string") == 0 ) {
      s2d->name_str = string_from_quoted_equality(*bkstr);
    } 
    else if( strstartcmp(*bkstr,"!static") == 0 ) {
      s2d->is_static = TRUE;
    }
    else if( strstartcmp(*bkstr,"length=") == 0 ) {
      s2d->length = number_from_quoted_equality(*bkstr);
      if( s2d->length < 0 || s2d->length > 512 ) {
	warn("sequence[%d] name length [%d] is ridiculous... setting to 20",s2d->number,s2d->length);
	s2d->length = 20;
      }
    }
    else if( strstartcmp(*bkstr,"direction=") ==0 ) {
      s2d->direction = string_from_quoted_equality(*bkstr);
    }
    else {
      warn("In reading sequence name line for sequence number %d, unable to interpret [%s]",s2d->number,*bkstr);
    }
  }

  if( s2d->direction == NULL )
    s2d->direction = stringalloc("DOWN");

  return TRUE;

}
  
  

# line 652 "display.dy"
boolean put_away_Aln2DisplayField_strpair(Aln2DisplayField * a2df,char * pair)
{
  char buffer[20];

  if( strstartcmp(pair,"string=") == 0) {
    if( a2df->string != NULL) {
      warn("Already a string or eval state ment in this Aln2DisplayField, overriding %s",a2df->string);
      ckfree(a2df->string);
    }
    a2df->string = string_from_quoted_equality(pair);
    a2df->is_string = TRUE;
    sprintf(buffer,"%d",(int)strlen(a2df->string));
    a2df->length = stringalloc(buffer);
    a2df->is_static = TRUE;
  }
  else if( strstartcmp(pair,"eval") == 0 ) {
    if( a2df->string != NULL) {
      warn("Already a string or eval state ment in this Aln2DisplayField, overriding %s",a2df->string);
      ckfree(a2df->string);
    }
    a2df->string = string_from_quoted_equality(pair);
  }
  else if( strstartcmp(pair,"substring") == 0 ) {
    if( a2df->string != NULL) {
      warn("Already a string or eval state ment in this Aln2DisplayField, overriding %s",a2df->string);
      ckfree(a2df->string);
    }
    a2df->string = string_from_quoted_equality(pair);
    a2df->is_sub    = TRUE;
    a2df->is_static = FALSE;
  }

  else if( strstartcmp(pair,"length=") == 0 ) {
    a2df->length = string_from_quoted_equality(pair);
  }
  else if( strstartcmp(pair,"direction=") == 0) {
    a2df->direction = string_from_quoted_equality(pair);
  }
  else if( strstartcmp(pair,"convert=") == 0 ) {
    a2df->convert = string_from_quoted_equality(pair);
  }
  else if( strstartcmp(pair,"!static") == 0) {
    a2df->is_static = TRUE;
  }
  else if( strstartcmp(pair,"!single") == 0) {
    a2df->is_single = TRUE;
  }
  else {
    warn("Unable to interpret [%s] as a valid Aln2DisplayField tag",pair);
    return FALSE;
  }
  
  return TRUE;
}

# line 707 "display.dy"
Aln2DisplayResource * read_Aln2DisplayResource_line(char * line)
{
  Aln2DisplayResource * out;
  char ** base;
  char ** bkstr;

  out = Aln2DisplayResource_alloc();
  if( out == NULL)
    return NULL;


  base = bkstr = breakstring(line,spacestr);


  for(bkstr++;*bkstr != NULL;bkstr++) {

    if( strstartcmp(*bkstr,"type=") == 0) {
     
      out->type = string_from_quoted_equality(*bkstr);
    }
    else if( strstartcmp(*bkstr,"name=") == 0 ) {
      out->name = string_from_quoted_equality(*bkstr);
    }
    else if( strstartcmp(*bkstr,"arg=") == 0) {
      out->arg = string_from_quoted_equality(*bkstr);
    }
    else {
      warn("Got strange tag [%s] in resource line..., while reading display... ignoring",*bkstr);
    }
  }
  ckfree(base);


  return out;
}

 


# line 746 "display.dy"
Index2Display * read_Index2Display_line(char * line)
{
  Index2Display * out;
  char ** bkstr;
  char ** base;

  if( strwhitestartcmp(line,"index",spacestr) != 0 ) {
    warn("Tried to read an Index2Display line with no index! Bad news!");
    return NULL;
  }

  
  base = bkstr = breakstring(line,spacestr);

  out = Index2Display_alloc();
  if( out == NULL )
    return NULL;

  bkstr++;
  
  out->name = stringalloc(*bkstr);

  for(bkstr++;*bkstr != NULL;bkstr++) {
    if( strstartcmp(*bkstr,"start") == 0 ) {
      if( out->eval != NULL ) {
	warn("you are replacing an index evaluation. Remember that each index can only have one length [or start if in sequence]");
	ckfree(out->eval);
      }
      out->eval = string_from_quoted_equality(*bkstr);
      out->is_start = TRUE;
    }

    else if( strstartcmp(*bkstr,"length") == 0 ) {
      if( out->eval != NULL ) {
	warn("you are replacing an index evaluation. Remember that each index can only have one length [or start if in sequence]");
	ckfree(out->eval);
      }
      out->eval = string_from_quoted_equality(*bkstr);
    }
    else {
      warn("Found an uninterpretable tag [%s] in index %s",out->name,*bkstr);
    }

  }

  ckfree(base);

  return out;
}




# line 799 "display.dy"
void show_Aln2Display(Aln2Display * a2d,FILE * ofp)
{
  register int i;
  
  fprintf(ofp,"Display %s\n",a2d->name);
  for(i=0;i<a2d->res_len;i++)
    show_Aln2DisplayResource(a2d->a2dr[i],ofp);
  for(i=0;i<a2d->len;i++)
    show_Sequence2Display(a2d->s2d[i],ofp);
}

# line 810 "display.dy"
void show_Aln2DisplayResource(Aln2DisplayResource * a2dr,FILE * ofp)
{
  fprintf(ofp,"Resource [%s] Type [%s] Argument [%s]\n",a2dr->type,a2dr->name,a2dr->arg);
}

# line 815 "display.dy"
void show_Sequence2Display(Sequence2Display * s2d,FILE * ofp)
{
  register int i;
  fprintf(ofp,"Sequence %d\n",s2d->number);
  for(i=0;i<s2d->len;i++)
    show_Label2Display(s2d->l2d[i],ofp);
}

# line 823 "display.dy"
void show_Label2Display(Label2Display * l2d,FILE * ofp)
{
  register int i;
  fprintf(ofp,"\tLabel %s\n",l2d->label);
  for(i=0;i<l2d->len;i++) {
    fprintf(ofp,"\t  Field %d",i);
    show_Aln2DisplayField(l2d->a2df[i],ofp);
  }
}

# line 833 "display.dy"
void show_Aln2DisplayField(Aln2DisplayField * a2df,FILE * ofp)
{
  fprintf(ofp," String = [%s], Length [%s], Direction [%s]\n",a2df->string,a2df->length,a2df->direction);
}


# line 772 "display.c"
/* Function:  hard_link_Aln2DisplayField(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [Aln2DisplayField *]
 *
 * Return [UNKN ]  Undocumented return value [Aln2DisplayField *]
 *
 */
Aln2DisplayField * hard_link_Aln2DisplayField(Aln2DisplayField * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a Aln2DisplayField object: passed a NULL object");    
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  Aln2DisplayField_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [Aln2DisplayField *]
 *
 */
Aln2DisplayField * Aln2DisplayField_alloc(void) 
{
    Aln2DisplayField * out; /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(Aln2DisplayField *) ckalloc (sizeof(Aln2DisplayField))) == NULL)    {  
      warn("Aln2DisplayField_alloc failed ");    
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->string = NULL;  
    out->length = NULL;  
    out->direction = NULL;   
    out->convert = NULL; 
    out->is_static = FALSE;  
    out->is_string = FALSE;  
    out->is_single = FALSE;  
    out->is_sub = FALSE; 


    return out;  
}    


/* Function:  free_Aln2DisplayField(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [Aln2DisplayField *]
 *
 * Return [UNKN ]  Undocumented return value [Aln2DisplayField *]
 *
 */
Aln2DisplayField * free_Aln2DisplayField(Aln2DisplayField * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a Aln2DisplayField obj. Should be trappable");  
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
    if( obj->string != NULL) 
      ckfree(obj->string);   
    if( obj->length != NULL) 
      ckfree(obj->length);   
    if( obj->direction != NULL)  
      ckfree(obj->direction);    
    if( obj->convert != NULL)    
      ckfree(obj->convert);  


    ckfree(obj); 
    return NULL; 
}    


/* Function:  hard_link_Index2Display(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [Index2Display *]
 *
 * Return [UNKN ]  Undocumented return value [Index2Display *]
 *
 */
Index2Display * hard_link_Index2Display(Index2Display * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a Index2Display object: passed a NULL object");   
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  Index2Display_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [Index2Display *]
 *
 */
Index2Display * Index2Display_alloc(void) 
{
    Index2Display * out;/* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(Index2Display *) ckalloc (sizeof(Index2Display))) == NULL)  {  
      warn("Index2Display_alloc failed ");   
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->name = NULL;    
    out->eval = NULL;    
    out->is_start = FALSE;   
    out->number = 0; 


    return out;  
}    


/* Function:  free_Index2Display(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [Index2Display *]
 *
 * Return [UNKN ]  Undocumented return value [Index2Display *]
 *
 */
Index2Display * free_Index2Display(Index2Display * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a Index2Display obj. Should be trappable"); 
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
    if( obj->eval != NULL)   
      ckfree(obj->eval);     


    ckfree(obj); 
    return NULL; 
}    


/* Function:  swap_Label2Display(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_Label2Display
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [Aln2DisplayField **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_Label2Display(Aln2DisplayField ** list,int i,int j)  
{
    Aln2DisplayField * temp; 
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_Label2Display(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_Label2Display which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [Aln2DisplayField **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_Label2Display(Aln2DisplayField ** list,int left,int right,int (*comp)(Aln2DisplayField * ,Aln2DisplayField * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_Label2Display(list,left,(left+right)/2);    
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_Label2Display (list,++last,i);  
      }  
    swap_Label2Display (list,left,last); 
    qsort_Label2Display(list,left,last-1,comp);  
    qsort_Label2Display(list,last+1,right,comp); 
}    


/* Function:  sort_Label2Display(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_Label2Display
 *
 *
 * Arg:         obj [UNKN ] Object containing list [Label2Display *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_Label2Display(Label2Display * obj,int (*comp)(Aln2DisplayField *, Aln2DisplayField *)) 
{
    qsort_Label2Display(obj->a2df,0,obj->len-1,comp);    
    return;  
}    


/* Function:  expand_Label2Display(obj,len)
 *
 * Descrip:    Really an internal function for add_Label2Display
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [Label2Display *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_Label2Display(Label2Display * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_Label2Display called with no need");  
      return TRUE;   
      }  


    if( (obj->a2df = (Aln2DisplayField ** ) ckrealloc (obj->a2df,sizeof(Aln2DisplayField *)*len)) == NULL)   {  
      warn("ckrealloc failed for expand_Label2Display, returning FALSE");    
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_Label2Display(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [Label2Display *]
 * Arg:        add [OWNER] Object to add to the list [Aln2DisplayField *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_Label2Display(Label2Display * obj,Aln2DisplayField * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_Label2Display(obj,obj->len + Label2DisplayLISTLENGTH) == FALSE) 
        return FALSE;    
      }  


    obj->a2df[obj->len++]=add;   
    return TRUE; 
}    


/* Function:  flush_Label2Display(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [Label2Display *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_Label2Display(Label2Display * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->a2df[i] != NULL)  {  
        free_Aln2DisplayField(obj->a2df[i]); 
        obj->a2df[i] = NULL; 
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  swap_ind_Label2Display(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_ind_Label2Display
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [Index2Display **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_ind_Label2Display(Index2Display ** list,int i,int j)  
{
    Index2Display * temp;    
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_ind_Label2Display(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_ind_Label2Display which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [Index2Display **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_ind_Label2Display(Index2Display ** list,int left,int right,int (*comp)(Index2Display * ,Index2Display * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_ind_Label2Display(list,left,(left+right)/2);    
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_ind_Label2Display (list,++last,i);  
      }  
    swap_ind_Label2Display (list,left,last); 
    qsort_ind_Label2Display(list,left,last-1,comp);  
    qsort_ind_Label2Display(list,last+1,right,comp); 
}    


/* Function:  sort_ind_Label2Display(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_ind_Label2Display
 *
 *
 * Arg:         obj [UNKN ] Object containing list [Label2Display *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_ind_Label2Display(Label2Display * obj,int (*comp)(Index2Display *, Index2Display *)) 
{
    qsort_ind_Label2Display(obj->i2d,0,obj->ind_len-1,comp); 
    return;  
}    


/* Function:  expand_ind_Label2Display(obj,len)
 *
 * Descrip:    Really an internal function for add_ind_Label2Display
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [Label2Display *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_ind_Label2Display(Label2Display * obj,int len) 
{


    if( obj->ind_maxlen > obj->ind_len )     {  
      warn("expand_Label2Displayind_ called with no need");  
      return TRUE;   
      }  


    if( (obj->i2d = (Index2Display ** ) ckrealloc (obj->i2d,sizeof(Index2Display *)*len)) == NULL)   {  
      warn("ckrealloc failed for expand_Label2Display, returning FALSE");    
      return FALSE;  
      }  
    obj->ind_maxlen = len;   
    return TRUE; 
}    


/* Function:  add_ind_Label2Display(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [Label2Display *]
 * Arg:        add [OWNER] Object to add to the list [Index2Display *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_ind_Label2Display(Label2Display * obj,Index2Display * add) 
{
    if( obj->ind_len >= obj->ind_maxlen) {  
      if( expand_ind_Label2Display(obj,obj->ind_len + Label2DisplayLISTLENGTH) == FALSE) 
        return FALSE;    
      }  


    obj->i2d[obj->ind_len++]=add;    
    return TRUE; 
}    


/* Function:  flush_ind_Label2Display(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [Label2Display *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_ind_Label2Display(Label2Display * obj) 
{
    int i;   


    for(i=0;i<obj->ind_len;i++)  { /*for i over list length*/ 
      if( obj->i2d[i] != NULL)   {  
        free_Index2Display(obj->i2d[i]); 
        obj->i2d[i] = NULL;  
        }  
      } /* end of for i over list length */ 


    obj->ind_len = 0;    
    return i;    
}    


/* Function:  Label2Display_alloc_std(void)
 *
 * Descrip:    Equivalent to Label2Display_alloc_len(Label2DisplayLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [Label2Display *]
 *
 */
Label2Display * Label2Display_alloc_std(void) 
{
    return Label2Display_alloc_len(Label2DisplayLISTLENGTH); 
}    


/* Function:  Label2Display_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [Label2Display *]
 *
 */
Label2Display * Label2Display_alloc_len(int len) 
{
    Label2Display * out;/* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = Label2Display_alloc()) == NULL)    
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->a2df = (Aln2DisplayField ** ) ckcalloc (len,sizeof(Aln2DisplayField *))) == NULL)   {  
      warn("Warning, ckcalloc failed in Label2Display_alloc_len");   
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    if((out->i2d = (Index2Display ** ) ckcalloc (len,sizeof(Index2Display *))) == NULL)  {  
      warn("Warning, ckcalloc failed in Label2Display_alloc_len");   
      return NULL;   
      }  
    out->ind_len = 0;    
    out->ind_maxlen = len;   


    return out;  
}    


/* Function:  hard_link_Label2Display(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [Label2Display *]
 *
 * Return [UNKN ]  Undocumented return value [Label2Display *]
 *
 */
Label2Display * hard_link_Label2Display(Label2Display * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a Label2Display object: passed a NULL object");   
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  Label2Display_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [Label2Display *]
 *
 */
Label2Display * Label2Display_alloc(void) 
{
    Label2Display * out;/* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(Label2Display *) ckalloc (sizeof(Label2Display))) == NULL)  {  
      warn("Label2Display_alloc failed ");   
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->label = NULL;   
    out->direction = NULL;   
    out->a2df = NULL;    
    out->len = out->maxlen = 0;  
    out->i2d = NULL; 
    out->ind_len = out->ind_maxlen = 0;  
    out->is_lone = FALSE;    


    return out;  
}    


/* Function:  free_Label2Display(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [Label2Display *]
 *
 * Return [UNKN ]  Undocumented return value [Label2Display *]
 *
 */
Label2Display * free_Label2Display(Label2Display * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a Label2Display obj. Should be trappable"); 
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
    if( obj->label != NULL)  
      ckfree(obj->label);    
    if( obj->direction != NULL)  
      ckfree(obj->direction);    
    if( obj->a2df != NULL)   {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->a2df[i] != NULL)    
          free_Aln2DisplayField(obj->a2df[i]);   
        }  
      ckfree(obj->a2df); 
      }  
    if( obj->i2d != NULL)    {  
      for(i=0;i<obj->ind_len;i++)    {  
        if( obj->i2d[i] != NULL) 
          free_Index2Display(obj->i2d[i]);   
        }  
      ckfree(obj->i2d);  
      }  


    ckfree(obj); 
    return NULL; 
}    


/* Function:  swap_Sequence2Display(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_Sequence2Display
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [Label2Display **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_Sequence2Display(Label2Display ** list,int i,int j)  
{
    Label2Display * temp;    
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_Sequence2Display(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_Sequence2Display which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [Label2Display **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_Sequence2Display(Label2Display ** list,int left,int right,int (*comp)(Label2Display * ,Label2Display * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_Sequence2Display(list,left,(left+right)/2); 
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_Sequence2Display (list,++last,i);   
      }  
    swap_Sequence2Display (list,left,last);  
    qsort_Sequence2Display(list,left,last-1,comp);   
    qsort_Sequence2Display(list,last+1,right,comp);  
}    


/* Function:  sort_Sequence2Display(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_Sequence2Display
 *
 *
 * Arg:         obj [UNKN ] Object containing list [Sequence2Display *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_Sequence2Display(Sequence2Display * obj,int (*comp)(Label2Display *, Label2Display *)) 
{
    qsort_Sequence2Display(obj->l2d,0,obj->len-1,comp);  
    return;  
}    


/* Function:  expand_Sequence2Display(obj,len)
 *
 * Descrip:    Really an internal function for add_Sequence2Display
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [Sequence2Display *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_Sequence2Display(Sequence2Display * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_Sequence2Display called with no need");   
      return TRUE;   
      }  


    if( (obj->l2d = (Label2Display ** ) ckrealloc (obj->l2d,sizeof(Label2Display *)*len)) == NULL)   {  
      warn("ckrealloc failed for expand_Sequence2Display, returning FALSE"); 
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_Sequence2Display(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [Sequence2Display *]
 * Arg:        add [OWNER] Object to add to the list [Label2Display *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_Sequence2Display(Sequence2Display * obj,Label2Display * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_Sequence2Display(obj,obj->len + Sequence2DisplayLISTLENGTH) == FALSE)   
        return FALSE;    
      }  


    obj->l2d[obj->len++]=add;    
    return TRUE; 
}    


/* Function:  flush_Sequence2Display(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [Sequence2Display *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_Sequence2Display(Sequence2Display * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->l2d[i] != NULL)   {  
        free_Label2Display(obj->l2d[i]); 
        obj->l2d[i] = NULL;  
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  swap_ind_Sequence2Display(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_ind_Sequence2Display
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [Index2Display **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_ind_Sequence2Display(Index2Display ** list,int i,int j)  
{
    Index2Display * temp;    
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_ind_Sequence2Display(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_ind_Sequence2Display which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [Index2Display **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_ind_Sequence2Display(Index2Display ** list,int left,int right,int (*comp)(Index2Display * ,Index2Display * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_ind_Sequence2Display(list,left,(left+right)/2); 
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_ind_Sequence2Display (list,++last,i);   
      }  
    swap_ind_Sequence2Display (list,left,last);  
    qsort_ind_Sequence2Display(list,left,last-1,comp);   
    qsort_ind_Sequence2Display(list,last+1,right,comp);  
}    


/* Function:  sort_ind_Sequence2Display(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_ind_Sequence2Display
 *
 *
 * Arg:         obj [UNKN ] Object containing list [Sequence2Display *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_ind_Sequence2Display(Sequence2Display * obj,int (*comp)(Index2Display *, Index2Display *)) 
{
    qsort_ind_Sequence2Display(obj->i2d,0,obj->ind_len-1,comp);  
    return;  
}    


/* Function:  expand_ind_Sequence2Display(obj,len)
 *
 * Descrip:    Really an internal function for add_ind_Sequence2Display
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [Sequence2Display *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_ind_Sequence2Display(Sequence2Display * obj,int len) 
{


    if( obj->ind_maxlen > obj->ind_len )     {  
      warn("expand_Sequence2Displayind_ called with no need");   
      return TRUE;   
      }  


    if( (obj->i2d = (Index2Display ** ) ckrealloc (obj->i2d,sizeof(Index2Display *)*len)) == NULL)   {  
      warn("ckrealloc failed for expand_Sequence2Display, returning FALSE"); 
      return FALSE;  
      }  
    obj->ind_maxlen = len;   
    return TRUE; 
}    


/* Function:  add_ind_Sequence2Display(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [Sequence2Display *]
 * Arg:        add [OWNER] Object to add to the list [Index2Display *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_ind_Sequence2Display(Sequence2Display * obj,Index2Display * add) 
{
    if( obj->ind_len >= obj->ind_maxlen) {  
      if( expand_ind_Sequence2Display(obj,obj->ind_len + Sequence2DisplayLISTLENGTH) == FALSE)   
        return FALSE;    
      }  


    obj->i2d[obj->ind_len++]=add;    
    return TRUE; 
}    


/* Function:  flush_ind_Sequence2Display(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [Sequence2Display *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_ind_Sequence2Display(Sequence2Display * obj) 
{
    int i;   


    for(i=0;i<obj->ind_len;i++)  { /*for i over list length*/ 
      if( obj->i2d[i] != NULL)   {  
        free_Index2Display(obj->i2d[i]); 
        obj->i2d[i] = NULL;  
        }  
      } /* end of for i over list length */ 


    obj->ind_len = 0;    
    return i;    
}    


/* Function:  Sequence2Display_alloc_std(void)
 *
 * Descrip:    Equivalent to Sequence2Display_alloc_len(Sequence2DisplayLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [Sequence2Display *]
 *
 */
Sequence2Display * Sequence2Display_alloc_std(void) 
{
    return Sequence2Display_alloc_len(Sequence2DisplayLISTLENGTH);   
}    


/* Function:  Sequence2Display_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [Sequence2Display *]
 *
 */
Sequence2Display * Sequence2Display_alloc_len(int len) 
{
    Sequence2Display * out; /* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = Sequence2Display_alloc()) == NULL) 
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->l2d = (Label2Display ** ) ckcalloc (len,sizeof(Label2Display *))) == NULL)  {  
      warn("Warning, ckcalloc failed in Sequence2Display_alloc_len");    
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    if((out->i2d = (Index2Display ** ) ckcalloc (len,sizeof(Index2Display *))) == NULL)  {  
      warn("Warning, ckcalloc failed in Sequence2Display_alloc_len");    
      return NULL;   
      }  
    out->ind_len = 0;    
    out->ind_maxlen = len;   


    return out;  
}    


/* Function:  hard_link_Sequence2Display(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [Sequence2Display *]
 *
 * Return [UNKN ]  Undocumented return value [Sequence2Display *]
 *
 */
Sequence2Display * hard_link_Sequence2Display(Sequence2Display * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a Sequence2Display object: passed a NULL object");    
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  Sequence2Display_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [Sequence2Display *]
 *
 */
Sequence2Display * Sequence2Display_alloc(void) 
{
    Sequence2Display * out; /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(Sequence2Display *) ckalloc (sizeof(Sequence2Display))) == NULL)    {  
      warn("Sequence2Display_alloc failed ");    
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->number = 0; 
    out->l2d = NULL; 
    out->len = out->maxlen = 0;  
    out->i2d = NULL; 
    out->ind_len = out->ind_maxlen = 0;  
    out->name_str = NULL;    
    out->length = 0; 
    out->is_static = FALSE;  
    out->direction = NULL;   


    return out;  
}    


/* Function:  free_Sequence2Display(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [Sequence2Display *]
 *
 * Return [UNKN ]  Undocumented return value [Sequence2Display *]
 *
 */
Sequence2Display * free_Sequence2Display(Sequence2Display * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a Sequence2Display obj. Should be trappable");  
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
    if( obj->l2d != NULL)    {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->l2d[i] != NULL) 
          free_Label2Display(obj->l2d[i]);   
        }  
      ckfree(obj->l2d);  
      }  
    if( obj->i2d != NULL)    {  
      for(i=0;i<obj->ind_len;i++)    {  
        if( obj->i2d[i] != NULL) 
          free_Index2Display(obj->i2d[i]);   
        }  
      ckfree(obj->i2d);  
      }  
    if( obj->name_str != NULL)   
      ckfree(obj->name_str);     
    if( obj->direction != NULL)  
      ckfree(obj->direction);    


    ckfree(obj); 
    return NULL; 
}    


/* Function:  hard_link_Aln2DisplayResource(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [Aln2DisplayResource *]
 *
 * Return [UNKN ]  Undocumented return value [Aln2DisplayResource *]
 *
 */
Aln2DisplayResource * hard_link_Aln2DisplayResource(Aln2DisplayResource * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a Aln2DisplayResource object: passed a NULL object"); 
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  Aln2DisplayResource_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [Aln2DisplayResource *]
 *
 */
Aln2DisplayResource * Aln2DisplayResource_alloc(void) 
{
    Aln2DisplayResource * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(Aln2DisplayResource *) ckalloc (sizeof(Aln2DisplayResource))) == NULL)  {  
      warn("Aln2DisplayResource_alloc failed "); 
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->type = NULL;    
    out->name = NULL;    
    out->arg = NULL; 


    return out;  
}    


/* Function:  free_Aln2DisplayResource(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [Aln2DisplayResource *]
 *
 * Return [UNKN ]  Undocumented return value [Aln2DisplayResource *]
 *
 */
Aln2DisplayResource * free_Aln2DisplayResource(Aln2DisplayResource * obj) 
{
    int return_early = 0;    


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a Aln2DisplayResource obj. Should be trappable");   
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
    if( obj->type != NULL)   
      ckfree(obj->type);     
    if( obj->name != NULL)   
      ckfree(obj->name);     
    if( obj->arg != NULL)    
      ckfree(obj->arg);  


    ckfree(obj); 
    return NULL; 
}    


/* Function:  swap_res_Aln2Display(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_res_Aln2Display
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [Aln2DisplayResource **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_res_Aln2Display(Aln2DisplayResource ** list,int i,int j)  
{
    Aln2DisplayResource * temp;  
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_res_Aln2Display(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_res_Aln2Display which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [Aln2DisplayResource **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_res_Aln2Display(Aln2DisplayResource ** list,int left,int right,int (*comp)(Aln2DisplayResource * ,Aln2DisplayResource * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_res_Aln2Display(list,left,(left+right)/2);  
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_res_Aln2Display (list,++last,i);    
      }  
    swap_res_Aln2Display (list,left,last);   
    qsort_res_Aln2Display(list,left,last-1,comp);    
    qsort_res_Aln2Display(list,last+1,right,comp);   
}    


/* Function:  sort_res_Aln2Display(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_res_Aln2Display
 *
 *
 * Arg:         obj [UNKN ] Object containing list [Aln2Display *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_res_Aln2Display(Aln2Display * obj,int (*comp)(Aln2DisplayResource *, Aln2DisplayResource *)) 
{
    qsort_res_Aln2Display(obj->a2dr,0,obj->res_len-1,comp);  
    return;  
}    


/* Function:  expand_res_Aln2Display(obj,len)
 *
 * Descrip:    Really an internal function for add_res_Aln2Display
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [Aln2Display *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_res_Aln2Display(Aln2Display * obj,int len) 
{


    if( obj->res_maxlen > obj->res_len )     {  
      warn("expand_Aln2Displayres_ called with no need");    
      return TRUE;   
      }  


    if( (obj->a2dr = (Aln2DisplayResource ** ) ckrealloc (obj->a2dr,sizeof(Aln2DisplayResource *)*len)) == NULL)     {  
      warn("ckrealloc failed for expand_Aln2Display, returning FALSE");  
      return FALSE;  
      }  
    obj->res_maxlen = len;   
    return TRUE; 
}    


/* Function:  add_res_Aln2Display(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [Aln2Display *]
 * Arg:        add [OWNER] Object to add to the list [Aln2DisplayResource *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_res_Aln2Display(Aln2Display * obj,Aln2DisplayResource * add) 
{
    if( obj->res_len >= obj->res_maxlen) {  
      if( expand_res_Aln2Display(obj,obj->res_len + Aln2DisplayLISTLENGTH) == FALSE) 
        return FALSE;    
      }  


    obj->a2dr[obj->res_len++]=add;   
    return TRUE; 
}    


/* Function:  flush_res_Aln2Display(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [Aln2Display *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_res_Aln2Display(Aln2Display * obj) 
{
    int i;   


    for(i=0;i<obj->res_len;i++)  { /*for i over list length*/ 
      if( obj->a2dr[i] != NULL)  {  
        free_Aln2DisplayResource(obj->a2dr[i]);  
        obj->a2dr[i] = NULL; 
        }  
      } /* end of for i over list length */ 


    obj->res_len = 0;    
    return i;    
}    


/* Function:  swap_Aln2Display(list,i,j)
 *
 * Descrip:    swap function: an internal for qsort_Aln2Display
 *             swaps two positions in the array
 *
 *
 * Arg:        list [UNKN ] List of structures to swap in [Sequence2Display **]
 * Arg:           i [UNKN ] swap position [int]
 * Arg:           j [UNKN ] swap position [int]
 *
 */
/* swap function for qsort function */ 
void swap_Aln2Display(Sequence2Display ** list,int i,int j)  
{
    Sequence2Display * temp; 
    temp=list[i];    
    list[i]=list[j]; 
    list[j]=temp;    
}    


/* Function:  qsort_Aln2Display(list,left,right,comp)
 *
 * Descrip:    qsort - lifted from K&R 
 *             sorts the array using quicksort
 *             Probably much better to call sort_Aln2Display which sorts from start to end
 *
 *
 * Arg:         list [UNKN ] List of structures to swap in [Sequence2Display **]
 * Arg:         left [UNKN ] left position [int]
 * Arg:        right [UNKN ] right position [int]
 * Arg:         comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void qsort_Aln2Display(Sequence2Display ** list,int left,int right,int (*comp)(Sequence2Display * ,Sequence2Display * )) 
{
    int i,last;  
    if( left >= right )  
      return;    


    swap_Aln2Display(list,left,(left+right)/2);  
    last = left; 
    for ( i=left+1; i <= right;i++)  {  
      if( (*comp)(list[i],list[left]) < 0)   
        swap_Aln2Display (list,++last,i);    
      }  
    swap_Aln2Display (list,left,last);   
    qsort_Aln2Display(list,left,last-1,comp);    
    qsort_Aln2Display(list,last+1,right,comp);   
}    


/* Function:  sort_Aln2Display(obj,comp)
 *
 * Descrip:    sorts from start to end using comp 
 *             sorts the array using quicksort by calling qsort_Aln2Display
 *
 *
 * Arg:         obj [UNKN ] Object containing list [Aln2Display *]
 * Arg:        comp [FUNCP] Function which returns -1 or 1 to sort on [int (*comp]
 *
 */
void sort_Aln2Display(Aln2Display * obj,int (*comp)(Sequence2Display *, Sequence2Display *)) 
{
    qsort_Aln2Display(obj->s2d,0,obj->len-1,comp);   
    return;  
}    


/* Function:  expand_Aln2Display(obj,len)
 *
 * Descrip:    Really an internal function for add_Aln2Display
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [Aln2Display *]
 * Arg:        len [UNKN ] Length to add one [int]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean expand_Aln2Display(Aln2Display * obj,int len) 
{


    if( obj->maxlen > obj->len )     {  
      warn("expand_Aln2Display called with no need");    
      return TRUE;   
      }  


    if( (obj->s2d = (Sequence2Display ** ) ckrealloc (obj->s2d,sizeof(Sequence2Display *)*len)) == NULL)     {  
      warn("ckrealloc failed for expand_Aln2Display, returning FALSE");  
      return FALSE;  
      }  
    obj->maxlen = len;   
    return TRUE; 
}    


/* Function:  add_Aln2Display(obj,add)
 *
 * Descrip:    Adds another object to the list. It will expand the list if necessary
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list [Aln2Display *]
 * Arg:        add [OWNER] Object to add to the list [Sequence2Display *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
/* will expand function if necessary */ 
boolean add_Aln2Display(Aln2Display * obj,Sequence2Display * add) 
{
    if( obj->len >= obj->maxlen) {  
      if( expand_Aln2Display(obj,obj->len + Aln2DisplayLISTLENGTH) == FALSE) 
        return FALSE;    
      }  


    obj->s2d[obj->len++]=add;    
    return TRUE; 
}    


/* Function:  flush_Aln2Display(obj)
 *
 * Descrip:    Frees the list elements, sets length to 0
 *             If you want to save some elements, use hard_link_xxx
 *             to protect them from being actually destroyed in the free
 *
 *
 * Arg:        obj [UNKN ] Object which contains the list  [Aln2Display *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int flush_Aln2Display(Aln2Display * obj) 
{
    int i;   


    for(i=0;i<obj->len;i++)  { /*for i over list length*/ 
      if( obj->s2d[i] != NULL)   {  
        free_Sequence2Display(obj->s2d[i]);  
        obj->s2d[i] = NULL;  
        }  
      } /* end of for i over list length */ 


    obj->len = 0;    
    return i;    
}    


/* Function:  Aln2Display_alloc_std(void)
 *
 * Descrip:    Equivalent to Aln2Display_alloc_len(Aln2DisplayLISTLENGTH)
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [Aln2Display *]
 *
 */
Aln2Display * Aln2Display_alloc_std(void) 
{
    return Aln2Display_alloc_len(Aln2DisplayLISTLENGTH); 
}    


/* Function:  Aln2Display_alloc_len(len)
 *
 * Descrip:    Allocates len length to all lists
 *
 *
 * Arg:        len [UNKN ] Length of lists to allocate [int]
 *
 * Return [UNKN ]  Undocumented return value [Aln2Display *]
 *
 */
Aln2Display * Aln2Display_alloc_len(int len) 
{
    Aln2Display * out;  /* out is exported at the end of function */ 


    /* Call alloc function: return NULL if NULL */ 
    /* Warning message alread in alloc function */ 
    if((out = Aln2Display_alloc()) == NULL)  
      return NULL;   


    /* Calling ckcalloc for list elements */ 
    if((out->a2dr = (Aln2DisplayResource ** ) ckcalloc (len,sizeof(Aln2DisplayResource *))) == NULL) {  
      warn("Warning, ckcalloc failed in Aln2Display_alloc_len"); 
      return NULL;   
      }  
    out->res_len = 0;    
    out->res_maxlen = len;   


    if((out->s2d = (Sequence2Display ** ) ckcalloc (len,sizeof(Sequence2Display *))) == NULL)    {  
      warn("Warning, ckcalloc failed in Aln2Display_alloc_len"); 
      return NULL;   
      }  
    out->len = 0;    
    out->maxlen = len;   


    return out;  
}    


/* Function:  hard_link_Aln2Display(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj [UNKN ] Object to be hard linked [Aln2Display *]
 *
 * Return [UNKN ]  Undocumented return value [Aln2Display *]
 *
 */
Aln2Display * hard_link_Aln2Display(Aln2Display * obj) 
{
    if( obj == NULL )    {  
      warn("Trying to hard link to a Aln2Display object: passed a NULL object"); 
      return NULL;   
      }  
    obj->dynamite_hard_link++;   
    return obj;  
}    


/* Function:  Aln2Display_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [Aln2Display *]
 *
 */
Aln2Display * Aln2Display_alloc(void) 
{
    Aln2Display * out;  /* out is exported at end of function */ 


    /* call ckalloc and see if NULL */ 
    if((out=(Aln2Display *) ckalloc (sizeof(Aln2Display))) == NULL)  {  
      warn("Aln2Display_alloc failed "); 
      return NULL;  /* calling function should respond! */ 
      }  
    out->dynamite_hard_link = 1; 
#ifdef PTHREAD   
    pthread_mutex_init(&(out->dynamite_mutex),NULL);     
#endif   
    out->name = NULL;    
    out->a2dr = NULL;    
    out->res_len = out->res_maxlen = 0;  
    out->s2d = NULL; 
    out->len = out->maxlen = 0;  
    out->other = NULL;   


    return out;  
}    


/* Function:  free_Aln2Display(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [Aln2Display *]
 *
 * Return [UNKN ]  Undocumented return value [Aln2Display *]
 *
 */
Aln2Display * free_Aln2Display(Aln2Display * obj) 
{
    int return_early = 0;    
    int i;   


    if( obj == NULL) {  
      warn("Attempting to free a NULL pointer to a Aln2Display obj. Should be trappable");   
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
    if( obj->a2dr != NULL)   {  
      for(i=0;i<obj->res_len;i++)    {  
        if( obj->a2dr[i] != NULL)    
          free_Aln2DisplayResource(obj->a2dr[i]);    
        }  
      ckfree(obj->a2dr); 
      }  
    if( obj->s2d != NULL)    {  
      for(i=0;i<obj->len;i++)    {  
        if( obj->s2d[i] != NULL) 
          free_Sequence2Display(obj->s2d[i]);    
        }  
      ckfree(obj->s2d);  
      }  
    if( obj->other != NULL)  
      free_Sequence2Display(obj->other);     


    ckfree(obj); 
    return NULL; 
}    



#ifdef _cplusplus
}
#endif
