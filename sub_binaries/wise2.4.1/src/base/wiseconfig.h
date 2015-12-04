/************************************************************************/
/* WiseTools, Release 1.0. This software is copyright                   */
/* Ewan Birney (c), 1995.                                               */
/* Anyone can distribute/compile and use this software                  */
/* free of charge as long as this notice remains intact                 */
/*                                                                      */
/* Any alterations should be cleared by Ewan Birney.                    */
/* Email birney@molbiol.ox.ac.uk                                        */
/* URL   http://www.molbiol.ox.ac.uk/www/users/birney/wise/topwise.html */
/*                                                                      */
/* This software is provided as is, and the author does not accept any  */
/* liability for its use or performance :)                              */
/************************************************************************/


#ifndef EWANCONFIG
#define EWANCONFIG

#include <stdio.h>

boolean read_set_config(void);

char ** filename_list(int * retval);

int config_is_key(char * key);
	/* if keyword is present in list then return 1 else return 0 */

char * config_string_from_key(char * key);
	/* if keyword is present then return string: if no keyword */
	/* will return NULL: however, string could be NULL: a NULL */
	/* return could mean either no keyword or NULL string      */
	/* associated with keyword. Test with config_is_key        */ 

char * config_single_from_key(char *  key);
	/* if keyword is present, returns the a string whith the first word */
	/* all packed up:  string has  been alloced... so free it if you like*/

int config_number_from_key(char * key,int * pass);
	/* will put the value of the first word associated with key */
	/* into *pass using atoi: will return 0 if success, and     */
	/* 1 if not */

int read_config_file(char * filename);
	/* reads in a config file */

/* internal functions: don't call 'em !!! */

char * read_config_string(char * base,char * string,FILE * ifp);
int index_from_keyword(char * key);

#endif /* of EWANCONFIG */

