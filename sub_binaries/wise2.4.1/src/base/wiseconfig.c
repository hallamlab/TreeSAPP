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

#include "wisebase.h"

#define MAXFIG 512
#define MAXLINE 512

typedef struct {
	char * keyword;
	char * string;
	char * filename; /* this is not alloc'd */
	} config;

static config list[MAXFIG];
static int maxnum=0;

static char * filenamelist[64];
static int noffiles=0;

char * get_usermailname(void)
{
	char * out;
	
	out = config_single_from_key("usermailname");

	if( out != NULL) return out;

	out = getenv("user");

	if( out != NULL )
		return stringalloc(out);

	out = getenv("USER");

	if( out != NULL)
		return stringalloc(out);

	return NULL;
}



boolean is_config_system(void)
{
	return (noffiles == 0 ? FALSE : TRUE);
}

char ** filename_list(int * retval)
{
	*retval=noffiles;
	return filenamelist;
}


boolean read_set_config(void)
{
	char * runner;
	

	if( (runner=getenv("WISESYSTEMFILE")) != NULL)
		read_config_file(runner);
	else	read_config_file("wise.cfg");
	
	if( (runner=getenv("WISEPERSONALFILE")) != NULL)
		read_config_file(runner);

	if( is_config_system() == FALSE)
		{
		  log_full_error(WARNING,0,"Unable to find config files, "
				 "going to try a variety of other filenames. Default has changed "
				 "to wise.cfg for all systems");
		read_config_file("wise.cfg");
		read_config_file(".wisecfg");
		if( is_config_system() == FALSE)
		  log_full_error(WARNING,0,"Unable to find any config files "
				 "will bug out now");
		return TRUE;
		}
	return TRUE;
}	

int index_from_keyword(char * key)
{
	register int i;
	
	if(noffiles == 0)
		{
		log_full_error(PEDANTIC,0,"Tried to pull keyword in config "
"but no file has been read");
		return -1;
		}

	for(i=0;i<maxnum;i++)
		if(strcmp(key,list[i].keyword) == 0)
			return i;
	return -1;
}

int config_is_key(char * key)
{
	if( index_from_keyword(key) != -1)
		return 1;
	else return 0;
}

char * config_string_from_key(char * key)
{
	register int i;
	
	if( (i=index_from_keyword(key)) == -1)
		return NULL;
	else return list[i].string;
}

char * config_single_from_key(char * key)
{
	char *  word;
	register int i;
	char buffer[MAXLINE];


	if(key == NULL)
		{
		log_full_error(WARNING,0,"Passed a NULL key into config_single_from_key...giving back nothing!");
		return NULL;
		}


	if( (i=index_from_keyword(key)) == -1)
		return NULL;

	if(list[i].string == NULL)
		{
		return NULL;
		}

	strcpy(buffer,list[i].string);

	word=strtok(buffer," \n\t\0");
	return (word != NULL ? stringalloc (word) : NULL);
}

int config_number_from_key(char * key,int * passed)
{
	char * runner;
	register int i;
	char buffer[MAXLINE];

	if( (i=index_from_keyword(key)) == -1)
		return 0;
	
	strcpy(buffer,list[i].string);
	if( (runner=strtok(buffer," \n\t")) == NULL)
		return 1;
	*passed=atoi(runner);
	return 0;
}

int read_config_file(char * filename)
{
	char buffer[MAXLINE];
	char * runner;
	FILE * ifp;
	int index=maxnum; /* maxnum declared as static above */

	if( (ifp=openfile(filename,"r")) == NULL) 
/* don't forget, openfile will append wiseconfig and personalconfig */
/* onto the filename */
		{

/* banal user checks for stupid users */

#ifdef CAREFUL

	if( looks_like_vms(filename))
	log_full_error(WARNING,0,"Filename %s looks like a VMS file "
"You may have compiled without #define UNIX",filename);
	
	if( looks_like_unix(filename))
	log_full_error(WARNING,0,"Filename %s looks like a unix file "
"You may have compiled wrongly",filename);
	
#endif
		return 1;
		}


	filenamelist[noffiles++]=stringalloc(filename);


	while( fgets(buffer,MAXLINE,ifp) != NULL)
		{
		if(buffer[0] == '!')
			continue;
		if( (runner=strtok(buffer," \n\t")) == NULL)
			continue;
		list[index].keyword=stringalloc(runner);
		list[index].filename=filenamelist[noffiles-1];
	/* the -1 on above line is 'cause we've post++ the noffiles */
	/* above */
		if( (runner=strtok(NULL,"\n\0")) == NULL)
			{
			list[index++].string=NULL;
			continue;
			}
		else	list[index++].string=read_config_string(NULL,runner,ifp);
	/* NB watch the postfix ++ on index above */

		}
	maxnum=index;
	fclose(ifp);
	return 0;
}

/* this bottem bit has to deal with recursive lines, UGH */

char * read_config_string(char * base,char * string,FILE * ifp)
{
	char * runner;
	char * backrunner;
	char buffer[MAXLINE];

	runner=string+strlen(string);
	runner--;
	if( *runner == '\n')
		{
		*runner='\0';
		runner--;
		}


	if(base != NULL)
		{
		runner=(char *) malloc (sizeof(char) *( strlen (base) + 
strlen (string) +1 ));
		if (runner == NULL)
	log_full_error(FATAL,0,"Fatal alloc error in reading config string, otta here");
		strcpy(runner,base);
		strcat(runner,string);
		}

	else	{
		runner=(char *) malloc (sizeof(char) * (strlen (string)+1));
		if (runner == NULL)
	log_full_error(FATAL,0,"Fatal alloc error in reading config string, otta here (2)");
		strcpy(runner,string);
		}
	backrunner=runner+strlen(runner);
	backrunner--;

	if( *backrunner == '\\')
		{
		*backrunner=' ';
		if( fgets(buffer,MAXLINE,ifp) != NULL)
			return read_config_string(runner,buffer,ifp);
		else	log_full_error(WARNING,0,"bad end backslash, encountered EOF");
		}
	return runner;
}

