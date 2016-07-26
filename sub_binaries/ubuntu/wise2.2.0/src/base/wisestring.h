#ifndef DYNAMITEwisestringHEADERFILE
#define DYNAMITEwisestringHEADERFILE
#ifdef _cplusplus
extern "C" {
#endif
#include "wisebase.h"



/**********************************************/
/* useful macros to put into fprintf lines... */
/* Makes sure you don't trash memory etc      */
/**********************************************/
#define CHECKSTRING(str) (str == NULL ? "NullString" : str)
#define CKS CHECKSTRING

/**********************************************/
/* useful standard strings for parsing        */
/* spacestr is general whitespace             */
/* breakstr is general non alpha/num          */
/* used alot in breakstring                   */
/**********************************************/
#define spacestr " \t\n\0"
#define breakstr "!\"#$%^&*()-+={}[]@';:?/.,\\|~` \n\t"

/**********************************************/
/* Not the nicest of macros. Stay away from it*/
/**********************************************/ 
#define NEXTWORD strtok(NULL,spacestr)



    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  get_number_from_slashed_string(qstr,slashstr)
 *
 * Descrip:    handy: pass a string like "xxx/yyy/zzz"    
 *             and a query like yyy, will return 2        
 *             (ie, maps positions in slashed string with 
 *             a number from start). Returns -1 if none   
 *
 *
 * Arg:            qstr [UNKN ] Undocumented argument [char *]
 * Arg:        slashstr [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_get_number_from_slashed_string(char * qstr,char * slashstr);
#define get_number_from_slashed_string Wise2_get_number_from_slashed_string


/* Function:  is_integer_string(string,val)
 *
 * Descrip:    checks that strings are ints or doubles    
 *             and then converts, storing value in val    
 *             if val == NULL, will not store (!)        
 *                                                         
 *             Does use sensible library functions        
 *             strtol...                          
 *
 *
 * Arg:        string [UNKN ] Undocumented argument [char *]
 * Arg:           val [UNKN ] Undocumented argument [int *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_is_integer_string(char * string,int * val);
#define is_integer_string Wise2_is_integer_string


/* Function:  is_double_string(string,val)
 *
 * Descrip:     checks that strings are doubles    
 *              and then converts, storing value in val    
 *               if val == NULL, will not store (!)        
 *                                                         
 *              Does use sensible library functions        
 *              strtod
 *
 *
 * Arg:        string [UNKN ] Undocumented argument [char *]
 * Arg:           val [UNKN ] Undocumented argument [double *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_is_double_string(char * string,double * val);
#define is_double_string Wise2_is_double_string


/* Function:  compress_space_around_punc(buffer,punc,space)
 *
 * Descrip: No Description
 *
 * Arg:        buffer [UNKN ] Undocumented argument [char *]
 * Arg:          punc [UNKN ] Undocumented argument [char *]
 * Arg:         space [UNKN ] Undocumented argument [char *]
 *
 */
void Wise2_compress_space_around_punc(char * buffer,char * punc,char * space);
#define compress_space_around_punc Wise2_compress_space_around_punc


/* Function:  striptoprint(line)
 *
 * Descrip:    useful strip functions to remove nasty chars
 *             does not allocate memory, simply uses       
 *             given memory, but returns the line pointer
 *             so you can use it in nested function calls                                
 *
 *
 * Arg:        line [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [char *]
 *
 */
char * Wise2_striptoprint(char * line);
#define striptoprint Wise2_striptoprint


/* Function:  stringalloc_next_quoted_string(buffer)
 *
 * Descrip:    takes str's of type <garbage> "xxxxx"   
 *             and gives back xxxxxx .             
 *             stringalloc'd piece so make sure you free it
 *
 *
 * Arg:        buffer [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [char *]
 *
 */
char * Wise2_stringalloc_next_quoted_string(char * buffer);
#define stringalloc_next_quoted_string Wise2_stringalloc_next_quoted_string


/* Function:  strwhitestartcmp(line,str,whitespace)
 *
 * Descrip:    sees if line starts with str, ignoring whitespace
 *
 *             returns 0 if they match, to look like strcmp
 *
 *
 * Arg:              line [UNKN ] Undocumented argument [char *]
 * Arg:               str [UNKN ] Undocumented argument [char *]
 * Arg:        whitespace [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_strwhitestartcmp(char * line,char * str,char * whitespace);
#define strwhitestartcmp Wise2_strwhitestartcmp


/* Function:  strwordcmp(buf,str,space)
 *
 * Descrip:    sees if buf matches str\s in perl regex, ie
 *             a word match
 *
 *             space defined \s
 *
 *             returns 0 if they match to look like strcmp
 *
 *
 * Arg:          buf [UNKN ] Undocumented argument [char *]
 * Arg:          str [UNKN ] Undocumented argument [char *]
 * Arg:        space [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_strwordcmp(char * buf,char * str,char * space);
#define strwordcmp Wise2_strwordcmp


/* Function:  strstartcmp(buf,str)
 *
 * Descrip:    sees if buf starts with str.
 *
 *             returns 0 if so, to mimic strcmp
 *
 *
 * Arg:        buf [UNKN ] Undocumented argument [char *]
 * Arg:        str [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_strstartcmp(char * buf,char * str);
#define strstartcmp Wise2_strstartcmp


/* Function:  print_numbered_line(num,ofp)
 *
 * Descrip:    prints lines like _1_________
 *
 *
 * Arg:        num [UNKN ] Undocumented argument [int]
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
void Wise2_print_numbered_line(int num,FILE * ofp);
#define print_numbered_line Wise2_print_numbered_line


/* Function:  print_line(ofp)
 *
 * Descrip:    prints _______________ (70 chars)
 *
 *
 * Arg:        ofp [UNKN ] Undocumented argument [FILE *]
 *
 */
void Wise2_print_line(FILE * ofp);
#define print_line Wise2_print_line


/* Function:  chop_newline(str)
 *
 * Descrip:    removes trailing newline if present
 *
 *
 * Arg:        str [UNKN ] Undocumented argument [char *]
 *
 */
void Wise2_chop_newline(char * str);
#define chop_newline Wise2_chop_newline


/* Function:  good_datastring_fromend(str)
 *
 * Descrip:    Tries to find the last 'database name' 
 *             type string from a string.
 *
 *             Does not allocate memory
 *
 *
 * Arg:        str [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [char *]
 *
 */
char * Wise2_good_datastring_fromend(char * str);
#define good_datastring_fromend Wise2_good_datastring_fromend


/* Function:  estrcasecmp(one,two)
 *
 * Descrip:    returns strcmp on the captilalised
 *             one and two bufferers (doesn't touch them).
 *
 *
 * Arg:        one [UNKN ] Undocumented argument [char *]
 * Arg:        two [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_estrcasecmp(char *  one,char *  two);
#define estrcasecmp Wise2_estrcasecmp


/* Function:  number_from_quoted_equality(s)
 *
 * Descrip:    supposedly gets a number from a
 *             string like xxx="12"
 *
 *
 * Arg:        s [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_number_from_quoted_equality(char * s);
#define number_from_quoted_equality Wise2_number_from_quoted_equality


/* Function:  number_from_equality(string)
 *
 * Descrip:    supposedly gets a number from
 *             xxxx=12
 *
 *
 * Arg:        string [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [int]
 *
 */
int Wise2_number_from_equality(char * string);
#define number_from_equality Wise2_number_from_equality


/* Function:  string_from_quoted_equality(string)
 *
 * Descrip:    gets the string from xxx="yyy". Returns
 *             yyy allocated, and messes around with string
 *
 *
 * Arg:        string [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [char *]
 *
 */
char * Wise2_string_from_quoted_equality(char * string);
#define string_from_quoted_equality Wise2_string_from_quoted_equality


/* Function:  string_from_charred_equality(string,quote)
 *
 * Descrip:    gets the string from xxx="yyy" ,where " comes
 *             from the quote argument Returns
 *             yyy allocated, and messes around with string
 *
 *
 * Arg:        string [UNKN ] Undocumented argument [char *]
 * Arg:         quote [UNKN ] Undocumented argument [char]
 *
 * Return [UNKN ]  Undocumented return value [char *]
 *
 */
char * Wise2_string_from_charred_equality(char * string,char quote);
#define string_from_charred_equality Wise2_string_from_charred_equality


/* Function:  breakstring(string,parsestr)
 *
 * Descrip:    A call to /breakstring_protect(string,parsestr,"\"")
 *
 *
 * Arg:          string [UNKN ] Undocumented argument [char *]
 * Arg:        parsestr [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [char **]
 *
 */
char ** Wise2_breakstring(char * string,char * parsestr);
#define breakstring Wise2_breakstring


/* Function:  breakstring_protect(string,parsestr,strpair)
 *
 * Descrip:    will parse out words in string using parse  
 *             as white space, like strtok, but strings    
 *             enclosed in characters from strpair are not 
 *             taken in parsed form.                       
 *             breakstring =                               
 *             breakstring_protect(string,parse,"\"");     
 *             hence will not break in double quotes       
 *                                                          
 *             unlike strtok they return char **           
 *             which is a list of char * of words          
 *             the last being NULL'd                       
 *                                                          
 *                                                          
 *             They returned an alloc'd char ** which you  
 *             are expected to free. Standard idiom is     
 *             base=brk=breakstring(buffer,spacestr)       
 *               (NB - spacestr #defin'd above )           
 *               ... do stuff using brk                     
 *               eg *brk = first word                       
 *                  *(++brk) = next word                    
 *                last word = NULL                          
 *              cleanup by ckfree(base)                     
 *
 *
 * Arg:          string [UNKN ] Undocumented argument [char *]
 * Arg:        parsestr [UNKN ] Undocumented argument [char *]
 * Arg:         strpair [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [char **]
 *
 */
char ** Wise2_breakstring_protect(char * string,char * parsestr,char * strpair);
#define breakstring_protect Wise2_breakstring_protect


/* Function:  strip_quote_chars(string,quote)
 *
 * Descrip:    removes chars in quote
 *
 *
 * Arg:        string [UNKN ] Undocumented argument [char *]
 * Arg:         quote [UNKN ] Undocumented argument [char *]
 *
 */
void Wise2_strip_quote_chars(char * string,char * quote);
#define strip_quote_chars Wise2_strip_quote_chars


/* Function:  padstring(buffer,string,maxlen)
 *
 * Descrip:    copies string into buffer, and if under maxlen,
 *             adds spaces. Does *not* put in '\0'
 *
 *
 * Arg:        buffer [UNKN ] Undocumented argument [char *]
 * Arg:        string [UNKN ] Undocumented argument [char *]
 * Arg:        maxlen [UNKN ] Undocumented argument [int]
 *
 */
void Wise2_padstring(char * buffer,char * string,int maxlen);
#define padstring Wise2_padstring


/* Function:  capitalise(word)
 *
 * Descrip:    toupper's each char in word
 *
 *
 * Arg:        word [UNKN ] Undocumented argument [char *]
 *
 */
void Wise2_capitalise(char * word);
#define capitalise Wise2_capitalise


/* Function:  show_line(line,max,*ofp)
 *
 * Descrip:    This shouws line putting a new line in every max
 *             chars, not minding word boundaries
 *
 *
 * Arg:        line [UNKN ] Undocumented argument [char *]
 * Arg:         max [UNKN ] Undocumented argument [int]
 * Arg:        *ofp [UNKN ] Undocumented argument [FILE]
 *
 */
void Wise2_show_line(char * line,int max,FILE *ofp);
#define show_line Wise2_show_line


/* Function:  show_text(line,max,*ofp)
 *
 * Descrip:    This shouws line putting a new line in every max
 *             chars, *minding* word boundaries
 *
 *
 * Arg:        line [UNKN ] Undocumented argument [char *]
 * Arg:         max [UNKN ] Undocumented argument [int]
 * Arg:        *ofp [UNKN ] Undocumented argument [FILE]
 *
 */
void Wise2_show_text(char * line,int max,FILE *ofp);
#define show_text Wise2_show_text


/* Function:  second_word_alloc(str,space)
 *
 * Descrip:    returns the second word alloc'd for
 *
 *             xxx yyyy
 *
 *             returns yyyy.
 *
 *
 * Arg:          str [UNKN ] Undocumented argument [char *]
 * Arg:        space [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [char *]
 *
 */
char * Wise2_second_word_alloc(char * str,char * space);
#define second_word_alloc Wise2_second_word_alloc


/* Function:  stringallocf(str,)
 *
 * Descrip:    Don't use this
 *
 *             sprintf's then allocs
 *
 *
 * Arg:        str [UNKN ] Undocumented argument [char *]
 * Arg:            [UNKN ] Undocumented argument [.]
 *
 * Return [UNKN ]  Undocumented return value [char *]
 *
 */
char * Wise2_stringallocf(char * str,...);
#define stringallocf Wise2_stringallocf


/* Function:  stringalloc(c)
 *
 * Descrip:    returns the allocated copy of c.
 *             Usually called strdup in other packages
 *
 *
 * Arg:        c [UNKN ] Undocumented argument [char *]
 *
 * Return [UNKN ]  Undocumented return value [char *]
 *
 */
char * Wise2_stringalloc(char * c);
#define stringalloc Wise2_stringalloc


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
boolean Wise2_only_whitespace(char * str,char * space);
#define only_whitespace Wise2_only_whitespace
boolean Wise2_looks_like_vms(const char * str);
#define looks_like_vms Wise2_looks_like_vms
boolean Wise2_looks_like_unix(const char * str);
#define looks_like_unix Wise2_looks_like_unix
char * Wise2_strend(char * bu,char * se);
#define strend Wise2_strend
char * Wise2_string_before_equality(char * string);
#define string_before_equality Wise2_string_before_equality
char * Wise2_sub_string(char * into,char * key,char * sub);
#define sub_string Wise2_sub_string

#ifdef _cplusplus
}
#endif

#endif
