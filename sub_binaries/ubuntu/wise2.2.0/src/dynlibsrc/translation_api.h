

/* Functions that create, manipulate or act on Translation
 *
 * Wise2_get_Protein_from_Translation
 * Wise2_hard_link_Translation
 * Wise2_Translation_alloc
 * Wise2_replace_start_Translation
 * Wise2_access_start_Translation
 * Wise2_replace_end_Translation
 * Wise2_access_end_Translation
 * Wise2_replace_parent_Translation
 * Wise2_access_parent_Translation
 * Wise2_replace_protein_Translation
 * Wise2_access_protein_Translation
 * Wise2_free_Translation [destructor]
 *
 */

/* API for object Translation */
/* Function:  Wise2_get_Protein_from_Translation(ts,ct)
 *
 * Descrip:    Gets the protein
 *
 *
 * Arg:        ts           translation [Wise2_Translation *]
 * Arg:        ct           codon table to use [Wise2_CodonTable *]
 *
 * Returns Protein sequence [Wise2_Protein *]
 *
 */
Wise2_Protein * Wise2_get_Protein_from_Translation( Wise2_Translation * ts,Wise2_CodonTable * ct);

/* Function:  Wise2_hard_link_Translation(obj)
 *
 * Descrip:    Bumps up the reference count of the object
 *             Meaning that multiple pointers can 'own' it
 *
 *
 * Arg:        obj          Object to be hard linked [Wise2_Translation *]
 *
 * Returns Undocumented return value [Wise2_Translation *]
 *
 */
Wise2_Translation * Wise2_hard_link_Translation( Wise2_Translation * obj);

/* Function:  Wise2_Translation_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Returns Undocumented return value [Wise2_Translation *]
 *
 */
Wise2_Translation * Wise2_Translation_alloc();

/* Function:  Wise2_replace_start_Translation(obj,start)
 *
 * Descrip:    Replace member variable start
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_Translation *]
 * Arg:        start        New value of the variable [int]
 *
 * Returns member variable start [boolean]
 *
 */
boolean Wise2_replace_start_Translation( Wise2_Translation * obj,int start);

/* Function:  Wise2_access_start_Translation(obj)
 *
 * Descrip:    Access member variable start
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_Translation *]
 *
 * Returns member variable start [int]
 *
 */
int Wise2_access_start_Translation( Wise2_Translation * obj);

/* Function:  Wise2_replace_end_Translation(obj,end)
 *
 * Descrip:    Replace member variable end
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_Translation *]
 * Arg:        end          New value of the variable [int]
 *
 * Returns member variable end [boolean]
 *
 */
boolean Wise2_replace_end_Translation( Wise2_Translation * obj,int end);

/* Function:  Wise2_access_end_Translation(obj)
 *
 * Descrip:    Access member variable end
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_Translation *]
 *
 * Returns member variable end [int]
 *
 */
int Wise2_access_end_Translation( Wise2_Translation * obj);

/* Function:  Wise2_replace_parent_Translation(obj,parent)
 *
 * Descrip:    Replace member variable parent
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_Translation *]
 * Arg:        parent       New value of the variable [Wise2_Transcript *]
 *
 * Returns member variable parent [boolean]
 *
 */
boolean Wise2_replace_parent_Translation( Wise2_Translation * obj,Wise2_Transcript * parent);

/* Function:  Wise2_access_parent_Translation(obj)
 *
 * Descrip:    Access member variable parent
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_Translation *]
 *
 * Returns member variable parent [Wise2_Transcript *]
 *
 */
Wise2_Transcript * Wise2_access_parent_Translation( Wise2_Translation * obj);

/* Function:  Wise2_replace_protein_Translation(obj,protein)
 *
 * Descrip:    Replace member variable protein
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_Translation *]
 * Arg:        protein      New value of the variable [Wise2_Protein *]
 *
 * Returns member variable protein [boolean]
 *
 */
boolean Wise2_replace_protein_Translation( Wise2_Translation * obj,Wise2_Protein * protein);

/* Function:  Wise2_access_protein_Translation(obj)
 *
 * Descrip:    Access member variable protein
 *             For use principly by API functions
 *
 *
 * Arg:        obj          Object holding the variable [Wise2_Translation *]
 *
 * Returns member variable protein [Wise2_Protein *]
 *
 */
Wise2_Protein * Wise2_access_protein_Translation( Wise2_Translation * obj);

/* This is the destructor function, ie, call this to free object*/
/* Function:  Wise2_free_Translation(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj          Object that is free'd [Wise2_Translation *]
 *
 * Returns Undocumented return value [Wise2_Translation *]
 *
 */
Wise2_Translation * Wise2_free_Translation( Wise2_Translation * obj);

