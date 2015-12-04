
#ifndef DATABASE_HEADER_FILE
#define DATABASE_HEADER_FILE

/*
 * Enum for database reload success/failure
 */


enum return_status {
  DB_RETURN_UNKNOWN = 143,
  DB_RETURN_OK,
  DB_RETURN_ERROR,
  DB_RETURN_END };

typedef int DB_Return_Type;

/*
 * Enum for search success/failure
 */

enum search_return {
  SEARCH_UNKNOWN = 152,
  SEARCH_ERROR,
  SEARCH_OK };

typedef int Search_Return_Type;


#endif
