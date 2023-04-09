#include <stdio.h>

int myreadline( char **buffer , int *len , FILE* f ) ;
char **mysplitline( int *nwords , char *line , int len ) ;







int myreadline( char **buffer , int *len , FILE* f ) {
  char c ;
  int nchars = 0 ;
  int KEEPREAD = 1 ;

  c = getc(f) ;
  KEEPREAD = ( c != '\n' && c != EOF ) ;
  while ( KEEPREAD ) {
    if ( nchars >= *len-1 ) {
      *len += 100 ;
      (*buffer) = (char*)realloc( *buffer, (*len)*sizeof(char) ) ;
      (*buffer)[*len-1] = '\0' ;
    }
    (*buffer)[nchars] = c ;
    nchars ++ ;
    c = getc(f) ;
    KEEPREAD = ( c != '\n' && c != EOF ) ;
  }

  if (c == EOF) return 0 ;
  else return 1 ;
}


char **mysplitline( int *nwords , char *line , int len ) {
  int maxnchar = 0 , i = 0 , nchar = 0 ;
  char **words = NULL ;

  if ( len > 0 && line != NULL ) {
    if( line[0] != ' ' && line[0] != '\t' ) *nwords = 1 ;
    for( i=1 ; i<len ; i++ ) {
      if( ( line[i] != ' ' && line[i] != '\t' ) && ( line[i-1] == ' ' || line[i-1] == '\t' ) ) {
	(*nwords) ++ ;
	if ( maxnchar < nchar ) maxnchar = nchar ;
	nchar = 0 ;
      }
      nchar ++ ;
    }
    words = (char**)calloc( *nwords , sizeof(char*) ) ;

    int wordlen = 0 , nw = 0 ;
    for( int i=0; i<len; i++ ) {
      if( line[i] != ' ' && line[i] != '\t' ) {
	for( wordlen = 0; ( (wordlen+i<len) && ( line[i+wordlen]!=' ' && line[i+wordlen]!='\t' ) ); wordlen++ ) ;
	words[nw] = (char*)calloc( wordlen+1 , sizeof(char) ) ;
	for( int j=i; j<i+wordlen; j++ ) words[nw][j-i] = line[j] ;
	words[nw][wordlen] = '\0' ;
	nw ++ ;
	i = i + wordlen - 1 ;
      }
    }
  }

  return words ;
}
