/*
 * TAPENADE Automatic Differentiation Engine
 * Copyright (C) 1999-2021 Inria
 * See the LICENSE.md file in the project root for more information.
 *
 */

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "adContext.h"
#include "adComplex.h"
#include "json.h"


typedef enum dbad_mode { TANGENT=1, ADJOINT=2 } dbad_mode_t;
static int dbad_mode, dbad_phase ;
static double dbad_ddeps = 1.e-6 ;
static double dbad_seed = 0.137 ;
static double dbad_currentSeed = 0.0 ;
static double* dbad_seedIndep = NULL;
static int     dbad_seedIndepCnt = 0;
static int     dbad_curSeedIndep = 0;
static double* dbad_seedDep = NULL;
static int     dbad_seedDepCnt = 0;
static int     dbad_curSeedDep = 0;
static double dbad_condensed_val, dbad_condensed_tgt, dbad_condensed_adj ;
const static char* jsonfile = "dbad.json"; //read eps/seed/... from this file (if present)


/***************************************************/
/*        JSON file reading
/**************************************************/
const char *read_file(const char *path) {

  FILE *file = fopen(path, "r");

  if (file == NULL) {
      //fprintf(stderr, "Expected file \"%s\" not found", path);
      return NULL;
  }

  fseek(file, 0, SEEK_END);
  long len = ftell(file);
  fseek(file, 0, SEEK_SET);
  char *buffer = malloc(len + 1);

  if (buffer == NULL) {
    fprintf(stderr, "Unable to allocate memory for file");
    fclose(file);
    return NULL;
  }

  fread(buffer, 1, len, file);
  buffer[len] = '\0';

  fclose(file);

  return (const char *)buffer;
}


const char* json_type_str(typed(json_element_type) type) {
    switch(type) {
        case JSON_ELEMENT_TYPE_NULL:
            return "null";
            break;
        case JSON_ELEMENT_TYPE_BOOLEAN:
            return "bool";
            break;
        case JSON_ELEMENT_TYPE_STRING:
            return "string";
            break;
        case JSON_ELEMENT_TYPE_NUMBER:
            return "number";
            break;
        case JSON_ELEMENT_TYPE_OBJECT:
            return "object";
            break;
        case JSON_ELEMENT_TYPE_ARRAY:
            return "array";
            break;
        default:
            fprintf(stderr, "unexpected json_element_type");
            return NULL;
            break;
    }
}


void json_elem_get_double(typed(json_element) elem, double* value, int* ok) {
    *ok = 0;
    if( elem.type!=JSON_ELEMENT_TYPE_NUMBER ) {
        fprintf(stderr, "unexpected type -->%s<-- when reading double\n", json_type_str(elem.type));
        return;
    }
    else {
        if( elem.value.as_number.type==JSON_NUMBER_TYPE_DOUBLE ) {
            *value = elem.value.as_number.value.as_double;
            *ok = 1;
        }
        else {//JSON_NUMBER_TYPE_LONG
            *value = (double) elem.value.as_number.value.as_long;
            *ok = 1;
        }
    }
}


void jsonobj_read_double(typed(json_object) *obj, const char* key, double* value, int* ok) {
    *ok = 0;
    result(json_element) key_result = json_object_find(obj, key);
    if( ! result_is_err(json_element)(&key_result)) {
        typed(json_element) elem = result_unwrap(json_element)(&key_result);
        json_elem_get_double(elem, value, ok);
    }
    else {
        // key not found or other error
    }
}


double* jsonobj_read_double_lst(typed(json_object) *obj, const char* key, int* cnt) {
    *cnt = 0;
    double* value_lst = NULL;
    double cur_value;
    int ok;

    result(json_element) key_result = json_object_find(obj, key);

    if( result_is_err(json_element)(&key_result)) {
        return value_lst;
    }
    //
    typed(json_element) elem = result_unwrap(json_element)(&key_result);
    if( elem.type!=JSON_ELEMENT_TYPE_ARRAY ) {
        fprintf(stderr, "unexpected type -->%s<-- when reading value list\n",
                json_type_str(elem.type));
    }
    else {
        *cnt = elem.value.as_array->count;
        value_lst = malloc(*cnt*sizeof(double));
        for(int j = 0; j <*cnt; ++j) {
            typed(json_element) cur_elem = elem.value.as_array->elements[j];
            json_elem_get_double(cur_elem, &cur_value, &ok);
            if( ok ) {
                value_lst[j] = cur_value;
            }
            else {
                //todo message
                free( value_lst );
                value_lst = NULL;
                break;
            }
        }
    }

    return value_lst;
}


void json_dbad_apply() {
    // epsilon/seed potentially may be provided via json file,
    // which takes precedence if present
    const char* json_str = read_file(jsonfile);
    double value;
    double* value_lst;
    int rc;
    if( json_str!=NULL ) {
        // parse file content to internal json format
        result(json_element) element_result = json_parse(json_str);
        free((void *)json_str); // dispose character buffer

        if (result_is_err(json_element)(&element_result)) {
            typed(json_error) error = result_unwrap_err(json_element)(&element_result);
            fprintf(stderr, "Error parsing JSON: %s\n", json_error_to_string(error));
            return;
        }
        else {
            if( dbad_phase==99 ) {
                printf("...start reading settings from file ***%s***\n", jsonfile);
            }

            // Extract the data
            typed(json_element) element = result_unwrap(json_element)(&element_result);
            typed(json_object) *obj = element.value.as_object;

            // check 'eps' (for divided differences)
            jsonobj_read_double(obj, "eps", &value, &rc);
            if( rc ) {
                dbad_ddeps = value;
                if( dbad_phase==99 ) {
                    printf("......setting dbad_ddeps=%15.10e\n", dbad_ddeps);
                }
            }

            // check 'seed' (used for randomised direction)
            jsonobj_read_double(obj, "seed", &value, &rc);
            if( rc ) {
                dbad_seed = value;
                if( dbad_phase==99 ) {
                    printf("......setting dbad_seed=%15.10e\n", dbad_seed);
                }
            }

            // seed for independents
            value_lst = jsonobj_read_double_lst(obj, "indep_seed", &rc);
            if( value_lst!=NULL ) {
                dbad_seedIndepCnt = rc;
                dbad_seedIndep = value_lst;
            }
            
            // seed for dependents
            value_lst = jsonobj_read_double_lst(obj, "dep_seed", &rc);
            if( value_lst!=NULL ) {
                dbad_seedDepCnt = rc;
                dbad_seedDep = value_lst;
            }

            // dispose memory
            json_free(&element);
        }
    }
}//json_dbad_apply
/***************************************************/
/*        end JSON stuff
/**************************************************/


double dbad_nextRandom() {
  dbad_currentSeed += dbad_seed ;
  if (dbad_currentSeed>=1.0) dbad_currentSeed-=1.0 ;
  /* Return a value in range [1.0 2.0[ */
  return dbad_currentSeed+1.0 ;
}


double dbad_nextIndep() {
    if( dbad_seedIndep==NULL ) {
        return dbad_nextRandom();
    }
    else {
        if( dbad_curSeedIndep<dbad_seedIndepCnt ) {
            return dbad_seedIndep[dbad_curSeedIndep++];
        }
        else {
            return 0;
        }
    }
}


double dbad_nextDep() {
    if( dbad_seedDep==NULL ) {
        return dbad_nextRandom();
    }
    else {
        if( dbad_curSeedDep<dbad_seedDepCnt ) {
            return dbad_seedDep[dbad_curSeedDep++];
        }
        else {
            return 0;
        }
    }
}


void dbad_set_phase(dbad_mode_t mode) {
    char* phase = getenv("DBAD_PHASE") ;
    if(mode==ADJOINT) {
        if (phase==NULL) {
            dbad_phase = 0 ;
        } else if (strcmp(phase,"99")==0) {
            dbad_phase = 99 ;
        } else {
            dbad_phase = 0 ;
        }
    }
    else {//TANGENT
        if (phase==NULL) {
            printf("Please set DBAD_PHASE environment variable to 1 (perturbed) or 2 (tangent)\n") ;
            exit(0) ;
        } else if (strcmp(phase,"2")==0) {
            dbad_phase = 2 ;
        } else if (strcmp(phase,"1")==0) {
            dbad_phase = 1 ;
        } else if (strcmp(phase,"99")==0) {
            dbad_phase = 99 ;
        } else {
            printf("DBAD_PHASE environment variable must be set to 1 or 2\n") ;
            exit(0) ;
        }
    }
}


void dbad_terminate() {
    if( dbad_seedIndep!=NULL ) {
        free( dbad_seedIndep );
        dbad_seedIndepCnt = 0;
        dbad_curSeedIndep = 0;
    }
    if( dbad_seedDep!=NULL ) {
        free( dbad_seedDep );
        dbad_seedDepCnt = 0;
        dbad_curSeedDep = 0;
    }
}


void adContextTgt_init(double epsilon, double seed) {
    // set phase
    dbad_set_phase(TANGENT);

    // default settings
    dbad_mode = 1 ;
    dbad_ddeps = epsilon ;
    dbad_seed = seed ;

    // potentially override with settings from json file
    json_dbad_apply();

    if (dbad_phase==2) {
        printf("Tangent code,  seed=%7.1e\n", dbad_seed) ;
        printf("=============================================\n") ;
        dbad_currentSeed = 0.0 ;
    } else if (dbad_phase==1) {
        printf("Perturbed run, seed=%7.1e, epsilon=%7.1e\n", dbad_seed, dbad_ddeps) ;
        printf("=============================================\n") ;
        dbad_currentSeed = 0.0 ;
    } else if (dbad_phase==99) {
        printf("INTERNAL INTERFACE TESTS, seed=%7.1e, epsilon=%7.1e\n", dbad_seed, dbad_ddeps) ;
        printf("=============================================\n") ;
    } else {
        printf("DBAD_PHASE environment variable must be set to 1 or 2\n") ;
        exit(0) ;
    }
}

void adContextTgt_initReal8(char* varname, double *indep, double *indepd) {
  *indepd = dbad_nextIndep() ;
  if (dbad_phase==1)
    *indep = (*indep)+dbad_ddeps*(*indepd) ;
  else if (dbad_phase==99)
    printf("initReal8 of %s: %24.16e //%24.16e\n", varname, *indep, *indepd) ;
}

void adContextTgt_initReal8Array(char* varname, double *indep, double *indepd, int length) {
  int i ;
  for (i=0 ; i<length ; ++i) {
    indepd[i] = dbad_nextIndep() ;
  }
  if (dbad_phase==1) {
    for (i=0 ; i<length ; ++i) {
      indep[i] = indep[i]+dbad_ddeps*indepd[i] ;
    }
  } else if (dbad_phase==99) {
    printf("initReal8Array of %s, length=%i:\n", varname, length) ;
    for (i=0 ; i<length ; ++i)
      printf("    %i:%24.16e //%24.16e",i,indep[i],indepd[i]) ;
    printf("\n") ;
  }
}

void adContextTgt_initReal4(char* varname, float *indep, float *indepd) {
  *indepd = (float)dbad_nextIndep() ;
  if (dbad_phase==1)
    *indep = (*indep)+dbad_ddeps*(*indepd) ;
  else if (dbad_phase==99)
    printf("initReal4 of %s: %24.16e //%24.16e\n", varname, *indep, *indepd) ;
}

void adContextTgt_initReal4Array(char* varname, float *indep, float *indepd, int length) {
  int i ;
  for (i=0 ; i<length ; ++i) {
    indepd[i] = (float)dbad_nextIndep() ;
  }
  if (dbad_phase==1) {
    for (i=0 ; i<length ; ++i) {
      indep[i] = indep[i]+dbad_ddeps*indepd[i] ;
    }
  } else if (dbad_phase==99) {
    printf("initReal4Array of %s, length=%i:\n", varname, length) ;
    for (i=0 ; i<length ; ++i)
      printf("    %i:%24.16e //%24.16e",i,indep[i],indepd[i]) ;
    printf("\n") ;
  }
}

void adContextTgt_initComplex16(char* varname, double complex *indep, double complex *indepd) {
  double rdot =  dbad_nextIndep() ;
  double idot =  dbad_nextIndep() ;
  *indepd = rdot + I*idot ;
  if (dbad_phase==1) {
    *indep = *indep + dbad_ddeps*(*indepd) ;
  } else if (dbad_phase==99)
    printf("initComplex16 of %s: %24.16e+i%24.16e //%24.16e+i%24.16e\n",
           varname, creal(*indep), cimag(*indep), creal(*indepd), cimag(*indepd)) ;
}

void adContextTgt_initComplex16Array(char* varname, double complex *indep, double complex *indepd, int length) {
  double rdot, idot ;
  int i ;
  for (i=0 ; i<length ; ++i) {
    rdot =  dbad_nextIndep() ;
    idot =  dbad_nextIndep() ;
    indepd[i] = rdot + I*idot ;
  }
  if (dbad_phase==1) {
    for (i=0 ; i<length ; ++i) {
      indep[i] = indep[i] + dbad_ddeps*indepd[i] ;
    }
  } else if (dbad_phase==99) {
    printf("initComplex16Array of %s, length=%i:\n", varname, length) ;
    for (i=0 ; i<length ; ++i)
      printf("    %i:%24.16e+i%24.16e //%24.16e+i%24.16e",
             i,creal(indep[i]),cimag(indep[i]),creal(indepd[i]),cimag(indepd[i])) ;
    printf("\n") ;
  }
}

void adContextTgt_initComplex8(char* varname, ccmplx *indep, ccmplx *indepd) {
  indepd->r = (float)dbad_nextIndep() ;
  indepd->i = (float)dbad_nextIndep() ;
  if (dbad_phase==1) {
    indep->r = indep->r + dbad_ddeps*indepd->r ;
    indep->i = indep->i + dbad_ddeps*indepd->i ;
  } else if (dbad_phase==99)
    printf("initComplex8 of %s: %24.16e+i%24.16e //%24.16e+i%24.16e\n",
           varname, indep->r, indep->i, indepd->r, indepd->i) ;
}

void adContextTgt_initComplex8Array(char* varname, ccmplx *indep, ccmplx *indepd, int length) {
  int i ;
  for (i=0 ; i<length ; ++i) {
    indepd[i].r = (float)dbad_nextIndep() ;
    indepd[i].i = (float)dbad_nextIndep() ;
  }
  if (dbad_phase==1) {
    for (i=0 ; i<length ; ++i) {
      indep[i].r = indep[i].r+dbad_ddeps*indepd[i].r ;
      indep[i].i = indep[i].i+dbad_ddeps*indepd[i].i ;
    }
  } else if (dbad_phase==99) {
    printf("initComplex8Array of %s, length=%i:\n", varname, length) ;
    for (i=0 ; i<length ; ++i)
      printf("    %i:%24.16e+i%24.16e //%24.16e+i%24.16e",
             i,indep[i].r,indep[i].i,indepd[i].r,indepd[i].i) ;
    printf("\n") ;
  }
}

void adContextTgt_startConclude() {
  dbad_currentSeed= 0.0 ;
  dbad_condensed_val = 0.0 ;
  dbad_condensed_tgt = 0.0 ;
}

void adContextTgt_concludeReal8(char* varname, double dep, double depd) {
    double depb = dbad_nextDep(); //dbad_nextRandom() ;
  dbad_condensed_val += depb*(dep) ;
  if (dbad_phase==2 || dbad_phase==1)
    dbad_condensed_tgt += depb*(depd) ;
  else if (dbad_phase==99)
    printf("concludeReal8 of %s [%24.16e *] %24.16e //%24.16e\n", varname, depb, dep, depd) ;
}

void adContextTgt_concludeReal8Array(char* varname, double *dep, double *depd, int length) {
  int i ;
  double depb ;
  if (dbad_phase==99) printf("concludeReal8Array of %s, length=%i:\n", varname, length) ;
  for (i=0 ; i<length ; ++i) {
      depb = dbad_nextDep(); //dbad_nextRandom() ;
    dbad_condensed_val += depb*dep[i] ;
    if (dbad_phase==2 || dbad_phase==1) {
       dbad_condensed_tgt += depb*depd[i] ;
    } else if (dbad_phase==99) {
      printf("    %i:[%24.16e *] %24.16e //%24.16e",i,depb,dep[i],depd[i]) ;
    }
  }
  if (dbad_phase==99) printf("\n") ;
}

void adContextTgt_concludeReal4(char* varname, float dep, float depd) {
    float depb = dbad_nextDep();//(float)dbad_nextRandom() ;
  dbad_condensed_val += depb*(dep) ;
  if (dbad_phase==2 || dbad_phase==1)
    dbad_condensed_tgt += depb*(depd) ;
  else if (dbad_phase==99)
    printf("concludeReal4 of %s [%24.16e *] %24.16e //%24.16e\n", varname, depb, dep, depd) ;
}

void adContextTgt_concludeReal4Array(char* varname, float *dep, float *depd, int length) {
  int i ;
  float depb ;
  if (dbad_phase==99) printf("concludeReal4Array of %s, length=%i:\n", varname, length) ;
  for (i=0 ; i<length ; ++i) {
      depb = dbad_nextDep();//(float)dbad_nextRandom() ;
    dbad_condensed_val += depb*dep[i] ;
    if (dbad_phase==2 || dbad_phase==1) {
       dbad_condensed_tgt += depb*depd[i] ;
    } else if (dbad_phase==99) {
      printf("    %i:[%24.16e *] %24.16e //%24.16e",i,depb,dep[i],depd[i]) ;
    }
  }
  if (dbad_phase==99) printf("\n") ;
}

void adContextTgt_concludeComplex16(char* varname, double complex dep, double complex depd) {
    double depbr = 1.0;//dbad_nextRandom() ;
    double depbi = 1.0;//dbad_nextRandom() ;
  dbad_condensed_val += depbr*creal(dep) + depbi*cimag(dep);
  if (dbad_phase==2 || dbad_phase==1)
    dbad_condensed_tgt += depbr*creal(depd) + depbi*cimag(depd) ;
  else if (dbad_phase==99)
    printf("concludeComplex16 of %s [%24.16e;%24.16e *] %24.16e+i%24.16e //%24.16e+i%24.16e\n",
           varname, depbr, depbi, creal(dep), cimag(dep), creal(depd), cimag(depd)) ;
}

void adContextTgt_concludeComplex16Array(char* varname, double complex *dep, double complex *depd, int length) {
  int i ;
  double depbr, depbi ;
  if (dbad_phase==99) printf("concludeComplex16Array of %s, length=%i:\n", varname, length) ;
  for (i=0 ; i<length ; ++i) {
    depbr = dbad_nextRandom() ;
    depbi = dbad_nextRandom() ;
    dbad_condensed_val += depbr*creal(dep[i]) + depbi*cimag(dep[i]);
    if (dbad_phase==2 || dbad_phase==1) {
      dbad_condensed_tgt += depbr*creal(depd[i]) + depbi*cimag(depd[i]) ;
    } else if (dbad_phase==99) {
      printf("    %i:[%24.16e;%24.16e *] %24.16e //%24.16e",
             i, depbr, depbi, creal(dep[i]), cimag(dep[i]), creal(depd[i]), cimag(depd[i])) ;
    }
  }
  if (dbad_phase==99) printf("\n") ;
}

void adContextTgt_concludeComplex8(char* varname, ccmplx *dep, ccmplx *depd) {
  float depbr = (float)dbad_nextRandom() ;
  float depbi = (float)dbad_nextRandom() ;
  dbad_condensed_val += depbr*(dep->r) + depbi*(dep->i) ;
  if (dbad_phase==2 || dbad_phase==1)
    dbad_condensed_tgt += depbr*(depd->r) + depbi*(depd->i) ;
  else if (dbad_phase==99)
    printf("concludeComplex8 of %s [%24.16e;%24.16e *] %24.16e+i%24.16e //%24.16e+i%24.16e\n",
           varname, depbr, depbi, dep->r, dep->i, depd->r, depd->i) ;
}

void adContextTgt_concludeComplex8Array(char* varname, ccmplx *dep, ccmplx *depd, int length) {
  int i ;
  float depbr, depbi ;
  if (dbad_phase==99) printf("concludeComplex8Array of %s, length=%i:\n", varname, length) ;
  for (i=0 ; i<length ; ++i) {
    depbr = (float)dbad_nextRandom() ;
    depbi = (float)dbad_nextRandom() ;
    dbad_condensed_val += depbr*(dep[i].r) + depbi*(dep[i].i) ;
    if (dbad_phase==2 || dbad_phase==1) {
      dbad_condensed_tgt += depbr*(depd[i].r) + depbi*(depd[i].i) ;
    } else if (dbad_phase==99) {
      printf("    %i:[%24.16e;%24.16e *] %24.16e+i%24.16e //%24.16e+i%24.16e",
             i, depbr, depbi, dep[i].r, dep[i].i, depd[i].r, depd[i].i) ;
    }
  }
  if (dbad_phase==99) printf("\n") ;
}

void adContextTgt_conclude() {
  if (dbad_phase==2) {
    printf("[seed:%7.1e] Condensed result : %24.16e\n", dbad_seed, dbad_condensed_val) ;
    printf("[seed:%7.1e] Condensed tangent: %24.16e\n", dbad_seed, dbad_condensed_tgt) ;
  } else if (dbad_phase==1) {
    printf("[seed:%7.1e] Condensed perturbed result : %24.16e (epsilon:%7.1e)\n",
           dbad_seed, dbad_condensed_val, dbad_ddeps) ;
    printf("[seed:%7.1e] Condensed perturbed tangent: %24.16e\n", dbad_seed, dbad_condensed_tgt) ;
  }

  // final cleanup
  dbad_terminate();
}

void adContextAdj_init(double seed) {
    // set phase
    dbad_set_phase(ADJOINT);

    // default settings
    dbad_mode = 0 ;
    dbad_seed = seed ;

    // potentially override with settings from json file
    json_dbad_apply();

//    char* phase = getenv("DBAD_PHASE") ;
    if(dbad_phase==99) {
        printf("INTERNAL INTERFACE TESTS, seed=%7.1e\n", seed) ;
    }
    printf("Adjoint code,  seed=%7.1e\n", seed) ;
    printf("===================================\n") ;
    dbad_currentSeed = 0.0 ;
}

void adContextAdj_initReal8(char* varname, double *dep, double *depb) {
    *depb = dbad_nextDep();//dbad_nextRandom() ;
  if (dbad_phase==99)
    printf("initReal8 of %s %24.16e\n", varname, *depb) ;
}

void adContextAdj_initReal8Array(char* varname, double *dep, double *depb, int length) {
  int i ;
  for (i=0 ; i<length ; ++i) {
      depb[i] = dbad_nextDep();//dbad_nextRandom() ;
  }
  if (dbad_phase==99) {
    printf("initReal8Array of %s, length=%i\n", varname, length) ;
    for (i=0 ; i<length ; ++i)
      printf("    %i:%24.16e", i, depb[i]) ;
    printf("\n") ;
  }
}

void adContextAdj_initReal4(char* varname, float *dep, float *depb) {
    *depb = (float)dbad_nextDep();//dbad_nextDep();//dbad_nextRandom() ;
  if (dbad_phase==99)
    printf("initReal4 of %s %24.16e\n", varname, *depb) ;
}

void adContextAdj_initReal4Array(char* varname, float *dep, float *depb, int length) {
  int i ;
  for (i=0 ; i<length ; ++i) {
    depb[i] = (float)dbad_nextDep();//dbad_nextRandom() ;
  }
  if (dbad_phase==99) {
    printf("initReal4Array of %s, length=%i\n", varname, length) ;
    for (i=0 ; i<length ; ++i)
      printf("    %i:%24.16e",i, depb[i]) ;
    printf("\n") ;
  }
}

void adContextAdj_initComplex16(char* varname, double complex *dep, double complex *depb) {
  double rbar =  dbad_nextDep();//dbad_nextRandom() ;
  double ibar =  dbad_nextDep();//dbad_nextRandom() ;
  *depb = rbar + I*ibar ;
  if (dbad_phase==99)
    printf("initComplex16 of %s %24.16e+i%24.16e\n", varname, creal(*depb), cimag(*depb)) ;
}

void adContextAdj_initComplex16Array(char* varname, double complex *dep, double complex *depb, int length) {
  double rbar, ibar ;
  int i ;
  for (i=0 ; i<length ; ++i) {
    rbar = dbad_nextDep();//dbad_nextRandom() ;
    ibar = dbad_nextDep();//dbad_nextRandom() ;
    depb[i] = rbar + I*ibar ;
  }
  if (dbad_phase==99) {
    printf("initComplex16Array of %s, length=%i\n", varname, length) ;
    for (i=0 ; i<length ; ++i)
      printf("    %i:%24.16e+i%24.16e",i, creal(depb[i]), cimag(depb[i])) ;
    printf("\n") ;
  }
}

void adContextAdj_initComplex8(char* varname, ccmplx *dep, ccmplx *depb) {
  depb->r = (float)dbad_nextDep();//dbad_nextRandom() ;
  depb->i = (float)dbad_nextDep();//dbad_nextRandom() ;
  if (dbad_phase==99)
    printf("initComplex8 of %s %24.16e+i%24.16e\n", varname, depb->r, depb->i) ;
}

void adContextAdj_initComplex8Array(char* varname, ccmplx *dep, ccmplx *depb, int length) {
  int i ;
  for (i=0 ; i<length ; ++i) {
    depb[i].r = (float)dbad_nextDep();//dbad_nextRandom() ;
    depb[i].i = (float)dbad_nextDep();//dbad_nextRandom() ;
  }
  if (dbad_phase==99) {
    printf("initComplex8Array of %s, length=%i\n", varname, length) ;
    for (i=0 ; i<length ; ++i)
      printf("    %i:%24.16e+i%24.16e", i, depb[i].r, depb[i].i) ;
    printf("\n") ;
  }
}

void adContextAdj_startConclude() {
  dbad_currentSeed= 0.0 ;
  dbad_condensed_adj = 0.0 ;
}

void adContextAdj_concludeReal8(char* varname, double dep, double depb) {
    double depd = dbad_nextIndep();//dbad_nextRandom() ;
  dbad_condensed_adj += depd*depb ;
  if (dbad_phase==99)
    printf("concludeReal8 of %s [%24.16e *]%24.16e\n", varname, depd, depb) ;
}

void adContextAdj_concludeReal8Array(char* varname, double *dep, double *depb, int length) {
  int i ;
  double depd ;
  if (dbad_phase==99) printf("concludeReal8Array of %s, length=%i:\n", varname, length) ;
  for (i=0 ; i<length ; ++i) {
    depd = dbad_nextIndep();//dbad_nextRandom() ;
    dbad_condensed_adj += depd*depb[i] ;
    if (dbad_phase==99) printf("    %i:[%24.16e *] %24.16e",i,depd,depb[i]) ;
  }
  if (dbad_phase==99) printf("\n") ;
}

void adContextAdj_concludeReal4(char* varname, float dep, float depb) {
  float depd = (float)dbad_nextIndep();//dbad_nextRandom() ;
  dbad_condensed_adj += depd*depb ;
  if (dbad_phase==99)
    printf("concludeReal4 of %s [%24.16e *]%24.16e\n", varname, depd, depb) ;
}

void adContextAdj_concludeReal4Array(char* varname, float *dep, float *depb, int length) {
  int i ;
  float depd ;
  if (dbad_phase==99) printf("concludeReal4Array of %s, length=%i:\n", varname, length) ;
  for (i=0 ; i<length ; ++i) {
    depd = (float)dbad_nextIndep();//dbad_nextRandom() ;
    dbad_condensed_adj += depd*depb[i] ;
    if (dbad_phase==99) printf("    %i:[%24.16e *] %24.16e",i,depd,depb[i]) ;
  }
  if (dbad_phase==99) printf("\n") ;
}

void adContextAdj_concludeComplex16(char* varname, double complex dep, double complex depb) {
  double depdr = dbad_nextIndep();//dbad_nextRandom() ;
  double depdi = dbad_nextIndep();//dbad_nextRandom() ;
  dbad_condensed_adj += depdr*creal(depb) + depdi*cimag(depb) ;
  if (dbad_phase==99)
    printf("concludeComplex16 of %s [%24.16e+i%24.16e *]%24.16e+i%24.16e\n", varname, depdr, depdi, creal(depb), cimag(depb)) ;
}

void adContextAdj_concludeComplex16Array(char* varname, double complex *dep, double complex *depb, int length) {
  int i ;
  double depdr, depdi ;
  if (dbad_phase==99) printf("concludeComplex16Array of %s, length=%i:\n", varname, length) ;
  for (i=0 ; i<length ; ++i) {
    depdr = dbad_nextIndep();//dbad_nextRandom() ;
    depdi =dbad_nextIndep();// dbad_nextRandom() ;
    dbad_condensed_adj += depdr*creal(depb[i]) + depdi*cimag(depb[i]) ;
    if (dbad_phase==99) printf("    %i:[%24.16e+i%24.16e *] %24.16e+i%24.16e",i,depdr,depdi,creal(depb[i]),cimag(depb[i])) ;
  }
  if (dbad_phase==99) printf("\n") ;
}

void adContextAdj_concludeComplex8(char* varname, ccmplx *dep, ccmplx *depb) {
  float depdr = (float)dbad_nextIndep();//dbad_nextRandom() ;
  float depdi = (float)dbad_nextIndep();//dbad_nextRandom() ;
  dbad_condensed_adj += depdr*depb->r + depdi*depb->i ;
  if (dbad_phase==99)
    printf("concludeComplex8 of %s [%24.16e+i%24.16e *]%24.16e+i%24.16e\n", varname, depdr, depdi, depb->r, depb->i) ;
}

void adContextAdj_concludeComplex8Array(char* varname, ccmplx *dep, ccmplx *depb, int length) {
  int i ;
  float depdr, depdi ;
  if (dbad_phase==99) printf("concludeComplex8Array of %s, length=%i:\n", varname, length) ;
  for (i=0 ; i<length ; ++i) {
    depdr = (float)dbad_nextIndep();//dbad_nextRandom() ;
    depdi = (float)dbad_nextIndep();//dbad_nextRandom() ;
    dbad_condensed_adj += depdr*depb[i].r + depdi*depb[i].i ;
    if (dbad_phase==99) printf("    %i:[%24.16e+i%24.16e *] %24.16e+i%24.16e",i,depdr,depdi,depb[i].r,depb[i].i) ;
  }
  if (dbad_phase==99) printf("\n") ;
}

void adContextAdj_conclude() {
  printf("[seed:%7.1e] Condensed adjoint: %24.16e\n", dbad_seed, dbad_condensed_adj) ;
  // final cleanup
  dbad_terminate();
}

//############## INTERFACE PROCEDURES CALLED FROM FORTRAN ################

void adcontexttgt_init_(double *epsilon, double *seed) {
  adContextTgt_init(*epsilon, *seed) ;
}

void adcontexttgt_initreal8_(char* varname, double *indep, double *indepd) {
  adContextTgt_initReal8(varname, indep, indepd) ;
}

void adcontexttgt_initreal8array_(char* varname, double *indep, double *indepd, int *length) {
  adContextTgt_initReal8Array(varname, indep, indepd, *length) ;
}

void adcontexttgt_initreal4_(char* varname, float *indep, float *indepd) {
  adContextTgt_initReal4(varname, indep, indepd) ;
}

void adcontexttgt_initreal4array_(char* varname, float *indep, float *indepd, int *length) {
  adContextTgt_initReal4Array(varname, indep, indepd, *length) ;
}

void adcontexttgt_initcomplex16_(char* varname, cdcmplx *indep, cdcmplx *indepd) {
  adContextTgt_initComplex16(varname, (double complex *)indep, (double complex *)indepd) ;
}

void adcontexttgt_initcomplex16array_(char* varname, cdcmplx *indep, cdcmplx *indepd, int *length) {
  adContextTgt_initComplex16Array(varname, (double complex *)indep, (double complex *)indepd, *length) ;
}

void adcontexttgt_initcomplex8_(char* varname, ccmplx *indep, ccmplx *indepd) {
  adContextTgt_initComplex8(varname, indep, indepd) ;
}

void adcontexttgt_initcomplex8array_(char* varname, ccmplx *indep, ccmplx *indepd, int *length) {
  adContextTgt_initComplex8Array(varname, indep, indepd, *length) ;
}

void adcontexttgt_startconclude_() {
  adContextTgt_startConclude() ;
}

void adcontexttgt_concludereal8_(char* varname, double *dep, double *depd) {
  if (dbad_phase==99)
      printf("concludereal8_ of %s: \n", varname);
  adContextTgt_concludeReal8(varname, *dep, *depd) ;
}

void adcontexttgt_concludereal8array_(char* varname, double *dep, double *depd, int *length) {
  if (dbad_phase==99)
      printf("concludereal8array_ of %s: \n", varname);
  adContextTgt_concludeReal8Array(varname, dep, depd, *length) ;
}

void adcontexttgt_concludereal4_(char* varname, float *dep, float *depd) {
  adContextTgt_concludeReal4(varname, *dep, *depd) ;
}

void adcontexttgt_concludereal4array_(char* varname, float *dep, float *depd, int *length) {
  adContextTgt_concludeReal4Array(varname, dep, depd, *length) ;
}

void adcontexttgt_concludecomplex16_(char* varname, cdcmplx *dep, cdcmplx *depd) {
  adContextTgt_concludeComplex16(varname, *((double complex *)dep), *((double complex *)depd)) ;
}

void adcontexttgt_concludecomplex16array_(char* varname, cdcmplx *dep, cdcmplx *depd, int *length) {
  adContextTgt_concludeComplex16Array(varname, (double complex *)dep, (double complex *)depd, *length) ;
}

void adcontexttgt_concludecomplex8_(char* varname, ccmplx *dep, ccmplx *depd) {
  if (dbad_phase==99)
      printf("concludecomplex8_ of %s: \n", varname);
  adContextTgt_concludeComplex8(varname, dep, depd) ;
}

void adcontexttgt_concludecomplex8array_(char* varname, ccmplx *dep, ccmplx *depd, int *length) {
  if (dbad_phase==99)
      printf("concludecomplex8array_ of %s: \n", varname);
  adContextTgt_concludeComplex8Array(varname, dep, depd, *length) ;
}

void adcontexttgt_conclude_() {
  adContextTgt_conclude() ;
}

void adcontextadj_init_(double *seed) {
  adContextAdj_init(*seed) ;
}

void adcontextadj_initreal8_(char* varname, double *dep, double *depb) {
  if (dbad_phase==99)
    printf("initreal8_ of %s \n", varname) ;
  adContextAdj_initReal8(varname, dep, depb) ;
}

void adcontextadj_initreal8array_(char* varname, double *dep, double *depb, int *length) {
  if (dbad_phase==99)
    printf("initreal8array_ of %s \n", varname) ;
  adContextAdj_initReal8Array(varname, dep, depb, *length) ;
}

void adcontextadj_initreal4_(char* varname, float *dep, float *depb) {
  adContextAdj_initReal4(varname, dep, depb) ;
}

void adcontextadj_initreal4array_(char* varname, float *dep, float *depb, int *length) {
  adContextAdj_initReal4Array(varname, dep, depb, *length) ;
}

void adcontextadj_initcomplex16_(char* varname, cdcmplx *dep, cdcmplx *depb) {
  adContextAdj_initComplex16(varname, (double complex *)dep, (double complex *)depb) ;
}

void adcontextadj_initcomplex16array_(char* varname, cdcmplx *dep, cdcmplx *depb, int *length) {
  adContextAdj_initComplex16Array(varname, (double complex *)dep, (double complex *)depb, *length) ;
}

void adcontextadj_initcomplex8_(char* varname, ccmplx *dep, ccmplx *depb) {
  adContextAdj_initComplex8(varname, dep, depb) ;
}

void adcontextadj_initcomplex8array_(char* varname, ccmplx *dep, ccmplx *depb, int *length) {
  adContextAdj_initComplex8Array(varname, dep, depb, *length) ;
}

void adcontextadj_startconclude_() {
  adContextAdj_startConclude() ;
}

void adcontextadj_concludereal8_(char* varname, double *dep, double *depb) {
  if (dbad_phase==99)
    printf("concludereal8_ of %s \n", varname) ;
  adContextAdj_concludeReal8(varname, *dep, *depb) ;
}

void adcontextadj_concludereal8array_(char* varname, double *dep, double *depb, int *length) {
  if (dbad_phase==99)
    printf("concludereal8array_ of %s \n", varname) ;
  adContextAdj_concludeReal8Array(varname, dep, depb, *length) ;
}

void adcontextadj_concludereal4_(char* varname, float *dep, float *depb) {
  if (dbad_phase==99)
    printf("concludereal4_ of %s \n", varname) ;
  adContextAdj_concludeReal4(varname, *dep, *depb) ;
}

void adcontextadj_concludereal4array_(char* varname, float *dep, float *depb, int *length) {
  if (dbad_phase==99)
    printf("concludereal4array_ of %s \n", varname) ;
  adContextAdj_concludeReal4Array(varname, dep, depb, *length) ;
}

void adcontextadj_concludecomplex16_(char* varname, cdcmplx *dep, cdcmplx *depb) {
  adContextAdj_concludeComplex16(varname, *((double complex *)dep), *((double complex *)depb)) ;
}

void adcontextadj_concludecomplex16array_(char* varname, cdcmplx *dep, cdcmplx *depb, int *length) {
  adContextAdj_concludeComplex16Array(varname, (double complex *)dep, (double complex *)depb, *length) ;
}

void adcontextadj_concludecomplex8_(char* varname, ccmplx *dep, ccmplx *depb) {
  adContextAdj_concludeComplex8(varname, dep, depb) ;
}

void adcontextadj_concludecomplex8array_(char* varname, ccmplx *dep, ccmplx *depb, int *length) {
  adContextAdj_concludeComplex8Array(varname, dep, depb, *length) ;
}

void adcontextadj_conclude_() {
  adContextAdj_conclude() ;
}
