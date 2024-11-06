/**
 * \file Needleman-Wunsch-recmemo.c
 * \brief recursive implementation with memoization of Needleman-Wunsch global alignment algorithm that computes the distance between two genetic sequences 
 * \version 0.1
 * \date 03/10/2022 
 * \author Jean-Louis Roch (Ensimag, Grenoble-INP - University Grenoble-Alpes) jean-louis.roch@grenoble-inp.fr
 *
 * Documentation: see Needleman-Wunsch-recmemo.h
 * Costs of basic base opertaions (SUBSTITUTION_COST, SUBSTITUTION_UNKNOWN_COST, INSERTION_COST) are
 * defined in Needleman-Wunsch-recmemo.h
 */


#include "Needleman-Wunsch-recmemo.h"
#include <stdio.h>  
#include <stdlib.h> 
#include <string.h> /* for strchr */
// #include <ctype.h> /* for toupper */

#include "characters_to_base.h" /* mapping from char to base */

/*****************************************************************************/
   
/* Context of the memoization : passed to all recursive calls */
/** \def NOT_YET_COMPUTED
 * \brief default value for memoization of minimal distance (defined as an impossible value for a distance, -1).
 */
#define NOT_YET_COMPUTED -1L 

/** \struct NW_MemoContext
 * \brief data for memoization of recursive Needleman-Wunsch algorithm 
*/
struct NW_MemoContext 
{
    char *X ; /*!< the longest genetic sequences */
    char *Y ; /*!< the shortest genetic sequences */
    size_t M; /*!< length of X */
    size_t N; /*!< length of Y,  N <= M */
    long **memo; /*!< memoization table to store memo[0..M][0..N] (including stopping conditions phi(M,j) and phi(i,N) */
} ;

/*
 *  static long EditDistance_NW_RecMemo(struct NW_MemoContext *c, size_t i, size_t j) 
 * \brief  EditDistance_NW_RecMemo :  Private (static)  recursive function with memoization \
 * direct implementation of Needleman-Wursch extended to manage FASTA sequences (cf TP description)
 * \param c : data passed for recursive calls that includes the memoization array 
 * \param i : starting position of the left sequence :  c->X[ i .. c->M ] 
 * \param j : starting position of the right sequence :  c->Y[ j .. c->N ] 
 */ 
static long EditDistance_NW_RecMemo(struct NW_MemoContext *c, size_t i, size_t j) 
/* compute and returns phi(i,j) using data in c -allocated and initialized by EditDistance_NW_Rec */
{
   if (c->memo[i][j] == NOT_YET_COMPUTED)
   {  
      long res ;
      char Xi = c->X[i] ;
      char Yj = c->Y[j] ;
      if (i == c->M) /* Reach end of X */
      {  if (j == c->N) res = 0;  /* Reach end of Y too */
         else res = (isBase(Yj) ? INSERTION_COST : 0) + EditDistance_NW_RecMemo(c, i, j+1) ;
      }
      else if (j == c->N) /* Reach end of Y but not end of X */
      {  res = (isBase(Xi) ? INSERTION_COST : 0) + EditDistance_NW_RecMemo(c, i+1, j) ;
      }
      else if (! isBase(Xi))  /* skip ccharacter in Xi that is not a base */
      {  ManageBaseError( Xi ) ;
         res = EditDistance_NW_RecMemo(c, i+1, j) ;
      }
      else if (! isBase(Yj))  /* skip ccharacter in Yj that is not a base */
      {  ManageBaseError( Yj ) ;
         res = EditDistance_NW_RecMemo(c, i, j+1) ;
      }
      else  
      {  /* Note that stopping conditions (i==M) and (j==N) are already stored in c->memo (cf EditDistance_NW_Rec) */ 
         long min = /* initialization  with cas 1*/
                   ( isUnknownBase(Xi) ?  SUBSTITUTION_UNKNOWN_COST 
                          : ( isSameBase(Xi, Yj) ? 0 : SUBSTITUTION_COST ) 
                   )
                   + EditDistance_NW_RecMemo(c, i+1, j+1) ; 
         { long cas2 = INSERTION_COST + EditDistance_NW_RecMemo(c, i+1, j) ;      
           if (cas2 < min) min = cas2 ;
         }
         { long cas3 = INSERTION_COST + EditDistance_NW_RecMemo(c, i, j+1) ;      
           if (cas3 < min) min = cas3 ; 
         }
         res = min ;
      }
       c->memo[i][j] = res ;
   }
   return c->memo[i][j] ;
}

/* EditDistance_NW_Rec :  is the main function to call, cf .h for specification 
 * It allocates and initailizes data (NW_MemoContext) for memoization and call the 
 * recursivefunction EditDistance_NW_RecMemo 
 * See .h file for documentation
 */
long EditDistance_NW_Rec(char* A, size_t lengthA, char* B, size_t lengthB)
{
   _init_base_match() ;
   struct NW_MemoContext ctx;
   if (lengthA >= lengthB) /* X is the longest sequence, Y the shortest */
   {  ctx.X = A ;
      ctx.M = lengthA ;
      ctx.Y = B ;
      ctx.N = lengthB ;
   }
   else
   {  ctx.X = B ;
      ctx.M = lengthB ;
      ctx.Y = A ;
      ctx.N = lengthA ;
   }
   size_t M = ctx.M ;
   size_t N = ctx.N ;
   {  /* Allocation and initialization of ctx.memo to NOT_YET_COMPUTED*/
      /* Note: memo is of size (N+1)*(M+1) but is stored as (M+1) distinct arrays each with (N+1) continuous elements 
       * It would have been possible to allocate only one big array memezone of (M+1)*(N+1) elements 
       * and then memo as an array of (M+1) pointers, the memo[i] being the address of memzone[i*(N+1)].
       */ 
      ctx.memo = (long **) malloc ( (M+1) * sizeof(long *)) ;
      if (ctx.memo == NULL) { perror("EditDistance_NW_Rec: malloc of ctx_memo" ); exit(EXIT_FAILURE); }
      for (int i=0; i <= M; ++i) 
      {  ctx.memo[i] = (long*) malloc( (N+1) * sizeof(long));
         if (ctx.memo[i] == NULL) { perror("EditDistance_NW_Rec: malloc of ctx_memo[i]" ); exit(EXIT_FAILURE); }
         for (int j=0; j<=N; ++j) ctx.memo[i][j] = NOT_YET_COMPUTED ;
      }
   }    
   
   /* Compute phi(0,0) = ctx.memo[0][0] by calling the recursive function EditDistance_NW_RecMemo */
   long res = EditDistance_NW_RecMemo( &ctx, 0, 0 ) ;
    
   { /* Deallocation of ctx.memo */
      for (int i=0; i <= M; ++i) free( ctx.memo[i] ) ;
      free( ctx.memo ) ;
   }
   return res ;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
                                             
                                             /* ITERATIF */
                                             
///////////////////////////////////////////////////////////////////////////////////////////////////////

struct NW_MemoContext_Iter
{
    char *X ; /*!< the longest genetic sequences */
    char *Y ; /*!< the shortest genetic sequences */
    size_t M; /*!< length of X */
    size_t N; /*!< length of Y,  N <= M */
    long *memo; /*!< memoization table to store memo[0..M] */
};

long EditDistance_NW_Iter(char* A, size_t lengthA, char* B, size_t lengthB)
{
   long res;

   _init_base_match();
   struct NW_MemoContext_Iter ctx;
   if (lengthA >= lengthB) /* X is the longest sequence, Y the shortest */
   {  ctx.X = A ;
      ctx.M = lengthA ;
      ctx.Y = B ;
      ctx.N = lengthB ;
   }
   else
   {  ctx.X = B ;
      ctx.M = lengthB ;
      ctx.Y = A ;
      ctx.N = lengthA ;
   }
   size_t M = ctx.M ;
   size_t N = ctx.N ;

   ctx.memo = (long*) malloc((N+1)*sizeof(long)) ;

   if (ctx.memo == NULL) { perror("EditDistance_NW_Rec: malloc of ctx_memo" ); exit(EXIT_FAILURE); }
   
   ctx.memo[N] = 0;
   
   for(int j=N-1; j >= 0; j-=1){
      ctx.memo[j] = (isBase(ctx.Y[j]) ? INSERTION_COST : 0)+ ctx.memo[j+1];
   }
   
   for(int i=M-1; i>=0; i-=1){
      long temp = ctx.memo[N];
      ctx.memo[N] = (isBase(ctx.X[i]) ? INSERTION_COST : 0) + temp;
      
      for (int j=N-1;j >= 0; j-=1){
         if (isBase(ctx.X[i]) == 0) {
            ManageBaseError(ctx.X[i]);
            temp = ctx.memo[j]; 
         } else if (isBase(ctx.Y[j]) == 0) {
            ManageBaseError(ctx.Y[j]);
            temp = ctx.memo[j]; 
            ctx.memo[j] = ctx.memo[j+1];
         } else {
            long sigma = 0;
            if isUnknownBase(ctx.X[i]){
               sigma = SUBSTITUTION_UNKNOWN_COST;
            } else if (isSameBase(ctx.X[i],ctx.Y[j]) == 0){
               sigma = SUBSTITUTION_COST;
            }
            long case_1 = sigma + temp;
            long case_2 = INSERTION_COST + ctx.memo[j];
            long case_3 = INSERTION_COST + ctx.memo[j+1];
            
            long min = case_1;
            if (case_2 < min){
               min = case_2;
            }
            if (case_3 < min){
               min = case_3;
            }
            temp = ctx.memo[j];
            ctx.memo[j]= min; 
         }
      }
   } 

   res = ctx.memo[0];
   free(ctx.memo);

   return res;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
                                             
                                             /* CACHE AWARE */
                                             
///////////////////////////////////////////////////////////////////////////////////////////////////////

struct NW_MemoContext_CA
{
   char *X ; /*!< the longest genetic sequences */
   char *Y ; /*!< the shortest genetic sequences */
   size_t M; /*!< length of X */
   size_t N; /*!< length of Y,  N <= M */
   long *memo_H; /*!< memoization table to store memo[0..N], valeurs des lignes matrice */
   long *memo_v; /*!< memoization table to store memo[0..K], valeurs colonnes blocs */
};

#define LINE_SIZE 64

long EditDistance_NW_CA(char* A, size_t lengthA, char* B, size_t lengthB, int Z)
{
   long res;

   _init_base_match();
   struct NW_MemoContext_CA ctx;
   if (lengthA >= lengthB) /* X is the longest sequence, Y the shortest */
   {  ctx.X = A ;
      ctx.M = lengthA ;
      ctx.Y = B ;
      ctx.N = lengthB ;
   }
   else
   {  ctx.X = B ;
      ctx.M = lengthB ;
      ctx.Y = A ;
      ctx.N = lengthA ;
   }
   size_t M = ctx.M ;
   size_t N = ctx.N ;

   /*Taille des blocs utilisés pour le cache aware*/
   int K = Z / LINE_SIZE;

   /*Allocation des mémoires*/
   ctx.memo_H = (long*) malloc((N+1)*sizeof(long)) ;
   if (ctx.memo_H == NULL) { perror("EditDistance_NW_Rec: malloc of ctx_memo" ); exit(EXIT_FAILURE); }

   ctx.memo_v = (long*) malloc((K)*sizeof(long)) ;
   if (ctx.memo_v == NULL) { perror("EditDistance_NW_Rec: malloc of ctx_memo" ); exit(EXIT_FAILURE); }

   long temp_bloc = 0;
   long temp_v = 0;

   for (int I=M; I>=0; I -= K){
      int i_end = (I-K+1 < 0) ? 0 : I-K+1;

      for (int J=N; J>=0; J-= K){
         int j_end = (J-K+1 < 0) ? 0 : J-K+1;
         long temp_diag_bloc = (I == M) ? NOT_YET_COMPUTED : ctx.memo_H[j_end];

         for (int i=I; i>=i_end; i-=1){
            int i_bloc = i - i_end;

            for (int j=J; j>=j_end; j-=1){
               if (i == M) { /* Cas de la lignes des blocs tout en bas */
                  if (j == N){
                     ctx.memo_H[j] = 0;
                  } else {
                     if (isBase(ctx.Y[j])){ //(j == N) ? ctx.memo_H[j] = 0 : ctx.memo_H[j] = (isBase(ctx.Y[j]) ? INSERTION_COST : 0) + (j == J ? ctx.memo_v[i_bloc] : ctx.memo_H[j+1]);
                        ctx.memo_H[j] = INSERTION_COST;
                     } else {
                        ctx.memo_H[j] = 0;
                     }

                     if (j == J){
                        ctx.memo_H[j]+= ctx.memo_v[i_bloc];
                     } else {
                        ctx.memo_H[j]+= ctx.memo_H[j+1];
                     }
                  }
               } else if (j == N) { /* Cas de la colonne des blocs de droite */
                  temp_bloc = ctx.memo_H[N];
                  ctx.memo_H[N] = (isBase(ctx.X[i]) ? INSERTION_COST : 0) + temp_bloc;
               } else if ( isBase(ctx.X[i]) == 0){
                  ManageBaseError(ctx.X[i]);
                  temp_bloc = ctx.memo_H[j];
               } else if (isBase(ctx.Y[j]) == 0){
                  ManageBaseError(ctx.Y[j]);
                  temp_bloc = ctx.memo_H[j];
                  ctx.memo_H[j] = (j == J ? ctx.memo_v[i_bloc] : ctx.memo_H[j+1]);
               } else {
                  long sigma = 0;
                  if (isUnknownBase(ctx.X[i])){
                     sigma = SUBSTITUTION_UNKNOWN_COST;
                  } else if (isSameBase(ctx.X[i],ctx.Y[j]) == 0){
                     sigma = SUBSTITUTION_COST;
                  }
                  long case_1 = sigma + (j == J ? temp_v : temp_bloc);
                  long case_2 = INSERTION_COST + ctx.memo_H[j];
                  long case_3 = INSERTION_COST + (j == J ? ctx.memo_v[i_bloc] : ctx.memo_H[j+1]);

                  long min = case_1;
                  if (case_2 < min){
                     min = case_2;
                  }
                  if (case_3 < min){
                     min = case_3;
                  }
                  temp_bloc = ctx.memo_H[j];
                  ctx.memo_H[j]= min; 
               }
            }

            temp_v = ctx.memo_v[i_bloc];
            ctx.memo_v[i_bloc] = ctx.memo_H[j_end];
         } 
         temp_v = temp_diag_bloc;
      }   
   }

   res = ctx.memo_H[0];

   /*Libération des mémoires*/
   free(ctx.memo_H);
   free(ctx.memo_v);

   return res;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
                                             
                                             /* CACHE OBLIVIOUS */
                                             
//////////////////////////////////////////////////////////////////////////////////////////////////////

struct NW_MemoContext_C0
{
   char *X ; /*!< the longest genetic sequences */
   char *Y ; /*!< the shortest genetic sequences */
   size_t M; /*!< length of X */
   size_t N; /*!< length of Y,  N <= M */
   long *memo_H; /*!< memoization table to store memo[0..N] */
   long *memo_v; /*!< memoization table to store memo[0..K]*/
};

#define seuil 200

void blockingRec(struct NW_MemoContext_CO *c, int i_begin, int i_end){
   if (i_begin-i_end > seuil){
      int i_middle = (i_begin+i_end)/2;
      blockingRec(c,i_begin,i_middle+1);
      blockingRec(c,i_middle,i_end);
   } else {
      /* calcul itératif d'un bloc */
      for (int i=i_begin; i>=i_end; i -=1){ 
         int i_bloc = i - i_end;

         /*calcul de la sous-colonne de droite de la colonne*/

      }

      for (int j=c->N; j>=0; j-=1){
         for (int i=i_begin; i>=i_end; k-=1){
            int i_bloc = i - i_end;

            /*calcul des autres sous-colonnes de la colonne*/
            
         }

      }
      
   }
} 


long EditDistance_NW_CO(char* A, size_t lengthA, char* B, size_t lengthB)
{
   long res;

   _init_base_match();
   struct NW_MemoContext_C0 ctx;
   if (lengthA >= lengthB) /* X is the longest sequence, Y the shortest */
   {  ctx.X = A ;
      ctx.M = lengthA ;
      ctx.Y = B ;
      ctx.N = lengthB ;
   }
   else
   {  ctx.X = B ;
      ctx.M = lengthB ;
      ctx.Y = A ;
      ctx.N = lengthA ;
   }
   size_t M = ctx.M ;
   size_t N = ctx.N ;

   /*Allocation des mémoires*/
   ctx.memo_H = (long*) malloc((N+1)*sizeof(long)) ;
   if (ctx.memo_H == NULL) { perror("EditDistance_NW_Rec: malloc of ctx_memo" ); exit(EXIT_FAILURE); }

   ctx.memo_v = (long*) malloc((S)*sizeof(long)) ;
   if (ctx.memo_v == NULL) { perror("EditDistance_NW_Rec: malloc of ctx_memo" ); exit(EXIT_FAILURE); }

   blockingRec(ctx,M,0);

   res = ctx.memo_H[0];

   /*Libération des mémoires*/
   free(ctx.memo_H);
   free(ctx.memo_v);

   return res;
}