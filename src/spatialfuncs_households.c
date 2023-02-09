#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Rdynload.h>
#include <Rdefines.h>




/*****************************************************************/
/*pi function optimized for iterating through a list of two types*/
/*finding ones that fulfill the distance requirement.            */
/*                                                               */
/*@param hh      the household column                            */
/*@param type    the type column                                 */
/*@param x       the x coordinate                                */
/*@param y       the y coordinate                                */
/*@param len     the length of the three data arrays             */
/*@param typeA   the "from" type                                 */
/*@param typeB   the "to" type                                   */
/*@param r_low   the low end of the range of values to look at   */
/*@param r       the sequence of upper distances to consider     */
/*@param len_r   the number of different Rs to consider          */
/*                   usually will be just the indicies           */
/*@param inds    the indices into the original array, helps boot */
/*@param rc      the array of values that we are going to return */
/*****************************************************************/
void get_pi_hh_typed (int *hh,
           int *type,
		   double *x,
		   double *y,
		   int *len,
		   int *typeA,
		   int *typeB,
		   double *r_low,
		   double *r,
		   int *len_r,
		   int *inds,
		   double *rc) {
  int i,j,k;
  long long num_cnt, denom_cnt; /*counters for those filling conditions*/
  double dist;

  /*repeat the calculation for all r*/
  for (i=0;i<*len_r;i++) {
    //Rprintf("%1.1f,%1.1f\n", r_low[i],r[i]); //DEBUG
    /*zero out counts*/
    num_cnt = 0;
    denom_cnt = 0;

    if (*typeA != -1) {

      for (j=0;j<*len;j++) {
	    if (type[j] != *typeA) continue;

    	for (k=0;k<*len;k++) {
    	  /*ignore pairs of the same observation*/
    	  if (inds[k]==inds[j]) continue;
    	  /*ignore pairs from the same household*/
    	  if (hh[k]==hh[j]) continue;
    	  
    	  dist = sqrt(pow(x[j]-x[k],2)+pow(y[j]-y[k],2));
    	  if ((dist<=r[i])  & (dist>=r_low[i])) denom_cnt++;
    
    	  if (type[k] != *typeB) continue;
    	  if ((dist<=r[i])  & (dist>=r_low[i])) num_cnt++;
    	}
      }

    } else {
      Rprintf("To be implemented\n");
      return;
    }
    //Rprintf("%d/%d\n",num_cnt,denom_cnt);//DEBUG
    rc[i] = (double)num_cnt/denom_cnt;
  }
}






/********************************************************************/
/*theta function optimized for iterating through a list of two types*/
/*finding ones that fulfill the distance requirement.               */
/*                                                                  */
/*@param hh      the household column                               */
/*@param type    the type column                                    */
/*@param x       the x coordinate                                   */
/*@param y       the y coordinate                                   */
/*@param len     the length of the three data arrays                */
/*@param typeA   the "from" type                                    */
/*@param typeB   the "to" type                                      */
/*@param r_low   the low end of the range of values to look at      */
/*@param r       the sequence of upper distances to consider        */
/*@param len_r   the number of different Rs to consider             */
/*                   usually will be just the indicies              */
/*@param inds    the indices into the original array, helps boot    */
/*@param rc      the array of values that we are going to return    */
/********************************************************************/
void get_theta_hh_typed (int *hh,
                      int *type,
                      double *x,
                      double *y,
                      int *len,
                      int *typeA,
                      int *typeB,
                      double *r_low,
                      double *r,
                      int *len_r,
                      int *inds,
                      double *rc) {
    
    int i,j,k;
    long long num_cnt, denom_cnt; /*counters for those filling conditions*/
double dist;

/*repeat the calculation for all r*/
for (i=0;i<*len_r;i++) {
    //Rprintf("%1.1f,%1.1f\n", r_low[i],r[i]); //DEBUG
    /*zero out counts*/
    num_cnt = 0;
    denom_cnt = 0;
    
    if (*typeA != -1) {
        
        for (j=0;j<*len;j++) {
            if (type[j] != *typeA) continue;
            
            for (k=0;k<*len;k++) {
                /*ignore pairs of the same observation*/
                if (inds[k]==inds[j]) continue;
                /*ignore pairs from the same household*/
                if (hh[k]==hh[j]) continue;
                
                dist = sqrt(pow(x[j]-x[k],2)+pow(y[j]-y[k],2));
                
                if ((dist<=r[i])  & (dist>=r_low[i])) {
                    if (type[k] == *typeB) {
                        num_cnt++;
                    } else {
                        denom_cnt++;
                    }
                }
            }
        }
        
    } else {
        Rprintf("To be implemented\n");
        return;
    }
    //Rprintf("%d/%d\n",num_cnt,denom_cnt);//DEBUG
    rc[i] = (double)num_cnt/denom_cnt;
}
}



/*****************************************************************/
/*pi function optimized for iterating through a list of two types*/
/*finding ones that fulfill the distance requirement.            */
/*                                                               */
/*@param type_1  the  first type column                          */
/*@param type_2  the second type column                          */
/*@param x       the x coordinate                                */
/*@param y       the y coordinate                                */
/*@param weight  the sample weight of the individual/cluster     */
/*@param len     the length of the three data arrays             */
/*@param typeA   the "from" type                                 */
/*@param typeB   the "to" type                                   */
/*@param r_low   the low end of the range of values to look at   */
/*@param r       the sequence of upper distances to consider     */
/*@param len_r   the number of different Rs to consider          */
/*                   usually will be just the indicies           */
/*@param inds    the indices into the original array, helps boot */
/*@param rc      the array of values that we are going to return */
/*****************************************************************/
void get_pi_hh_typed_wts (int *hh,
                       int *type,
                       double *x,
                       double *y,
                       double *weight,
                       int *len,
                       int *typeA,
                       int *typeB,
                       double *r_low,
                       double *r,
                       int *len_r,
                       int *inds,
                       double *rc) {
    
    int i,j,k;
    int num_cnt, denom_cnt; /*counters for those filling conditions*/
    double dist;
    
    /*repeat the calculation for all r*/
    for (i=0;i<*len_r;i++) {
        
        /*zero out counts*/
        num_cnt = 0;
        denom_cnt = 0;
        
        if (*typeA != -1) {
            
            for (j=0;j<*len;j++) {
                if (type[j] != *typeA) continue;
                
                for (k=0;k<*len;k++) {
                    /*ignore pairs of the same type*/
                    if (inds[k]==inds[j]) continue;
                    /*ignore pairs from the same household*/
                    if (hh[k]==hh[j]) continue;
                    
                    dist = sqrt(pow(x[j]-x[k],2)+pow(y[j]-y[k],2));
                    if ((dist<=r[i])  & (dist>=r_low[i])) denom_cnt = denom_cnt + weight[k]*weight[j];
                    
                    if (type[k] != *typeB) continue;
                    if ((dist<=r[i])  & (dist>=r_low[i])) num_cnt = num_cnt + weight[k]*weight[j];
                }
            }
            
        } else {
            Rprintf("To be implemented\n");
            return;
        }
        //Rprintf("%d/%d\n",num_cnt,denom_cnt);//DEBUG
        rc[i] = (double)num_cnt/denom_cnt;
    }
}




/*****************************************************************/
/* tau function for typed data                                   */
/* interates through a list of two types finding ones            */
/*that fulfill the distance requirement.                         */
/* this version takes in coordinates instead of distances        */
/*                                                               */
/*@param hh      the household column                            */
/*@param type_1  the  first type column                          */
/*@param type_2  the second type column                          */
/*@param x       the x coordinate                                */
/*@param y       the y coordinate                                */
/*@param len     the length of the three data arrays             */
/*@param typeA   the "from" type                                 */
/*@param typeB   the "to" type                                   */
/*@param r_low   the low end of the range of values to look at   */
/*@param r       the sequence of upper distances to consider     */
/*@param len_r   the number of different Rs to consider          */
/*                   usually will be just the indicies           */
/*@param inds    the indices into the original array, helps boot */
/*@param comp_type the comparison type                           */
/*@param rc      the array of values that we are going to return */
/*****************************************************************/
void get_tau_hh_typed (int *hh,
            int *type,
		    double *x,
		    double *y,
		    int *len,
		    int *typeA,
		    int *typeB,
		    double *r_low,
		    double *r,
		    int *len_r,
		    int *inds,
		    int *comp_type,
		    double *rc) {

    int i = 0;
    double divisor;
    double tmp_r_low = 0;
    double tmp_r = DBL_MAX;
    int tmp_len_r = 1;

    if (*comp_type==0) {
      /*get the divisor in the pi function*/
      get_pi_hh_typed(hh, type,x,y,len,typeA,typeB,&tmp_r_low, &tmp_r,
		   &tmp_len_r,inds,&divisor);

      /*get the main pi function*/
      get_pi_hh_typed(hh, type,x,y,len,typeA,typeB,r_low,r,len_r,inds,rc);
    } else if (*comp_type==1) {
      /*get the divisor in the pi function*/
      get_theta_hh_typed(hh, type,x,y,len,typeA,typeB,&tmp_r_low, &tmp_r,
		   &tmp_len_r,inds,&divisor);

      /*get the main pi function*/
      get_theta_hh_typed(hh, type, x, y, len, typeA, typeB, r_low, r, len_r,inds, rc);
    }

    for (i = 0; i < *len_r; i++) {
      if (rc[i] == divisor) { //deals with the Inf/Inf problem
    	if (rc[i] != 0) {
    	  rc[i] = 1;
    	} else {
    	  rc[i] = NAN;
    	}
      } else  {
    	rc[i] = rc[i]/divisor;
      }
    }
}





/*****************************************************************/
/* tau function for typed data                                   */
/* interates through a list of two types finding ones            */
/*that fulfill the distance requirement.                         */
/* this version takes in coordinates instead of distances        */
/*                                                               */
/*@param type_1  the  first type column                          */
/*@param type_2  the second type column                          */
/*@param x       the x coordinate                                */
/*@param y       the y coordinate                                */
/*@param weight  the sample weight of the individual/cluster     */
/*@param len     the length of the three data arrays             */
/*@param typeA   the "from" type                                 */
/*@param typeB   the "to" type                                   */
/*@param r_low   the low end of the range of values to look at   */
/*@param r       the sequence of upper distances to consider     */
/*@param len_r   the number of different Rs to consider          */
/*                   usually will be just the indicies           */
/*@param inds    the indices into the original array, helps boot */
/*@param rc      the array of values that we are going to return */
/*****************************************************************/
void get_tau_hh_typed_wts (int *hh,
                        int *type,
                        double *x,
                        double *y,
                        double *weight,
                        int *len,
                        int *typeA,
                        int *typeB,
                        double *r_low,
                        double *r,
                        int *len_r,
                        int *inds,
                        double *rc) {
    int i = 0;
    double divisor;
    double tmp_r_low = 0;
    double tmp_r = DBL_MAX;
    int tmp_len_r = 1;
    
    /*get the divisor in the pi function*/
    get_pi_hh_typed_wts(hh,type,x,y,weight,len,typeA,typeB,&tmp_r_low,&tmp_r,
                     &tmp_len_r,inds,&divisor);
    
    /*get the main pi function*/
    get_pi_hh_typed_wts(hh,type,x,y,weight,len,typeA,typeB,r_low,r,len_r,inds,rc);
    
    for (i = 0; i < *len_r; i++) {
        rc[i] = rc[i]/divisor;
    }
}


