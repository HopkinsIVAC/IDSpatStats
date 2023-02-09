#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Rdynload.h>
#include <Rdefines.h>
#include <stdbool.h>



/*****************************************************************/
/*pi function optimized for iterating through a list of two types*/
/*finding ones that fulfill the distance requirement.            */
/*                                                               */
/*@param type_1  the first type column                           */
/*@param type_2  the second type column                          */
/*@param x       the x coordinate                                */
/*@param y       the y coordinate                                */
/*@param weight  the sample weight of the individual/cluster     */
/*@param alpha   the weighted number of individuals/cluster      */
/*@param delta   the weighted number of non-vacc indivs/cluster  */
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
void get_pi_clustsurvey_window (double *p,
                                    double *x,
                                    double *y,
                                    int *urban,
                                    int *s,
                                    double *alpha,
                                    double *delta,
                                    int *len,
                                    double *r_low,
                                    double *r,
                                    int *len_r,
                                    int *inds,
                                    double *rc) {
    
    int i,j,k;
    int skip_r, skip_phi;  /*counters for those filling conditions*/
    double sum_numerator_i, sum_denom_i, sum_delta_j, sum_alpha_j;
    double dist;
    int dist_window_j, dist_window_i;
    bool correct_dist;
    
    //repeat the calculation for all r
    for (k=0;k<*len_r;k++) {
        sum_numerator_i = 0;
        sum_denom_i = 0;
        skip_r = 1;
        
        for (i=0;i<*len;i++) {
            // Skip cluster if no individuals with the characteristic
            if (p[i] == 0) continue;
            
            sum_delta_j = 0;
            sum_alpha_j = 0;
            skip_phi = 1;
            
            if(urban[i]==1) {
                dist_window_i = 2;
            } else {
                dist_window_i = 5;
            }
            
            // Cycle through paired clusters with cluster j
            //  - need to add s and r to each cluster within the correct distance to the sums
            //  - need to sum the sum_typeAj; 
            for (j=0;j<*len;j++) {
                
                // Check distance of clusters
                dist = sqrt(pow(x[i]-x[j],2) + pow(y[i]-y[j],2));
                
                if(urban[j]==1) {
                    dist_window_j = 2;
                } else {
                    dist_window_j = 5;
                }
                
                correct_dist = (dist<=(r[k] + dist_window_j + dist_window_i)) & (dist>=(r_low[k] - dist_window_j - dist_window_i));
                
                if (correct_dist) {
                    
                    if (inds[j] == inds[i]) {
                        sum_delta_j = sum_delta_j + delta[j] - 1;
                        sum_alpha_j = sum_alpha_j + alpha[j] - 1;
                    } 
                    else if (inds[j] != inds[i]) {
                        sum_delta_j = sum_delta_j + delta[j];
                        sum_alpha_j = sum_alpha_j + alpha[j];
                    }
                    
                    if (sum_alpha_j > 0) {
                        skip_phi = 0;
                        skip_r = 0;
                    }
                }
            }
            
            if (skip_phi == 1) continue;
            
            sum_numerator_i = sum_numerator_i + (delta[i] * sum_delta_j);
            sum_denom_i = sum_denom_i + (delta[i] * sum_alpha_j);
        }
        
        //Skip if no phi values for that distance of r
        if (skip_r == 1){
            rc[k] = NA_REAL;
        }else {
            rc[k] = (double) sum_numerator_i / sum_denom_i;
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
void get_tau_clustsurvey_window (double *p,
                                 double *x,
                                 double *y,
                                 int *urban,
                                 int *s,
                                 double *weight,
                                 double *delta,
                                 double *alpha,
                                 int *len,
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
    get_pi_clustsurvey_window(p,x,y,urban,s,delta,alpha,len,&tmp_r_low,&tmp_r,&tmp_len_r,inds,&divisor);
    
    /*get the main pi function*/
    get_pi_clustsurvey_window(p,x,y,urban,s,delta,alpha,len,r_low,r,len_r,inds,rc);
    
    for (i = 0; i < *len_r; i++) {
        rc[i] = rc[i]/divisor;
    }
}






/*****************************************************************/
/*pi function optimized for iterating through a list of two types*/
/*finding ones that fulfill the distance requirement.            */
/*                                                               */
/*@param type_1  the first type column                           */
/*@param type_2  the second type column                          */
/*@param x       the x coordinate                                */
/*@param y       the y coordinate                                */
/*@param weight  the sample weight of the individual/cluster     */
/*@param alpha   the weighted number of individuals/cluster      */
/*@param delta   the weighted number of non-vacc indivs/cluster  */
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
void get_pi_clustsurvey_wts_window (double *p,
                             double *x,
                             double *y,
                             int *urban,
                             int *s,
                             double *weight,
                             double *alpha,
                             double *delta,
                             int *len,
                             double *r_low,
                             double *r,
                             int *len_r,
                             int *inds,
                             double *rc) {
    
    int i,j,k;
    int skip_r, skip_phi;  /*counters for those filling conditions*/
double sum_numerator_i, sum_denom_i, sum_delta_j, sum_alpha_j;
double dist;
int dist_window_j, dist_window_i;
bool correct_dist;

//repeat the calculation for all r
for (k=0;k<*len_r;k++) {
    sum_numerator_i = 0;
    sum_denom_i = 0;
    skip_r = 1;
    
    for (i=0;i<*len;i++) {
        // Skip cluster if no individuals with the characteristic
        if (p[i] == 0) continue;
        
        sum_delta_j = 0;
        sum_alpha_j = 0;
        skip_phi = 1;
        
        if(urban[i]==1) {
            dist_window_i = 2;
        } else {
            dist_window_i = 5;
        }
        
        // Cycle through paired clusters with cluster j
        //  - need to add s and r to each cluster within the correct distance to the sums
        //  - need to sum the sum_typeAj; 
        for (j=0;j<*len;j++) {
            // Check distance of clusters
            dist = sqrt(pow(x[i]-x[j],2) + pow(y[i]-y[j],2));
            
            if(urban[j]==1) {
                dist_window_j = 2;
            } else {
                dist_window_j = 5;
            }
            
            if (k==0) {
                correct_dist = (dist<=r[k]) & (dist>=r_low[k]);
            } else {
                correct_dist = (dist<=(r[k] + dist_window_j + dist_window_i)) & (dist>=(r_low[k] - dist_window_j - dist_window_i));
            }
            
            if (correct_dist) {

                if (inds[j] == inds[i]) {
                    sum_delta_j = sum_delta_j + delta[j] - weight[j];
                    sum_alpha_j = sum_alpha_j + alpha[j] - weight[j];
                } 
                else if (inds[j] != inds[i]) {
                    sum_delta_j = sum_delta_j + delta[j];
                    sum_alpha_j = sum_alpha_j + alpha[j];
                }
                
                if (sum_alpha_j > 0) {
                    skip_phi = 0;
                    skip_r = 0;
                }
            }
        }
        
        if (skip_phi == 1) continue;
        
        sum_numerator_i = sum_numerator_i + (delta[i] * sum_delta_j);
        sum_denom_i = sum_denom_i + (delta[i] * sum_alpha_j);
    }
    
    //Skip if no phi values for that distance of r
    if (skip_r == 1){
        rc[k] = NA_REAL;
    }else {
        rc[k] = (double) sum_numerator_i / sum_denom_i;
    }
}
}






/*****************************************************************/
/*pi function optimized for iterating through a list of two types*/
/*finding ones that fulfill the distance requirement.            */
/*                                                               */
/*@param type_1  the first type column                           */
/*@param type_2  the second type column                          */
/*@param x       the x coordinate                                */
/*@param y       the y coordinate                                */
/*@param weight  the sample weight of the individual/cluster     */
/*@param alpha   the weighted number of individuals/cluster      */
/*@param delta   the weighted number of non-vacc indivs/cluster  */
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
void get_pi_clustsurvey_hh_wts_window (double *p,
                                double *x,
                                double *y,
                                int *urban,
                                int *s,
                                double *weight,
                                double *alpha,
                                double *delta,
                                double *rho1,
                                double *rho2,
                                int *len,
                                double *r_low,
                                double *r,
                                int *len_r,
                                int *inds,
                                double *rc) {
    
    int i,j,k;
    int skip_r, skip_phi;  /*counters for those filling conditions*/
double sum_numerator_i, sum_denom_i, sum_delta_j, sum_alpha_j;
double dist;
int dist_window_j, dist_window_i;
bool correct_dist;

//repeat the calculation for all r
for (k=0;k<*len_r;k++) {
    sum_numerator_i = 0;
    sum_denom_i = 0;
    skip_r = 1;
    
    for (i=0;i<*len;i++) {
        // Skip cluster if no individuals with the characteristic
        if (p[i] == 0) continue;
        
        sum_delta_j = 0;
        sum_alpha_j = 0;
        skip_phi = 1;
        
        if(urban[i]==1) {
            dist_window_i = 2;
        } else {
            dist_window_i = 5;
        }
        
        // Cycle through paired clusters with cluster j
        //  - need to add s and r to each cluster within the correct distance to the sums
        //  - need to sum the sum_typeAj; 
        for (j=0;j<*len;j++) {
            // Check distance of clusters
            dist = sqrt(pow(x[i]-x[j],2) + pow(y[i]-y[j],2));
            
            if(urban[j]==1) {
                dist_window_j = 2;
            } else {
                dist_window_j = 5;
            }
            
            if (k==0) {
                correct_dist = (dist<=r[k]) & (dist>=r_low[k]);
            } else {
                correct_dist = (dist<=(r[k] + dist_window_j + dist_window_i)) & (dist>=(r_low[k] - dist_window_j - dist_window_i));
            }
            
            if (correct_dist) {
                
                if (inds[i] != inds[j]) {
                    sum_delta_j = sum_delta_j + delta[j];
                    sum_alpha_j = sum_alpha_j + alpha[j];
                } 
                else if (inds[i] == inds[j]) {
                    sum_delta_j = sum_delta_j + (delta[j]-weight[j])*(1-rho1[j]);
                    sum_alpha_j = sum_alpha_j + (alpha[j]-weight[j])*(1-rho2[j]);
                }
                
                if (sum_alpha_j>0) {
                    skip_phi = 0;
                    skip_r = 0;
                }
            }
        }
        
        if (skip_phi == 1) continue;
        
        sum_numerator_i = sum_numerator_i + (delta[i] * sum_delta_j);
        sum_denom_i = sum_denom_i + (delta[i] * sum_alpha_j);
    }
    
    //Skip if no phi values for that distance of r
    if (skip_r == 1){
        rc[k] = NA_REAL;
    }else {
        rc[k] = (double) sum_numerator_i / sum_denom_i;
    }
}
}





/*****************************************************************/
/*pi function optimized for iterating through a list of two types*/
/*finding ones that fulfill the distance requirement.            */
/*                                                               */
/*@param type_1  the first type column                           */
/*@param type_2  the second type column                          */
/*@param x       the x coordinate                                */
/*@param y       the y coordinate                                */
/*@param weight  the sample weight of the individual/cluster     */
/*@param alpha   the weighted number of individuals/cluster      */
/*@param delta   the weighted number of non-vacc indivs/cluster  */
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
void get_pi_clustsurvey_hh_wts_window2 (double *p,
                                       double *x,
                                       double *y,
                                       int *urban,
                                       int *s,
                                       double *weight,
                                       double *alpha,
                                       double *delta,
                                       double *rho1,
                                       double *rho2,
                                       int *len,
                                       double *r_low,
                                       double *r,
                                       int *len_r,
                                       int *inds,
                                       double *rc) {
    
    int i,j,k;
    int skip_r, skip_phi;  /*counters for those filling conditions*/
double sum_numerator_i, sum_denom_i, sum_delta_j, sum_alpha_j;
double dist;
int dist_window_j, dist_window_i;
bool correct_dist;

//repeat the calculation for all r
for (k=0;k<*len_r;k++) {
    sum_numerator_i = 0;
    sum_denom_i = 0;
    skip_r = 1;
    
    for (i=0;i<*len;i++) {
        // Skip cluster if no individuals with the characteristic
        if (p[i] == 0) continue;
        
        sum_delta_j = 0;
        sum_alpha_j = 0;
        skip_phi = 1;
        
        if(urban[i]==1) {
            dist_window_i = 2;
        } else {
            dist_window_i = 5;
        }
        
        // Cycle through paired clusters with cluster j
        //  - need to add s and r to each cluster within the correct distance to the sums
        //  - need to sum the sum_typeAj; 
        for (j=0;j<*len;j++) {
            // Check distance of clusters
            dist = sqrt(pow(x[i]-x[j],2) + pow(y[i]-y[j],2));
            
            if(urban[j]==1) {
                dist_window_j = 2;
            } else {
                dist_window_j = 5;
            }
            
            if (k==0) {
                correct_dist = (dist<=r[k]) & (dist>=r_low[k]);
            } else {
                correct_dist = (dist<=(r[k] + dist_window_j + dist_window_i)) & (dist>=(r_low[k] - dist_window_j - dist_window_i));
            }
            
            if (correct_dist) {
                
   
                if ((inds[i] == inds[j]) & (k==0)) {
                    sum_delta_j = sum_delta_j + (delta[j]-weight[j])*(1-rho1[j]);
                    sum_alpha_j = sum_alpha_j + (alpha[j]-weight[j])*(1-rho2[j]);
                } else {
                    sum_delta_j = sum_delta_j + delta[j];
                    sum_alpha_j = sum_alpha_j + alpha[j];
                }
                
                if (sum_alpha_j>0) {
                    skip_phi = 0;
                    skip_r = 0;
                }
            }
        }
        
        if (skip_phi == 1) continue;
        
        sum_numerator_i = sum_numerator_i + (delta[i] * sum_delta_j);
        sum_denom_i = sum_denom_i + (delta[i] * sum_alpha_j);
    }
    
    //Skip if no phi values for that distance of r
    if (skip_r == 1){
        rc[k] = NA_REAL;
    }else {
        rc[k] = (double) sum_numerator_i / sum_denom_i;
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
void get_tau_clustsurvey_wts_window (double *p,
                              double *x,
                              double *y,
                              int *urban,
                              int *s,
                              double *weight,
                              double *delta,
                              double *alpha,
                              int *len,
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
    get_pi_clustsurvey_wts_window(p,x,y,urban,s,weight,delta,alpha,len,&tmp_r_low,&tmp_r,&tmp_len_r,inds,&divisor);
    
    /*get the main pi function*/
    get_pi_clustsurvey_wts_window(p,x,y,urban,s,weight,delta,alpha,len,r_low,r,len_r,inds,rc);
    
    for (i = 0; i < *len_r; i++) {
        rc[i] = rc[i]/divisor;
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
void get_tau_clustsurvey_hh_wts_window (double *p,
                                 double *x,
                                 double *y,
                                 int *urban,
                                 int *s,
                                 double *weight,
                                 double *delta,
                                 double *alpha,
                                 double *rho1,
                                 double *rho2,
                                 int *len,
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
    get_pi_clustsurvey_hh_wts_window(p,x,y,urban,s,weight,delta,alpha,rho1,rho2,len,&tmp_r_low,&tmp_r,&tmp_len_r,inds,&divisor);
    
    /*get the main pi function*/
    get_pi_clustsurvey_hh_wts_window(p,x,y,urban,s,weight,delta,alpha,rho1,rho2,len,r_low,r,len_r,inds,rc);
    
    for (i = 0; i < *len_r; i++) {
        rc[i] = rc[i]/divisor;
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
void get_tau_clustsurvey_hh_wts_window2 (double *p,
                                        double *x,
                                        double *y,
                                        int *urban,
                                        int *s,
                                        double *weight,
                                        double *delta,
                                        double *alpha,
                                        double *rho1,
                                        double *rho2,
                                        int *len,
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
    get_pi_clustsurvey_hh_wts_window2(p,x,y,urban,s,weight,delta,alpha,rho1,rho2,len,&tmp_r_low,&tmp_r,&tmp_len_r,inds,&divisor);
    
    /*get the main pi function*/
    get_pi_clustsurvey_hh_wts_window2(p,x,y,urban,s,weight,delta,alpha,rho1,rho2,len,r_low,r,len_r,inds,rc);
    
    for (i = 0; i < *len_r; i++) {
        rc[i] = rc[i]/divisor;
    }
}

