#include <Rcpp.h>
#include <iostream>
#include <fstream>
#include <string>
#include <iterator>
#include <cmath>
#include <math.h>
#include <numeric>
#include <algorithm>
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <map>


using namespace std;
using namespace Rcpp;



struct get_log
{
    double operator()(double c) const
    {
        return log(c);
    }
};
struct get_lgamma
{
    double operator()(double c) const
    {
        return lgamma(c+1);
    }
};

bool IsNAN(double i) { return (std::isnan(i)); }




//[[Rcpp::export]]
vector<double> calProb_each(NumericMatrix decompose_mtx2,
    vector<double> ek_count,
    vector<vector<double>> pk_dist, 
    vector<int> N_list,
    double recruits_num,
    vector<vector<double>> num_subtract) {
    
    double d_ave = inner_product(pk_dist[0].begin(), pk_dist[0].end(), pk_dist[1].begin(), 0.0); //average degree of population
    double uni_size = accumulate(ek_count.begin(),ek_count.end(),0.0);
    double dup_num = recruits_num - uni_size;

    double d_bar = inner_product(pk_dist[0].begin(), pk_dist[0].end(), ek_count.begin(), 0.0)/uni_size;//average degree of recruits
    vector<double>numerator_log, denom_log;
    
    for (int N : N_list) {
        double denom_log_i=0.0;
        for (int i = 0; i < num_subtract[0].size(); i++) {
            denom_log_i = denom_log_i + 
                num_subtract[1][i]*log(N*d_ave- d_bar * num_subtract[0][i]);
        }
        denom_log.push_back(denom_log_i);
        
        double numerator_log_i = 0.0;
        for (int ek_loc = 0; ek_loc < ek_count.size(); ek_loc++) {
            if (ek_count[ek_loc] > 0) {
                //Nk_prod_log, same for 0 or duplicates
                vector<double> Nk(ek_count[ek_loc]); // vector with 'ek_count[ek_loc]' ints.
                iota(begin(Nk), end(Nk), N * pk_dist[1][ek_loc] - ek_count[ek_loc] + 1);
                vector<double>Nk_log(Nk.size());
                transform(Nk.begin(), Nk.end(), Nk_log.begin(), get_log());
                numerator_log_i = numerator_log_i + accumulate(Nk_log.begin(), Nk_log.end(), 0.0) + ek_count[ek_loc] * log(pk_dist[0][ek_loc]); // ek_prod_log, here 0 duplicates only
                Nk.clear();
                Nk_log.clear();
            }
        }
        numerator_log.push_back(numerator_log_i);
       
    }
    //numerator_log-denom_log
    transform(numerator_log.begin(), numerator_log.end(),
        denom_log.begin(), numerator_log.begin(), minus<double>());
    if (dup_num == 0) {
        vector<double>ek_log(ek_count.size());
        transform(ek_count.begin(), ek_count.end(), ek_log.begin(), get_lgamma());
        double sum = accumulate(ek_log.begin(), ek_log.end(), 0.0);
        double times = lgamma(uni_size + 1) - sum;
        for_each(numerator_log.begin(), numerator_log.end(), [&](double& final) {final =exp(final+ times); });
        replace_if(numerator_log.begin(), numerator_log.end(), IsNAN, 0.0);
    }
    if (dup_num > 0) {
        //list combinations of dups
        vector<double> ek_count2;
        vector<double> ek_count2_lgamma(ek_count.size());
        vector<double> ek_diff(ek_count.size(), 0.0);
        vector<double> dups_log;
        NumericVector decom2;
        double p_dup=0.0;
        double p_temp=0.0;
        for (int i = 0; i < decompose_mtx2.ncol(); i++) {
            decom2 = decompose_mtx2.column(i);
            ek_count2 = as<vector<double>>(decom2);
            sort(ek_count2.begin(), ek_count2.end());
            do {
                p_temp=0.0;
                transform(ek_count.begin(), ek_count.end(), ek_count2.begin(),
                    ek_diff.begin(), minus<double>());
                if (*min_element(ek_diff.begin(), ek_diff.end()) < 0)
                    continue;
                //duplicated degrees
                for (int loc = 0; loc < ek_count2.size(); loc++)
                    p_temp += ek_count2[loc] * log(pk_dist[0][loc]);
        
        //count repeated counts
        vector<double>order;
        map<double, double> order_matter;
        for_each(ek_count2.begin(), ek_count2.end(), [&order_matter](double val) { order_matter[val]++; });
        for (auto p : order_matter) {
            order.push_back(p.second );
        }
        vector<double>order_lgamma(order.size());
        transform(order.begin(), order.end(), order_lgamma.begin(), get_lgamma());


               /* //count repeated counts
               vector<double> rep;
               vector<double> nums(ek_count2.begin(),ek_count2.end());
               sort(nums.begin(), nums.end());
               for(auto it = begin(nums); it != end(nums); ) {  
                   int dups = count(it, end(nums), *it);
                   if ( dups > 1 ){
                       rep.push_back(dups);
                   }
                   for(auto last = *it;*++it == last;);
               }*/

                transform(ek_count2.begin(), ek_count2.end(), ek_count2_lgamma.begin(), get_lgamma());
                transform(ek_diff.begin(), ek_diff.end(), ek_diff.begin(), get_lgamma());

                p_temp +=//- accumulate(order_lgamma.begin(), order_lgamma.end(), 0.0)+-lgamma(1+accumulate(begin(rep),end(rep),1,multiplies<double>()))+ 
                    lgamma(dup_num + 1) - accumulate(ek_count2_lgamma.begin(), ek_count2_lgamma.end(), 0.0) +
                    lgamma(uni_size - dup_num + 1) - accumulate(ek_diff.begin(), ek_diff.end(), 0.0);
                p_dup += exp(p_temp);
            } while (next_permutation(ek_count2.begin(), ek_count2.end()));
        }
            for_each(numerator_log.begin(), numerator_log.end(), [&](double& final) {final =p_dup*exp(final); });
            replace_if(numerator_log.begin(), numerator_log.end(), IsNAN, 0.0);
    }

    
    return numerator_log;
}
//[[Rcpp::export]]
NumericVector calProb_dup(NumericMatrix pk_dist_mtx, NumericVector N_list_vec, 
    NumericMatrix decompose_mtx,NumericMatrix decompose_mtx2, double recruits_num, NumericMatrix num_subtract_mtx) {
    //read in pk_dist
    vector<vector<double>> pk_dist;
    vector<double> temp;
    NumericVector col;
    col = pk_dist_mtx.column(0);
    temp = as<vector<double>>(col);
    pk_dist.push_back(temp);
    col = pk_dist_mtx.column(1);
    temp = as<vector<double>>(col);
    pk_dist.push_back(temp);
    //read in N_list
    vector<int> N_list = as<vector<int>> (N_list_vec);
    //read in num_subtract    
    vector<vector<double>> num_subtract;
    col = num_subtract_mtx.column(0);
    temp = as<vector<double>>(col);
    num_subtract.push_back(temp);
    col = num_subtract_mtx.column(1);
    temp = as<vector<double>>(col);
    num_subtract.push_back(temp);


    //save results
    vector<double>p_list(N_list.size(),0);
    vector<double>p;
    //list combinations
    vector<double> ek_count;
    NumericVector decom;

    for (int i = 0; i < decompose_mtx.ncol(); i++){
        decom=decompose_mtx.column(i);
        ek_count=as<vector<double>>(decom);
        sort(ek_count.begin(), ek_count.end());
        do{
            p=calProb_each(decompose_mtx2, ek_count, pk_dist, N_list, recruits_num,num_subtract);
            
            transform(p_list.begin(), p_list.end(),p.begin(), p_list.begin(), plus<double>());
        }while(next_permutation(ek_count.begin(), ek_count.end()));
    }
 
  
   NumericVector out=wrap(p_list);
   return out;
 
}


