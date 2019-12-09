#pragma once
#include <array>

#include "custom_vec.hpp"

//TODO : vecteur de mot clé pour pouvoir comparer un ensemble d'attribut (indiv, pop, habitat) => fonction de comparaison
double calc_Q0(data_plane_vec_c const &genotype);
double calc_Q1(data_plane_vec_c const &genotype);
double calc_Q2(data_plane_vec_c const &genotype);
//Fst = Q1, Q2
double calc_Fstat(double Qx, double Qy);

struct result_qstat
{
    double Q0_intra_ind{0};
    double Q1_intra_pop{0};
    double Q2_inter_pop{0};
};

//Other way to calc
result_qstat calc_qstat(data_plane_vec_c const &genotype);

std::vector<double> calc_qr(data_plane_vec_c const &genotype, int dist_max);

template<typename value>
double var(std::vector<value> vec);
template<typename value>
double mean(std::vector<value> vec);
template<typename value>
double cov_X_Y(std::vector<value> X_vec, std::vector<value> Y_vec);
