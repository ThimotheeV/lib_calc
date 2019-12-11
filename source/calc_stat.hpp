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
result_qstat calc_qstat_loc_by_loc(data_plane_vec_c const &genotype, int locus);
result_qstat calc_qstat_all_loc(data_plane_vec_c const &genotype);
std::vector<double> calc_qr_loc_by_loc(data_plane_vec_c const &genotype, int dist_max, int locus);
std::vector<double> calc_qr_all_loc(data_plane_vec_c const &genotype, int dist_max);
std::vector<double> Fst_qstat_loc_by_loc(std::vector<double> const &calc_qr_loc_by_loc, int locus);
std::vector<double> Fst_qstat_all_loc(std::vector<double> const &calc_qr_all_loc);

template <typename value>
double mean(std::vector<value> X_vec);
template <typename value>
double var(std::vector<value> X_vec, double meanX);
template <typename value1, typename value2>
double cov_X_Y(std::vector<value1> X_vec, double meanX, std::vector<value2> Y_vec, double meanY);
template <typename value1, typename value2>
std::array<double, 2> linear_regres_X_Y(std::vector<value1> X_vec, std::vector<value2> Y_vec);