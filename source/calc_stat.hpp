#pragma once
#include <array>

#include "custom_vec.hpp"

//TODO : calculer tout ensemble
//TODO : vecteur de mot clé pour pouvoir comparer un ensemble d'attribut (indiv, pop, habitat) => fonction de comparaison
std::array<int, 2> calc_Q0(data_plane_vec_c const &genotype);
std::array<int, 2> calc_Q1(data_plane_vec_c const &genotype);
std::array<int, 2> calc_Q2(data_plane_vec_c const &genotype);
std::array<int, 2> calc_Fst(std::array<int, 2> Q1, std::array<int, 2> Q2);

struct result_qstat
{
    std::array<int, 2> Q0{0,0};
    std::array<int, 2> Q1{0,0};
    std::array<int, 2> Q2{0,0};
};

//Other way to calc
result_qstat calc_qstat(data_plane_vec_c const &genotype);