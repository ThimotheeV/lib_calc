#pragma once

#include <string>
#include <vector>
#include <array>
#include <cmath>

#include "arg_pars.hpp"
namespace gss
{
    std::string const read_file(std::string const &filename);
    std::string read_write_cmdline(int argc, char **argv);
}
template <std::size_t ploidy>
struct genepop_input_c
{
    genepop_input_c(){};
    //10 class => 11 limits
    genepop_input_c(std::string path_to_genepop_file, int nbr_dist_class = 10, std::string path_to_chr_map_file = "", int nbr_chr_dist_class = 0);

    std::array<int, ploidy> trim_locus(std::string locus);
    void calc_dist_class_btw_deme(int nbr_dist_class);
    void calc_dist_btw_loc(std::vector<std::tuple<int, std::string, double, int>> const &map_vec, int chr, int nbr_chr_dist_class);

    std::string Header;
    std::vector<std::string> Locus_name;
    //Chr<loc1, loc2>
    std::vector<std::vector<std::vector<double>>> Dist_btw_loc;
    int Nbr_chr_dist_class;
    //Chr<loc1, loc2>
    std::vector<std::vector<std::vector<int>>> Dist_class_btw_loc;

    std::vector<std::vector<std::string>> Pop_name;
    std::vector<std::vector<double>> Dist_btw_deme;
    int Nbr_dist_class;
    std::vector<std::vector<int>> Dist_class_btw_deme;

    std::vector<std::vector<std::string>> Indiv_name;
    //Pop<Indiv<Name, Locus>
    std::vector<std::vector<std::vector<std::array<int, ploidy>>>> Genotype;
};

//simple class to read file and chose
struct selector_input_c
{
    selector_input_c() = default;
    selector_input_c(std::string path_to_file);

    std::string Input_name;

    bool Hobs{false};
    bool Hexp{false};
    bool Nb_allele{false};
    bool Var{false};
    bool MGW{false};

    bool SFS{false};
    int min_gene_for_SFS{-1};

    bool F_stat{false};
    bool Q_stat{false};

    bool Qr{false};
    bool Ar{false};
    bool Er{false};

    bool Eta{false};

    bool Missing_data{false};
    int Nbr_class{1};
};

struct result_c
{
    result_c(int nbr_class);

    double Hobs_mean{-1};
    double Hobs_var{-1};
    double Hexp_mean{-1};
    double Hexp_var{-1};
    double Nb_allele_mean{-1};
    double Nb_allele_var{-1};
    double Var_mean{-1};
    double Var_var{-1};
    double MGW_mean{-1};
    double MGW_var{-1};

    double Fis_mean{-1};
    double Fis_var{-1};
    double Fst_mean{-1};
    double Fst_var{-1};

    double Qwi_mean{-1};
    double Qwi_var{-1};
    double Qwd_mean{-1};
    double Qwd_var{-1};
    double Qbd_mean{-1};
    double Qbd_var{-1};

    std::vector<double> Qr;

    std::array<double, 3> Ar_reg{-1};
    std::array<double, 3> Er_reg{-1};
};

#include "input.tpp"