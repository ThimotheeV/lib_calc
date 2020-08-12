#pragma once

#include <string>
#include <vector>
#include <array>

std::string remove_spaces_underscores(std::string str);
std::string remove_spaces_in_range(std::string str, int pos_beg, int pos_end);
std::string lower_case(std::string str);
std::vector<std::string> sep_by_char(std::string const &str, char sep);
std::vector<std::string> trim_unix_windows_file_by_line(std::string const &str);
std::vector<std::vector<std::string>> trim_str_vec_by_string(std::vector<std::string> const &vec_str, std::string sep);
std::string const read_file(std::string const &filename);

template <std::size_t ploidy>
struct genepop_input_c
{
    genepop_input_c(){};
    //10 class => 11 limits
    genepop_input_c(std::string path_to_file, int nbr_class = 10);

    std::array<int, ploidy> trim_locus(std::string locus);
    void calc_dist_class_btw_deme(int nbr_class);

    std::string Header;
    std::vector<std::string> Locus_name;
    std::vector<std::vector<std::string>> Pop_name;
    std::vector<std::vector<double>> Dist_btw_deme;
    int Nbr_dist_class;
    std::vector<std::vector<int>> Dist_class_btw_deme;
    std::vector<std::vector<std::string>> Indiv_name;
    //Pop<Indiv<Name, Locus>
    std::vector<std::vector<std::vector<std::array<int, ploidy>>>> Genotype;
};

#include "input.tpp" 