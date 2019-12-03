#pragma once

#include <string>
#include <vector>
#include <array>

std::string trim_spaces_underscores(std::string str);
std::string lower_case(std::string str);
std::vector<std::string> trim_by_char(std::string const &str, char sep);
std::vector<std::string> trim_unix_windows_file_by_line(std::string const &str);
std::vector<std::vector<std::string>> trim_vec_by_string(std::vector<std::string> const &vec_str, std::string sep);
std::string const read_file(std::string const &filename);

template <std::size_t ploidy>
struct genepop_input_c
{
    genepop_input_c(){};
    genepop_input_c(std::string path_to_file);

    std::array<int, ploidy> trim_locus(std::string locus);

    //Header + Locus name
    std::vector<std::string> Locus_name;
    std::vector<std::string> Pop_name;
    std::vector<std::vector<std::string>> Indiv_name;
    //Pop<Indiv<Name, Locus>
    std::vector<std::vector<std::vector<std::array<int, ploidy>>>> Genotype;
};