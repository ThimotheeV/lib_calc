#pragma once
#include <string>
#include <map>

#include "input.hpp"

void output_stat_files(selector_input_c const &selec, result_c const &result, std::string path_to_output_file);
void output_eta_stat_files(std::vector<std::array<double, 5>> result, std::string path_to_output_file);
void output_exp_regr_eta_stat_files(std::vector<std::array<double, 5>> result);
void output_sfs_stat_files(std::map<int, double> const &result, std::string path_to_output_file);

namespace gss
{
    template <typename values>
    void print_output(std::string path_to_settings_file, std::vector<values> vec_value, std::string open_file_mode);
}
#include "output_file.tpp"