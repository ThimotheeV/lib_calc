#pragma once
#include <string>

#include "input.hpp"

void output_stat_files(selector_input_c const &selec, result_c const &result);
void output_eta_stat_files(std::vector<std::array<double, 3>> result);

template <typename values>
void print_output(std::string path_to_file, std::vector<values> vec_value, std::string open_file_mode);

#include "output_file.tpp"