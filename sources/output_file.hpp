#pragma once
#include <string>

struct output_file_c
{
    //for test
    output_file_c() = default;

    output_file_c(std::string const &directory);

    std::string Directory;
};