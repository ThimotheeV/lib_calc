#include <iostream>
#include <iomanip>
#include <fstream>
#include <cassert>

#include "exe.hpp"

int main(int argc, char *argv[])
{
    std::string file_selector;
    if (argc == 1)
    {
        file_selector = "./selector.txt";
    }
    else
    {
        file_selector = std::string(argv[1]); 
    }

    selector_input_c selector(file_selector);

    genepop_input_c<2> input(selector.Input_name, selector.Nbr_class);
    data_plane_vec_c data_plane_vec(input);

    result_c result = run(selector, data_plane_vec);

    output_stat_files(selector, result);
}