#include <iostream>
#include <iomanip>
#include <fstream>
#include <cassert>

#include "run.hpp"

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

    genepop_input_c<2> input(selector.Data_filename, selector.Nbr_geo_dist_class, selector.Genetic_map_name, selector.Nbr_chr_dist_class);
    data_plane_vec_c data_plane_vec(input);

    run(selector, data_plane_vec);
}