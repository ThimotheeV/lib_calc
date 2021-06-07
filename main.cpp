#include <iostream>
#include <iomanip>
#include <fstream>
#include <cassert>
#include <chrono>

#include "run.hpp"

std::string const datestring = __DATE__;
std::string const timestring = __TIME__;
std::string const version = "0.9 (Built on " + datestring + " at " + timestring + ")";

int main(int argc, char *argv[])
{
    //Allow to approximately calculate run-time duration
    auto debut = std::chrono::high_resolution_clock::now();

    std::cout << "************Begin run************" << std::endl;
    std::string file_selector;

    if (argc == 1)
    {
        file_selector = "GSumStatSettings.txt";
    }
    else
    {
        file_selector = std::string(argv[1]);
    }

    selector_input_c selector(file_selector);

    if (selector.Ploidy == 1)
    {
        genepop_input_c<1> input(selector.Data_filename, selector.Geo_dist_class_nbr, selector.Genetic_map_name, selector.Chr_dist_class_nbr);
        data_plane_vec_c data_plane_vec(input);
        run(selector, data_plane_vec);
    }
    else
    {
        genepop_input_c<2> input(selector.Data_filename, selector.Geo_dist_class_nbr, selector.Genetic_map_name, selector.Chr_dist_class_nbr);
        data_plane_vec_c data_plane_vec(input);
        run(selector, data_plane_vec);
    }

    std::cout << "\nTotal execution time  is ";
    auto fin = std::chrono::high_resolution_clock::now();
    auto time = std::chrono::duration_cast<std::chrono::nanoseconds>(fin - debut).count() * 0.000000001;
    std::cout << time << " seconds" << std::endl;

    std::cout << "\nNormal ending of GSumStat.\n"
              << std::endl;
}