#include <iostream>

#include "calc_stat.hpp"
#include "output_file.hpp"

int main(int argc, char *argv[])
{
    std::string filename(argv[1]);
    genepop_input_c<1> input(filename);
    data_plane_vec_c data_plane_vec(input);
    auto result = calc_qstat_all_loc(data_plane_vec);

    std::cout << "Q1 : " << result.Q1_intra_pop << std::endl;
    return 0;
}