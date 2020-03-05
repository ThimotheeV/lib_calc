#include <iostream>

#include "calc_stat.hpp"
#include "output_file.hpp"

//TODO : Test dif opti float/double jeux de données de lézard avec 13 000 SNP

int main(int argc, char *argv[])
{
    std::string filename(argv[1]);
    genepop_input_c<2> input(filename);
    data_plane_vec_c data_plane_vec(input);

    auto ar = ar_by_pair(data_plane_vec);

    for (auto frac : ar)
    {
        std::cout<<frac.at(0)<<" : "<<frac.at(1)<<std::endl;
    }
}