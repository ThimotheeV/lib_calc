#include <algorithm>
#include <sstream>
#include <fstream>
#include <iterator>
#include <functional>
#include <iostream>
#include <cmath>

#include "input.hpp"

std::string const gss::read_file(std::string const &filename)
{
    //Put filename in a stream
    //ifstream constructor
    std::ifstream file_str(filename);
    if (!file_str)
    {
        throw std::logic_error("( Could not open file : " + filename + ". I exit. )");
    }
    else
    {
        //Test if the stream is good
        if (file_str.bad())
        {
            throw std::logic_error("( Could not open file : " + filename + ". I exit. )");
        }
        //Don't skip space
        file_str >> std::noskipws;
        //Convert the stream in a string
        const std::string &fstr{std::istream_iterator<char>{file_str}, {}};
        return fstr;
    }
}

/**
 * @brief Parser file to param 
 * 
 */

std::string gss::read_write_cmdline(int argc, char **argv)
{
    std::string result;
    if (argc > 1)
    {
        for (int arg = 1; arg < argc - 1; ++arg)
        {
            result += argv[arg];
            result += "\n";
        }
        result += argv[argc - 1];
    }
    //Put filename in a stream
    std::ofstream file_str("cmdline_settings.txt");
    if (!file_str)
    {
        throw std::logic_error("( Could not open file 'cmdline_settings.txt'. I exit. )");
    }
    else
    {
        file_str << argv[0] << "\n";
        file_str << result;
        file_str.close();
    }

    return result;
}

selector_input_c::selector_input_c(std::string path_to_settings_file)
{
    const auto &file_str = gss::read_file(path_to_settings_file);
    const auto &str_vec = gss::slice_unix_windows_file_by_line(file_str);

    for (auto const &line : str_vec)
    {
        //Clean and parse key=value in a vector[key, value]
        const auto &clean_line = gss::remove_spaces_tab_in_range(line, 0, static_cast<int>(line.size()));
        auto line_pair_key_value = gss::slice_by_char(clean_line, '=');
        const auto &stat_name = gss::remove_underscores(gss::str_tolower(line_pair_key_value[0]));

        if (stat_name == "datafilename")
        {
            Data_filename = line_pair_key_value[1];
            continue;
        }
        if (stat_name == "geneticmapname")
        {
            Genetic_map_name = line_pair_key_value[1];
            continue;
        }
        if (stat_name == "hobs")
        {
            Hobs = gss::convert_str_bool(line_pair_key_value[1]);
            continue;
        }
        if (stat_name == "hexp")
        {
            Hexp = gss::convert_str_bool(line_pair_key_value[1]);
            continue;
        }
        if (stat_name == "nballele")
        {
            Nb_allele = gss::convert_str_bool(line_pair_key_value[1]);
            continue;
        }
        if (stat_name == "var")
        {
            Var = gss::convert_str_bool(line_pair_key_value[1]);
            continue;
        }
        if (stat_name == "mgw")
        {
            MGW = gss::convert_str_bool(line_pair_key_value[1]);
            continue;
        }
        if (stat_name == "sfs")
        {
            SFS = gss::convert_str_bool(line_pair_key_value[1]);
            continue;
        }
        if (stat_name == "mingene")
        {
            //Needed to read integer in scientific format
            min_gene_for_SFS = static_cast<int>(std::stod(line_pair_key_value[1]));
            continue;
        }
        if (stat_name == "fstat")
        {
            F_stat = gss::convert_str_bool(line_pair_key_value[1]);
            continue;
        }
        if (stat_name == "qstat")
        {
            Q_stat = gss::convert_str_bool(line_pair_key_value[1]);
            continue;
        }
        if (stat_name == "qr")
        {
            Qr = gss::convert_str_bool(line_pair_key_value[1]);
            continue;
        }
        if (stat_name == "ar")
        {
            Ar = gss::convert_str_bool(line_pair_key_value[1]);
            continue;
        }
        if (stat_name == "er")
        {
            Er = gss::convert_str_bool(line_pair_key_value[1]);
            continue;
        }
        if (stat_name == "eta")
        {
            Eta = gss::convert_str_bool(line_pair_key_value[1]);
            continue;
        }
        if (stat_name == "missingdata")
        {
            Missing_data = gss::convert_str_bool(line_pair_key_value[1]);
            continue;
        }
        if (stat_name == "geodistclassnbr")
        {
            Geo_dist_class_nbr = static_cast<int>(std::stod(line_pair_key_value[1]));
            continue;
        }
        if (stat_name == "chrdistclassnbr")
        {
            Chr_dist_class_nbr = static_cast<int>(std::stod(line_pair_key_value[1]));
            continue;
        }
        throw std::invalid_argument("( Unknown keyworld : " + line + " from cmdline_settings.txt or glib_settings.txt. I exit. )");
    }
}

result_c::result_c(int nbr_geo_dist_class)
{
    Qr = std::vector<double>(nbr_geo_dist_class, -1);
}