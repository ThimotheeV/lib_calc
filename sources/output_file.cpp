#include <iostream>
#include <iomanip>
#include <fstream>

#include "output_file.hpp"
#include "data_plane_vec.hpp"

//TODO : data_plane_vec_c => flat_data_vec
//Coherance of display is handle by previous func => verification of calculus possibility
void output_stat_files(selector_input_c const &selec, result_c const &result)
{
    std::vector<std::string> head;
    head.reserve(20 + selec.Nbr_class + 2 + 2);

    if (selec.Hobs)
    {
        head.emplace_back("Hobs_mean");
        head.emplace_back("Hobs_var");
    }
    if (selec.Hexp)
    {
        head.emplace_back("Hexp_mean");
        head.emplace_back("Hexp_var");
    }
    if (selec.Nb_allele)
    {
        head.emplace_back("Nb_allele_mean");
        head.emplace_back("Nb_allele_var");
    }
    if (selec.Var)
    {
        head.emplace_back("Var_mean");
        head.emplace_back("Var_var");
    }
    if (selec.MGW)
    {
        head.emplace_back("MGW_mean");
        head.emplace_back("MGW_var");
    }

    if (selec.F_stat)
    {
        head.emplace_back("Fis_mean");
        head.emplace_back("Fis_var");
        head.emplace_back("Fst_mean");
        head.emplace_back("Fst_var");
    }

    if (selec.Q_stat)
    {
        if (result.Qwi_mean != -1)
        {
            head.emplace_back("Qwi_mean");
            head.emplace_back("Qwi_var");
        }

        head.emplace_back("Qwd_mean");
        head.emplace_back("Qwd_var");

        if (result.Qbd_mean != -1)
        {
            head.emplace_back("Qbd_mean");
            head.emplace_back("Qbd_var");
        }
    }

    if (selec.Qr)
    {
        for (int dist = 0; dist < selec.Nbr_class; ++dist)
        {
            std::string temp = "Q" + std::to_string(dist);
            head.emplace_back(temp);
        }
    }

    if (selec.Ar)
    {
        head.emplace_back("Ar_slop");
        head.emplace_back("Ar_intersec");
    }

    if (selec.Er)
    {
        head.emplace_back("Er_slop");
        head.emplace_back("Er_intersec");
    }
    
    head.shrink_to_fit();

    print_output("./Stats.txt", head, "over");

    //+++++++++++++++++++++++++++++++++++++++++++++++//

    //for output
    std::vector<double> stats_run;
    stats_run.reserve(20 + selec.Nbr_class + 2 + 2);

    if (selec.Hobs)
    {
        stats_run.emplace_back(result.Hobs_mean);
        stats_run.emplace_back(result.Hobs_var);
    }
    if (selec.Hexp)
    {
        stats_run.emplace_back(result.Hexp_mean);
        stats_run.emplace_back(result.Hexp_var);
    }
    if (selec.Nb_allele)
    {
        stats_run.emplace_back(result.Nb_allele_mean);
        stats_run.emplace_back(result.Nb_allele_var);
    }
    if (selec.Var)
    {
        stats_run.emplace_back(result.Var_mean);
        stats_run.emplace_back(result.Var_var);
    }
    if (selec.MGW)
    {
        stats_run.emplace_back(result.MGW_mean);
        stats_run.emplace_back(result.MGW_var);
    }

    if (selec.F_stat)
    {
        stats_run.emplace_back(result.Fis_mean);
        stats_run.emplace_back(result.Fis_var);
        stats_run.emplace_back(result.Fst_mean);
        stats_run.emplace_back(result.Fst_var);
    }

    if (selec.Q_stat)
    {
        if (result.Qwi_mean != -1)
        {
            stats_run.emplace_back(result.Qwi_mean);
            stats_run.emplace_back(result.Qwi_var);
        }

        stats_run.emplace_back(result.Qwd_mean);
        stats_run.emplace_back(result.Qwd_var);

        if (result.Qbd_mean != -1)
        {
            stats_run.emplace_back(result.Qbd_mean);
            stats_run.emplace_back(result.Qbd_var);
        }
    }

    if (selec.Qr)
    {
        for (int dist = 0; dist < selec.Nbr_class; ++dist)
        {
            stats_run.emplace_back(result.Qr.at(dist));
        }
    }

    if (selec.Ar)
    {
        stats_run.emplace_back(result.Ar_reg.at(0));
        stats_run.emplace_back(result.Ar_reg.at(1));
    }

    if (selec.Er)
    {
        stats_run.emplace_back(result.Er_reg.at(0));
        stats_run.emplace_back(result.Er_reg.at(1));
    }

    stats_run.shrink_to_fit();
    print_output("./Stats.txt", stats_run, "app");
}
