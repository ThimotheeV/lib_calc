#include <iostream>
#include <iomanip>
#include <fstream>
#include <algorithm>

#include "output_file.hpp"
#include "data_plane_vec.hpp"

//TODO : data_plane_vec_c => flat_data_vec
//Coherance of display is handle by previous func => verification of calculus possibility
void output_stat_files(selector_input_c const &selec, result_c const &result)
{
    std::vector<std::string> head;
    head.reserve(20 + selec.Geo_dist_class_nbr + 2 + 2);

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
        for (int dist = 0; dist < selec.Geo_dist_class_nbr; ++dist)
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

    gss::print_output("./Stats.txt", head, "over");

    //+++++++++++++++++++++++++++++++++++++++++++++++//

    //for output
    std::vector<double> stats_run;
    stats_run.reserve(20 + selec.Geo_dist_class_nbr + 2 + 2);

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
        for (int dist = 0; dist < selec.Geo_dist_class_nbr; ++dist)
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
    gss::print_output("./Stats.txt", stats_run, "app");
}

void output_eta_stat_files(std::vector<std::array<double, 5>> result)
{
    //<pair of deme * pair of locus,<dist-deme, dist-locus, value eta>>
    std::sort(result.begin(), result.end(),
              [](auto const &a, auto const &b)
              {
                  if (a.at(0) < b.at(0))
                  {
                      return true;
                  }
                  else
                  {
                      if (a.at(0) == b.at(0))
                      {
                          return a.at(1) < b.at(1);
                      }
                      else
                      {
                          return false;
                      }
                  }
              });

    std::vector<std::string> head;
    head.reserve(6);

    head.emplace_back("Geo_dist_btw_deme");
    head.emplace_back("Dist_btw_locus_pb");
    head.emplace_back("Sum_Phi");
    head.emplace_back("Sum_Q1_join");
    head.emplace_back("Sum_Eta_denum");
    head.emplace_back("Num_pt");

    gss::print_output("./Stats_dl.txt", head, "over");

    //+++++++++++++++++++++++++++++++++++++++++++++++//

    double dist_geo = result[0].at(0);
    double dist_chr = result[0].at(1);
    double phi = 0;
    double q1_join = 0;
    double eta_denum = 0;
    double num_pt = 0;
    for (auto const &values : result)
    {
        if ((dist_geo == values.at(0)) && (dist_chr == values.at(1)))
        {
            phi += values.at(2);
            q1_join += values.at(3);
            eta_denum += values.at(4);
            num_pt += 1;
        }
        else
        {
            gss::print_output<double>("./Stats_dl.txt", {dist_geo, dist_chr, phi, q1_join, eta_denum, num_pt},
                                      "app");
            dist_geo = values.at(0);
            dist_chr = values.at(1);
            num_pt = 0;
            phi = 0;
            q1_join = 0;
            eta_denum = 0;
        }
    }
}

#include "calc_stat_dl.hpp"
void output_exp_regr_eta_stat_files(std::vector<std::array<double, 5>> result)
{
    //<pair of deme * pair of locus,<dist-deme, dist-locus, value eta>>
    //sort by dist-locus and by dist-deme
    std::sort(result.begin(), result.end(),
              [](auto const &a, auto const &b)
              {
                  if (a.at(1) < b.at(1))
                  {
                      return true;
                  }
                  else
                  {
                      if (a.at(1) == b.at(1))
                      {
                          return a.at(0) < b.at(0);
                      }
                      else
                      {
                          return false;
                      }
                  }
              });

    std::vector<std::string> head;
    head.reserve(4);

    head.emplace_back("Dist_btw_locus_pb");
    head.emplace_back("a");
    head.emplace_back("b");
    head.emplace_back("b_g");

    gss::print_output("./eta_exp_regr.txt", head, "over");

    auto eta_by_chr_dist = std::vector<std::array<double, 3>>{};
    double num_pt = 0;
    double dist_chr = result[0].at(1);
    double dist_geo = result[0].at(0);
    double eta = 0;

    for (auto const &values : result)
    {
        if (dist_chr == values.at(1))
        {
            if (dist_geo == values.at(0))
            {
                eta += (values.at(2) - values.at(3)) / values.at(4);
                num_pt += 1;
            }
            else
            {
                eta_by_chr_dist.push_back({values.at(0), eta, num_pt});
                eta = 0;
                num_pt = 0;
                dist_geo = values.at(0);
            }
        }
        else
        {
            // std::cout<<"--------------------- Dist_btw_locus_pb : "<<dist_chr<<" ---------------------"<<std::endl;
            auto temp = exp_regr(eta_by_chr_dist);
            gss::print_output<double>("./eta_exp_regr.txt", {dist_chr, temp.at(0), temp.at(1), temp.at(2)}, "app");
            dist_chr = values.at(1);
            dist_geo = values.at(0);
            eta = 0;
            num_pt = 0;
            eta_by_chr_dist.resize(0);
        }
    }
}

void output_sfs_stat_files(std::map<int, double> const &result)
{
    std::vector<std::string> head;
    head.reserve(2);

    head.emplace_back("Loci_count");
    head.emplace_back("Allele_count_class");

    gss::print_output("./SFS.txt", head, "over");

    //+++++++++++++++++++++++++++++++++++++++++++++++//

    for (auto const &values : result)
    {
        gss::print_output<double>("./SFS.txt", {static_cast<double>(values.first), values.second}, "app");
    }
}