//TODO : Handle empty line
//TODO : Vérifier que le nbr de locus est le même pour tout les indivs
//Constructor of genepop_input_c
#include <iostream>
#include <algorithm>

template <std::size_t ploidy>
genepop_input_c<ploidy>::genepop_input_c(std::string path_to_genepop_file, int nbr_dist_class, std::string path_to_chr_map_file, int nbr_chr_dist_class)
{
    auto const file_str = gss::str_tolower(gss::read_file(path_to_genepop_file));
    auto const file_str_vec = gss::slice_unix_windows_file_by_line(file_str);
    //TODO : Utiliser un vecteur de séparateur
    std::vector<std::vector<std::string>> const deme_vec = gss::slice_str_vec_by_string(file_str_vec, "pop");
    //Special case of header => first line + locus name => in deme_vec[0]
    Locus_name.reserve(deme_vec[0].size() - 1);
    Header = deme_vec[0][0];
    for (int locus = 1; locus < deme_vec[0].size(); ++locus)
    {
        auto const temp_str_vec = gss::slice_by_char(deme_vec[0][locus], ',');
        //No need to copy header in memory if temp_str_vec haven't more than 1 string in it
        Locus_name.reserve(Locus_name.size() + temp_str_vec.size() - 1);
        Locus_name.insert(Locus_name.end(), temp_str_vec.begin(), temp_str_vec.end());
    }
    Locus_name.shrink_to_fit();

    for (int locus = 0; locus < Locus_name.size(); ++locus)
    {
        Locus_name[locus] = gss::remove_spaces_tab_underscores(Locus_name[locus]);
    }

    //Trim deme by indiv and locus
    Genotype.resize(deme_vec.size() - 1);
    Indiv_name.resize(deme_vec.size() - 1);
    Pop_name.resize(deme_vec.size() - 1);

    //deme level, first was header
    for (int nbr_deme = 1; nbr_deme < deme_vec.size(); ++nbr_deme)
    {
        //indiv level
        std::vector<std::vector<std::array<int, ploidy>>> temp_deme;
        //deme_vec[nbr_deme].size() => nbr of indiv by deme
        temp_deme.resize(deme_vec[nbr_deme].size());
        Indiv_name[nbr_deme - 1].resize(deme_vec[nbr_deme].size());

        for (int indiv = 0; indiv < deme_vec[nbr_deme].size(); ++indiv)
        {
            auto const temp_vec = gss::slice_by_char(deme_vec[nbr_deme][indiv], ',');
            //Store name for each indiv
            //handle empty name
            if (temp_vec.size() > 1)
            {
                Indiv_name[nbr_deme - 1][indiv] = temp_vec[0];
            }
            else
            {
                Indiv_name[nbr_deme - 1][indiv] = "NaN";
            }
            //Store loci for each indiv
            auto locus_vec = gss::slice_by_char(*temp_vec.rbegin(), ' ');
            //locus level
            auto temp_indiv = std::vector<std::array<int, ploidy>>(locus_vec.size());
            for (int locus = 0; locus < locus_vec.size(); ++locus)
            {
                if (locus_vec[locus].size() > 0)
                {
                    temp_indiv[locus] = trim_locus(locus_vec[locus]);
                }
            }
            temp_deme[indiv] = temp_indiv;
        }
        Pop_name[nbr_deme - 1] = gss::slice_by_char(Indiv_name[nbr_deme - 1].back(), ' ');
        Genotype[nbr_deme - 1] = temp_deme;
    }

    if (nbr_dist_class == 0)
    {
        nbr_dist_class = Pop_name.size();
    }

    calc_dist_class_btw_deme(nbr_dist_class);

    //Handle chr dist

    //No genetic map
    if (path_to_chr_map_file.size() == 0)
    {
        Dist_btw_loc = std::vector<std::vector<std::vector<double>>>{};
        Nbr_chr_dist_class = 0;
        Dist_class_btw_loc = std::vector<std::vector<std::vector<int>>>{};
    }
    else
    {
        auto const map_file = gss::str_tolower(gss::read_file(path_to_chr_map_file));
        auto const map_file_vec = gss::slice_unix_windows_file_by_line(map_file);
        if (map_file_vec.size() != Locus_name.size())
        {
            throw std::logic_error("Number of locus != in genepop and .map file. I quit.");
        }

        std::vector<std::vector<std::string>> crude_map_vec(map_file_vec.size());
        for (int i = 0; i < crude_map_vec.size(); ++i)
        {
            crude_map_vec[i] = gss::slice_by_char(map_file_vec[i], '\t');
            if (crude_map_vec[i].size() != 4)
            {
                throw std::logic_error("In .map file, at line " + std::to_string(i + 1) + " 4 columns asked, only  " + std::to_string(crude_map_vec[i].size()) + " provide. I quit.");
            }
        }
        //Sort chr
        std::sort(crude_map_vec.begin(), crude_map_vec.end(),
                  [](auto &str_vec1, auto &str_vec2)
                  {
                      return std::stoi(str_vec1.at(0)) < std::stoi(str_vec2.at(0));
                  });

        //Num chr, name locus, dist in cM and pos in chr
        std::vector<std::tuple<int, std::string, double, int>> refine_map_vec(map_file_vec.size());

        int chr = 0;
        int prev_chr = std::stoi(crude_map_vec[0].at(0));
        for (int i = 0; i < crude_map_vec.size(); ++i)
        {
            int actual_chr = std::stoi(crude_map_vec[i].at(0));
            if (prev_chr != actual_chr)
            {
                ++chr;
            }
            prev_chr = actual_chr;
            std::get<0>(refine_map_vec[i]) = chr;
            std::get<1>(refine_map_vec[i]) = gss::remove_spaces_tab_underscores(crude_map_vec[i].at(1));
            std::get<2>(refine_map_vec[i]) = std::stod(crude_map_vec[i].at(2));
            std::get<3>(refine_map_vec[i]) = std::stoi(crude_map_vec[i].at(3));
        }

        calc_dist_btw_loc(refine_map_vec, chr, nbr_chr_dist_class);
    }
}

//TODO : Données manquante peuvent etre 00XX, 000XXX, 0 sinon erreur dans quelle deme, indiv, locus !
template <std::size_t ploidy>
std::array<int, ploidy> genepop_input_c<ploidy>::trim_locus(std::string str)
{
    std::array<int, ploidy> locus;

    //Diploid
    if (ploidy > 1)
    {
        if (str.size() > 3)
        {
            int middle = str.size() / 2;
            locus.at(0) = std::stoi(str.substr(0, middle));
            locus.at(1) = std::stoi(str.substr(middle, str.size() - 1));
        }
        else
        {
            locus.at(0) = std::stoi(str);
            locus.at(1) = 0;
        }
    }
    //Haploid
    else
    {
        locus.at(0) = std::stoi(str);
    }

    //TODO : A ameliorer
    if ((locus.at(0) < 0) || (locus.at(0) > 999))
    {
        throw std::logic_error("");
    }

    return locus;
}

template <std::size_t ploidy>
void genepop_input_c<ploidy>::calc_dist_class_btw_deme(int nbr_dist_class)
{
    Nbr_dist_class = nbr_dist_class;
    Dist_btw_deme = std::vector<std::vector<double>>(Pop_name.size(), std::vector<double>(Pop_name.size()));
    std::vector<std::array<double, 2>> deme_place(Pop_name.size());
    auto deme_place_itr = deme_place.begin();

    for (auto deme : Pop_name)
    {
        //keep the last coord for each deme
        try
        {
            deme_place_itr->at(0) = std::stof(deme[0]);
            deme_place_itr->at(1) = std::stof(deme[1]);
        }
        catch (const std::exception &e)
        {
        }
        ++deme_place_itr;
    }

    double max_dist = 0;
    for (auto deme1 = 0; deme1 < Dist_btw_deme.size(); ++deme1)
    {
        for (auto deme2 = 0; deme2 < Dist_btw_deme.size(); ++deme2)
        {
            //euclidien dist
            Dist_btw_deme[deme1][deme2] = sqrt(pow(deme_place[deme1].at(0) - deme_place[deme2].at(0), 2) + pow(deme_place[deme1].at(1) - deme_place[deme2].at(1), 2));
            if (Dist_btw_deme[deme1][deme2] > max_dist)
            {
                max_dist = Dist_btw_deme[deme1][deme2];
            }
        }
    }

    double dist_btw_class = max_dist / nbr_dist_class;
    Dist_class_btw_deme = std::vector<std::vector<int>>(Pop_name.size(), std::vector<int>(Pop_name.size()));
    {
        for (auto deme1 = 0; deme1 < Dist_btw_deme.size(); ++deme1)
        {
            for (auto deme2 = 0; deme2 < Dist_btw_deme.size(); ++deme2)
            {
                //if Dist_btw_deme  =  0 class 0 else floor(Dist_btw_deme/dist_btw_class) to be [class_limit_min, class_limit_max[
                //
                Dist_class_btw_deme[deme1][deme2] = (Dist_btw_deme[deme1][deme2] ? ceil(Dist_btw_deme[deme1][deme2] / dist_btw_class) - 1 : 0);
                //Handle trouble with value who are close to max_dist and who can be ceil at class + 1
                if (Dist_class_btw_deme[deme1][deme2] == nbr_dist_class)
                {
                    --Dist_class_btw_deme[deme1][deme2];
                }
            }
        }
    }
}

//WARNING : Use pos not cM
template <std::size_t ploidy>
void genepop_input_c<ploidy>::calc_dist_btw_loc(std::vector<std::tuple<int, std::string, double, int>> const &map_vec, int chr, int nbr_chr_dist_class)
{
    //number of chr and not max chr num
    chr += 1;
    //Organise assign correct locus index to locus in map
    std::vector<std::vector<int>> dist_from_chr_ori(chr);
    for (int locus = 0; locus < Locus_name.size(); ++locus)
    {
        std::string name = Locus_name[locus];
        bool not_found_name = true;
        int i = 0;
        while (not_found_name || (i < map_vec.size()))
        {
            if (name == std::get<1>(map_vec.at(i)))
            {
                not_found_name = false;
                auto dist = std::get<3>(map_vec.at(i));
                dist_from_chr_ori[std::get<0>(map_vec.at(i))].emplace_back(dist);
            }
            ++i;
        }
        if (not_found_name)
        {
            throw std::logic_error("Locus name : " + name + " in genepop file don't found in .map file. I quit.");
        }
    }

    //Compute Dist_btw_loc and max_dist
    Dist_btw_loc.resize(chr);
    double max_dist = 0;
    for (auto chr_num = 0; chr_num < chr; ++chr_num)
    {
        auto nbr_loc_per_chr = dist_from_chr_ori[chr_num].size();
        Dist_btw_loc[chr_num].resize(nbr_loc_per_chr);
        for (auto loc1 = 0; loc1 < nbr_loc_per_chr; ++loc1)
        {
            Dist_btw_loc[chr_num][loc1].resize(nbr_loc_per_chr);
            for (auto loc2 = 0; loc2 < nbr_loc_per_chr; ++loc2)
            {
                auto dist = abs(dist_from_chr_ori[chr_num][loc1] - dist_from_chr_ori[chr_num][loc2]);
                Dist_btw_loc[chr_num][loc1][loc2] = dist;
                if (dist > max_dist)
                {
                    max_dist = dist;
                }
            }
        }
    }

    //Compute Dist_class_btw_loc
    if (nbr_chr_dist_class > 0)
    {
        Nbr_chr_dist_class = nbr_chr_dist_class;
        Dist_class_btw_loc.resize(chr);
        double class_length = max_dist / Nbr_chr_dist_class;

        for (auto chr_num = 0; chr_num < chr; ++chr_num)
        {
            auto nbr_loc_per_chr = dist_from_chr_ori[chr_num].size();
            Dist_class_btw_loc[chr_num].resize(nbr_loc_per_chr);
            for (auto loc1 = 0; loc1 < nbr_loc_per_chr; ++loc1)
            {
                Dist_class_btw_loc[chr_num][loc1].resize(nbr_loc_per_chr);
                for (auto loc2 = 0; loc2 < nbr_loc_per_chr; ++loc2)
                {
                    //if Dist_btw_loc  =  0 class 0 else floor(Dist_btw_loc/class_length) to be [class_limit_min, class_limit_max[
                    Dist_class_btw_loc[chr_num][loc1][loc2] = (Dist_btw_loc[chr_num][loc1][loc2] ? ceil(Dist_btw_loc[chr_num][loc1][loc2] / class_length) - 1 : 0);
                    //Handle trouble with value who are close to max_dist and who can be ceil at class + 1
                    if (Dist_class_btw_loc[chr_num][loc1][loc2] == Nbr_chr_dist_class)
                    {
                        --Dist_class_btw_loc[chr_num][loc1][loc2];
                    }
                }
            }
        }
    }
    else
    {
        Nbr_chr_dist_class = 0;
    }
}