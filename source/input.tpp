//TODO : Handle empty line
//TODO : Vérifier que le nbr de locus est le même pour tout les indivs
//Constructor of genepop_input_c
template <std::size_t ploidy>
genepop_input_c<ploidy>::genepop_input_c(std::string path_to_file, int nbr_class)
{
    auto const file_str = lower_case(read_file(path_to_file));
    auto const file_str_vec = trim_unix_windows_file_by_line(file_str);
    //TODO : Utiliser un vecteur de séparateur
    std::vector<std::vector<std::string>> const pop_vec = trim_str_vec_by_string(file_str_vec, "pop");

    //Special case of header => first line + locus name => in pop_vec[0]
    Locus_name.reserve(pop_vec[0].size() - 1);
    Header = pop_vec[0][0];
    for (int locus = 1; locus < pop_vec[0].size(); ++locus)
    {
        auto const temp_str_vec = sep_by_char(pop_vec[0][locus], ',');
        //No need to copy header in memory if temp_str_vec haven't more than 1 string in it
        Locus_name.reserve(Locus_name.size() + temp_str_vec.size() - 1);
        Locus_name.insert(Locus_name.end(), temp_str_vec.begin(), temp_str_vec.end());
    }

    for (int locus = 0; locus < Locus_name.size(); ++locus)
    {
        Locus_name[locus] = remove_spaces_underscores(Locus_name[locus]);
    }

    //Trim pop by indiv and locus
    Genotype.resize(pop_vec.size() - 1);
    Indiv_name.resize(pop_vec.size() - 1);
    Pop_name.resize(pop_vec.size() - 1);

    //pop level, first was header
    for (int nbr_pop = 1; nbr_pop < pop_vec.size(); ++nbr_pop)
    {
        //indiv level
        std::vector<std::vector<std::array<int, ploidy>>> temp_pop;
        //pop_vec[nbr_pop].size() => nbr of indiv by pop
        temp_pop.resize(pop_vec[nbr_pop].size());
        Indiv_name[nbr_pop - 1].resize(pop_vec[nbr_pop].size());

        for (int indiv = 0; indiv < pop_vec[nbr_pop].size(); ++indiv)
        {
            auto const temp_vec = sep_by_char(pop_vec[nbr_pop][indiv], ',');
            //Store name for each indiv
            Indiv_name[nbr_pop - 1][indiv] = temp_vec[0];

            //Store loci for each indiv
            auto locus_vec = sep_by_char(temp_vec[1], ' ');

            //locus level
            auto temp_indiv = std::vector<std::array<int, ploidy>>(locus_vec.size());
            for (int locus = 0; locus < locus_vec.size(); ++locus)
            {
                if (locus_vec[locus].size() > 0)
                {
                    temp_indiv[locus] = trim_locus(locus_vec[locus]);
                }
            }
            temp_pop[indiv] = temp_indiv;
        }
        Pop_name[nbr_pop - 1] = sep_by_char(Indiv_name[nbr_pop - 1].back(), ' ');
        Genotype[nbr_pop - 1] = temp_pop;
    }

    if (nbr_class == 0)
    {
        nbr_class = Pop_name.size();
    }

    calc_dist_class_btw_pop(nbr_class);
}

//TODO : Données manquante peuvent etre 00XX, 000XXX, 0 sinon erreur dans quelle pop, indiv, locus !
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
void genepop_input_c<ploidy>::calc_dist_class_btw_pop(int nbr_class)
{
    Nbr_dist_class = nbr_class;
    Dist_btw_pop = std::vector<std::vector<double>>(Pop_name.size(), std::vector<double>(Pop_name.size()));
    std::vector<std::array<double, 2>> pop_place(Pop_name.size());
    auto pop_place_itr = pop_place.begin();

    for (auto pop : Pop_name)
    {
        //keep the last coord for each pop
        try
        {
            pop_place_itr->at(0) = std::stof(pop[0]);
            pop_place_itr->at(1) = std::stof(pop[1]);
        }
        catch (const std::exception &e)
        {
        }
        ++pop_place_itr;
    }

    double max_dist = 0;
    for (auto pop1 = 0; pop1 < Dist_btw_pop.size(); ++pop1)
    {
        for (auto pop2 = 0; pop2 < Dist_btw_pop.size(); ++pop2)
        {
            //euclidien dist
            Dist_btw_pop[pop1][pop2] = sqrt(pow(pop_place[pop1].at(0) - pop_place[pop2].at(0), 2) + pow(pop_place[pop1].at(1) - pop_place[pop2].at(1), 2));
            if (Dist_btw_pop[pop1][pop2] > max_dist)
            {
                max_dist = Dist_btw_pop[pop1][pop2];
            }
        }
    }

    double dist_btw_class = max_dist / nbr_class;
    Dist_class_btw_pop = std::vector<std::vector<int>>(Pop_name.size(), std::vector<int>(Pop_name.size()));
    {
        for (auto pop1 = 0; pop1 < Dist_btw_pop.size(); ++pop1)
        {
            for (auto pop2 = 0; pop2 < Dist_btw_pop.size(); ++pop2)
            {
                //if Dist_btw_pop  =  0 class 0 else floor(Dist_btw_pop/dist_btw_class) to be [class_limit_min, class_limit_max[
                //
                Dist_class_btw_pop[pop1][pop2] = (Dist_btw_pop[pop1][pop2] ? ceil(Dist_btw_pop[pop1][pop2] / dist_btw_class) - 1 : 0);
                //Handle trouble with value who are close to max_dist and who can be ceil at class + 1
                if (Dist_class_btw_pop[pop1][pop2] == nbr_class)
                {
                    --Dist_class_btw_pop[pop1][pop2];
                }
            }
        }
    }
}