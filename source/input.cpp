#include <algorithm>
#include <sstream>
#include <fstream>
#include <iterator>
#include <functional>
#include <iostream>
#include <cmath>

#include "input.hpp"

//TODO : A refactoriser
std::vector<std::string> trim_by_char(std::string const &str, char sep)
{
    std::vector<std::string> result;
    result.reserve(str.size());

    int beg = 0;
    int pos = beg;

    while (pos < str.size())
    {
        while (str[pos] == sep)
        {
            ++pos;
        }

        beg = pos;

        while (str[pos] != sep)
        {
            ++pos;
        }

        result.push_back(str.substr(beg, pos - beg));
    }

    result.shrink_to_fit();
    return result;
}

//TODO : A refactoriser
std::vector<std::string> trim_unix_windows_file_by_line(std::string const &str)
{
    std::vector<std::string> result;
    result.reserve(str.size());

    int beg = 0;
    int pos = beg;
    bool adver = 1;

    while (pos < str.size())
    {
        while ((str[pos] == '\n') || (str[pos] == '\r'))
        {
            if (adver && (str[pos] == '\r'))
            {
                std::cerr << "File in WINDOWS format" << std::endl;
                adver = !adver;
            }
            ++pos;
        }

        beg = pos;

        while ((str[pos] != '\n') && (str[pos] != '\r'))
        {
            ++pos;
        }

        result.push_back(str.substr(beg, pos - beg));
    }

    return result;
}

std::vector<std::vector<std::string>> trim_str_vec_by_string(std::vector<std::string> const &vec_str, std::string sep)
{
    std::vector<std::vector<std::string>> result;
    result.reserve(vec_str.size());

    int pos = 0;
    result.push_back(std::vector<std::string>());
    result[pos].reserve(vec_str.size());

    for (int i = 0; i < vec_str.size(); ++i)
    {
        if (vec_str[i] != sep)
        {
            result[pos].push_back(vec_str[i]);
        }
        else
        {
            result[pos].shrink_to_fit();
            ++pos;
            result.push_back(std::vector<std::string>());
            result[pos].reserve(vec_str.size());
        }
    }
    result.shrink_to_fit();
    return result;
}

std::string remove_spaces_underscores(std::string str)
{
    str = remove_spaces_in_range(str, 0, str.size());

    auto underscore = std::find(str.begin(), str.end(), '_');

    while (underscore != str.end())
    {
        str.erase(underscore);
        underscore = std::find(str.begin(), str.end(), '_');
    }

    return str;
}

std::string remove_spaces_in_range(std::string str, int pos_beg, int pos_end)
{
    auto beg_itr = str.begin();
    auto space = std::find(beg_itr + pos_beg, beg_itr + pos_end, ' ');

    while ((space != str.end()) && (pos_end > 0))
    {
        str.erase(space);
        --pos_end;
        space = std::find(space, beg_itr + pos_end, ' ');
    }

    return str;
}

std::string lower_case(std::string str)
{
    for (auto &chara : str)
    {
        chara = std::tolower(chara);
    }

    return str;
}

//TODO : fonction qui permet d'aller cherhcer le prochain caractère interessant
//string.find()

std::string const read_file(std::string const &filename)
{
    //Put filename in a stream
    std::ifstream file_str(filename);
    if (!file_str)
    {
        //TODO : mettre le nom du fichier en cerr !
        throw std::logic_error("File not open");
    }
    else
    {
        //Don't skip space
        file_str >> std::noskipws;
        //Convert the stream in a string
        const std::string &fstr{std::istream_iterator<char>{file_str}, {}};
        //Test if the stream is good
        if (file_str.bad())
        {
            throw std::logic_error("File io error");
        }
        return fstr;
    }
}

//TODO : Handle empty line
//TODO : A adapter pour l'antislash windows
//TODO : Vérifier que le nbr de locus est le même pour tout les indivs
//Constructor of genepop_input_c
template <std::size_t ploidy>
genepop_input_c<ploidy>::genepop_input_c(std::string path_to_file)
{
    auto const file_str = lower_case(read_file(path_to_file));
    auto const file_str_vec = trim_unix_windows_file_by_line(file_str);
    //TODO : Utiliser un vecteur de séparateur
    std::vector<std::vector<std::string>> const pop_vec = trim_str_vec_by_string(file_str_vec, "pop");

    // std::cout<<pop_vec[0].size()<<std::endl;

    //Special case of header => first line + locus name => in pop_vec[0]
    Locus_name.reserve(pop_vec[0].size()-1);
    Header = pop_vec[0][0];
    for (int locus = 1; locus < pop_vec[0].size(); ++locus)
    {
        auto const temp_str_vec = trim_by_char(pop_vec[0][locus], ',');
        //No need to copy header in memory if temp_str_vec haven't more than 1 string in it
        Locus_name.reserve(Locus_name.size() + temp_str_vec.size() - 1);
        Locus_name.insert(Locus_name.end(), temp_str_vec.begin(), temp_str_vec.end());
    }

    for (int locus = 1; locus < Locus_name.size(); ++locus)
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
        temp_pop.reserve(pop_vec[nbr_pop].size());
        Indiv_name[nbr_pop - 1].reserve(pop_vec[nbr_pop].size());

        for (auto const &indiv : pop_vec[nbr_pop])
        {
            auto const temp_vec = trim_by_char(indiv, ',');
            //Remove space before and after text
            auto tmp_vec = remove_spaces_in_range(temp_vec[0], 0, temp_vec[0].find_first_not_of(' '));
            tmp_vec = remove_spaces_in_range(tmp_vec, tmp_vec.find_last_not_of(' '), tmp_vec.size());
            Indiv_name[nbr_pop - 1].push_back(tmp_vec);

            tmp_vec = remove_spaces_in_range(temp_vec[1], 0, temp_vec[1].find_first_not_of(' '));
            auto const locus_vec = trim_by_char(tmp_vec, ' ');
            //locus level
            if (locus_vec.size() != Locus_name.size())
            {
                std::string str = "Pop " + std::to_string(nbr_pop - 1) + " have " + std::to_string(locus_vec.size()) + " locus but file " + path_to_file + " describ " + std::to_string(Locus_name.size()) + " locus.";
                throw std::logic_error(str);
            }
            auto temp_indiv = std::vector<std::array<int, ploidy>>(locus_vec.size());
            for (int locus = 0; locus < locus_vec.size(); ++locus)
            {
                temp_indiv[locus] = trim_locus(locus_vec[locus]);
            }
            temp_pop.emplace_back(std::move(temp_indiv));
        }
        Pop_name[nbr_pop - 1] = trim_by_char(*(Indiv_name[nbr_pop - 1].end() - 1), ' ');
        Genotype[nbr_pop - 1] = std::move(temp_pop);
    }

    calc_dist_btw_pop();
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
void genepop_input_c<ploidy>::calc_dist_btw_pop()
{
    Dist_btw_pop = std::vector<std::vector<float>>(Pop_name.size(), std::vector<float>(Pop_name.size()));
    std::vector<std::array<int, 2>> pop_place(Pop_name.size());
    auto pop_place_itr = pop_place.begin();

    for (auto pop : Pop_name)
    {
        try
        {
            pop_place_itr->at(0) = std::stoi(pop[0]);
            pop_place_itr->at(1) = std::stoi(pop[1]);
        }
        catch (const std::exception &e)
        {
            std::cerr << pop[0] << " is not a numeric value, can't be use to calculate distance between population. Reset matrix of distance between pop." << std::endl;
            Dist_btw_pop = std::vector<std::vector<float>>{};
            return;
        }

        ++pop_place_itr;
    }

    for (auto pop1 = 0; pop1 < Dist_btw_pop.size(); ++pop1)
    {
        for (auto pop2 = 0; pop2 < Dist_btw_pop.size(); ++pop2)
        {
            //euclidien dist
            Dist_btw_pop[pop1][pop2] = sqrt(pow(pop_place[pop1].at(0) - pop_place[pop2].at(0), 2) + pow(pop_place[pop1].at(1) - pop_place[pop2].at(1), 2));
        }
    }
}
// explicit instantiations (need for separate template declaration on .hpp and .cpp)
template struct genepop_input_c<1>;
template struct genepop_input_c<2>;