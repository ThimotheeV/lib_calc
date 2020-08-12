#include <algorithm>
#include <sstream>
#include <fstream>
#include <iterator>
#include <functional>
#include <iostream>
#include <cmath>

#include "input.hpp"

//TODO : A refactoriser
std::vector<std::string> sep_by_char(std::string const &str, char sep)
{
    std::vector<std::string> result;
    result.reserve(str.size());

    int beg = 0;
    int pos = 0;

    while (pos < str.size())
    {
        //Ignored char
        while ((pos < str.size()) && (str[pos] == sep))
        {
            ++pos;
        }

        beg = pos;

        while ((pos < str.size()) && (str[pos] != sep))
        {
            ++pos;
        }

        //Handle sep as caracter endline
        if (beg != pos)
        {
            result.push_back(str.substr(beg, pos - beg));
        }
    }

    result.shrink_to_fit();
    return result;
}

//Remove empty line
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
                std::cout << "File in WINDOWS format" << std::endl;
                adver = !adver;
            }
            ++pos;
        }

        beg = pos;

        while ((str[pos] != '\n') && (str[pos] != '\r') && (pos < str.size()))
        {
            ++pos;
        }
        if (beg != pos)
        {
            result.push_back(str.substr(beg, pos - beg));
        }
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