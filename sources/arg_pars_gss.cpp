#include <algorithm>
#include <iostream>

#include "arg_pars_gss.hpp"

//TRIMMING

std::string gss::remove_underscores(std::string str)
{
    auto underscore = std::find(str.begin(), str.end(), '_');

    while (underscore != str.end())
    {
        str.erase(underscore);
        underscore = std::find(str.begin(), str.end(), '_');
    }

    return str;
}

// Remove if char is space, tabulation or "_"
std::string gss::remove_spaces_tab_underscores(std::string str)
{
    str = gss::remove_spaces_tab_in_range(str, 0, str.size());

    return gss::remove_underscores(str);
}

std::string gss::remove_spaces_tab_in_range(std::string const &str, std::size_t pos_beg, std::size_t pos_end)
{
    std::string result;
    result.reserve(str.size());
    result += str.substr(0, pos_beg);

    std::size_t beg = pos_beg;
    std::size_t pos = beg;

    while (pos < str.size() && pos < pos_end)
    {
        while ((str[pos] == '\t') || (str[pos] == ' '))
        {
            ++pos;
        }

        beg = pos;

        while (((str[pos] != '\t') && (str[pos] != ' ')) && (pos < str.size()))
        {
            ++pos;
        }
        //Handle sep as last element before end of line
        if ((str[beg] != '\n') && (str[beg] != '\r') && (beg != pos))
        {
            result += str.substr(beg, pos - beg);
        }
    }

    result.shrink_to_fit();
    return result;
}

std::string gss::remove_comma_semicolons_in_range(std::string str, int pos_beg, int pos_end)
{
    auto beg_itr = str.begin();
    auto space = std::find(beg_itr + pos_beg, beg_itr + pos_end, ',');
    auto tab = std::find(beg_itr + pos_beg, beg_itr + pos_end, ';');
    if (space > tab)
    {
        space = tab;
    }

    while ((space != str.end()) && (pos_end > 0))
    {
        str.erase(space);
        --pos_end;
        space = std::find(space, beg_itr + pos_end, ',');
        tab = std::find(beg_itr + pos_beg, beg_itr + pos_end, ';');
        if (space > tab)
        {
            space = tab;
        }
    }

    return str;
}

std::string gss::replace_spaces_tab_in_range_by_comma(std::string str, int pos_beg, int pos_end)
{
    auto beg_itr = str.begin();
    auto space = std::find(beg_itr + pos_beg, beg_itr + pos_end, ' ');
    auto tab = std::find(beg_itr + pos_beg, beg_itr + pos_end, '\t');
    if (space > tab)
    {
        space = tab;
    }

    while ((space != str.end()) && (pos_end > 0))
    {
        str.replace(space, space + 1, ",");
        //--pos_end;
        space = std::find(space, beg_itr + pos_end, ' ');
        tab = std::find(beg_itr + pos_beg, beg_itr + pos_end, '\t');
        if (space > tab)
        {
            space = tab;
        }
    }

    return str;
}

std::string gss::replace_spaces_tab_underscores_by_comma(std::string str)
{
    str = gss::replace_spaces_tab_in_range_by_comma(str, 0, static_cast<int>(str.size()));

    auto underscore = std::find(str.begin(), str.end(), '_');

    while (underscore != str.end())
    {
        str.replace(underscore, underscore + 1, ",");
        underscore = std::find(str.begin(), str.end(), '_');
    }

    return str;
}

std::string gss::str_tolower(std::string str)
{
    for (auto &chara : str)
    {
        chara = std::tolower(chara);
    }

    return str;
}

std::vector<std::string> gss::slice_by_char(std::string const &str, char sep)
{
    std::vector<std::string> result;
    result.reserve(str.size());

    std::size_t beg = 0;
    std::size_t pos = beg;

    while (pos < str.size())
    {
        while (str[pos] == sep)
        {
            ++pos;
        }

        beg = pos;

        while ((str[pos] != sep) && (pos < str.size()))
        {
            ++pos;
        }
        //Handle sep as last element before end of line
        if ((str[beg] != '\n') && (str[beg] != '\r') && (beg != pos))
        {
            auto temp = str.substr(beg, pos - beg);
            result.push_back(temp);
        }
    }

    result.shrink_to_fit();
    return result;
}

std::vector<std::string> gss::slice_by_tab_or_space(std::string const &str)
{
    std::vector<std::string> result;
    result.reserve(str.size());

    std::size_t beg = 0;
    std::size_t pos = beg;

    while (pos < str.size())
    {
        while ((str[pos] == '\t') || (str[pos] == ' '))
        {
            ++pos;
        }

        beg = pos;

        while (((str[pos] != '\t') && (str[pos] != ' ')) && (pos < str.size()))
        {
            ++pos;
        }
        //Handle sep as last element before end of line
        if ((str[beg] != '\n') && (str[beg] != '\r') && (beg != pos))
        {
            auto temp = str.substr(beg, pos - beg);
            result.push_back(temp);
        }
    }

    result.shrink_to_fit();
    return result;
}

std::vector<std::string> gss::slice_by_comma_semicolon(std::string const &str)
{
    std::vector<std::string> result;
    result.reserve(str.size());

    std::size_t beg = 0;
    std::size_t pos = beg;

    while (pos < str.size())
    {
        while (str[pos] == ',' || str[pos] == ';')
        {
            ++pos;
        }

        beg = pos;

        while ((str[pos] != ',') && (str[pos] != ';') && (pos < str.size()))
        {
            ++pos;
        }
        result.push_back(str.substr(beg, pos - beg));
    }

    result.shrink_to_fit();
    return result;
}

std::vector<std::string> gss::slice_unix_windows_file_by_line(std::string const &str)
{
    std::vector<std::string> result;
    result.reserve(str.size());

    std::size_t beg = 0;
    std::size_t pos = 0;

    while (pos < str.size())
    {
        //Between line && delete space in line
        while ((pos < str.size()) && ((str[pos] == '\n') || (str[pos] == '\r') || (str[pos] == ' ')))
        {
            ++pos;
        }

        beg = pos;
        //Inside line
        while ((pos < str.size()) && (str[pos] != '\n') && (str[pos] != '\r'))
        {
            ++pos;
        }
        // Ignore comment line and empty line
        if ((str[beg] != '%') && (str[beg] != '#') && (str[beg] != '/') && (str[beg] != '@') && (pos != beg))
        {
            result.push_back(str.substr(beg, pos - beg));
        }
    }

    result.shrink_to_fit();
    return result;
}

std::vector<std::vector<std::string>> gss::slice_str_vec_by_string(std::vector<std::string> const &vec_str, std::string const &sep)
{
    std::vector<std::vector<std::string>> result;
    result.reserve(vec_str.size());

    std::size_t pos = 0;
    result.emplace_back();
    result[pos].reserve(vec_str.size());

    for (std::size_t i = 0; i < vec_str.size(); ++i)
    {
        if (vec_str[i] != sep)
        {
            result[pos].push_back(vec_str[i]);
        }
        else
        {
            result[pos].shrink_to_fit();
            ++pos;
            result.emplace_back();
            result[pos].reserve(vec_str.size());
        }
    }
    result.shrink_to_fit();
    return result;
}

//TODO : an be template
std::vector<int> gss::int_vector_parse_by_comma(std::string const &str)
{
    std::vector<std::string> temp_result = gss::slice_by_char(str, ',');

    std::vector<int> result(temp_result.size());
    auto result_itr = result.begin();
    for (auto const &value : temp_result)
    {
        *result_itr = std::stoi(value);
        ++result_itr;
    }
    return result;
}

std::vector<double> gss::double_vector_parse_by_comma(std::string const &str)
{
    std::vector<std::string> temp_result = gss::slice_by_char(str, ',');

    std::vector<double> result(temp_result.size());
    auto result_itr = result.begin();
    for (auto const &value : temp_result)
    {
        *result_itr = std::stod(value);
        ++result_itr;
    }
    return result;
}

std::vector<std::string> gss::str_vector_parse_by_comma(std::string const &str)
{
    std::vector<std::string> result = gss::slice_by_char(str, ',');

    return result;
}

std::vector<std::string> gss::str_vector_parse_by_comma_semicolon(std::string const &str)
{
    std::vector<std::string> result = gss::slice_by_comma_semicolon(str);

    return result;
}

bool gss::convert_str_bool(std::string const &str)
{

    if (gss::str_tolower(str) == "true" || gss::str_tolower(str) == "t" || gss::str_tolower(str) == "yes" || gss::str_tolower(str) == "y")
    {
        return true;
    }
    else
    {
        if (gss::str_tolower(str) == "false" || gss::str_tolower(str) == "f" || gss::str_tolower(str) == "no" || gss::str_tolower(str) == "n")
        {
            return false;
        }
        else
        {
            throw std::logic_error("( In convert_str_bool() : " + str + " isn't a boolean value. I exit. )");
        }
    }
}

bool gss::str_has_only_digits(std::string const &str)
{
    return str.find_first_not_of("0123456789") == std::string::npos;
}
