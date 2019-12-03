#include <iostream>
#include <iomanip>
#include <fstream>

#include "output_file.hpp"

output_file_c::output_file_c(std::string const &directory)
{
    Directory = directory;

    if (Directory.empty())
    {
        return;
    }

    std::ofstream test(Directory + "/test.txt");
    if (test.is_open())
    {
        test << "test\n";
        test.close();
    }
    else
    {
        throw std::invalid_argument("Unable to open test.txt");
    }
}
