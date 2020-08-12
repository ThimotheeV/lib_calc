#include <iostream>
#include <iomanip>
#include <fstream>

#include "output_file.hpp"

//TODO : Contruire en modulaire à partir de la liste de stat donné en entrée
//=> Chaque module gère une stat (ou un ensemble de stat) particulier
//Comment faire pour garder le fichier en écriture ouvert.

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
