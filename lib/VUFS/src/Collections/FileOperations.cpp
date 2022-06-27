#include "FileOperations.h"
#include <iostream>


int FileOperations::getFileSize(const string& filename)
{
    ifstream fin;
    fin.open(filename.c_str(), ios_base::binary | ios_base::in);

    if (!fin.good() || fin.eof() || !fin.is_open()) { 
        return 0; 
    }

    fin.seekg(0, ios_base::beg);
    ifstream::pos_type begin_pos = fin.tellg();
    fin.seekg(0, ios_base::end);

    return static_cast<int>(fin.tellg() - begin_pos);
}

void FileOperations::remove_both_comment_and_blank_lines_from_ascii_file(const std::string& filename, const char& characterToRemoveLine, string& asciiFileInString)
{
    ifstream fin(filename.c_str());
    if (!fin) 
    {
        std::cerr << "fail to open: " << filename << "." << endl;
        return;
    }

    std::string line;
    const char LF('\n');
    while (std::getline(fin, line, LF))
    {
        if(!line.empty() && line.front() != characterToRemoveLine)
        {
            line.push_back(LF); // std::getline() function discards LF
            asciiFileInString.append(line);
        }
    }
}
