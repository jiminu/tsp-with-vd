#ifndef _FILEOPERATIONS_H_
#define _FILEOPERATIONS_H_

#include <string>
#include <fstream>
using namespace std;

class FileOperations
{
public:
    static int getFileSize(const string& filename);
    static void remove_both_comment_and_blank_lines_from_ascii_file(const std::string& filename, const char& characterToRemoveLine, string& asciiFileInString);
};

#endif


