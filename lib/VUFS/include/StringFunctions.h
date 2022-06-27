#ifndef _STRINGFUNCTIONS_H_
#define _STRINGFUNCTIONS_H_

#include <string>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <vector>

using namespace std;


const string STR_SEPS = " \t\r\n";


class StringFunctions
{
public:

    static string     subString( const string& aString, const unsigned int& pos, const int& length );
    static string     subString( const string& aString, int& pos );
    static string     strToken( const string& aString, int& pos );
    
    static string     strToLower( const string& aString );
    static string     strToUpper( const string& aString );

    static string     strTrim( const string& aString );
    static string     strTrimLeft( const string& aString );
    static string     strTrimRight( const string& aString );

    static bool       isStringLetterInAlphabet( const string& aStringLetter );
    static bool       isStringLetterInNumber( const string& aStringLetter );

    static string     getExtensionFromFileName( const string& fullFileName );
    static string     getFileNameWithoutPathAndExtension( const string& fullFileName );
    static string     getFileNameWithoutPath( const string& fullFileName );
    static string     getPathFromFileName( const string& fullFileName );

    static string     getBlankStringForFixedField( const int& fieldSize, const string& targetRec );

    static string     convertIntegerToString( const int& targetInt );
    static string     convertDoubleToString( const double& targetDouble );
    static string     convertDoubleToString( const double& targetDouble, const int& decimalSize );

    static string     convertIntegerToString(const int& targetInt, const int& size);
    static string     convertDoubleToString(const double& targetDouble, const int& size, const int& precision);

    static void       updateStartPositionWithoutWhiteSpace( const string& strLine, int& startPosition );
    static int        getPositionStartWithNumber( const string& strLine );

    // August 24, 2020 by Joonghyun
    static int        tokens_from_string(const string& inputString, const string& delimiters, vector<string>& tokenVector);

    //  by Y. Cho
    static void       pop_space_at_front(           string& str );
    static void       pop_white_space_at_front(     string& str );
    static string     pop_alphabetic_char_at_front( string& str );
    static string     pop_numeric_char_at_front(    string& str );
};

#endif
