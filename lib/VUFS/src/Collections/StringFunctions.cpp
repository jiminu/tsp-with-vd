#include "StringFunctions.h"
#include <cmath>
#include <cstdio>
using namespace std;

string StringFunctions::subString( const string& aString, const unsigned int& pos, const int& length )
{
    if ( pos < aString.length() ) {
        return aString.substr( pos, length );
    }
    else {
        return string();
    }
}



string StringFunctions::subString( const string& aString, int& pos )
{
    int startPosition = pos;
    int endPosition   = 0;
    

    StringFunctions::updateStartPositionWithoutWhiteSpace( aString, startPosition );
    endPosition = aString.find_first_of( STR_SEPS, startPosition );
    pos = endPosition;
    
    return StringFunctions::subString( aString, startPosition, endPosition-startPosition ).c_str();
}



string StringFunctions::strToken( const string& aString, int& pos )
{
    int startPosition = pos;
    int endPosition;

    int sourceLength  = aString.length();
    for ( ; startPosition<sourceLength; startPosition++ )  {
        if ( aString[startPosition] != ' ' && aString[startPosition] != '\t' )
            break;
    }

    endPosition = aString.find_first_of( STR_SEPS, startPosition );
    pos = endPosition;
    
    return aString.substr( startPosition, endPosition-startPosition );
}



string StringFunctions::strToLower( const string& aString )
{
    int lengthOfString = aString.length();

	string tempString = aString;
	
    for( int i=0; i<lengthOfString; i++ ) {
        tempString[i] = tolower(tempString[i]);
    }

    return tempString;    
}



string StringFunctions::strToUpper( const string& aString )
{
    int lengthOfString = aString.length();

	string tempString = aString;
	
    for( int i=0; i<lengthOfString; i++ ) {
        tempString[i] = toupper(tempString[i]);
    }

    return tempString;    
}



string StringFunctions::strTrim( const string& aString )
{
    return strTrimLeft( strTrimRight( aString ) );
}



string StringFunctions::strTrimLeft( const string& aString )
{
    string tempStr = aString;

    int n = aString.find_first_not_of(" \n\f\t\v\r");
    tempStr.replace(0, n, "");
    
    return tempStr;
}



string StringFunctions::strTrimRight( const string& aString )
{
    string tempStr = aString;
    int n = aString.find_last_not_of(" \n\f\t\v\r");   
    tempStr.replace(n+1, aString.length()-n, "");

    return tempStr;
}



bool StringFunctions::isStringLetterInAlphabet( const string& aStringLetter )
{
    string tempStr = strToUpper( aStringLetter.substr(0, 1) );

    if( (int)tempStr[0] >= 'A' && (int)tempStr[0] <= 'Z' ) {
        return true;
    }
    else
        return false;
}



bool StringFunctions::isStringLetterInNumber( const string& aStringLetter )
{
    string tempStr = strToUpper( aStringLetter.substr(0, 1) );

    if( (int)tempStr[0] >= '0' && (int)tempStr[0] <= '9' ) {
        return true;
    }
    else
        return false;
}



string StringFunctions::getExtensionFromFileName( const string& fullFileName )
{
    int posOfExtEnd   = fullFileName.length();
    int posOfExtStart = fullFileName.rfind( ".", posOfExtEnd )+1;
    int lengthOfExt   = posOfExtEnd - posOfExtStart;

    return fullFileName.substr( posOfExtStart, lengthOfExt );
}



string StringFunctions::getFileNameWithoutPathAndExtension( const string& fullFileName )
{
    int posOfExtEnd   = fullFileName.length();
    int posOfExtStart = fullFileName.rfind( ".", posOfExtEnd )+1;
    int lengthOfExt   = posOfExtEnd - posOfExtStart;

    string nameOfExt  = fullFileName.substr( posOfExtStart, lengthOfExt );

    #if defined(WIN32) || defined(WIN64) || defined(_WIN32) || defined(_WIN64)
		int posOfFileNameStart         = fullFileName.rfind( "\\", posOfExtEnd )+1;   // FOR WINDOWS SYSTEM
	#else
		int posOfFileNameStart         = fullFileName.rfind( "/", posOfExtEnd )+1;	  // FOR LINUX SYSTEM
	#endif

    int legnthOfFileNameWithoutExt = posOfExtStart - posOfFileNameStart-1;
    
    return fullFileName.substr( posOfFileNameStart, legnthOfFileNameWithoutExt );
}



string StringFunctions::getFileNameWithoutPath( const string& fullFileName )
{
    int posOfExtEnd   = fullFileName.length();
    int posOfExtStart = fullFileName.rfind( ".", posOfExtEnd )+1;
    int lengthOfExt   = posOfExtEnd - posOfExtStart;

    string nameOfExt  = fullFileName.substr( posOfExtStart, lengthOfExt );

    #if defined(WIN32) || defined(WIN64) || defined(_WIN32) || defined(_WIN64)
		int posOfFileNameStart         = fullFileName.rfind( "\\", posOfExtEnd )+1;   // FOR WINDOWS SYSTEM
	#else
		int posOfFileNameStart         = fullFileName.rfind( "/", posOfExtEnd )+1;	  // FOR LINUX SYSTEM
	#endif

    int legnthOfFileNameWithExt    = posOfExtEnd - posOfFileNameStart;
    
    return fullFileName.substr( posOfFileNameStart, legnthOfFileNameWithExt );
}



string StringFunctions::getPathFromFileName( const string& fullFileName )
{
    int posOfExtEnd   = fullFileName.length();

    #if defined(WIN32) || defined(WIN64) || defined(_WIN32) || defined(_WIN64)
        int posOfSep  = fullFileName.rfind( "\\", posOfExtEnd );   // FOR WINDOWS SYSTEM
    #else
        int posOfSep  = fullFileName.rfind( "/", posOfExtEnd );	  // FOR LINUX SYSTEM
    #endif

    string pathOfFile = fullFileName.substr(0, posOfSep + 1);

    return pathOfFile;
}



string StringFunctions::getBlankStringForFixedField( const int& fieldSize, const string& targetRec )
{
    int    numOfSpaceOnDemand = fieldSize - targetRec.length();
    string blankString;
    
    int i_space = 0;
    for ( i_space=0; i_space<numOfSpaceOnDemand; i_space++ ) {
        blankString = blankString + " ";
    }
    
    return blankString;
}



string StringFunctions::convertIntegerToString( const int& targetInt )
{
    char   tempBuff[32];
    string tempStr;

    // sprintf( tempBuff, "%d", targetInt );
#if defined(WIN32) || defined(WIN64) || defined(_WIN32) || defined(_WIN64)
	sprintf_s( tempBuff, "%d", targetInt );
#else
	sprintf( tempBuff, "%d", targetInt );
#endif
    tempStr = tempBuff;

    return tempStr;
}



string StringFunctions::convertDoubleToString( const double& targetDouble )
{
    char   tempBuff[32];
    string tempStr;

    //sprintf( tempBuff, "%f", targetDouble );
#if defined(WIN32) || defined(WIN64) || defined(_WIN32) || defined(_WIN64)
    sprintf_s(tempBuff, "%f", targetDouble);
#else
    sprintf(tempBuff, "%f", targetDouble);
#endif
    tempStr = tempBuff;

    return tempStr;    
}



string StringFunctions::convertDoubleToString( const double& targetDouble, const int& decimalSize )
{
    char   tempBuff[32];
    string tempStr;

    double roundoff;
    if (decimalSize > 0) {
        roundoff = (targetDouble >= 0.0) ? 0.5: -0.5;
        for (int i = 0; i < decimalSize; ++i) {
            roundoff /= 10.;
        }
    }

    //sprintf_s(tempBuff, "%f", targetDouble + roundoff);
#if defined(WIN32) || defined(WIN64) || defined(_WIN32) || defined(_WIN64)
    sprintf_s( tempBuff, "%f", targetDouble + roundoff);
#else
    sprintf(tempBuff, "%f", targetDouble + roundoff);
#endif

    tempStr = tempBuff;

    int posOfDot = tempStr.find( ".", 0 );

    tempStr = tempStr.substr( 0, posOfDot+decimalSize+1 );

    return tempStr;
}



string StringFunctions::convertIntegerToString(const int& targetInt, const int& size)
{
    char   tempBuff[32];
    string targetString;

    //sprintf(tempBuff, "%d", targetInt);
#if defined(WIN32) || defined(WIN64) || defined(_WIN32) || defined(_WIN64)
    sprintf_s(tempBuff, "%d", targetInt);
#else
    sprintf(tempBuff, "%d", targetInt);
#endif

    targetString = tempBuff;

    string blankSpace = StringFunctions::getBlankStringForFixedField(size, targetString);

    return (blankSpace + targetString);
}



string StringFunctions::convertDoubleToString(const double& targetDouble, const int& size, const int& precision)
{
    char   tempBuff[32];
    string targetString;

    int digit = size - (precision + 1);
    double digitOfTarget = log10(abs(targetDouble));

    string  sprintfField("%.");

    if (targetDouble == 0.0) {
#if defined(WIN32) || defined(WIN64) || defined(_WIN32) || defined(_WIN64)
        sprintf_s(tempBuff, "%d", precision);
#else
        sprintf(tempBuff, "%d", precision);
#endif
        sprintfField += tempBuff;
        sprintfField += "f";
    }
    else {
        if (digitOfTarget > digit || digitOfTarget + precision < 0) {
            int sizeEFormat = 7; 
            if (precision + sizeEFormat < size) {
#if defined(WIN32) || defined(WIN64) || defined(_WIN32) || defined(_WIN64)
                sprintf_s(tempBuff, "%d", precision);
#else
                sprintf(tempBuff, "%d", precision);
#endif
            }
            else {
#if defined(WIN32) || defined(WIN64) || defined(_WIN32) || defined(_WIN64)
                sprintf_s(tempBuff, "%d", size - sizeEFormat);
#else
                sprintf(tempBuff, "%d", size - sizeEFormat);
#endif
            }
            sprintfField += tempBuff;
            sprintfField += "E";

        }
        else {
#if defined(WIN32) || defined(WIN64) || defined(_WIN32) || defined(_WIN64)
            sprintf_s(tempBuff, "%d", precision);
#else
            sprintf(tempBuff, "%d", precision);
#endif
            sprintfField += tempBuff;
            sprintfField += "f";
        }
    }
#if defined(WIN32) || defined(WIN64) || defined(_WIN32) || defined(_WIN64)
    sprintf_s(tempBuff, sprintfField.c_str(), targetDouble);
#else
    sprintf(tempBuff, sprintfField.c_str(), targetDouble);
#endif
    targetString = tempBuff;

    string blankSpace = StringFunctions::getBlankStringForFixedField(size, targetString);

    return (blankSpace + targetString);
}



void StringFunctions::updateStartPositionWithoutWhiteSpace( const string& strLine, int& startPosition )
{
    int sourceLength = strLine.length();
    for ( ; startPosition<sourceLength; startPosition++ )  {
        if ( strLine[startPosition] != ' ' && strLine[startPosition] != '\t' )
            break;
    }
}



int StringFunctions::getPositionStartWithNumber( const string& strLine )
{
    int sourceLength = strLine.length();
    int i_strLine    = 0;

    for ( ; i_strLine<sourceLength; i_strLine++ )  {

        if( (int)strLine[i_strLine] >= '0' && (int)strLine[i_strLine] <= '9' )
            break;
    }
    
    return i_strLine;
}

int StringFunctions::tokens_from_string(const string& inputString, const string& delimiters, vector<string>& tokenVector)
{
    std::size_t prev = 0, pos;
    while ((pos = inputString.find_first_of(delimiters, prev)) != std::string::npos)
    {
        if (pos > prev)
            tokenVector.push_back(inputString.substr(prev, pos - prev));
        prev = pos + 1;
    }
    if (prev < inputString.length())
        tokenVector.push_back(inputString.substr(prev, std::string::npos));

    return tokenVector.size();
}


//  by Y. Cho

void     StringFunctions::pop_space_at_front(string& str )
{
    string::size_type off   = 0;
    string::size_type count = str.size();

    for ( string::const_iterator i_pos=str.begin(); i_pos!=str.end(); ++i_pos ) {
        if ( *i_pos == ' ' ) {
            ++off;
            --count;
        }
        else {
            break;
        }
    }

    str = str.substr(off, count);
}


void     StringFunctions::pop_white_space_at_front(string& str)
{
    string::size_type off = 0;
    string::size_type count = str.size();

    for (string::const_iterator i_pos = str.begin(); i_pos != str.end(); ++i_pos) {
        if (*i_pos == ' ' || *i_pos == '\t' || *i_pos == '\n') {
            ++off;
            --count;
        }
        else {
            break;
        }
    }

    str = str.substr(off, count);
}


string     StringFunctions::pop_alphabetic_char_at_front( string& str )
{
    pop_white_space_at_front(str);

    string::size_type count = 0;
    for ( string::const_iterator i_pos=str.begin(); i_pos!=str.end(); ++i_pos ) {
        if ( isalpha( *i_pos ) || *i_pos == ' ' ) {
            ++count;
        }
        else {
            break;
        }
    }

    string alphabeticChars = str.substr( 0, count );
    str = str.substr( count, str.size()-count );

    return alphabeticChars;
}


string     StringFunctions::pop_numeric_char_at_front(    string& str )
{
    if (str.empty()) {
        return string();
    }

    pop_white_space_at_front(str);

    string::size_type count = 0;
    string::const_iterator i_pos=str.begin();
    if ( isdigit( *i_pos ) || *i_pos == '-' ) {
        ++count;
        ++i_pos;

        for ( ; i_pos!=str.end(); ++i_pos ) {
            if ( isdigit( *i_pos ) || *i_pos == '.' ) {
                ++count;
            }
            else {
                break;
            }
        }
    }

    string numericChars = str.substr( 0, count );
    str = str.substr( count, str.size()-count );

    return numericChars;
}


