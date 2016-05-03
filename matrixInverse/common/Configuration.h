//Created by Erman Gurses on 2/8/16.
//Modified by Rutuja Patil on 3/28/2016.
/*****************************************************************************
 * Configuration.h header file
 *
 * See top of Configuration.cpp for the default command line parameters
 * that this class can handle.  We could just expose the addParamInt, etc.
 * if we wanted to enable users to extend the command-line options that can
 * be handled.
 *
 * Assumptions:
 *      Max number of options is 50.
 *
 * Usage:
 *      // Parse the command-line inputs.
 *      Configuration config(argc,argv);
 *
 *      // Get any particular command-line input.
 *      int bx = config.getInt("bx");
 *
 *      // Generate an LDAP string to output configuration.
 *      string outString = config.toLDAPString();
 *
 ****************************************************************************/
#ifndef CONFIGURATION_H
#define CONFIGURATION_H

#define MAX_NUM_OPTIONS 50

#include <string>
#include <map>
#include <set>
#include <list>
#include <getopt.h>

class Configuration
{
  public:
    Configuration();
    ~Configuration();
    int getInt(std::string fieldname);
    bool getBool(std::string fieldname);
    std::string getString(std::string fieldname);
    void parse(int argc, char* argv[]);
    std::string toLDAPString();

  protected:
  
    void addParamInt(std::string fieldname, char single_char,
                     int init_val, std::string help_str);
    void addParamBool(std::string fieldname, char single_char,
                             bool init_val, std::string help_str);    
    void processArgVal(std::string fieldname, char * arg_val);

    void addParamString(std::string fieldname,char single_char,
                   std::string init_val,std::string help_str);
     int parseInt( char* string );

  private:

    std::string mExecName;
    std::string mUsageString;
    std::string mCommandChars;          // For use with long_options 
    int mCount; // count of objects
    std::list<int> mCharOrder;
    std::map<int,std::string> mCharToFieldname;
    std::set<int> mKnownChars;
    std::map<std::string,int> mIntMap;
    std::map<std::string,bool> mBoolMap;
    std::map<std::string,std::string> mStringMap;

    struct option mLongOptions[MAX_NUM_OPTIONS];

};
#endif         
