//Created by Erman Gurses on 2/8/16.
//Modified by Rutuja Patil on 3/8/2016.
/*****************************************************************************
 * Configuration.cpp Implementation File
 *
 ****************************************************************************/
#include <iostream>
#include <string>
#include <map>
#include <stdlib.h>
#include <sstream>
#include <cstring>

#include "Configuration.h"

/*****************************************************************************
 *
 * parseInt() helper function
 *
 ****************************************************************************/
int Configuration::parseInt(char* str){
    return (int) strtol( str, NULL, 10 );

}

/*****************************************************************************
 *
 * Constructor that parses command line.
 *
 ****************************************************************************/
Configuration::Configuration() : mCount(0) {   

    mUsageString = "usage: %s\n";
    
    // Initializing command-line argument names, init value, and help strings.
    // Block size
    addParamInt("Nx",'k',50,"--Nx <problem-size> for x-axis in elements (N)");
    addParamInt("Ny",'l',50,"--Ny <problem-size> for y-axis in elements (N)");
    addParamInt("Nz",'m',50,"--Nz <problem-size> for z-axis in elements (N)");
    addParamInt("bx",'x',1,"--bx <block-size> for x-axis in elements");
    addParamInt("by",'y',1,"--by <block-size> for y-axis in elements");
    addParamInt("num_threads",'p',1,"-p <num_threads>, number of cores");
    addParamBool("n",'n', false, "-n do not print time");
    addParamBool("v", 'v', false, "-v verify output\n");
  
    // Terminating the mLongOptions array.
    mLongOptions[mCount].name = 0;
    mLongOptions[mCount].has_arg = 0;
    mLongOptions[mCount].flag = 0;
    mLongOptions[mCount].val = 0;

}
/*****************************************************************************
 *
 * Destructor
 *
 ****************************************************************************/
Configuration::~Configuration(){

 for (int  i = 0; i < mCount; i++){
  delete mLongOptions[i].name;
 }
} 
/******************************************************************************
*
* Method to parse() commandLine
*
* *******************************************************************************/
void Configuration::parse(int argc,char* argv[]){
    int c;
    mExecName = std::string(argv[0]);     //executable is always here
    do {
        int option_index = 0;
        c = getopt_long(argc, argv, mCommandChars.c_str(), mLongOptions,
                        &option_index);
         //If this is a known character then process it.
        if (mKnownChars.find(c)!=mKnownChars.end()) {
            processArgVal(mCharToFieldname[c], optarg);
        }
    } while (c!=-1);

}

/*****************************************************************************
 *
 * toLDAPString() Method
 *
 ****************************************************************************/
std::string Configuration::toLDAPString()
{ 
    std::stringstream ldap;
    ldap << "Executable:" << mExecName << ",";
    bool first = true;
    for (std::list<int>::iterator iter=mCharOrder.begin();
            iter!=mCharOrder.end(); iter++) {
        if (!first) { ldap << ","; }
        else { first = false; }
        ldap << mCharToFieldname[*iter] << ":" 
           << mIntMap[mCharToFieldname[*iter]];

    }
    ldap << ",";
    return  ldap.str();
}

/*****************************************************************************
 *
 * Methods that add fields to configuration.
 *
 ****************************************************************************/
void Configuration::addParamInt(std::string fieldname, char single_char, 
        int init_val, std::string help_str) {
    mIntMap[fieldname] = init_val;
    
    mCharToFieldname[single_char] = fieldname;
    mKnownChars.insert(single_char);
    mCharOrder.push_back(single_char);
    
    mUsageString += help_str + "\n";
    
    mCommandChars += std::string(1,single_char) + ":";

    // Set up the long options data structure.
    // Assuming all arguments are optional for now.
    char* temp = new char[fieldname.length()+1];
    std::strcpy (temp, fieldname.c_str());
    mLongOptions[mCount].name = temp;
    // required_argument indicates need for int value after flag
    mLongOptions[mCount].has_arg = required_argument;
    mLongOptions[mCount].flag = 0;
    mLongOptions[mCount].val = single_char;
    mCount++;

    
}

/*****************************************************************************
 *
 * Methods that add fields to configuration.
 *
 ****************************************************************************/
void Configuration::addParamBool(std::string fieldname, char single_char,
        bool init_val, std::string help_str) {
    mBoolMap[fieldname] = init_val;
    mCharToFieldname[single_char] = fieldname;
    mKnownChars.insert(single_char);
    mUsageString += help_str + "\n";
    mCommandChars += std::string(1,single_char);
 
    // Set up the long options data structure.
    // Assuming all arguments are optional for now.
    char* temp = new char[fieldname.length()+1];
    std::strcpy (temp, fieldname.c_str());
    mLongOptions[mCount].name = temp;
    // required_argument indicates need for int value after flag
    mLongOptions[mCount].has_arg = no_argument;
    mLongOptions[mCount].flag = 0;
    mLongOptions[mCount].val = single_char;
    mCount++;
}

/*************************************************************************
 *
 * Method that add fields to configuration
 *
 * ***********************************************************************/
void Configuration::addParamString(std::string fieldname,char single_char,
      std::string init_val,std::string help_str){

      mStringMap[fieldname] = init_val;
      mCharToFieldname[single_char] = fieldname;
      
      mKnownChars.insert(single_char);
      mCharOrder.push_back(single_char);

      mUsageString += help_str + "\n";

      mCommandChars += std::string(1,single_char) + ":";

    // Set up the long options data structure.
    // Assuming all arguments are optional for now.
    char* temp = new char[fieldname.length()+1];
    std::strcpy (temp, fieldname.c_str());
    mLongOptions[mCount].name = temp;
    // required_argument indicates need for string value after flag
      mLongOptions[mCount].has_arg = required_argument;
      mLongOptions[mCount].flag = 0;
      mLongOptions[mCount].val = single_char;
      mCount++;
    
 }  

/*****************************************************************************
 *
 * Method that processes parsed field values.
 *
 ****************************************************************************/
void Configuration::processArgVal(std::string fieldname, char* arg_val) {

    // use the fieldname to locate which map we are looking in
    if (mIntMap.find(fieldname)!=mIntMap.end()) {
        mIntMap[fieldname] = parseInt( arg_val );
    } else if(mBoolMap.find(fieldname) != mBoolMap.end()){
    	mBoolMap[fieldname] = 1;
    } else if(mStringMap.find(fieldname) != mStringMap.end()){
        mStringMap[fieldname] =  std::string( arg_val );
    }else {
        std::cerr << "Unknown fieldname: " << fieldname << std::endl;
        exit(-1);
    }
}
/*****************************************************************************
 *
 * Method that returns int type fieldValue given fieldname.
 *
 ****************************************************************************/

int Configuration::getInt(std::string fieldname){

    if (mIntMap.find(fieldname) != mIntMap.end()){
        return mIntMap[fieldname];
    } else {
        std::cerr << "Unknown fieldname: " << fieldname << std::endl;
        exit(-1);
    }
}

/*****************************************************************************
 *
 * Method that returns bool type fieldValue given fieldname.
 *
 ****************************************************************************/

bool Configuration::getBool(std::string fieldname){

    if (mBoolMap.find(fieldname) != mBoolMap.end()){
        return mBoolMap[fieldname];
    } else {
        std::cerr << "Unknown fieldname: " << fieldname << std::endl;
        exit(-1);
    }
}


/**************************************************************************
*
* Method that returns String type filedValue given fieldname.
*
**************************************************************************/

std::string Configuration::getString(std::string fieldname){

   if(mStringMap.find(fieldname) != mStringMap.end()){
      return mStringMap[fieldname];
   }else{
     std::cerr << "Unknown fieldname: " << fieldname << std::endl;
     exit(-1);

   }

}     
