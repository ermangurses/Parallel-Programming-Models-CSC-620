//Created by Erman Gurses on 2/8/16.
/*****************************************************************************
 * Measurements.h header file
 *
 * This class can handle measurement parameters that are paired as
 *   (FieldName,FieldValue)
 * The type of the pairs are:
 *   (string, int)
 *   (string, float)
 *   (string, string)
 *              
 * Usage:
 *   //Set (FieldName,FieldValue) pairs
 *     Measurement measurement;
 *     // (string, int) pair
 *     measurement.setField("GlobalToLocalMemory",13434);
 *
 *     // (string, float) pair
 *     measurement.setField("DataTransferFromDevice",14.3423f);
 *
 *     // (string, string) pair
 *     measurement.setField("ExecutionStatus","SUCCESS");
 *
 *   //Get particular FieldValue given FieldName *
 *     int GlobalToLocalMemory = 
                               measurement.getField("GlobalToLocalMemory");
 *     float DataTransferFromDevice = 
                               measurement.getField("DataTransferFromDevice");
 *     string ExecutionStatus = 
                               measurement.getField("ExecutionStatus");
 *
 *   //Generate an LDAP string to output configuration.
 *      string outString = measurement.toLDAPString();
 ****************************************************************************/

#ifndef MEASUREMENTS_H
#define MEASUREMENTS_H

#include <iostream>
#include <string>
#include <map>
#include <sstream>

class Measurements{
  public:
      Measurements();
      int getFieldInt(std::string fieldName);
      float getFieldFloat(std::string fieldName);
      double getFieldDouble(std::string fieldName);      
      std::string getFieldString(std::string fieldName);
      void setField(std::string fieldName, int fieldVal);
      void setField(std::string fieldName, float fieldVal);
      void setField(std::string fieldName, double fieldVal);            
      void setField(std::string fieldName, std::string fieldVal);
      std::string toLDAPString();
    private:
      std::stringstream ldap;
      std::map<std::string, int>::iterator iterInt;
      std::map<std::string, float>::iterator iterFloat;
      std::map<std::string, double>::iterator iterDouble;            
      std::map<std::string, std::string>::iterator iterString;
      std::map<std::string, int> IntMap;
      std::map<std::string, float> FloatMap;
      std::map<std::string, double> DoubleMap;       
      std::map<std::string, std::string> StringMap;

};
#endif
