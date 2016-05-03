//Created by Erman Gurses on 2/8/16.
/*****************************************************************************
 * Measurements.cpp Implementation File
 *
 ****************************************************************************/
#include <iostream>
#include <string>
#include <map>
#include <stdlib.h>
#include <iomanip>

#include "Measurements.h"

/*****************************************************************************
 *
 * Default Constructor Measurements()
 *
 ****************************************************************************/
Measurements::Measurements(){

}
/*****************************************************************************
 *
 * setField Method for int type fieldVal
 *
 ****************************************************************************/
void Measurements::setField(std::string fieldName, int fieldVal){
  IntMap[fieldName]=fieldVal;
}
/*****************************************************************************
 *
 * setField Method for float type fieldVal
 *
 ****************************************************************************/
void Measurements::setField(std::string fieldName, float fieldVal){
  FloatMap[fieldName]=fieldVal;
}
/*****************************************************************************
 *
 * setField Function for double type fieldVal
 *
 ****************************************************************************/
void Measurements::setField(std::string fieldName, double fieldVal){
  DoubleMap[fieldName]=fieldVal;
}
/*****************************************************************************
 *
 * setField Method for string type fieldVal
 *
 ****************************************************************************/
void Measurements::setField(std::string fieldName, std::string fieldVal){
  StringMap[fieldName]=fieldVal;
}
/*****************************************************************************
 *
 * Method that returns float type field value given field name.
 *
 ****************************************************************************/
int Measurements::getFieldInt( std::string fieldName){

  if (IntMap.find(fieldName) != IntMap.end()){
      return IntMap[fieldName];
  } else {
      std::cerr << "Unknown fieldName: " << fieldName << std::endl;
      exit(-1);
  }
}
/*****************************************************************************
 *
 * Method that returns float type field value given field name.
 *
 ****************************************************************************/
float Measurements::getFieldFloat(std::string fieldName){

  if (FloatMap.find(fieldName) != FloatMap.end()){
    return FloatMap[fieldName];
  } else {
    std::cerr << "Unknown fieldName: " << fieldName << std::endl;
    exit(-1);
  }
}
/*****************************************************************************
 *
 * Function that returns double type field value given field name.
 *
 ****************************************************************************/
double Measurements::getFieldDouble(std::string fieldName){

  if (DoubleMap.find(fieldName) != DoubleMap.end()){
    return DoubleMap[fieldName];
  } else {
    std::cerr << "Unknown fieldName: " << fieldName << std::endl;
    exit(-1);
  }
}
/*****************************************************************************
 *
 * Method that returns string type field value given field name.
 *
 ****************************************************************************/
std::string Measurements::getFieldString(std::string fieldName){

  if (StringMap.find(fieldName) != StringMap.end()){
    return StringMap[fieldName];
  } else {
    std::cerr << "Unknown fieldName: " << fieldName << std::endl;
    exit(-1);
  }
}
/*****************************************************************************
 *
 * toLDAPString() Method creates ldap string and returns it.
 *
 ****************************************************************************/
std::string Measurements::toLDAPString(){

  for (iterInt = IntMap.begin(); iterInt!= IntMap.end(); iterInt++){
    ldap << iterInt->first;
    ldap <<":";
    ldap << iterInt->second;
    ldap <<",";
  }
  for (iterFloat = FloatMap.begin(); iterFloat!= FloatMap.end(); iterFloat++){
    ldap << iterFloat->first;
    ldap << ":";
    ldap << std::fixed <<  iterFloat->second;
    ldap <<",";
  }
  for (iterDouble = DoubleMap.begin(); iterDouble!= DoubleMap.end(); 
                                                               iterDouble++){
    ldap << iterDouble->first;
    ldap << ":";
    ldap << std::fixed <<  iterDouble->second;
    ldap <<",";
  }
  for (iterString = StringMap.begin(); iterString!= StringMap.end();
     iterString++){

    ldap << iterString->first;
    ldap <<":";
    ldap <<  iterString->second;
    ldap <<",";
  }
  return ldap.str();
}
