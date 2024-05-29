#include "optionsdata.h"

OptionsData::OptionsData()
{

}

void OptionsData::setProgramName(QString programName){
    this->programName = programName;
}

QString OptionsData::getProgramName(){
    return this->programName;
}

void OptionsData::setHitlistSize(int hitlistSize){
    this->hitlistSize = hitlistSize;
}

int OptionsData::getHitlistSize(){
    return this->hitlistSize;
}

void OptionsData::setDatabaseName(QString databaseName){
    this->databaseName = databaseName;
}

QString OptionsData::getDatabaseName(){
    return this->databaseName;
}

void OptionsData::setIsParsed(bool isParsed){
    this->isParsed = isParsed;
}

bool OptionsData::getIsParsed(){
    return this->isParsed;
}

void OptionsData::setEvalue(double evalue){
    this->evalue = evalue;
}

double OptionsData::getEvalue(){
    return this->evalue;
}

void OptionsData::setPenalty(double penalty){
    this->penalty = penalty;
}

double OptionsData::getPenalty(){
    return this->penalty;
}

void OptionsData::setReward(double reward){
    this->reward = reward;
}

double OptionsData::getReward(){
    return this->reward;
}
