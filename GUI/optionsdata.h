#ifndef OPTIONSDATA_H
#define OPTIONSDATA_H

#include <QString>

class OptionsData
{
public:
    OptionsData();
    void setProgramName(QString programName);
    QString getProgramName();
    void setHitlistSize(int hitlistSize);
    int getHitlistSize();
    void setDatabaseName(QString databaseName);
    QString getDatabaseName();
    void setIsParsed(bool isParsed);
    bool getIsParsed();
    void setEvalue(double evalue);
    double getEvalue();
    void setPenalty(double penalty);
    double getPenalty();
    void setReward(double reward);
    double getReward();

private:
    QString programName;
    int hitlistSize;
    QString databaseName;
    bool isParsed;
    double evalue;
    double penalty;
    double reward;
};

#endif // OPTIONSDATA_H
