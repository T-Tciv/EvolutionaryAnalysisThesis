#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <algorithm>
#include <map>
#include <tuple>
#include <random>
#include <set>

using namespace std;

struct SamplePartData
{
    unsigned int sampleCount;
    unsigned int sampleBadStatusCount;
    unsigned int sampleNotNCcount;

    vector<unsigned int> queryGapsCounts;
    vector<unsigned int> resultGapsCounts;
    vector<unsigned int> mismatchesCounts;

    vector<unsigned int> lengths;
    vector< pair<unsigned int, unsigned int> > sites;
};

class AvData 
{
private:
    unsigned int _count;

    unsigned int _left_length;
    unsigned int _site_length;
    unsigned int _right_length;

    unsigned int _left_query_gap_count;
    unsigned int _left_result_gap_count;
    unsigned int _left_mismatch_count;
    unsigned int _site_query_gap_count;
    unsigned int _site_result_gap_count;
    unsigned int _site_mismatch_count;
    unsigned int _right_query_gap_count;
    unsigned int _right_result_gap_count;
    unsigned int _right_mismatch_count;

    double _left_query_gap_density;
    double _left_result_gap_density;
    double _left_mismatch_density;
    double _site_query_gap_density;
    double _site_result_gap_density;
    double _site_mismatch_density;
    double _right_query_gap_density;
    double _right_result_gap_density;
    double _right_mismatch_density;

    double _adjacent_query_gap_density;
    double _adjacent_result_gap_density;
    double _adjacent_mismatch_density;

public:
    AvData() {
        _count = 0;
        _left_length = 0;
        _site_length = 0;
        _right_length = 0;
        _left_query_gap_count = 0;
        _left_result_gap_count = 0;
        _left_mismatch_count = 0;
        _site_query_gap_count = 0;
        _site_result_gap_count = 0;
        _site_mismatch_count = 0;
        _right_query_gap_count = 0;
        _right_result_gap_count = 0;
        _right_mismatch_count = 0;
        _left_query_gap_density = 0.0;
        _left_result_gap_density = 0.0;
        _left_mismatch_density = 0.0;
        _site_query_gap_density = 0.0;
        _site_result_gap_density = 0.0;
        _site_mismatch_density = 0.0;
        _right_query_gap_density = 0.0;
        _right_result_gap_density = 0.0;
        _right_mismatch_density = 0.0;

        _adjacent_query_gap_density = 0.0;
        _adjacent_result_gap_density = 0.0;
        _adjacent_mismatch_density = 0.0;
    }

    void increaseData(unsigned int left_length, unsigned int site_length, unsigned int right_length, 
                 unsigned int left_query_gap_count, unsigned int left_result_gap_count, unsigned int left_mismatch_count,
                 unsigned int site_query_gap_count, unsigned int site_result_gap_count, unsigned int site_mismatch_count,
                 unsigned int right_query_gap_count, unsigned int right_result_gap_count, unsigned int right_mismatch_count) {
        ++_count;
        _left_length += left_length;
        _right_length += right_length;
        _site_length += site_length;
        _left_query_gap_count += left_query_gap_count;
        _left_result_gap_count += left_result_gap_count;
        _left_mismatch_count += left_mismatch_count;
        _site_query_gap_count += site_query_gap_count;
        _site_result_gap_count += site_result_gap_count;
        _site_mismatch_count += _site_mismatch_count;
        _right_query_gap_count += right_query_gap_count;
        _right_result_gap_count += right_result_gap_count;
        _right_mismatch_count += right_mismatch_count;

        _left_query_gap_density += (left_length == 0) ? 0 : (1.0 * left_query_gap_count / left_length);
        _left_result_gap_density += (left_length == 0) ? 0 : (1.0 * left_result_gap_count / left_length);
        _left_mismatch_density += (left_length == 0) ? 0 : (1.0 * left_mismatch_count / left_length);
        _site_query_gap_density += (site_length == 0) ? 0 : (1.0 * site_query_gap_count / site_length);
        _site_result_gap_density += (site_length == 0) ? 0 : (1.0 * site_result_gap_count / site_length);
        _site_mismatch_density += (site_length == 0) ? 0 : (1.0 * site_mismatch_count / site_length);
        _right_query_gap_density += (right_length == 0) ? 0 : (1.0 * right_query_gap_count / right_length);
        _right_result_gap_density += (right_length == 0) ? 0 : (1.0 * right_result_gap_count / right_length);
        _right_mismatch_density += (right_length == 0) ? 0 : (1.0 * right_mismatch_count / right_length);

        unsigned int adjacent_length = left_length + right_length;
        _adjacent_query_gap_density = (adjacent_length == 0) ? 0 : (1.0 * (left_query_gap_count + right_query_gap_count) / adjacent_length);
        _adjacent_result_gap_density = (adjacent_length == 0) ? 0 : (1.0 * (left_result_gap_count + right_result_gap_count) / adjacent_length);
        _adjacent_mismatch_density = (adjacent_length == 0) ? 0 : (1.0 * (left_mismatch_count + right_mismatch_count) / adjacent_length);
    }

    string getPrintableAvData() {
        stringstream string_out;
        string delim_string = "\t";

        string_out << _count << delim_string;
        string_out << _left_length * 1.0 / _count << delim_string;
        string_out << _left_query_gap_count * 1.0 / _count << delim_string << _left_query_gap_density / _count << delim_string;
        string_out << _left_result_gap_count * 1.0 / _count << delim_string << _left_result_gap_density / _count << delim_string;
        string_out << _left_mismatch_count * 1.0 / _count << delim_string << _left_mismatch_density / _count << delim_string;
        string_out << _site_length * 1.0 / _count << delim_string;
        string_out << _site_query_gap_count * 1.0 / _count << delim_string << _site_query_gap_density / _count << delim_string;
        string_out << _site_result_gap_count * 1.0 / _count << delim_string << _site_result_gap_density / _count << delim_string;
        string_out << _site_mismatch_count * 1.0 / _count << delim_string << _site_mismatch_density / _count << delim_string;
        string_out << _right_length * 1.0 / _count << delim_string;
        string_out << _right_query_gap_count * 1.0 / _count << delim_string << _right_query_gap_density / _count << delim_string;
        string_out << _right_result_gap_count * 1.0 / _count << delim_string << _right_result_gap_density / _count << delim_string;
        string_out << _right_mismatch_count * 1.0 / _count << delim_string << _right_mismatch_density / _count << delim_string;
        string_out << endl;

        return string_out.str();
    }

    string getPrintableTotalAvData() {
        stringstream string_out;
        string delim_string = "\t";

        string_out << _count << delim_string;

        string_out << _left_length * 1.0 / _count << delim_string;
        string_out << (_left_query_gap_count + _left_result_gap_count + _left_mismatch_count) * 1.0 / _count << delim_string;
        string_out << (_left_query_gap_density + _left_result_gap_density + _left_mismatch_density) / _count << delim_string;

        string_out << _site_length * 1.0 / _count << delim_string;
        string_out << (_site_query_gap_count + _site_result_gap_count + _site_mismatch_count) * 1.0 / _count << delim_string;
        string_out << (_site_query_gap_density + _site_result_gap_density + _site_mismatch_density) / _count << delim_string;

        string_out << _right_length * 1.0 / _count << delim_string;
        string_out << (_right_query_gap_count + _right_result_gap_count + _right_mismatch_count) * 1.0 / _count << delim_string;
        string_out << (_right_query_gap_density + _right_result_gap_density + _right_mismatch_density) / _count << delim_string;

        string_out << endl;

        return string_out.str();
    }

    vector<double> getDensitiesQuotients() {
        double totalLeftDensity = _left_query_gap_density + _left_result_gap_density + _left_mismatch_density;
        double totalRightDensity = _right_query_gap_density + _right_result_gap_density + _right_mismatch_density;
        double totalSiteDensity = _site_query_gap_density + _site_result_gap_density + _site_mismatch_density;
        
        vector<double> quotients =             {
                                                totalSiteDensity * 1.0 / _count,
                                                totalSiteDensity / totalLeftDensity,
                                                _site_query_gap_density / _left_query_gap_density,
                                                _site_result_gap_density / _left_result_gap_density,
                                                _site_mismatch_density / _left_mismatch_density,
                                                totalSiteDensity / totalRightDensity,
                                                _site_query_gap_density / _right_query_gap_density,
                                                _site_result_gap_density / _right_result_gap_density,
                                                _site_mismatch_density / _right_mismatch_density,
                                                totalSiteDensity / (_adjacent_query_gap_density + _adjacent_result_gap_density + _adjacent_mismatch_density),
                                                _site_query_gap_density / _adjacent_query_gap_density,
                                                _site_result_gap_density / _adjacent_result_gap_density,
                                                _site_mismatch_density / _adjacent_mismatch_density,
                                               };

        return quotients;
    }
};

enum class SpecialCharacter : char
{ Gap = '-', SeqEnd = '+', FileExt = '.', FileSystemDiv = '/', FileNameDiv = '_', TableDelim = '\t' };

enum class SequenceStatus : int
{ Ok = 0, OneSide = 1, OnlySite = 2, NoSite = 3, Several = 4};

string removeSpaces(string line) 
{ 
    auto noSpaceEnd = remove(line.begin(), line.end(), ' ');
    line.erase(noSpaceEnd, line.end());
    return line; 
}

vector< pair<unsigned int, unsigned int> > getBindingSitesPositions(string sequence)
{
    unsigned int position = 0;
    unsigned int length = sequence.size();
    vector< pair<unsigned int, unsigned int> > positions;
    unsigned int currentLeft = 0;
    unsigned int currentRight = length - 1;

    while (position < length) {
        while (position < length && (!isupper(sequence[position]))) {
            ++position;
        }

        currentLeft = position;

        while (position < length && (!islower(sequence[position]))) {
            if(isupper(sequence[position])) {
                currentRight = position;
            }
            ++position;
        }

        if (currentLeft <= currentRight) {
            positions.push_back( pair<unsigned int, unsigned int>(currentLeft, currentRight) );
        }      
    }

    return positions;
}

string getIdFromFastaLine(string fastaLine)
{
    // Пример Fasta-строки:
    // > S2260;+; AC:X00855:324 TRRD;Hs:DHR ID: BS_E2F_DP1_NS00001
    // Из неё в этой функции нужно получить S2260
    // Но в ней может и не быть точек с запятой и знаков, база заполнена по-разному

    int firstSemicolonPos = fastaLine.find_first_of(';');
    int firstSpacePos = fastaLine.find_first_of(' ');
    int secondSpacePos = fastaLine.substr(firstSpacePos + 1).find_first_of(' ') + firstSpacePos + 1;

    int idLength = min(firstSemicolonPos, secondSpacePos) - firstSpacePos - 1;
    string id = fastaLine.substr(firstSpacePos + 1, idLength);

    return id;
}

void processAligments(string id, string result_acc, string query, string result, string strand, int hit_number, string evalue, 
                        ofstream &out, AvData & avData, SamplePartData & samplePartData) 
{
    unsigned int length = query.size();

    vector< pair<unsigned int, unsigned int> > sites =  getBindingSitesPositions(query);

    stringstream string_out;
    char delim_string = static_cast<char>(SpecialCharacter::TableDelim);
    string_out << id << delim_string << result_acc << delim_string << strand << delim_string << evalue << delim_string;
    int seq_status;
    pair<unsigned int, unsigned int> binding_site_positions;

    if (sites.size() == 0) {
        seq_status = static_cast<int>( SequenceStatus::NoSite );
        ++samplePartData.sampleBadStatusCount;
        string_out << seq_status << delim_string;
        string_out << "-" << delim_string << "-" << delim_string << "-" << delim_string;
        string_out << "-" << delim_string << "-" << delim_string << "-" << delim_string;
        string_out << "-" << delim_string << "-" << delim_string << "-" << delim_string << hit_number << endl;

        out << string_out.str();
        return;
    }

    if (sites.size() > 1) {
        #if 1
        binding_site_positions = {sites[0].first, sites[sites.size() - 1].second};
        #else
        seq_status = static_cast<int>( SequenceStatus::Several );
        ++samplePartData.sampleBadStatusCount;
        string_out << seq_status << delim_string;
        string_out << "-" << delim_string << "-" << delim_string << "-" << delim_string;
        string_out << "-" << delim_string << "-" << delim_string << "-" << delim_string;
        string_out << "-" << delim_string << "-" << delim_string << "-" << delim_string << hit_number << endl;

        out << string_out.str();
        return;
        #endif
    }

    binding_site_positions = sites[0];

    samplePartData.lengths.push_back(length);
    samplePartData.sites.push_back(binding_site_positions);

    unsigned int binding_site_left = binding_site_positions.first;
    unsigned int binding_site_right = binding_site_positions.second;

    unsigned int left_length = binding_site_left;
    unsigned int site_length = binding_site_right - binding_site_left + 1;
    unsigned int right_length = length - binding_site_right - 1;

    if (left_length == 0 && right_length == 0) {
        seq_status = static_cast<int>( SequenceStatus::OnlySite );
    } else if (left_length == 0 || right_length == 0) {
        seq_status = static_cast<int>( SequenceStatus::OneSide );
    } else {
        seq_status = static_cast<int>( SequenceStatus::Ok );
    }

    char gap = static_cast<char>( SpecialCharacter::Gap );
    char query_nucl;
    char result_nucl;
    unsigned int left_query_gap_count = 0;
    unsigned int left_result_gap_count = 0;
    unsigned int left_mismatch_count = 0;
    unsigned int site_query_gap_count = 0;
    unsigned int site_result_gap_count = 0;
    unsigned int site_mismatch_count = 0;
    unsigned int right_query_gap_count = 0;
    unsigned int right_result_gap_count = 0;
    unsigned int right_mismatch_count = 0;

    if (binding_site_left >= length) {
        left_length = 0;
    }

    if (binding_site_right >= length) {
        right_length = 0;
    }

    for (unsigned int i {0}; i < binding_site_left; i++) { 
        query_nucl = toupper(query[i]);
        result_nucl = toupper(result[i]);

        if (query_nucl == gap) {
            ++left_query_gap_count;
        }
        if (result_nucl == gap) {
            ++left_result_gap_count;
        }
        if (query_nucl != gap && result_nucl != gap && query_nucl != result_nucl) {
            ++left_mismatch_count;
        }
    }

    for (unsigned int i {binding_site_left}; i <= binding_site_right; i++) { 
        query_nucl = toupper(query[i]);
        result_nucl = toupper(result[i]);

        if (query_nucl == gap) {
            ++site_query_gap_count;
        }
        if (result_nucl == gap) {
            ++site_result_gap_count;
        }
        if (query_nucl != gap && result_nucl != gap && query_nucl != result_nucl) {
            ++site_mismatch_count;
        }
    }

    for (unsigned int i {binding_site_right + 1}; i < length; i++) { 
        query_nucl = toupper(query[i]);
        result_nucl = toupper(result[i]);

        if (query_nucl == gap) {
            ++right_query_gap_count;
        }
        if (result_nucl == gap) {
            ++right_result_gap_count;
        }
        if (query_nucl != gap && result_nucl != gap && query_nucl != result_nucl) {
            ++right_mismatch_count;
        }
    }

    samplePartData.mismatchesCounts.push_back(left_mismatch_count + site_mismatch_count + right_mismatch_count);
    samplePartData.queryGapsCounts.push_back(left_query_gap_count + site_query_gap_count + right_query_gap_count);
    samplePartData.resultGapsCounts.push_back(left_result_gap_count + site_result_gap_count + right_result_gap_count);

    double left_query_gap_density = 1.0 * left_query_gap_count / left_length;
    double left_result_gap_density = 1.0 * left_result_gap_count / left_length;
    double left_mismatch_density = 1.0 * left_mismatch_count / left_length;
    double site_query_gap_density = 1.0 * site_query_gap_count / site_length;
    double site_result_gap_density = 1.0 * site_result_gap_count / site_length;
    double site_mismatch_density = 1.0 * site_mismatch_count / site_length;
    double right_query_gap_density = 1.0 * right_query_gap_count / right_length;
    double right_result_gap_density = 1.0 * right_result_gap_count / right_length;
    double right_mismatch_density = 1.0 * right_mismatch_count / right_length;

    string_out << seq_status << delim_string;
    string_out << left_length << delim_string;
    string_out << left_query_gap_count << delim_string << left_query_gap_density << delim_string;
    string_out << left_result_gap_count << delim_string << left_result_gap_density << delim_string;
    string_out << left_mismatch_count << delim_string << left_mismatch_density << delim_string;
    string_out << site_length << delim_string;
    string_out << site_query_gap_count << delim_string << site_query_gap_density << delim_string;
    string_out << site_result_gap_count << delim_string << site_result_gap_density << delim_string;
    string_out << site_mismatch_count << delim_string << site_mismatch_density << delim_string;
    string_out << right_length << delim_string;
    string_out << right_query_gap_count << delim_string << right_query_gap_density << delim_string;
    string_out << right_result_gap_count << delim_string << right_result_gap_density << delim_string;
    string_out << right_mismatch_count << delim_string << right_mismatch_density << delim_string;
    string_out << hit_number << endl;

    out << string_out.str();

    if (result_acc.find("NC_") != 0) {
        ++samplePartData.sampleNotNCcount;
        cerr << "Not full chromosome | AC: " << result_acc << endl;
        return;
    }

    //cout << id << " " << result_acc << endl;

    ++samplePartData.sampleCount;
    avData.increaseData(left_length, site_length, right_length, 
                 left_query_gap_count, left_result_gap_count, left_mismatch_count,
                 site_query_gap_count, site_result_gap_count, site_mismatch_count,
                 right_query_gap_count, right_result_gap_count, right_mismatch_count);

}

vector<double> getSampleDensitiesOfErrorsSet(set<unsigned int> sampleErrorsSet, 
                                                            vector<unsigned int> lengths, 
                                                            vector< pair<unsigned int, unsigned int> > sites)
{
    double leftDensitySum = 0.0;
    double siteDensitySum = 0.0;
    double rightDensitySum = 0.0;
    double adjacentDensitySum = 0.0;

    unsigned int sequenceStartShift = 0; // Сдвиг позиций ошибок в общей цепочке из всех последовательностей выборки
    unsigned int currentSequenceNumber = 0;

    // Переменные для хранения позиций начала и конца участка-сайта в текущей последовательности
    unsigned int siteStartPos = sites[0].first;
    unsigned int siteEndPos = sites[0].second;

    // Количество ошибок рассматриваемого в данный момент типа (циклы для всех трёх типов идут последовательно)
    // в конкретной последовательности (т.е. зануляются при переходе к следующей)
    unsigned int leftErrorsCount = 0;
    unsigned int siteErrorsCount = 0;
    unsigned int rightErrorsCount = 0;

    // Булева переменная, чтобы определять, нужно ли при переходе к новой последовательности сохранять данные о плотностях в старой
    // (а это может и не потребоваться, если мы пропускаем какую-то последовательность, потому что в ней нет ошибок)
    // выставляется в true при выходе из цикла while, в котором и происходит переход от одной последовательности к другой
    bool isPrevSequenceProcessed = false;

    // Идём по множеству ошибок нужного типа, параллельно сдвигаясь по вектору последовательностей 
    for (auto it = sampleErrorsSet.begin(); it != sampleErrorsSet.end(); ++it) {
        unsigned int currentErrorPos = *it;
        // (переходим к следущей последовательности, если позиция текущей ошибки больше, чем позиция конца текущей последовательности)
        while (sequenceStartShift + lengths[currentSequenceNumber] < currentErrorPos) {
            // Если мы обработали последовательность, нужно посчитать в ней плотности
            if (isPrevSequenceProcessed) { 
                // Если левой части нет, то плотность ошибок в ней принимается нулевой
                unsigned int leftLength = siteStartPos;
                double leftDensity = static_cast<double>(leftErrorsCount) / leftLength;
                leftDensitySum += (leftLength == 0) ? 0 : leftDensity;
                // Последовательности без сайта просто не обрабатываются, поэтому здесь случай с делением на ноль деле не нужен
                unsigned int siteLength = siteEndPos - siteStartPos + 1;
                double siteDensity = static_cast<double>(siteErrorsCount) / siteLength;
                siteDensitySum += siteDensity;
                // Для правой, как для левой
                unsigned int rightLength = lengths[currentSequenceNumber] - siteEndPos - 1;
                double rightDensity = static_cast<double>(rightErrorsCount) / rightLength;
                rightDensitySum += (rightLength == 0) ? 0 : rightDensity;

                // Считаем для по двум фланкирующим участкам разом
                unsigned int adjacentLength = leftLength + rightLength;
                double adjacentDensity = (adjacentLength == 0) 
                                         ? 
                                         0 
                                         : 
                                         ((static_cast<double>(leftErrorsCount) + static_cast<double>(rightErrorsCount)) / adjacentLength);
                adjacentDensitySum += adjacentDensity;

                isPrevSequenceProcessed = false;
                
                #if 0
                cout << "seq:" << currentSequenceNumber << " dens: " << leftDensity << " " << siteDensity << " " << rightDensity << endl;
                #endif
                #if 0
                cout << "seq:" << currentSequenceNumber << " errcnt: " << leftErrorsCount << " " << siteErrorsCount << " " << rightErrorsCount << endl;
                #endif
            }

            sequenceStartShift += lengths[currentSequenceNumber];
            ++currentSequenceNumber;
            siteStartPos = sites[currentSequenceNumber].first;
            siteEndPos = sites[currentSequenceNumber].second;
            #if 0
            cout << "seq:" << currentSequenceNumber << "start: " << sequenceStartShift << " siteStart: " << siteStartPos << " siteEnd: " << siteEndPos << endl;
            #endif
            leftErrorsCount = 0;
            siteErrorsCount = 0;
            rightErrorsCount = 0;
        }

        isPrevSequenceProcessed = true;

        // В зависимости от положения позиции ошибки относительно сайта в этой последовательности (последовательность = left+site+right)
        // увеличиваем число ошибок данного типа в текущей последовательности
        if (currentErrorPos < sequenceStartShift + siteStartPos) {
            #if 0
            cout << "seq:" << currentSequenceNumber << " (left) " << currentErrorPos << " < " << sequenceStartShift << " + " << siteStartPos << endl;
            #endif
            ++leftErrorsCount;
        } else if (currentErrorPos > sequenceStartShift + siteEndPos) {
            #if 0
            cout << "seq:" << currentSequenceNumber << " (right) " << currentErrorPos << " > " << sequenceStartShift << " + " << siteEndPos << endl;
            #endif
            ++rightErrorsCount;
        } else {
            #if 0
            cout << "seq:" << currentSequenceNumber << " (site) " << currentErrorPos << endl;
            #endif
            ++siteErrorsCount;
        }
    }

    // Тот же подсчёт плотностей, что и в цикле, но здесь он нужен, потому что иначе последняя ошибка не будет учтена
    // TODO: вынести в отдельную функцию, что ли..........
    unsigned int leftLength = siteStartPos;
    double leftDensity = static_cast<double>(leftErrorsCount) / leftLength;
    leftDensitySum += (leftLength == 0) ? 0 : leftDensity;
    unsigned int siteLength = siteEndPos - siteStartPos + 1;
    double siteDensity = static_cast<double>(siteErrorsCount) / siteLength;
    siteDensitySum += siteDensity;
    unsigned int rightLength = lengths[currentSequenceNumber] - siteEndPos - 1;
    double rightDensity = static_cast<double>(rightErrorsCount) / rightLength;
    rightDensitySum += (rightLength == 0) ? 0 : rightDensity;
    unsigned int adjacentLength = leftLength + rightLength;
    double adjacentDensity = (adjacentLength == 0) 
        ? 
        0 
        : 
        ((static_cast<double>(leftErrorsCount) + static_cast<double>(rightErrorsCount)) / adjacentLength);
    adjacentDensitySum += adjacentDensity;

    #if 0
    cout << "seq:" << currentSequenceNumber << " errcnt: " << leftErrorsCount << " " << siteErrorsCount << " " << rightErrorsCount << endl;
    #endif

    vector<double> densities = {leftDensitySum, siteDensitySum, rightDensitySum, adjacentDensitySum};
    return densities;
}

double getRightMonteCarloValue(vector<double> quotients, double realQuotient)
{
    /*
    lower_bound -- Returns an iterator pointing to the first element in the range [first, last) that is not less than value
    uper_bound -- Returns an iterator pointing to the first element in the range [first, last) that is greater than value
    */
    sort(quotients.begin(), quotients.end());
    return (quotients.size() - (lower_bound(quotients.begin(), quotients.end(), realQuotient) - quotients.begin())) * 100.0 / quotients.size();
}

double getLeftMonteCarloValue(vector<double> quotients, double realQuotient)
{
    /*
    lower_bound -- Returns an iterator pointing to the first element in the range [first, last) that is not less than value
    uper_bound -- Returns an iterator pointing to the first element in the range [first, last) that is greater than value
    */
    sort(quotients.begin(), quotients.end());
    return (upper_bound(quotients.begin(), quotients.end(), realQuotient) - quotients.begin()) * 100.0 / quotients.size();  
    //TODO: убрать сортировку
}

vector<double> getSequenceDensitiesOfErrorsSet(set<unsigned int> seqErrorsSet, 
                                               unsigned int length, 
                                               pair<unsigned int, unsigned int> site)
{
    double leftCount = 0.0;
    double siteCount = 0.0;
    double rightCount = 0.0;

    for (auto it = seqErrorsSet.begin(); it != seqErrorsSet.end(); ++it) {
        unsigned int currentErrorPos = *it;

        if (currentErrorPos < site.first) {
            ++leftCount;
        } else if (currentErrorPos > site.second) {
            ++rightCount;
        } else {
            ++siteCount;
        }
    }

    double leftDensity = (site.first == 0) ? 0 : (leftCount * 1.0 / site.first);
    double siteDensity = siteCount * 1.0 / (site.second - site.first + 1);
    double rightDensity = (site.second == (length - 1)) ? 0 : (rightCount * 1.0 / (length - 1 - site.second));   

    double adjacentDensity = ((site.first == 0) && (site.second == (length - 1)))
                             ?
                             (0)
                             :
                             ((leftCount + rightCount) * 1.0 / (site.first + length - 1 - site.second));

    vector<double> densities = {leftDensity, siteDensity, rightDensity, adjacentDensity};
    return densities;
}

vector<double> runSingleSequenceMonteCarloExperiment (vector<unsigned int> queryGapsCounts, 
                                                      vector<unsigned int> resultGapsCounts, 
                                                      vector<unsigned int> mismatchesCounts, 
                                                      vector<unsigned int> lengths, 
                                                      vector< pair<unsigned int, unsigned int> > sites
)
{
    vector<double> average_densities_quotients;

    unsigned int sequencesCount = lengths.size();

    double totalLeftQueryGapsDensity = 0;
    double totalLeftResultGapsDensity = 0;
    double totalLeftMismatchesDensity = 0;

    double totalRightQueryGapsDensity = 0;
    double totalRightResultGapsDensity = 0;
    double totalRightMismatchesDensity = 0;

    double totalSiteQueryGapsDensity = 0;
    double totalSiteResultGapsDensity = 0;
    double totalSiteMismatchesDensity = 0;

    double totalAdjacentQueryGapsDensity = 0;
    double totalAdjacentResultGapsDensity = 0;
    double totalAdjacentMismatchesDensity = 0;

    random_device rnd_device; // a seed source for the random number engine
    mt19937 mersenne_engine{ rnd_device() }; // mersenne_twister_engine seeded with rd()

    // Для каждой последовательности:
    // 1. Генерируем позиции ошибок в последовательности
    // 2. Считаем плотности распределений ошибок в последовательности
    // 3. Суммируем
    //cout << "cycle start" << endl;
    for (unsigned int i = 0; i < sequencesCount; ++i) {
        //cout << "generation start" << endl;
        uniform_int_distribution<unsigned int> gen_error_pos{ 0, lengths[i] - 1 };
        //cout << "sets filling" << endl;

        set<unsigned int> seqQueryGapsSet;
        set<unsigned int> seqResultGapsSet;
        set<unsigned int> seqMismatchesSet;

        while (seqQueryGapsSet.size() != queryGapsCounts[i]) {
            seqQueryGapsSet.insert(gen_error_pos(mersenne_engine));
        }

        while (seqResultGapsSet.size() != resultGapsCounts[i]) {
            unsigned int error_pos = gen_error_pos(mersenne_engine);
            if (seqQueryGapsSet.count(error_pos) == 0) {
                seqResultGapsSet.insert(error_pos);
            }
        }

        while (seqMismatchesSet.size() != mismatchesCounts[i]) {
            unsigned int error_pos = gen_error_pos(mersenne_engine);
            if (seqQueryGapsSet.count(error_pos) == 0 && seqResultGapsSet.count(error_pos) == 0) {
                seqMismatchesSet.insert(error_pos);
            }
        }

        //cout << "generation end" << endl;

        vector<double> sequenceQueryGapsSetDensities = getSequenceDensitiesOfErrorsSet(seqQueryGapsSet, lengths[i], sites[i]);
        vector<double> sequenceResultGapsSetDensities = getSequenceDensitiesOfErrorsSet(seqResultGapsSet, lengths[i], sites[i]);
        vector<double> sequenceMismatchesSetDensities = getSequenceDensitiesOfErrorsSet(seqMismatchesSet, lengths[i], sites[i]);

        totalLeftQueryGapsDensity += sequenceQueryGapsSetDensities[0];
        totalSiteQueryGapsDensity += sequenceQueryGapsSetDensities[1];
        totalRightQueryGapsDensity += sequenceQueryGapsSetDensities[2];
        totalAdjacentQueryGapsDensity += sequenceQueryGapsSetDensities[3];

        totalLeftResultGapsDensity += sequenceResultGapsSetDensities[0];
        totalSiteResultGapsDensity += sequenceResultGapsSetDensities[1];
        totalRightResultGapsDensity += sequenceResultGapsSetDensities[2];
        totalAdjacentResultGapsDensity += sequenceResultGapsSetDensities[3];

        totalLeftMismatchesDensity += sequenceMismatchesSetDensities[0];
        totalSiteMismatchesDensity += sequenceMismatchesSetDensities[1];
        totalRightMismatchesDensity += sequenceMismatchesSetDensities[2];
        totalAdjacentMismatchesDensity += sequenceMismatchesSetDensities[3];
    }

    double totalSiteDensity = totalSiteQueryGapsDensity + totalSiteResultGapsDensity + totalSiteMismatchesDensity;
    double totalLeftDensity = totalLeftQueryGapsDensity + totalLeftResultGapsDensity + totalLeftMismatchesDensity;
    double totalRightDensity = totalRightQueryGapsDensity + totalRightResultGapsDensity + totalRightMismatchesDensity;
    double totalAdjacentDensity = totalAdjacentQueryGapsDensity + totalAdjacentResultGapsDensity + totalAdjacentMismatchesDensity;

    average_densities_quotients =          {
                                            totalSiteDensity * 1.0 / sequencesCount,
                                            totalSiteDensity / totalLeftDensity,
                                            totalSiteQueryGapsDensity / totalLeftQueryGapsDensity,
                                            totalSiteResultGapsDensity / totalLeftResultGapsDensity,
                                            totalSiteMismatchesDensity / totalLeftMismatchesDensity,
                                            totalSiteDensity / totalRightDensity,
                                            totalSiteQueryGapsDensity / totalRightQueryGapsDensity,
                                            totalSiteResultGapsDensity / totalRightResultGapsDensity,
                                            totalSiteMismatchesDensity / totalRightMismatchesDensity,  
                                            totalSiteDensity / totalAdjacentDensity,
                                            totalSiteQueryGapsDensity / totalAdjacentQueryGapsDensity,
                                            totalSiteResultGapsDensity / totalAdjacentResultGapsDensity,
                                            totalSiteMismatchesDensity / totalAdjacentMismatchesDensity                                           
                                            };

    return average_densities_quotients;
}    

vector<double> runSingleMonteCarloExperiment (unsigned int queryGapsCount, unsigned int resultGapsCount, unsigned int mismatchesCount, 
                                    vector<unsigned int> lengths, 
                                    vector< pair<unsigned int, unsigned int> > sites
)
{
    vector<double> average_densities_quotients;

    // 1. Раскидываем разные виды ошибок по выборке: 
    // чтобы не решать задачу "как раскидать ошибки по каждой последовательности так, чтобы сумма количеств ошибок равнялась общему числу ошибок",
    // сложим длины последовательностей, как если бы выборка была единой нуклеотидной цепочкой, раскидаем по ней ошибки,
    // а потом "разрежем" выборку обратно на последовательности
    unsigned int sampleLength = 0;
    for (auto& currentSequenceLength : lengths) {
        sampleLength += currentSequenceLength;
    }

    // Генерируем сами позиции ошибок, сохраняя их в соответсвующие множества
    random_device rnd_device; // a seed source for the random number engine
    mt19937 mersenne_engine{ rnd_device() }; // mersenne_twister_engine seeded with rd()
    uniform_int_distribution<unsigned int> gen_error_pos{ 0, sampleLength - 1 };

    set<unsigned int> sampleQueryGapsSet;
    set<unsigned int> sampleResultGapsSet;
    set<unsigned int> sampleMismatchesSet;

    #if 1
    // Каждая позиция в последовательности может либо не содержать ошибку, либо содержать ошибку ровно одного вида (вставка, делеция, несовпадение)
    // поэтому в каждом следующем цикле добавлется проверка на наличие сгенерированной позиции в множестве ошибок предыдущего типа
    while (sampleQueryGapsSet.size() != queryGapsCount) {
        sampleQueryGapsSet.insert(gen_error_pos(mersenne_engine));
    }

    while (sampleResultGapsSet.size() != resultGapsCount) {
        unsigned int error_pos = gen_error_pos(mersenne_engine);
        if (sampleQueryGapsSet.count(error_pos) == 0) {
            sampleResultGapsSet.insert(error_pos);
        }
    }

    while (sampleMismatchesSet.size() != mismatchesCount) {
        unsigned int error_pos = gen_error_pos(mersenne_engine);
        if (sampleQueryGapsSet.count(error_pos) == 0 && sampleResultGapsSet.count(error_pos) == 0) {
            sampleMismatchesSet.insert(error_pos);
        }
    }
    #else
    sampleQueryGapsSet = {3, 4, 12};
    sampleResultGapsSet = {14, 16, 26, 30};
    sampleMismatchesSet = {0, 5, 19, 31, 32};
    #endif

    // Код для вывода сгенерированных позиций
    # if 0
    cout << "Позиции ошибок (по типам):" << endl;
    for (auto it = sampleQueryGapsSet.begin(); it != sampleQueryGapsSet.end(); ++it) {
        cout << *it << " ";
    }
    cout << endl;
    for (auto it = sampleResultGapsSet.begin(); it != sampleResultGapsSet.end(); ++it) {
        cout << *it << " ";
    }
    cout << endl;
        for (auto it = sampleMismatchesSet.begin(); it != sampleMismatchesSet.end(); ++it) {
        cout << *it << " ";
    }
    cout << endl;
    #endif

    // 2. Разбиваем ошибки по последовательностям и считаем плотности
    // (в вызываемой функции мы для выбранного типа ошибки обходим соответсвующее множество позиций, переходя во внутреннем цикле к следующей последовательности,
    // когда позиция очередной ошибки выходит за границы текущей последовательности) 

    // !!!!!!!!!!!!!!!!!!! ЭТО СУММЫ ПЛОТНОСТЕЙ, ДЛЯ ПОЛУЧЕНИЯ СРЕДНЕЙ НУЖНО ЕЩЁ НА КОЛИЧЕСТВО ПОСЛЕДОВАТЕЛЬНОСТЕЙ ПОДЕЛИТЬ
    vector<double> sampleQueryGapsSetDensities = getSampleDensitiesOfErrorsSet(sampleQueryGapsSet, lengths, sites);
    vector<double> sampleResultGapsSetDensities = getSampleDensitiesOfErrorsSet(sampleResultGapsSet, lengths, sites);
    vector<double> sampleMismatchesSetDensities = getSampleDensitiesOfErrorsSet(sampleMismatchesSet, lengths, sites);

    // Код для выводы плотностей
    # if 0
    cout << "Плотности распределения ошибок (по типам):" << endl;
    for (auto it = sampleQueryGapsSetDensities.begin(); it != sampleQueryGapsSetDensities.end(); ++it) {
        cout << *it << " ";
    }
    cout << endl;
    for (auto it = sampleResultGapsSetDensities.begin(); it != sampleResultGapsSetDensities.end(); ++it) {
        cout << *it << " ";
    }
    cout << endl;
        for (auto it =  sampleMismatchesSetDensities.begin(); it !=  sampleMismatchesSetDensities.end(); ++it) {
        cout << *it << " ";
    }
    cout << endl;
    #endif

    // 3. Получаем общие для всех ошибок плотности по участкам, потому что я только сейчас вспомнила, что мне не нужны отдельные плотности
    double totalSiteDensity = sampleQueryGapsSetDensities[1] + sampleResultGapsSetDensities[1] + sampleMismatchesSetDensities[1];
    double totalLeftDensity = sampleQueryGapsSetDensities[0] + sampleResultGapsSetDensities[0] + sampleMismatchesSetDensities[0];
    double totalRightDensity = sampleQueryGapsSetDensities[2] + sampleResultGapsSetDensities[2] + sampleMismatchesSetDensities[2];
    double totalAdjacentDensity = sampleQueryGapsSetDensities[3] + sampleResultGapsSetDensities[3] + sampleMismatchesSetDensities[3];

    average_densities_quotients =          {
                                            totalSiteDensity * 1.0 / lengths.size(),
                                            totalSiteDensity / totalLeftDensity,
                                            sampleQueryGapsSetDensities[1] / sampleQueryGapsSetDensities[0],
                                            sampleResultGapsSetDensities[1] / sampleResultGapsSetDensities[0],
                                            sampleMismatchesSetDensities[1] / sampleMismatchesSetDensities[0],
                                            totalSiteDensity / totalRightDensity,
                                            sampleQueryGapsSetDensities[1]/ sampleQueryGapsSetDensities[2],
                                            sampleResultGapsSetDensities[1] / sampleResultGapsSetDensities[2],
                                            sampleMismatchesSetDensities[1]/ sampleMismatchesSetDensities[2],  
                                            totalSiteDensity / totalAdjacentDensity,
                                            sampleQueryGapsSetDensities[1] / sampleQueryGapsSetDensities[3],
                                            sampleResultGapsSetDensities[1] / sampleResultGapsSetDensities[3],
                                            sampleMismatchesSetDensities[1] / sampleMismatchesSetDensities[3]                                          
                                            };
    return average_densities_quotients;
}

unsigned int sumAllErrors(vector<unsigned int> errors)
{
    unsigned int sum = 0;
    for (unsigned int i = 0; i < errors.size(); ++i) {
        sum += errors[i];
    }

    return sum;
}

void runMonteCarloMethod (vector<unsigned int> queryGaps, vector<unsigned int> resultGaps, vector<unsigned int> mismatches, 
                                    vector<unsigned int> lengths, 
                                    vector< pair<unsigned int, unsigned int> > sites,
                                    vector<double> realQuotients,
                                    ofstream &out
)
{
    unsigned int queryGapsCount = sumAllErrors(queryGaps);
    unsigned int resultGapsCount = sumAllErrors(resultGaps);
    unsigned int mismatchesCount = sumAllErrors(mismatches);

    vector< vector<double> > quotientsVectors = {{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}};
    for (unsigned int i = 0; i < 10000; ++i) {
        //cout << "single start" << endl;
        vector<double> average_densities_quotients = 
            runSingleMonteCarloExperiment(queryGapsCount, resultGapsCount, mismatchesCount, lengths, sites);
            //runSingleSequenceMonteCarloExperiment(queryGaps, resultGaps, mismatches, lengths, sites);

        //cout << "single end" << endl;
        for (unsigned int j = 0; j < average_densities_quotients.size(); ++j) {
            quotientsVectors[j].push_back(average_densities_quotients[j]);
        }
    }

    cout << 1 << endl;

    #if 0
    for (auto it = siteToRightQuotients.begin(); it != siteToRightQuotients.end(); ++it) {
        cout << *it << " ";
    }
    #endif

    char delim_string = static_cast<char>(SpecialCharacter::TableDelim);

    for (unsigned int i = 0; i < quotientsVectors.size(); ++i) {
        out << realQuotients[i] << delim_string;
        out << getLeftMonteCarloValue(quotientsVectors[i], realQuotients[i]) << delim_string;
        out << getRightMonteCarloValue(quotientsVectors[i], realQuotients[i]) << delim_string;
    }

    out << endl;  
/*
    cout << "(left/site) " << realQuotients[0] << " left count:" << getLeftMonteCarloCount(leftToSiteQuotients, realQuotients[0]);
    cout << " right count:" << getRightMonteCarloCount(leftToSiteQuotients, realQuotients[0]) << endl;
    
    cout << "(right/site) " << realQuotients[1] << " left count:" << getLeftMonteCarloCount(rightToSiteQuotients, realQuotients[1]);
    cout << " right count:" << getRightMonteCarloCount(rightToSiteQuotients, realQuotients[1]) << endl;
*/
}

vector<double> processSamplePart(string dataFileName, SamplePartData & samplePartData, ofstream &sampleDensitiesOutFile)
{
    ifstream dataFile(dataFileName);

    string inFileLine;
    string outFileNameStart = "stats_";
    AvData avData;

    if (dataFile.is_open()) {
        // Получаем строку с именем файла
        getline(dataFile, inFileLine);
        outFileNameStart += inFileLine.substr(inFileLine.find_first_of('/') + 1, inFileLine.find_last_of('/') - inFileLine.find_first_of('/') - 1);
        string outFileName = outFileNameStart + ".txt";
        ofstream outFile;
        outFile.open(outFileName);

        // Заводим локальные переменные для сбора статистики
        int count;
        string queryId;
        string seqId;
        string query;
        string sequence;
        string evalue;
        string strand;
        bool isBestResult = true;

        // Проходимся по файлу для сбора статистики
        while (getline(dataFile, inFileLine)) {
            if (inFileLine[0] == '>') {
                queryId = getIdFromFastaLine(inFileLine);
                isBestResult = true;
            } else if (isBestResult) {
                if (inFileLine.find("_Count: ") == 0) {
                    count = stoi(removeSpaces(inFileLine.substr(inFileLine.find("_Count:") + 7)));
                } else if (inFileLine.find("ID:") == 0) {
                    seqId = removeSpaces(inFileLine.substr(inFileLine.find("ID:") + 3));
                } else if (inFileLine.find("_Blast aligment") == 0) {
                    getline(dataFile, query);
                    getline(dataFile, sequence);
                    processAligments(queryId, seqId, query, sequence, strand, 1, evalue, outFile, avData, samplePartData);
                    isBestResult = false;
                } else if (inFileLine.find("STRAND = ") == 0) {
                    strand = removeSpaces(inFileLine.substr(inFileLine.find("STRAND = ") + 9));
                } else if (inFileLine.find("_Evalue: ") == 0) {
                    evalue = removeSpaces(inFileLine.substr(inFileLine.find("_Evalue: ") + 9));
                }
            }
        }
        dataFile.close();
        outFile.close();
    }

    sampleDensitiesOutFile << avData.getPrintableTotalAvData() << endl; 
    return avData.getDensitiesQuotients();
}

int main (int argc, char** argv) {
    if (argc != 3) {
        cout << "Ошибка: Программе требуется 2 аргумента - идентификатор выборки, имя файла вывода" << endl;
        return -1;
    }

    #if 0
    //aa-tT--CG-T-agC {{4, 10}, {14, 14}}
    //--TaTTTaT-Tg--  {{2, 2}, {4, 6}, {8, 10}}
    vector< pair<unsigned int, unsigned int> > meow = getBindingSitesPositions("--TaTTTaT-Tg--");
    for (auto it = meow.begin(); it != meow.end(); ++it) {
        pair<unsigned int, unsigned int> site = *it;
        cout << site.first << ":" << site.second << endl;
    }
    cout << endl;
    #endif

    char tableDelimeter = static_cast<char>(SpecialCharacter::TableDelim);

    SamplePartData hs_hs {0, 0, 0, {}, {}, {}, {}, {}};
    SamplePartData hs_mm {0, 0, 0, {}, {}, {}, {}, {}};
    SamplePartData mm_hs {0, 0, 0, {}, {}, {}, {}, {}};
    SamplePartData mm_mm {0, 0, 0, {}, {}, {}, {}, {}};

    string sampleId = argv[1];
    string dataFileName = argv[2];

    ofstream sampleMonteCarloOutFile;
    sampleMonteCarloOutFile.open("montecarlo_" + sampleId + ".txt");

    ofstream sampleDensitiesOutFile;
    sampleDensitiesOutFile.open("densities_" + sampleId + ".txt");

    cout << "hshs" << endl;
    vector<double> hshsQuotients = processSamplePart(sampleId + "/Hs_" + sampleId + "_Hs/" + dataFileName, hs_hs, sampleDensitiesOutFile);
    cout << "mc" << endl;
    runMonteCarloMethod (hs_hs.queryGapsCounts, hs_hs.resultGapsCounts, hs_hs.mismatchesCounts, 
                                    hs_hs.lengths, 
                                    hs_hs.sites,
                                    hshsQuotients,
                                    sampleMonteCarloOutFile);
    cout << endl;
    cout << "hsmm" << endl;
    vector<double> hsmmQuotients = processSamplePart(sampleId + "/Hs_" + sampleId + "_Mm/" + dataFileName, hs_mm, sampleDensitiesOutFile);
    runMonteCarloMethod (hs_mm.queryGapsCounts, hs_mm.resultGapsCounts, hs_mm.mismatchesCounts, 
                                    hs_mm.lengths, 
                                    hs_mm.sites,
                                    hsmmQuotients,
                                    sampleMonteCarloOutFile);
    cout << endl;
    cout << "mmhs" << endl;
    vector<double> mmhsQuotients = processSamplePart(sampleId + "/Mm_" + sampleId + "_Hs/" + dataFileName, mm_hs, sampleDensitiesOutFile);
    runMonteCarloMethod (mm_hs.queryGapsCounts,  mm_hs.resultGapsCounts,  mm_hs.mismatchesCounts, 
                                    mm_hs.lengths, 
                                    mm_hs.sites,
                                    mmhsQuotients,
                                    sampleMonteCarloOutFile);
    cout << endl;
    cout << "mmmm" << endl;
    vector<double> mmmmQuotients = processSamplePart(sampleId + "/Mm_" + sampleId + "_Mm/" + dataFileName, mm_mm, sampleDensitiesOutFile);
    runMonteCarloMethod (mm_mm.queryGapsCounts, mm_mm.resultGapsCounts, mm_mm.mismatchesCounts, 
                                    mm_mm.lengths, 
                                    mm_mm.sites,
                                    mmmmQuotients,
                                    sampleMonteCarloOutFile);

    sampleDensitiesOutFile.close();
    return 0;
}