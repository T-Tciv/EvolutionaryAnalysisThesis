#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>

using namespace std;

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
        while (position < length && !isalpha(sequence[position])) {
            ++position;
        }

        while (position < length && islower(sequence[position])) {
            ++position;
        }

        currentLeft = position;

        while (position < length && isupper(sequence[position])) {
            ++position;
        }

        currentRight = position - 1;

        if (currentLeft <= currentRight) {
            positions.push_back( pair<unsigned int, unsigned int>(currentLeft, currentRight) );
        }      
    }

    return positions;
}

int main (int argc, char** argv) {
    string line;
    string sampleId;
    string prevField = "";
    string infoLine = "";
    string dataLine = "";
    ifstream trrd_in("in.txt");
    ofstream fasta_out;
    fasta_out.open("input.fst");

    if (trrd_in.is_open() && fasta_out.is_open()) {
        while (getline(trrd_in, line)) {
            line = removeSpaces(line);
            string line_type = line.substr(0, 2);

            if (line_type == "NM") {
                if (infoLine.size() > 0) {
                    infoLine += "ID: " + sampleId;
                    fasta_out << infoLine << endl;
                    fasta_out << dataLine << endl;
                }  

                infoLine = "> ";
                infoLine += line.substr(2) += " ";
            } else if (line_type == "BF") {
                infoLine += "AC:" + line.substr(2) += " ";
            } else if (line_type == "TR") {
                infoLine += line.substr(2) += " ";
            } else if (line_type == "SQ") {
                dataLine = line.substr(2);
            } else if (line_type == "ID") {
                if (infoLine.size() > 0) {
                    infoLine += "ID: " + sampleId;
                    fasta_out << infoLine << endl;
                    fasta_out << dataLine << endl;
                    infoLine = "";
                }

                sampleId = line.substr(2);
            } else if (line_type == "AC" || prevField == "DT") {
                sampleId += "_" + line.substr(2);
            } else if (line_type == "AS" && prevField == "DT") {
                //
            }

            prevField = line_type;
        }

        infoLine += "ID: " + sampleId;
        fasta_out << infoLine << endl;
        fasta_out << dataLine << endl;
        
        trrd_in.close();
        fasta_out.close();
    }

    cout << "parser: input.fst done" << endl;

    return 0;
}