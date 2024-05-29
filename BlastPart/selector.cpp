#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>

using namespace std;

int main (int argc, char** argv) {
    if (argc != 5) {
        return -1;
    } 

    string tf = argv[1];
    string org = argv[2];
    string fasta_input_file_name = argv[3];
    string output_file_name = argv[4];

    string line;
    ifstream fasta_in(fasta_input_file_name);
    ofstream sample_out;
    sample_out.open(output_file_name);

    bool isSelected = false;
    string trrdSubstring = "TRRD;";
    string organism;
    string dataLine;

    if (fasta_in.is_open() && sample_out.is_open()) {
        while (getline(fasta_in, line)) {
            if (line[0] == '>') {
                unsigned int trrdPos = line.find(trrdSubstring);
                if (trrdPos + trrdSubstring.size() > line.size()) {
                    cout << line << endl;
                }
                organism = line.substr(trrdPos + trrdSubstring.size(), 2);
                unsigned int tfIdPos = line.find(tf);

                if (organism == org && tfIdPos < line.size()) {
                    isSelected = true;
                    sample_out << line << endl;
                }
            } else if (isSelected) {
                sample_out << line << endl;
                isSelected = false;
            }
        }
        
        fasta_in.close();
        sample_out.close();
    } else {
        cout << "[selector: " << output_file_name << "] cannot open file(s)" << endl;
        return false;
    }

    cout << "[selector: " << output_file_name << "] done" << endl;

    return 0;
}