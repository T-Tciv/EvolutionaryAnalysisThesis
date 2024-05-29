#include <ncbi_pch.hpp>
#include <corelib/ncbiapp.hpp>
#include <corelib/ncbienv.hpp>
#include <corelib/ncbiargs.hpp>
#include <objmgr/object_manager.hpp>
#include <objects/seqalign/Seq_align_set.hpp>
#include <algo/blast/api/sseqloc.hpp>
#include <algo/blast/api/local_blast.hpp>
#include <algo/blast/api/uniform_search.hpp>
#include <algo/blast/api/blast_types.hpp>
#include <algo/blast/api/blast_aux.hpp>
#include <algo/blast/api/objmgr_query_data.hpp>
#include <algo/blast/api/blast_options.hpp>
#include <algo/blast/api/blast_options_handle.hpp>
#include <algo/blast/api/blast_nucl_options.hpp>
#include <algo/blast/api/remote_blast.hpp>
#include <algo/blast/blastinput/blast_input.hpp>
#include <algo/blast/blastinput/blast_fasta_input.hpp>

#include <corelib/ncbistd.hpp>
#include <objects/seq/Bioseq.hpp>
#include <objects/seqloc/Seq_id.hpp>
#include <objects/seqloc/Seq_loc.hpp>
#include <objects/seqloc/Seq_interval.hpp>
#include <objects/seq/Seq_inst.hpp>
#include <objmgr/scope.hpp>
#include <objmgr/seq_vector.hpp>
#include <objmgr/seqdesc_ci.hpp>
#include <objmgr/feat_ci.hpp>
#include <objmgr/align_ci.hpp>
#include <objtools/data_loaders/genbank/gbloader.hpp>

#include <objects/seqalign/Seq_align.hpp>

#include <algo/align/nw/nw_band_aligner.hpp>
#include <algo/align/nw/mm_aligner.hpp>
#include <algo/align/nw/nw_formatter.hpp>
#include <objects/general/Object_id.hpp>

#include <iostream>
#include <fstream>
#include <string>

// объявляет использование пространства имён NCBI
USING_NCBI_SCOPE;
// объявляет использование пространства имён BLAST
USING_SCOPE(blast);
USING_SCOPE(objects);

enum class SpecialCharacter : char
{ Gap = '-', SeqEnd = '+', FileExt = '.', FileSystemDiv = '/', FileNameDiv = '_' };

enum class SequenceStatus : int
{ Ok = 0, OneSide = 1, OnlySite = 2, NoSite = 3};

// Новое приложение пишется путем наследования класса от CNcbiApplication и реализации его методов Init(), Run() и Exit()
// Управление выполнением новому приложению передается путем вызова метода AppMain() объекта приложения
// Метод AppMain() вызывает методы Init(), Run() и Exit()
class BlastApp : public CNcbiApplication
{
private:
    virtual void Init(void);
    virtual int  Run(void);
    virtual void Exit(void);

    string file_name_start;

    string exampleFunction();
    pair<unsigned int, unsigned int> GetBindingSitePositions(string query);
    vector<string> ProcessQueriesData(string fasta_file_name);
    void SetDefaultNucleotideOptions(string program, CBlastNucleotideOptionsHandle* opts_handle);
    void SetOptions(CRef<CBlastOptionsHandle> opts_handle);
    string GetGappedQuery(string query, vector<unsigned int> gaps);
    string GetDescriptionByGi(TGi gi, int start_pos, int start_diff, int end_diff,
                           unsigned int align_length, int left_indent, 
                           int right_indent, string strand, ofstream &out);
    string GetComplementarySequence(string sequence);
    string GetDataByGi(TGi gi, unsigned int start_pos, unsigned int align_length, int left_indent, int right_indent, string strand);
    string GetGappedResult(TGi gi, int start_pos, unsigned int align_length,
                           int left_indent, 
                           int right_indent, string strand,
                           vector<unsigned int> gaps);
    vector<unsigned int> GetMismatches(string query, string result);
    void PrintAligment(string query, string result);
    void PrintWindowResults(string query, string result);
    void PrintCommonDensityResult(string id, string result_acc, string query, string result, string strand, int hit_number, double evalue, ofstream &out);
    void PrintDensityResult(string id, string result_acc, string query, string result, string strand, int hit_number, double evalue, ofstream &out);

    void PrintBlastSettings(CRef<CBlastOptionsHandle> opts_handle);
};

void BlastApp::Init(void)
{
    // Создание указателя на объект класса описаний аргументов коммандной строки
    // Обычно требуется создать объект CArgDescriptions перед любой другой инициализацией в функции Init()
    // В противном случае CNcbiApplication создает описание по умолчанию
    unique_ptr<CArgDescriptions> arg_desc(new CArgDescriptions);

    // Задаём описание программы
    arg_desc->SetUsageContext(GetArguments().GetProgramBasename(), "BLAST program");

    // Добавление обязательного аргумента "program" и ограничение его (выбор одного из видов BLAST)
    arg_desc->AddKey("program", "NucleotideProgramName", "One of blastn, megablast", CArgDescriptions::eString);
    arg_desc->SetConstraint("program", &(*new CArgAllow_Strings,"blastn", "megablast"));

    // Добавление обязательного аргумента с именем входного файла
    arg_desc->AddKey("infile", "Queryfile", "A FASTA file with the query", CArgDescriptions::eInputFile);

    // Добавление аргумента с указанием пути к директории для выходных файлов (по умолчанию -- текущая директория)
    arg_desc->AddDefaultKey("outDir", "outDir", "Output directory", CArgDescriptions::eString, "");

    // Добавление аргумента с ограничениями через entrez
    arg_desc->AddDefaultKey("entrez", "EntrezQuery", "This is the string with entrez query", CArgDescriptions::eString, "");

    // Добавляем отступы от концов найденной последовательности-совпадения
    arg_desc->AddDefaultKey("leftIndent", "leftIndent", "Indent from the start of alignment", CArgDescriptions::eInteger, "0");
    arg_desc->AddDefaultKey("rightIndent", "rightIndent", "Indent from the end of alignment", CArgDescriptions::eInteger, "0");

    // Выбор базы данных
    // refseq_genomes
    arg_desc->AddDefaultKey("db", "DataBase", "This is the name of the database", CArgDescriptions::eString, "nt");

    // "Expect" (E) value - параметр, определяющий ожидаемый уровень совпадения последовательностей
    arg_desc->AddDefaultKey("evalue", "evalue", "E-value threshold for saving hits", CArgDescriptions::eDouble, "0.05");

    // Установка значений penalty и reward
    arg_desc->AddDefaultKey("penalty", "penalty", "Penalty score for a mismatch", CArgDescriptions::eInteger, "-3");
    arg_desc->AddDefaultKey("reward", "reward", "Reward score for a match", CArgDescriptions::eInteger, "2");

    // Установка цены за открытие и продолжение gap
    arg_desc->AddDefaultKey("gapOpenCost", "gapOpenCost", "Gap opening cost", CArgDescriptions::eInteger, "5");
    arg_desc->AddDefaultKey("gapExtCost", "gapExtCost", "Gap extension cost", CArgDescriptions::eInteger, "2");

    // Максимальное количество выводимых совпадений
    arg_desc->AddDefaultKey("hitsize", "hitsize", "Hitlist size", CArgDescriptions::eInteger, "1");

    // "Освобождение" unique_ptr arg_desc от "ответственности" за его объект 
    // и передача указателя на этот объект в SetupArgDescriptions()
    // т.е. установка описаний аргументов в данном приложении
    SetupArgDescriptions(arg_desc.release());
}

pair<unsigned int, unsigned int> BlastApp::GetBindingSitePositions(string query) {
    unsigned int length = query.size();

    unsigned int left = 0;
    unsigned int right = length - 1;

    for (left = 0; left < length; left++) {
        if (isupper(query[left])) {
            break;
        }
    }

    for (unsigned int i {left}; i < length; i++) {
        if (islower(query[i])) {
            right = i - 1;
            break;
        }
    }

    for (unsigned int i {right + 1}; i < length; i++) {
        if (isupper(query[i])) {
            left = length;
            right = 0;
            break;
        }
    }

    return pair<unsigned int, unsigned int>(left, right);
}

void BlastApp::PrintBlastSettings(CRef<CBlastOptionsHandle> opts_handle) {
    const CArgs& args = GetArgs();

    stringstream args_stream;
    string infile = args["infile"].AsString();

    args_stream << "prog:" << args["program"].AsString() << " ";
    args_stream << "infile:" << infile << " ";
    args_stream << "entrez:" << args["entrez"].AsString() << " ";
    args_stream << "left_indent:" << args["leftIndent"].AsInteger() << " ";
    args_stream << "right_indent:" << args["rightIndent"].AsString() << " ";
    args_stream << "db:" << args["db"].AsString() << " ";

    args_stream << "evalue:" << opts_handle->GetEvalueThreshold() << " ";
    args_stream << "hitsize:" << opts_handle->GetHitlistSize() << " ";

    if (CBlastNucleotideOptionsHandle* nucl_options_handle = dynamic_cast<CBlastNucleotideOptionsHandle*>(&*opts_handle)) {
        args_stream << "word_size:" << nucl_options_handle->GetWordSize() << " ";
        args_stream << "penalty:" << nucl_options_handle->GetMismatchPenalty() << " ";
        args_stream << "reward:" << nucl_options_handle->GetMatchReward() << " ";
        args_stream << "gap_op_cost:" << nucl_options_handle->GetGapOpeningCost() << " ";
        args_stream << "gap_ext_cost:" << nucl_options_handle->GetGapExtensionCost() << " ";
        args_stream << "mask_for_lt:" << nucl_options_handle->GetMaskAtHash() << " ";
        args_stream << "low_complex_dust_filter_enable:" << nucl_options_handle->GetDustFiltering() << " ";
        args_stream << "complex_adj_mode:" << nucl_options_handle->GetComplexityAdjMode() << " ";
        args_stream << "is_repeat_filtering:" << nucl_options_handle->GetRepeatFiltering() << " ";
        args_stream << "matrix_name:" << nucl_options_handle->GetMatrixName() << " ";
    }
    
    args_stream << endl;
    
    string args_line = args_stream.str();

    ofstream out;
    string settings_file_name = file_name_start + "_settings.txt";
    out.open(settings_file_name);
    if (out.is_open()) {
        out << "-----------------" << endl;
        out << args_line << endl;
        cout << args_line << endl;
        out << "-----------------" << endl;
    }

    out.close();    
}

void BlastApp::SetDefaultNucleotideOptions(string program, CBlastNucleotideOptionsHandle* nucl_options_handle) {
    if (program == "blastn") {
        nucl_options_handle->SetTraditionalBlastnDefaults();
        nucl_options_handle->SetEvalueThreshold(0.05);
    }

    if (program == "megablast") {
        nucl_options_handle->SetTraditionalMegablastDefaults();
    }
}

void BlastApp::SetOptions(CRef<CBlastOptionsHandle> opts_handle)
{
    // Получение аргументов, установленных в приложение с помощью SetupArgDescriptions()
    const CArgs& args = GetArgs();

    // Приведение к типу нуклеотидных опций
    // (если тип не совпадает с ожидаемым, dynamic_cast вернёт nullptr) 
    if (CBlastNucleotideOptionsHandle* nucl_options_handle = dynamic_cast<CBlastNucleotideOptionsHandle*>(&*opts_handle)) {
        SetDefaultNucleotideOptions(args["program"].AsString(), nucl_options_handle);
        
        if (args["evalue"].AsDouble())
            nucl_options_handle->SetEvalueThreshold(args["evalue"].AsDouble());
        
        if (args["hitsize"].AsInteger() && args["hitsize"].AsInteger() > 0)
            nucl_options_handle->SetHitlistSize(args["hitsize"].AsInteger());

        if (args["reward"].AsInteger())
            nucl_options_handle->SetMatchReward(args["reward"].AsInteger());
            
        if (args["penalty"].AsInteger())
            nucl_options_handle->SetMismatchPenalty(args["penalty"].AsInteger());

        if (args["gapOpenCost"].AsInteger())
            nucl_options_handle->SetGapOpeningCost(args["gapOpenCost"].AsInteger());

        if (args["gapExtCost"].AsInteger())
            nucl_options_handle->SetGapExtensionCost(args["gapExtCost"].AsInteger());
    }

    PrintBlastSettings(opts_handle);
}

vector<string> BlastApp::ProcessQueriesData(string fasta_file_name) {
    string line;
    string query = "";
    vector<string> queries_data;
 
    ifstream in(fasta_file_name);
    if (in.is_open())
    {
        while (getline(in, line))
        {
            if (line[0] == '>') {
                queries_data.push_back(line);
                break;
            } 
            //cout << 1 << endl;
        }

        while (getline(in, line))
        {
            if (line[0] == '>') {
                queries_data.push_back(query);
                query = "";
                queries_data.push_back(line);
            } else {
                if (isspace(line[line.size() - 1])) {
                    line = line.substr(0, line.size() - 1);
                } 
                
                query += line;
            }
        }

        queries_data.push_back(query);
    }
    in.close(); 

    return queries_data;
}

string BlastApp::GetGappedQuery(string query, vector<unsigned int> gaps)
{
    string gapped_result_sequence = query;
    unsigned int gaps_count = 0;
    sort(gaps.begin(), gaps.end()); 

    for (unsigned int gap : gaps) {
        gapped_result_sequence.insert(gap + gaps_count, 1, static_cast<char>( SpecialCharacter::Gap ));
        ++gaps_count;
    }

    return gapped_result_sequence;
}

string BlastApp::GetComplementarySequence(string sequence) {
    string result;

    for (char nucl : sequence) {
        char newNucl;

        switch (nucl) {
            case 'a': newNucl = 't';
            break;
            case 't': newNucl = 'a';
            break;
            case 'c': newNucl = 'g';
            break;
            case 'g': newNucl = 'c';
            break;
            case 'A': newNucl = 'T';
            break;
            case 'T': newNucl = 'A';
            break;
            case 'C': newNucl = 'G';
            break;
            case 'G': newNucl = 'C';
            break;
            default: newNucl = nucl;
            break;
        }

        result += newNucl;
    }

    reverse(result.begin(), result.end());

    return result;
}

string BlastApp::GetDescriptionByGi(TGi gi, int start_pos, int start_diff, int end_diff, 
                           unsigned int align_length, int left_indent, 
                           int right_indent, string strand, ofstream &out)
{
    stringstream sequenceDataStream;
    // Создаём object manager
    CRef<CObjectManager> object_manager = CObjectManager::GetInstance();

    // Создаём GenBank data loader и связываем его с OM.
    // GenBank loader автоматически помечается default loader
    CGBDataLoader::RegisterInObjectManager(*object_manager);

    CScope scope(*object_manager);
    //  Добавляем дефолтные data loaders (т.е. GenBank data loader) в эту scope
    scope.AddDefaults();

    CSeq_id seq_id(CSeq_id::E_Choice::e_Gi , gi);

    // Получаем Bioseq handle для seq_id.
    CBioseq_Handle bioseq_handle = scope.GetBioseqHandle(seq_id);

    if ( !bioseq_handle ) {
        ERR_POST(Fatal << "Bioseq not found, with GI=" << gi);
    }

    out << "-----------------" << endl;
    out << "GI" << gi <<  endl;

    // Получаем последовательность, используя CSeqVector.
    CSeqVector seq_vector = bioseq_handle.GetSeqVector(CBioseq_Handle::eCoding_Iupac);

    // const CBioseq_set::TDescr & description = bioseq_handle.GetDescr();
    list description = bioseq_handle.GetDescr().Get();

    sequenceDataStream << "DESCRIPTION: " << endl;

    for (CRef<CSeqdesc> descr_elem : description) {
        if (!!descr_elem->IsTitle()) {
            sequenceDataStream << MSerial_AsnText << descr_elem->GetTitle() << endl;
        } else if (!!descr_elem->IsName()) {
            sequenceDataStream << MSerial_AsnText << descr_elem->GetName() << endl;
        }

        out << MSerial_AsnText << descr_elem <<  endl;
    }

    string label = scope.GetLabel(seq_id);
    sequenceDataStream << "ID: " << label << endl;
    sequenceDataStream << "GI: " << gi << endl;

    unsigned int seq_length = seq_vector.size();
    sequenceDataStream << "LENGTH = " << seq_length << endl;
    sequenceDataStream << "STRAND = " << strand << endl;

    if (strand == "MINUS") {
        int temp = start_diff;
        start_diff = end_diff;
        end_diff = temp;
    }

    sequenceDataStream << "ALIGMENT START = " << start_pos << "+" << start_diff << endl;
    sequenceDataStream << "ALIGMENT END = " << start_pos + align_length - 1 << "-" << end_diff << endl;

    int left = max((int) start_pos - left_indent, 0);
    int right = min((int) start_pos + (int) align_length + right_indent, (int) seq_length);

    sequenceDataStream << "WITH INDENTS: START = " << left << ", END = " << right - 1 << endl;

    string desc = sequenceDataStream.str();

    cout << desc << endl;
    out << "-----------------" << endl;

    return label;
}

string BlastApp::GetDataByGi(TGi gi, unsigned int start_pos, unsigned int align_length, int left_indent, int right_indent, string strand)
{
    // Создаём object manager
    CRef<CObjectManager> object_manager = CObjectManager::GetInstance();

    // Создаём GenBank data loader и связываем его с OM.
    // GenBank loader автоматически помечается default loader
    CGBDataLoader::RegisterInObjectManager(*object_manager);

    CScope scope(*object_manager);
    //  Добавляем дефолтные data loaders (т.е. GenBank data loader) в эту scope
    scope.AddDefaults();

    CSeq_id seq_id(CSeq_id::E_Choice::e_Gi , gi);

    // Получаем Bioseq handle для seq_id.
    CBioseq_Handle bioseq_handle = scope.GetBioseqHandle(seq_id);

    if ( !bioseq_handle ) {
        ERR_POST(Fatal << "Bioseq not found, with GI=" << gi);
    }

    // Получаем последовательность, используя CSeqVector.
    CSeqVector seq_vector = bioseq_handle.GetSeqVector(CBioseq_Handle::eCoding_Iupac);

    unsigned int seq_length = seq_vector.size();

    int left = max((int) start_pos - left_indent, 0);
    int right = min((int) start_pos + (int) align_length + right_indent, (int) seq_length);

    string sequence;
    seq_vector.GetSeqData(left, right, sequence);
    return sequence;
}

string BlastApp::GetGappedResult(TGi gi, int start_pos, unsigned int align_length,
                           int left_indent, 
                           int right_indent, string strand,
                           vector<unsigned int> gaps)
{
    string gapped_result_sequence = "";
    int start_shift = 0;

    if (start_pos < 0) {
        align_length += start_pos;
        start_shift = -start_pos;
        start_pos = 0;
    }

    gapped_result_sequence.append(GetDataByGi(gi, static_cast<unsigned int> (start_pos), align_length, left_indent, right_indent, strand));

    unsigned int gaps_count = 0;
    sort(gaps.begin(), gaps.end()); 

    for (unsigned int gap : gaps) {
        gapped_result_sequence.insert(gap + gaps_count, 1, static_cast<char>( SpecialCharacter::Gap ));
        ++gaps_count;
    }

    if (strand == "MINUS") {
        gapped_result_sequence = GetComplementarySequence(gapped_result_sequence); 
    }

    // В случае, если на последний символ найденной последовательности приходится не последний символ последовательности-запроса,
    // мы для выраванивания в выводе дозаполняем последовательность специальными символами
    if (static_cast<int> (align_length) - static_cast<int> (gapped_result_sequence.size()) > 0) {
        unsigned int end_shift = align_length - gapped_result_sequence.size();

        gapped_result_sequence.append(static_cast<unsigned int>(end_shift), static_cast<char>( SpecialCharacter::SeqEnd ));
    }

    // В случае, если на нулевой символ найденной последовательности приходится ненулевой символ последовательности-запроса,
    // мы для выраванивания в выводе дозаполняем последовательность специальными символами
    gapped_result_sequence.insert(0, static_cast<unsigned int>(start_shift), static_cast<char>( SpecialCharacter::SeqEnd ));
/*
    if (strand == "MINUS") {
        gapped_result_sequence = GetComplementarySequence(gapped_result_sequence); 
    }
*/
    return gapped_result_sequence;
}

void BlastApp::PrintAligment(string query, string result) {
    unsigned int length = query.size();
    unsigned int row_length = 60;
    unsigned int rows_number = length / row_length;
    bool isNotDivisible = (length % row_length != 0);

    for (unsigned int i {0}; i < rows_number; i++) { 
        cout << query.substr(i * row_length, row_length) << endl;
        cout << result.substr(i * row_length, row_length) << "\n" << endl;
    }

    if (isNotDivisible) {
        cout << query.substr(rows_number * row_length) << endl;
        cout << result.substr(rows_number * row_length) << endl;
    }
}

void BlastApp::PrintWindowResults(string query, string result) {
    unsigned int length = query.size();
    const unsigned int step = 10; 

    pair<unsigned int, unsigned int> binding_site_positions = GetBindingSitePositions(query);

    unsigned int binding_site_left = binding_site_positions.first;
    unsigned int binding_site_right = binding_site_positions.second;

    cout << "Binding site left: " << binding_site_left << endl;
    cout << "Binding site right: " << binding_site_right << endl;

    PrintAligment(query, result);

    int mismatches_window[step] {};
    int query_gaps_window[step] {};
    int result_gaps_window[step] {};

    int mismatches = 0;
    int query_gaps = 0;
    int result_gaps = 0;

    vector<int> mismatches_count;
    vector<int> query_gaps_count;
    vector<int> result_gaps_count;

    char gap = static_cast<char>( SpecialCharacter::Gap );

    for (unsigned int i {0}; i < step; i++) { 
        if (query[i] == gap) {
            query_gaps_window[i % step] = 1;
        }

        if (result[i] == gap) {
            result_gaps_window[i % step] = 1;
        }

        if (query[i] != gap && result[i] != gap && toupper(query[i]) != toupper(result[i])) {
           mismatches_window[i % step] = 1;
        } 
    }

    for (unsigned int i {step}; i < length; i++) { 
        for (unsigned int j {0}; j < step; j++) {
            mismatches += mismatches_window[j];
            query_gaps += query_gaps_window[j];
            result_gaps += result_gaps_window[j];
        }

        mismatches_count.push_back(mismatches);
        query_gaps_count.push_back(query_gaps);
        result_gaps_count.push_back(result_gaps);

        mismatches = 0;
        query_gaps = 0;
        result_gaps = 0;

        mismatches_window[i % step] = 0;
        query_gaps_window[i % step] = 0;
        result_gaps_window[i % step] = 0;

        if (query[i] == gap) {
            query_gaps_window[i % step] = 1;
        }

        if (result[i] == gap) {
            result_gaps_window[i % step] = 1;
        }

        if (query[i] != gap && result[i] != gap && toupper(query[i]) != toupper(result[i])) {
           mismatches_window[i % step] = 1;
        }  
    }

    for (unsigned int j {0}; j < step; j++) {
        mismatches += mismatches_window[j];
        query_gaps += query_gaps_window[j];
        result_gaps += result_gaps_window[j];
    }

    mismatches_count.push_back(mismatches);
    query_gaps_count.push_back(query_gaps);
    result_gaps_count.push_back(result_gaps);

    for (unsigned int i {step - 1}; i < length; i++) { 
        cout << query[i] << " ";
    }

    cout << endl;

    for (unsigned int i {step - 1}; i < length; i++) { 
        if (i <= binding_site_left - step  || i >= binding_site_right + step) {
            cout << "_ ";
        } else if (i >= binding_site_left && i <= binding_site_right) {
            cout << "* ";
        } else {
            cout << "- ";
        }
    }

    cout << endl;

    cout << "Mismatches: " << endl;

    for (int count : mismatches_count) { 
        cout << count << " ";
    }

    cout << endl;
    cout << "Query gaps: " << endl;

    for (int count : query_gaps_count) { 
        cout << count << " ";
    }

    cout << endl;
    cout << "Result gaps: " << endl;

    for (int count : result_gaps_count) { 
        cout << count << " ";
    }

    cout << endl;
}

void BlastApp::PrintCommonDensityResult(string id, string result_acc, string query, string result, string strand, int hit_number, double evalue, ofstream &out) {
    PrintAligment(query, result);

    unsigned int length = query.size();

    pair<unsigned int, unsigned int> binding_site_positions = GetBindingSitePositions(query);

    stringstream string_out;
    string delim_string = "\t";
    string_out << id << delim_string << result_acc << delim_string << strand << delim_string << evalue << delim_string;
    int seq_status;

    unsigned int binding_site_left = binding_site_positions.first;
    unsigned int binding_site_right = binding_site_positions.second;

    if (binding_site_left > binding_site_right) {   
        seq_status = static_cast<int>( SequenceStatus::NoSite );
        string_out << seq_status << delim_string;
        string_out << "-" << delim_string << "-" << delim_string << "-" << delim_string;
        string_out << "-" << delim_string << "-" << delim_string << "-" << delim_string;
        string_out << "-" << delim_string << "-" << delim_string << "-" << delim_string << hit_number << endl;

        out << string_out.str();
        cout << string_out.str();
        return;
    }

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
    unsigned int left_count = 0;
    unsigned int site_count = 0;
    unsigned int right_count = 0;

    if (binding_site_left >= length) {
        left_length = 0;
    }

    if (binding_site_right >= length) {
        right_length = 0;
    }

    for (unsigned int i {0}; i < binding_site_left; i++) { 
        query_nucl = toupper(query[i]);
        result_nucl = toupper(result[i]);

        if (query_nucl == gap || result_nucl == gap || query_nucl != result_nucl) {
            ++left_count;
        }
    }

    for (unsigned int i {binding_site_left}; i <= binding_site_right; i++) { 
        query_nucl = toupper(query[i]);
        result_nucl = toupper(result[i]);

        if (query_nucl == gap || result_nucl == gap || query_nucl != result_nucl) {
            ++site_count;
        }
    }

    for (unsigned int i {binding_site_right + 1}; i < length; i++) { 
        query_nucl = toupper(query[i]);
        result_nucl = toupper(result[i]);

        if (query_nucl == gap || result_nucl == gap || query_nucl != result_nucl) {
            ++right_count;
        }
    }

    double left_density = 1.0 * left_count / left_length;
    double site_density = 1.0 * site_count / site_length;
    double right_density = 1.0 * right_count / right_length;

    string_out << seq_status << delim_string;
    string_out << left_count << delim_string << left_length << delim_string << left_density << delim_string;
    string_out << site_count << delim_string << site_length << delim_string << site_density << delim_string;
    string_out << right_count << delim_string << right_length << delim_string << right_density << delim_string << hit_number << endl;

    out << string_out.str();
    cout << string_out.str();
}

void BlastApp::PrintDensityResult(string id, string result_acc, string query, string result, string strand, int hit_number, double evalue, ofstream &out) {
    PrintAligment(query, result);

    unsigned int length = query.size();

    pair<unsigned int, unsigned int> binding_site_positions = GetBindingSitePositions(query);

    stringstream string_out;
    string delim_string = "\t";
    string_out << id << delim_string << result_acc << delim_string << strand << delim_string << evalue << delim_string;
    int seq_status;

    unsigned int binding_site_left = binding_site_positions.first;
    unsigned int binding_site_right = binding_site_positions.second;

    if (binding_site_left > binding_site_right) {   
        seq_status = static_cast<int>( SequenceStatus::NoSite );
        string_out << seq_status << delim_string;
        string_out << "-" << delim_string << "-" << delim_string << "-" << delim_string;
        string_out << "-" << delim_string << "-" << delim_string << "-" << delim_string;
        string_out << "-" << delim_string << "-" << delim_string << "-" << delim_string << hit_number << endl;

        out << string_out.str();
        cout << string_out.str();
        return;
    }

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
    cout << string_out.str();
}

string BlastApp::exampleFunction() {
    CNWAligner aligner("TTCATCTCTAAATCTCTCTCATATATATCG", 30, "TTCGATCTCTTCTCCAGATAAATCG", 25);
    //aligner.SetWg(1);
    //aligner.SetWs(-2);
    //aligner.SetWm(1);
    //aligner.SetWms(-2);
    cout << aligner.Run();
    return aligner.GetTranscriptString(); 
}

int BlastApp::Run(void)
{
    // Получение аргументов, установленных в приложение с помощью SetupArgDescriptions()
    const CArgs& args = GetArgs();

    cout << exampleFunction();

    // Получаем шаблон для имени выходных файлов по имени входного (удаляем информацию о директории, удаляем расширение)
    // Добавляем данные о пути к директории для выходных файлов 
    string input_file_name = args["infile"].AsString();
    string output_dir = args["outDir"].AsString();
    char file_system_divider = static_cast<char>( SpecialCharacter::FileSystemDiv );   // В Unix это '/'
    char file_extension_symbol = static_cast<char>( SpecialCharacter::FileExt ); // Т.е. просто точка
    char file_name_divider = static_cast<char>( SpecialCharacter::FileNameDiv );     // Я использую '_'
    if (output_dir.size() > 0 && output_dir[output_dir.size() - 1] != file_system_divider) {
        output_dir += file_system_divider;
    }

    unsigned int input_file_name_dot_pos = input_file_name.find(file_extension_symbol);
    unsigned int last_slash_pos = input_file_name.find_last_of(file_system_divider);
    if (last_slash_pos < input_file_name.size() && last_slash_pos != string::npos) {
        file_name_start = input_file_name.substr(last_slash_pos + 1, input_file_name_dot_pos - last_slash_pos - 1) + file_name_divider + args["hitsize"].AsString();
    } else {
        file_name_start = input_file_name.substr(0, input_file_name_dot_pos) + file_name_divider + args["hitsize"].AsString();
    }

    file_name_start = output_dir + file_name_start;

    cout << "File name start: " << file_name_start << endl;
    vector<string> queries_data = ProcessQueriesData(args["infile"].AsString());

    // Получение ENUM-значения по названию типа BLAST
    EProgram program = ProgramNameToEnum(args["program"].AsString());

    // Обработка и проверка опций
    CRef<CBlastOptionsHandle> options(CBlastOptionsFactory::Create(program, CBlastOptions::eRemote));
    SetOptions(options);
    options->Validate();

    // Созадние объект-менеджера (класс для загрузки данных)
    // CRef<> - автоматически удалит объект на выходе
    // Пока CRef<> существует GetInstance() будет возвращать один и тот же объект
    CRef<CObjectManager> objmgr = CObjectManager::GetInstance();
    if (!objmgr) {
        //throw std::runtime_error("Could not initialize object manager");
        ERR_POST(Fatal << "Could not initialize object manager");
    }

    // Последовательность не белковая
    SDataLoaderConfig data_loader_config(false);
    bool parseFastaFirstLine = false;
    CBlastInputSourceConfig iconfig(data_loader_config, objects::eNa_strand_other, false, parseFastaFirstLine);
    CBlastFastaInputSource fasta_input(args["infile"].AsInputFile(), iconfig);

    CScope scope(*objmgr);

    CBlastInput blast_input(&fasta_input);

    TSeqLocVector query_loc_vector = blast_input.GetAllSeqLocs(scope);

    CRef<IQueryFactory> query_factory(new CObjMgr_QueryFactory(query_loc_vector));

    // Класс для установки в объект класса CRemoteBlast баз данных для поиска
    // Параметры: имя базы данных, тип молекул, entrez query
    // args["entrez"].AsString()
    // Mus musculus[Organism]
    // Homo sapiens[Organism]
    const CSearchDatabase database(args["db"].AsString(), CSearchDatabase::eBlastDbIsNucleotide, args["entrez"].AsString());
  
    // Класс CRemoteBlast используется для отправки запросов, обрабатываемых на серверах NCBI
    // Параметры: последовательность-запрос, опции для запроса, используемые бд (с помощью CSearchDatabase)
    CRemoteBlast remote_blast_worker(query_factory, options, database);

    // SubmitSync() - метод, который должен вызвать CRemoteBlast после построения объекта
    // После его вызова метод GetResultSet() будет возвращать CSearchResultSet
    // Возвращает true, если запрос был отправлен, false иначе

    if (remote_blast_worker.SubmitSync() == false)
         throw std::runtime_error("No results returned by SubmitSync");

    // Вывод идентификатора запроса (RID), выданного системой SPLITD
    // По нему можно получить предыдущие результаты (не старее 36 часов)
    cout << "RID: " << remote_blast_worker.GetRID() << '\n';

    // Класс CSearchResultSet — это контейнер для объектов CSearchResults
    // по одному для каждого запроса, отправленного в поиск
    CSearchResultSet sequences_list = *remote_blast_worker.GetResultSet();

    ofstream out_full;
    ofstream out_blast;
    ofstream meta_blast;
    ofstream meta_seq;

    string output_full_file_name = file_name_start + "_full.txt";
    out_full.open(output_full_file_name);

    string output_blast_file_name = file_name_start + "_blast.txt";
    out_blast.open(output_blast_file_name);

    string meta_blast_file_name = file_name_start + "_meta_blast.txt";
    meta_blast.open(meta_blast_file_name);

    string meta_seq_file_name = file_name_start + "_meta_seq.txt";
    meta_seq.open(meta_seq_file_name);

    unsigned int sequences_count = sequences_list.GetNumResults();

    for (unsigned int j {0}; j < sequences_count; j++) {

        if (meta_blast.is_open()) {
            meta_blast << MSerial_AsnText << sequences_list[j].GetSeqAlign();
        } else {
            cerr << "File " << meta_blast_file_name << " is not open" << endl;
            return -1;
        }

        list data = sequences_list[j].GetSeqAlign()->Get();

        // FASTA-описание текущей последовательности 
        string description = queries_data[2 * j];
        //unsigned int query_accession_start_pos = description.find("AC:") + 3;
        //unsigned int query_accession_end_pos = query_accession_start_pos + description.substr(query_accession_start_pos).find(':');
        //unsigned int query_align_end_pos = query_accession_end_pos + description.substr(query_accession_end_pos).find(' ');
        //string query_accession = description.substr(query_accession_start_pos, query_accession_end_pos - query_accession_start_pos);
        //string query_align_pos = description.substr(query_accession_end_pos + 2, query_align_end_pos - query_accession_end_pos - 1);
        unsigned int trdd_id_start = description.find(" ") + 1;
        unsigned int trrd_id_len = description.find(';') - trdd_id_start;
        string trrd_id = description.substr(trdd_id_start, trrd_id_len);

        cout << '\n' << description << endl;
        cout << "_Count: " << data.size() << endl;
        // Сама последовательность (цепочка нуклеотидов) лежит за ним
        string query = queries_data[2 * j + 1];

        int hit_number = 1;

        for (CRef<CSeq_align> r_align : data) {
            TGi gi = r_align->GetSeq_id(1).GetGi(); 
            double evalue;
            r_align->GetNamedScore(CSeq_align::EScoreType::eScore_EValue, evalue);
            cout << "_Evalue: " << evalue << endl;

            string strand = "PLUS";

            if (ENa_strand::eNa_strand_minus == r_align->GetSeqStrand(0)) {
                strand = "MINUS";
            } else if (ENa_strand::eNa_strand_plus != r_align->GetSeqStrand(0)) {
                ERR_POST(Fatal << "Could not define strand");
            }

            int left_indent = args["leftIndent"].AsInteger();
            int right_indent = args["rightIndent"].AsInteger();

            // Позиции выравнивания
            unsigned int start_query = r_align->GetSeqStart(0);
            unsigned int end_query = r_align->GetSeqStop(0);
            unsigned int start_result = r_align->GetSeqStart(1);

            // Длина исходной последовательности-запроса
            unsigned int query_length = query.size();
            // Длина найденной схожей части
            unsigned int align_length = r_align->GetAlignLength();

            vector<unsigned int> deletions_positions;
            vector<CSeq_align::SIndel> deletions = r_align->GetIndels(1);

            for (CSeq_align::SIndel deletion : deletions) {
                for (unsigned int gap_shift {0}; gap_shift < deletion.length; gap_shift++)
                {
                    deletions_positions.push_back(deletion.product_pos);
                }         
            }

            vector<unsigned int> insertions_positions;
            vector<CSeq_align::SIndel> insertions = r_align->GetIndels(0);
            unsigned int insertions_shift = start_query;
            if (strand == "MINUS") {
                insertions_shift = query_length - end_query - 1;

                int temp = left_indent;
                left_indent = right_indent;
                right_indent = temp;
            }             

            for (CSeq_align::SIndel insertion : insertions) {
                for (unsigned int gap_shift {0}; gap_shift < insertion.length; gap_shift++) {
                    //insertions_positions.push_back(insertion.product_pos);
                    insertions_positions.push_back(insertion.genomic_pos - start_result + insertions_shift + left_indent);
                }
                
            }

            // Расчёт сдвигов начала и конца выравнивания
            int start_diff = static_cast<int>(start_query);
            int end_diff = static_cast<int>(query_length) - static_cast<int>(end_query) - 1;

            // Удаления и вставки (пропуски в исходной и найденной последовательностях соответственно) 
            unsigned int gaps_query_count = deletions_positions.size();
            unsigned int gaps_result_count = insertions_positions.size();

            string gapped_query_sequence = GetGappedQuery(query, deletions_positions);

            int real_start_pos = static_cast<int>(start_result) - start_diff;

            if (strand == "MINUS") {
                real_start_pos = static_cast<int>(start_result) - end_diff;
            } 

            string result_acc;

            if (meta_seq.is_open()) {
                result_acc = GetDescriptionByGi(gi, real_start_pos, start_diff, end_diff, query_length + gaps_query_count - gaps_result_count,
                                       left_indent, right_indent, strand, meta_seq);
            } else {
                cerr << "File " << meta_seq_file_name << " is not open" << endl;
                return -1;
            }

            string gapped_result_sequence = GetGappedResult(gi, real_start_pos, query_length + gaps_query_count - gaps_result_count, 
                                                            left_indent, right_indent, strand, 
                                                            insertions_positions);

            cout << "_Full aligment" << endl; 
            PrintAligment(gapped_query_sequence, gapped_result_sequence);

            cout << "_Blast aligment" << endl; 
            string query_align = gapped_query_sequence.substr(start_diff, align_length);
            string result_align = gapped_result_sequence.substr(start_diff, align_length);
            PrintAligment(query_align, result_align);
/*
            if (out_full.is_open()) {
                PrintDensityResult(trrd_id, result_acc, gapped_query_sequence, gapped_result_sequence, strand, hit_number, evalue, out_full);
            } else {
                cerr << "File " << output_full_file_name << " is not open" << endl;
                return -1;
            }

            if (out_blast.is_open()) {
                string query_align = gapped_query_sequence.substr(start_diff, align_length);
                string result_align = gapped_result_sequence.substr(start_diff, align_length);

                PrintDensityResult(trrd_id, result_acc, query_align, result_align, strand, hit_number, evalue, out_blast);
            } else {
                cerr << "File " << output_blast_file_name << " is not open" << endl;
                return -1;
            }
*/            
            ++hit_number;
        }
    }

    out_full.close();
    out_blast.close();
    meta_blast.close();
    meta_seq.close();

    return 0;
}

void BlastApp::Exit(void)
{
    // Вроде освобождать нечего
    cout << "Exit" << endl;
}

#ifndef SKIP_DOXYGEN_PROCESSING
int NcbiSys_main(int argc, ncbi::TXChar* argv[])
{
    return BlastApp().AppMain(argc, argv);
}
#endif /* SKIP_DOXYGEN_PROCESSING */
