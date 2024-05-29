#include "blastlogic.h"

BlastLogic::BlastLogic()
{
    this->inputFileName = inputFileName;

    this->leftIndent = 0;
    this->rightIndent = 0;

    this->programName = "blastn";
    this->hitlistSize = 1;
    this->databaseName = "nt";
    this->isParsed = false;

    this->evalue = 0.000000005;
    this->penalty = -2;
    this->reward = 1;
    this->entrezQuery = "";
}

void BlastLogic::setInpuFileName(string inputFileName) {
    this->inputFileName = inputFileName;
}

void BlastLogic::setIndents(int leftIndent, int rightIndent) {
    this->leftIndent = leftIndent;
    this->rightIndent = rightIndent;
}

void BlastLogic::setBlastSettings(string programName, int hitlistSize, string databaseName, bool isParsed){
    this->programName = programName;
    this->hitlistSize = hitlistSize;
    this->databaseName = databaseName;
    this->isParsed = isParsed;
}

void BlastLogic::setEvalue(double evalue){
    this->evalue = evalue;
}

void BlastLogic::setPenalty(double penalty){
    this->penalty = penalty;
}

void BlastLogic::setReward(double reward){
    this->reward = reward;
}

void BlastLogic::setEntrezQuery(string entrezQuery){
    this->entrezQuery = entrezQuery;
}


vector<string> BlastLogic::makeBlast()
{
    // Получение ENUM-значения по названию типа BLAST
    EProgram program = ProgramNameToEnum(programName);

    // Обработка и проверка опций
    CRef<CBlastOptionsHandle> opts(CBlastOptionsFactory::Create(program, CBlastOptions::eRemote));

    opts->SetEvalueThreshold(0.5);

    opts->SetHitlistSize(hitlistSize);

    if (CBlastNucleotideOptionsHandle* nucl_handle =
        dynamic_cast<CBlastNucleotideOptionsHandle*>(&*opts)) {
        nucl_handle->SetMatchReward(0);
        nucl_handle->SetMismatchPenalty(0);
    }

    opts->Validate();  // Can throw CBlastException::eInvalidOptions for invalid option.

    // Созадние объект-менеджера (класс для загрузки данных)
    // CRef<> - автоматически удалит объект на выходе
    // Пока CRef<> существует GetInstance() будет возвращать один и тот же объект
    CRef<CObjectManager> objmgr = CObjectManager::GetInstance();
    if (!objmgr) {
         throw std::runtime_error("Could not initialize object manager");
    }

    // Определяем, является ли последовательность белковой (двойное отрицание для неявного приведения к bool)
    const bool is_protein = !!Blast_QueryIsProtein(opts->GetOptions().GetProgramType());
    SDataLoaderConfig dlconfig(is_protein);
    CBlastInputSourceConfig iconfig(dlconfig, objects::eNa_strand_other, false, isParsed);

    std::ifstream input_file(inputFileName);

    CBlastFastaInputSource fasta_input(input_file, iconfig);

    CScope scope(*objmgr);

    CBlastInput blast_input(&fasta_input);

    TSeqLocVector query_loc = blast_input.GetAllSeqLocs(scope);

    CRef<IQueryFactory> query_factory(new CObjMgr_QueryFactory(query_loc));

    // Определение типа молекул в бд по типу BLAST
    // Это не то же самое, что is_protein, т.к. тип последовательности-запроса и тип бд могут не совпадать
    CSearchDatabase::EMoleculeType db_molecule_type = (program == eBlastp || program == eBlastx || program == eRPSBlast || program == eRPSTblastn)
    ? CSearchDatabase::eBlastDbIsProtein : CSearchDatabase::eBlastDbIsNucleotide;

    // Класс для установки в объект класса CRemoteBlast баз данных для поиска
    // Параметры: имя базы данных, тип молекул, entrez query
    const CSearchDatabase target_db(databaseName, db_molecule_type);

    // Класс CRemoteBlast используется для отправки запросов, обрабатываемых на серверах NCBI
    // Параметры: последовательность-запрос, опции для запроса, используемые бд (с помощью CSearchDatabase)
    CRemoteBlast blaster(query_factory, opts, target_db);

    // SubmitSync() - метод, который должен вызвать CRemoteBlast после построения объекта
    // После его вызова метод GetResultSet() будет возвращать CSearchResultSet
    // Возвращает true, если запрос был отправлен, false иначе
    bool status = blaster.SubmitSync();

    if (status == false)
         throw std::runtime_error("No results returned by SubmitSync");

    // Вывод идентификатора запроса (RID), выданного системой SPLITD
    // По нему можно получить предыдущие результаты (не старее 36 часов)
    cerr << "RID: " << blaster.GetRID() << '\n';

    // Класс CSearchResultSet — это контейнер для объектов CSearchResults
    // по одному для каждого запроса, отправленного в поиск
    CSearchResultSet results = *blaster.GetResultSet();

    vector<string> blastResults;

    for (unsigned int i = 0; i < results.GetNumResults(); i++) {
         list data = results[i].GetSeqAlign()->Get();

         for (CRef<CSeq_align> r_align : data) {
             TGi gi = r_align->GetSeq_id(1).GetGi();
             unsigned int start_pos = r_align->GetSeqStart(1);
             unsigned int align_length = r_align->GetAlignLength();
             string strand = "PLUS";

             if (ENa_strand::eNa_strand_minus == r_align->GetSeqStrand(0)) {
                 strand = "MINUS";
             } else if (ENa_strand::eNa_strand_plus != r_align->GetSeqStrand(0)) {
                 ERR_POST(Fatal << "Could not define strand");
             }

             string sequenceData = PrintSequence(gi, start_pos, align_length, leftIndent, rightIndent, strand);
             blastResults.push_back(sequenceData);
         }
    }

    return blastResults;
}

string BlastLogic::PrintSequence(TGi gi, unsigned int start_pos, unsigned int align_length, int left_indent, int right_indent, string strand) {
    stringstream sequenceDataStream;
    // Создаём object manager
    CRef<CObjectManager> object_manager = CObjectManager::GetInstance();

    // Создаём GenBank data loader и связываем его с OM.
    // GenBank loader автоматически помечается default loader
    CGBDataLoader::RegisterInObjectManager(*object_manager);

    CScope scope(*object_manager);
    //  Добавляем дефолтные data loaders (т.е. GenBank data loader) в эту scope
    scope.AddDefaults();

    CSeq_id seq_id;
    seq_id.SetGi(gi);

    // Получаем Bioseq handle для seq_id.
    CBioseq_Handle bioseq_handle = scope.GetBioseqHandle(seq_id);

    if ( !bioseq_handle ) {
        ERR_POST(Fatal << "Bioseq not found, with GI=" << gi);
    }

    // Получаем последовательность, используя CSeqVector.
    CSeqVector seq_vect = bioseq_handle.GetSeqVector(CBioseq_Handle::eCoding_Iupac);

    // const CBioseq_set::TDescr & description = bioseq_handle.GetDescr();
    list description = bioseq_handle.GetDescr().Get();

    sequenceDataStream << "TITLE: ";

    for (CRef<CSeqdesc> descr_elem : description) {
        if (!!descr_elem->IsTitle()) {
            sequenceDataStream << MSerial_AsnText << descr_elem->GetTitle();
            break;
        }
    }

    sequenceDataStream <<  endl;

    // * With GenBank loader this request should not load the whole Bioseq.
    CBioseq_Handle::TId ids = scope.GetIds(seq_id);

    sequenceDataStream << "IDs: ";
    for (const auto& id : ids) {
        if ( &id != &ids.front() ) {
            sequenceDataStream << " + ";
        }
        sequenceDataStream << id.AsString();
    }
    sequenceDataStream << endl;

    unsigned int seq_length = seq_vect.size();
    sequenceDataStream << "LENGTH = " << seq_length << ", STRAND = " << strand << endl;
    sequenceDataStream << "ALIGMENT START = " << start_pos << ", ALIGMENT END = " << start_pos + align_length - 1 << endl;

    int left = max((int) start_pos - left_indent, 0);
    int right = min((int) start_pos + (int) align_length + right_indent, (int) seq_length);

    sequenceDataStream << "WITH INDENTS: START = " << left << ", END = " << right - 1 << endl;

    string str;
    seq_vect.GetSeqData(left, right, str);
    sequenceDataStream << str << endl;

    return sequenceDataStream.str();
}
