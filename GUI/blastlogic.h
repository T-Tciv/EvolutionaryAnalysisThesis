#ifndef BLASTLOGIC_H
#define BLASTLOGIC_H

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

USING_NCBI_SCOPE;
USING_SCOPE(blast);
USING_SCOPE(objects);

class BlastLogic
{
public:
    BlastLogic();
    vector<string> makeBlast();
    string PrintSequence(TGi gi, unsigned int start_pos, unsigned int align_length, int left_indent, int right_indent, string strand);

    void setInpuFileName(string inputFileName);

    void setIndents(int leftIndent, int rightIndent);

    void setBlastSettings(string programName, int hitlistSize, string databaseName, bool isParsed);

    void setEvalue(double evalue);
    void setPenalty(double penalty);
    void setReward(double reward);
    void setEntrezQuery(string entrezQuery);

private:
    string inputFileName;

    int leftIndent;
    int rightIndent;

    string programName;
    int hitlistSize;
    string databaseName;
    bool isParsed;

    double evalue;
    double penalty;
    double reward;
    string entrezQuery;
};

#endif // BLASTLOGIC_H
