// ***************************************************************************
// bamtools_sort.h (c) 2010 Derek Barnett, Erik Garrison
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 7 April 2011 (DB)
// ---------------------------------------------------------------------------
// Sorts a BAM file
// ***************************************************************************

// Adapted from the original BamTools sort tool for portculis

#pragma once

#include <string>
#include <vector>

#include <api/BamMultiReader.h>

using std::string;
using std::vector;

namespace BamTools {

// defaults
//
// ** These defaults should be tweaked & 'optimized' per testing ** //
//
//    I say 'optimized' because each system will naturally perform
//    differently.  We will attempt to determine a sensible
//    compromise that should perform well on average.
const unsigned int SORT_DEFAULT_MAX_BUFFER_COUNT  = 500000;  // max numberOfAlignments for buffer
const unsigned int SORT_DEFAULT_MAX_BUFFER_MEMORY = 1024;    // Mb
    
// ---------------------------------------------
// SortSettings implementation

struct SortSettings {

    // flags
    bool HasInputBamFilename;
    bool HasMaxBufferCount;
    bool HasMaxBufferMemory;
    bool HasOutputBamFilename;
    bool IsSortingByName;

    // filenames
    string InputBamFilename;
    string OutputBamFilename;

    // parameters
    unsigned int MaxBufferCount;
    unsigned int MaxBufferMemory;

    // constructor
    SortSettings(void)
        : HasInputBamFilename(false)
        , HasMaxBufferCount(false)
        , HasMaxBufferMemory(false)
        , HasOutputBamFilename(false)
        , IsSortingByName(false)
        , InputBamFilename("")
        , OutputBamFilename("")
        , MaxBufferCount(SORT_DEFAULT_MAX_BUFFER_COUNT)
        , MaxBufferMemory(SORT_DEFAULT_MAX_BUFFER_MEMORY)
    { }
};

    
class SortTool {
  
    // ctor & dtor
    public:
        SortTool(SortSettings* settings);
        ~SortTool(void) { }
        
    // 'public' interface
    public:
        bool Run(void);
        
    // internal methods
    protected:
        bool CreateSortedTempFile(vector<BamAlignment>& buffer);
        bool GenerateSortedRuns(void);
        bool MergeSortedRuns(void);
        bool WriteTempFile(const vector<BamAlignment>& buffer, const string& tempFilename);
        void SortBuffer(vector<BamAlignment>& buffer);
        
    // data members
    private:
        SortSettings* m_settings;
        string m_tempFilenameStub;
        int m_numberOfRuns;
        string m_headerText;
        RefVector m_references;
        vector<string> m_tempFilenames;
};


} // namespace BamTools

