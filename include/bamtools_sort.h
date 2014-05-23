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

#include <cstdio>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
using namespace std;

#include <api/SamConstants.h>
#include <api/BamMultiReader.h>
#include <api/BamWriter.h>
#include <api/algorithms/Sort.h>
using namespace BamTools;
using namespace BamTools::Algorithms;

namespace portculis {
namespace bamtools {

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
    bool hasInputBamFilename;
    bool hasMaxBufferCount;
    bool hasMaxBufferMemory;
    bool hasOutputBamFilename;
    bool isSortingByName;

    // filenames
    string inputBamFilename;
    string outputBamFilename;

    // parameters
    unsigned int maxBufferCount;
    unsigned int maxBufferMemory;

    // constructor
    SortSettings(void)
        : hasInputBamFilename(false)
        , hasMaxBufferCount(false)
        , hasMaxBufferMemory(false)
        , hasOutputBamFilename(false)
        , isSortingByName(false)
        , inputBamFilename("")
        , outputBamFilename("")
        , maxBufferCount(SORT_DEFAULT_MAX_BUFFER_COUNT)
        , maxBufferMemory(SORT_DEFAULT_MAX_BUFFER_MEMORY)
    { }
};

    
class SortTool {
  
// ctor & dtor
public:
        SortTool(SortSettings* settings) 
        : m_settings(settings)
        , m_numberOfRuns(0) 
    { 
        // set filename stub depending on inputfile path
        // that way multiple sort runs don't trip on each other's temp files
        if (m_settings) {
            size_t extensionFound = m_settings->inputBamFilename.find(".bam");
            if ( extensionFound != string::npos )
                m_tempFilenameStub = m_settings->inputBamFilename.substr(0,extensionFound);
            m_tempFilenameStub.append(".sort.temp.");
        }
    }
    
    ~SortTool(void) { }
        
// 'public' interface
public:
    bool run(void) {
 
        // this does a single pass, chunking up the input file into smaller sorted temp files, 
        // then write out using BamMultiReader to handle merging

        if ( generateSortedRuns() )
            return mergeSortedRuns();
        else 
            return false;
    } 

// internal methods
protected:
    bool createSortedTempFile(vector<BamAlignment>& buffer) {
 
        // do sorting
        sortBuffer(buffer);

        // write sorted contents to temp file, store success/fail
        stringstream tempStr;
        tempStr << m_tempFilenameStub << m_numberOfRuns;
        bool success = writeTempFile( buffer, tempStr.str() );

        // save temp filename for merging later
        m_tempFilenames.push_back(tempStr.str());

        // clear buffer contents & update run counter
        buffer.clear();
        ++m_numberOfRuns;

        // return success/fail of writing to temp file
        // TODO: a failure returned here is not actually caught and handled anywhere
        return success;
    }
    
    bool generateSortedRuns(void) {
    
        // open input BAM file
        BamReader reader;
        if ( !reader.Open(m_settings->inputBamFilename) ) {
            cerr << "bamtools sort ERROR: could not open " << m_settings->inputBamFilename
                 << " for reading... Aborting." << endl;
            return false;
        }

        // get basic data that will be shared by all temp/output files 
        SamHeader header = reader.GetHeader();
        if ( !header.HasVersion() )
            header.Version = Constants::SAM_CURRENT_VERSION;
        header.SortOrder = ( m_settings->isSortingByName
                           ? Constants::SAM_HD_SORTORDER_QUERYNAME
                           : Constants::SAM_HD_SORTORDER_COORDINATE );
        m_headerText = header.ToString();
        m_references = reader.GetReferenceData();

        // set up alignments buffer
        BamAlignment al;
        vector<BamAlignment> buffer;
        buffer.reserve( (size_t)(m_settings->maxBufferCount*1.1) );
        bool bufferFull = false;

        // if sorting by name, we need to generate full char data
        // so can't use GetNextAlignmentCore()
        if ( m_settings->isSortingByName ) {

            // iterate through file
            while ( reader.GetNextAlignment(al)) {

                // check buffer's usage
                bufferFull = ( buffer.size() >= m_settings->maxBufferCount );

                // store alignments until buffer is "full"
                if ( !bufferFull )
                    buffer.push_back(al);

                // if buffer is "full"
                else {
                    // so create a sorted temp file with current buffer contents
                    // then push "al" into fresh buffer
                    createSortedTempFile(buffer);
                    buffer.push_back(al);
                }
            }
        }

        // sorting by position, can take advantage of GNACore() speedup
        else {

            // iterate through file
            while ( reader.GetNextAlignmentCore(al) ) {

                // check buffer's usage
                bufferFull = ( buffer.size() >= m_settings->maxBufferCount );

                // store alignments until buffer is "full"
                if ( !bufferFull )
                    buffer.push_back(al);

                // if buffer is "full"
                else {
                    // create a sorted temp file with current buffer contents
                    // then push "al" into fresh buffer
                    createSortedTempFile(buffer);
                    buffer.push_back(al);
                }
            }
        }

        // handle any leftover buffer contents
        if ( !buffer.empty() )
            createSortedTempFile(buffer);

        // close reader & return success
        reader.Close();
        return true;
    }
    
    
    bool mergeSortedRuns(void) {
  
        // open up multi reader for all of our temp files
        // this might get broken up if we do a multi-pass system later ??
        BamMultiReader multiReader;
        if ( !multiReader.Open(m_tempFilenames) ) {
            cerr << "bamtools sort ERROR: could not open BamMultiReader for merging temp files... Aborting."
                 << endl;
            return false;
        }

        // open writer for our completely sorted output BAM file
        BamWriter mergedWriter;
        if ( !mergedWriter.Open(m_settings->outputBamFilename, m_headerText, m_references) ) {
            cerr << "bamtools sort ERROR: could not open " << m_settings->outputBamFilename
                 << " for writing... Aborting." << endl;
            multiReader.Close();
            return false;
        }

        // while data available in temp files
        BamAlignment al;
        while ( multiReader.GetNextAlignmentCore(al) )
            mergedWriter.SaveAlignment(al);

        // close files
        multiReader.Close();
        mergedWriter.Close();

        // delete all temp files
        vector<string>::const_iterator tempIter = m_tempFilenames.begin();
        vector<string>::const_iterator tempEnd  = m_tempFilenames.end();
        for ( ; tempIter != tempEnd; ++tempIter ) {
            const string& tempFilename = (*tempIter);
            remove(tempFilename.c_str());
        }

        // return success
        return true;
    }
    
    bool writeTempFile(const vector<BamAlignment>& buffer, const string& tempFilename) {
        // open temp file for writing
        BamWriter tempWriter;
        if ( !tempWriter.Open(tempFilename, m_headerText, m_references) ) {
            cerr << "bamtools sort ERROR: could not open " << tempFilename
                 << " for writing." << endl;
            return false;
        }

        // write data
        vector<BamAlignment>::const_iterator buffIter = buffer.begin();
        vector<BamAlignment>::const_iterator buffEnd  = buffer.end();
        for ( ; buffIter != buffEnd; ++buffIter )  {
            const BamAlignment& al = (*buffIter);
            tempWriter.SaveAlignment(al);
        }

        // close temp file & return success
        tempWriter.Close();
        return true;
    }
    
    void sortBuffer(vector<BamAlignment>& buffer) {
 
        // ** add further custom sort options later ?? **

        // sort buffer by desired method
        if ( m_settings->isSortingByName )
            std::stable_sort( buffer.begin(), buffer.end(), Sort::ByName() );
        else
            std::stable_sort( buffer.begin(), buffer.end(), Sort::ByPosition() );
    }

// data members
private:
    SortSettings* m_settings;
    string m_tempFilenameStub;
    int m_numberOfRuns;
    string m_headerText;
    RefVector m_references;
    vector<string> m_tempFilenames;
};

}
}

