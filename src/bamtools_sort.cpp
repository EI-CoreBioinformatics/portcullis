// ***************************************************************************
// bamtools_sort.cpp (c) 2010 Derek Barnett, Erik Garrison
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 27 March 2012 (DB)
// ---------------------------------------------------------------------------
// Sorts an input BAM file
// ***************************************************************************


#include "bamtools_sort.h"

#include <api/SamConstants.h>
#include <api/BamMultiReader.h>
#include <api/BamWriter.h>
#include <api/algorithms/Sort.h>
using namespace BamTools;
using namespace BamTools::Algorithms;

#include <cstdio>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
using namespace std;

namespace BamTools {
  

    
} // namespace BamTools

// constructor
SortTool::SortTool(SortSettings* settings) 
    : m_settings(settings)
    , m_numberOfRuns(0) 
{ 
    // set filename stub depending on inputfile path
    // that way multiple sort runs don't trip on each other's temp files
    if (m_settings) {
        size_t extensionFound = m_settings->InputBamFilename.find(".bam");
        if ( extensionFound != string::npos )
            m_tempFilenameStub = m_settings->InputBamFilename.substr(0,extensionFound);
        m_tempFilenameStub.append(".sort.temp.");
    }
}

// generates mutiple sorted temp BAM files from single unsorted BAM file
bool SortTool::GenerateSortedRuns(void) {
    
    // open input BAM file
    BamReader reader;
    if ( !reader.Open(m_settings->InputBamFilename) ) {
        cerr << "bamtools sort ERROR: could not open " << m_settings->InputBamFilename
             << " for reading... Aborting." << endl;
        return false;
    }
    
    // get basic data that will be shared by all temp/output files 
    SamHeader header = reader.GetHeader();
    if ( !header.HasVersion() )
        header.Version = Constants::SAM_CURRENT_VERSION;
    header.SortOrder = ( m_settings->IsSortingByName
                       ? Constants::SAM_HD_SORTORDER_QUERYNAME
                       : Constants::SAM_HD_SORTORDER_COORDINATE );
    m_headerText = header.ToString();
    m_references = reader.GetReferenceData();
    
    // set up alignments buffer
    BamAlignment al;
    vector<BamAlignment> buffer;
    buffer.reserve( (size_t)(m_settings->MaxBufferCount*1.1) );
    bool bufferFull = false;

    // if sorting by name, we need to generate full char data
    // so can't use GetNextAlignmentCore()
    if ( m_settings->IsSortingByName ) {

        // iterate through file
        while ( reader.GetNextAlignment(al)) {

            // check buffer's usage
            bufferFull = ( buffer.size() >= m_settings->MaxBufferCount );

            // store alignments until buffer is "full"
            if ( !bufferFull )
                buffer.push_back(al);

            // if buffer is "full"
            else {
                // so create a sorted temp file with current buffer contents
                // then push "al" into fresh buffer
                CreateSortedTempFile(buffer);
                buffer.push_back(al);
            }
        }
    }

    // sorting by position, can take advantage of GNACore() speedup
    else {

        // iterate through file
        while ( reader.GetNextAlignmentCore(al) ) {

            // check buffer's usage
            bufferFull = ( buffer.size() >= m_settings->MaxBufferCount );

            // store alignments until buffer is "full"
            if ( !bufferFull )
                buffer.push_back(al);

            // if buffer is "full"
            else {
                // create a sorted temp file with current buffer contents
                // then push "al" into fresh buffer
                CreateSortedTempFile(buffer);
                buffer.push_back(al);
            }
        }
    }

    // handle any leftover buffer contents
    if ( !buffer.empty() )
        CreateSortedTempFile(buffer);
    
    // close reader & return success
    reader.Close();
    return true;
}

bool SortTool::CreateSortedTempFile(vector<BamAlignment>& buffer) {
 
    // do sorting
    SortBuffer(buffer);
  
    // write sorted contents to temp file, store success/fail
    stringstream tempStr;
    tempStr << m_tempFilenameStub << m_numberOfRuns;
    bool success = WriteTempFile( buffer, tempStr.str() );
    
    // save temp filename for merging later
    m_tempFilenames.push_back(tempStr.str());
    
    // clear buffer contents & update run counter
    buffer.clear();
    ++m_numberOfRuns;
    
    // return success/fail of writing to temp file
    // TODO: a failure returned here is not actually caught and handled anywhere
    return success;
}

// merges sorted temp BAM files into single sorted output BAM file
bool SortTool::MergeSortedRuns(void) {
  
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
    if ( !mergedWriter.Open(m_settings->OutputBamFilename, m_headerText, m_references) ) {
        cerr << "bamtools sort ERROR: could not open " << m_settings->OutputBamFilename
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

bool SortTool::Run(void) {
 
    // this does a single pass, chunking up the input file into smaller sorted temp files, 
    // then write out using BamMultiReader to handle merging
    
    if ( GenerateSortedRuns() )
        return MergeSortedRuns();
    else 
        return false;
} 
    
void SortTool::SortBuffer(vector<BamAlignment>& buffer) {
 
    // ** add further custom sort options later ?? **
    
    // sort buffer by desired method
    if ( m_settings->IsSortingByName )
        std::stable_sort( buffer.begin(), buffer.end(), Sort::ByName() );
    else
        std::stable_sort( buffer.begin(), buffer.end(), Sort::ByPosition() );
}
    
bool SortTool::WriteTempFile(const vector<BamAlignment>& buffer,
                                              const string& tempFilename)
{
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
