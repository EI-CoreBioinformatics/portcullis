// ***************************************************************************
// Adapted from bamtools_coverage.cpp/h (c) 2010 Derek Barnett, Erik Garrison
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 24 July 2013
// ---------------------------------------------------------------------------
// Prints coverage data for a single BAM file 
// ***************************************************************************

#pragma once

#include <string>
#include <vector>
using std::string;
using std::pair;
using std::vector;

#include <boost/exception/all.hpp>
#include <boost/foreach.hpp>
#include <boost/shared_ptr.hpp>
using boost::shared_ptr;

#include <api/BamReader.h>
#include <api/BamAlignment.h>
#include <utils/bamtools_pileup_engine.h>
using namespace BamTools;

namespace portculis {
namespace bamtools {
    
typedef boost::error_info<struct BamToolsCoverageError,string> BamToolsCoverageErrorInfo;
struct BamToolsCoverageException: virtual boost::exception, virtual std::exception { };
    
typedef vector<uint32_t> CoverageLevels;
typedef vector<shared_ptr<CoverageLevels> > GenomeCoverageLevels;

class CoverageVisitor : public PileupVisitor {
  
    public:
        CoverageVisitor(shared_ptr<GenomeCoverageLevels> _genomeCoverageLevels, bool _strandSpecific)
            : PileupVisitor(),
            genomeCoverageLevels(_genomeCoverageLevels), strandSpecific(_strandSpecific)
        { }
            
        virtual ~CoverageVisitor() { }
  
    // PileupVisitor interface implementation
    public:
	
        // populates coverage info
        void Visit(const PileupPosition& pileupData) {
            
            //cout << "CVG: " << pileupData.Position << "/" << coverageLevels->size() << "," << pileupData.PileupAlignments.size() << endl;
            
            // Only count alignments piled up on the strand of interest
            uint32_t posCount = 0;
            uint32_t negCount = 0;
            
            posCount = pileupData.PileupAlignments.size();
            // Figure out a way to check for positive or negative counts
            BOOST_FOREACH(PileupAlignment pa, pileupData.PileupAlignments) {
                if (pa.Alignment.IsReverseStrand()) {
                    negCount++;
                }
                else {
                    posCount++;
                }
            }
            
            if (strandSpecific) {
            
                int32_t p1 = pileupData.Position * 2;
                int32_t p2 = p1+1;
                
                (*(*genomeCoverageLevels)[pileupData.RefId])[p1] = posCount;
                (*(*genomeCoverageLevels)[pileupData.RefId])[p2] = negCount;                
            }
            else {
                (*(*genomeCoverageLevels)[pileupData.RefId])[pileupData.Position] = posCount + negCount;
            }
        }
        
        shared_ptr<GenomeCoverageLevels> GetCoverageLevels() const {
            return genomeCoverageLevels;
        }

        
    private:
        shared_ptr<GenomeCoverageLevels> genomeCoverageLevels;
        bool strandSpecific;        
};


// ---------------------------------------------
// CoverageToolPrivate implementation

class CoverageTool {
  
// ctor & dtor
public:
    CoverageTool(bool _strandSpecific)
        : strandSpecific(_strandSpecific)
    { 
    }

    virtual ~CoverageTool(void) {         
    }
    
// interface
public:
    bool pileup(string inputFile) {
        
        //open our BAM reader
        BamReader reader;
        if ( !reader.Open(inputFile) ) {
            BOOST_THROW_EXCEPTION(BamToolsCoverageException() << BamToolsCoverageErrorInfo(string(
                    "bamtools coverage ERROR: could not open input BAM file: ") + inputFile));               
        }

        // retrieve references
        references = reader.GetReferenceData();
        
        // Generate the data
        shared_ptr<GenomeCoverageLevels> data(new GenomeCoverageLevels(references.size()));    
            
        for(size_t i = 0; i < references.size(); i++) {
            (*data)[i] = (strandSpecific ? 
                shared_ptr<CoverageLevels>(new CoverageLevels(references[i].RefLength * 2)) :
                shared_ptr<CoverageLevels>(new CoverageLevels(references[i].RefLength)));
        }
        
        // set up our output 'visitor'
        CoverageVisitor* cv = new CoverageVisitor(data, strandSpecific);

        // set up pileup engine with 'visitor'
        PileupEngine pileup;
        pileup.AddVisitor(cv);

        // process input data
        BamAlignment al;    
        while ( reader.GetNextAlignment(al) ) 
            pileup.AddAlignment(al);
        pileup.Flush();

        // clean up 
        reader.Close();
        delete cv;
        cv = 0;
        
        // return success
        return true;
    }
    
    shared_ptr<GenomeCoverageLevels> GetCoverageLevels() const {
        return genomeCoverageLevels;
    }
    
    void save(string outputFile) {
        
        std::ofstream outfile(outputFile.c_str());
        if(!outfile.is_open()) {
            BOOST_THROW_EXCEPTION(BamToolsCoverageException() << BamToolsCoverageErrorInfo(string(
                    "Could not open coverage output file for writing: ") + outputFile));
        }

        for(size_t i = 0; i < genomeCoverageLevels->size(); i++) {
            
            outfile << ">" << i << endl;

            shared_ptr<CoverageLevels> cvgLvls = (*genomeCoverageLevels)[i];

            outfile << (*cvgLvls)[0];
            for(size_t j = 1; j < cvgLvls->size(); j++) {
                outfile << " " << (*cvgLvls)[j];
            }

            outfile << endl;           
        }
        
        outfile.close();
    }

// data members
private: 
    shared_ptr<GenomeCoverageLevels> genomeCoverageLevels; 
    RefVector references;
    bool strandSpecific;
};  

}
}