//  ********************************************************************
//  This file is part of Portcullis.
//
//  Portcullis is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  Portcullis is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with Portcullis.  If not, see <http://www.gnu.org/licenses/>.
//  *******************************************************************

#pragma once

#include <fstream>
#include <string>
#include <iostream>
#include <vector>
using std::boolalpha;
using std::ifstream;
using std::string;
using std::vector;

#include <api/BamAlignment.h>

#include <boost/exception/all.hpp>
#include <boost/timer/timer.hpp>
#include <boost/filesystem.hpp>

#include "junction_system.hpp"
using boost::timer::auto_cpu_timer;
using boost::lexical_cast;
using boost::filesystem::absolute;
using boost::filesystem::copy_file;
using boost::filesystem::remove;
using boost::filesystem::exists;
using boost::filesystem::create_symlink;
using boost::filesystem::create_directory;
using boost::filesystem::symbolic_link_exists;
namespace po = boost::program_options;


namespace portcullis {
    
typedef boost::error_info<struct BamFilterError,string> BamFilterErrorInfo;
struct BamFilterException: virtual boost::exception, virtual std::exception { };

enum ClipMode {
    HARD,
    SOFT,
    COMPLETE
};

static string clipToString(ClipMode cm) {
    
    switch(cm) {
        case HARD:
            return "HARD";
        case SOFT:
            return "SOFT";
        case COMPLETE:
            return "COMPLETE";
    }

    return "COMPLETE";
}

static ClipMode clipFromString(string cm) {
    
    if (boost::iequals(cm, "HARD")) {
        return HARD;
    }
    else if (boost::iequals(cm, "SOFT")) {
        return SOFT;
    }
    else if (boost::iequals(cm, "COMPLETE")) {
        return COMPLETE;
    }
    
    BOOST_THROW_EXCEPTION(BamFilterException() << BamFilterErrorInfo(string(
                    "Unrecognised clip mode: ") + cm));
}

class BamFilter {

private:
    
    string junctionFile;
    string bamFile;
    string outputBam;
    StrandSpecific strandSpecific;
    ClipMode clipMode;
    bool saveMSRs;
    bool verbose;
    
public:
    
    BamFilter(const string& _junctionFile, const string& _bamFile, const string& _outputBam, bool _verbose) {
        junctionFile = _junctionFile;
        bamFile = _bamFile;
        outputBam = _outputBam;
        verbose = _verbose;
        strandSpecific = StrandSpecific::UNSTRANDED;
        clipMode = ClipMode::HARD;
        saveMSRs = false;
        
        // Test if provided genome exists
        if (!exists(junctionFile)) {
            BOOST_THROW_EXCEPTION(BamFilterException() << BamFilterErrorInfo(string(
                        "Could not find junction file at: ") + junctionFile));
        }
        
        // Test if provided genome exists
        if (!exists(bamFile)) {
            BOOST_THROW_EXCEPTION(BamFilterException() << BamFilterErrorInfo(string(
                        "Could not find BAM file at: ") + bamFile));
        }
    }
    
    virtual ~BamFilter() {
    }
    
    
protected:
    
    /**
     * Checks a given alignment to see if it exists in the given junction system
     * @param al Alignment to check
     * @param refs References
     * @param js The junction system containing good junctions to keep
     * @return Whether or not the alignment contains a junction found in the junction system
     */
    bool containsJunctionInSystem(BamAlignment& al, RefVector& refs, JunctionSystem& js) {
        
        int32_t refId = al.RefID;
        string refName = refs[refId].RefName;
        int32_t refLength = refs[refId].RefLength;
        int32_t lStart = al.Position;        
        int32_t lEnd = lStart;
        int32_t rStart = lStart;
        int32_t rEnd = lStart;

        for(size_t i = 0; i < al.CigarData.size(); i++) {

            CigarOp op = al.CigarData[i];
            if (op.Type == Constants::BAM_CIGAR_REFSKIP_CHAR) {
                rStart = lEnd + op.Length;
                
                // Create the intron
                shared_ptr<Intron> location(new Intron(RefSeq(refId, refName, refLength), lEnd, rStart, 
                    strandFromBool(al.IsReverseStrand())));                        

                if (js.getJunction(*location) != nullptr) {
                    // Break out of the loop leaving write set to true
                    return true;
                }                    
            }
            else if (BamUtils::opFollowsReference(op.Type)) {
                lEnd += op.Length;                
            }
        }
        
        return false;
    }
    
    shared_ptr<BamAlignment> clipMSR(BamAlignment& al, RefVector& refs, JunctionSystem& js, bool& allBad) {
        
        int32_t refId = al.RefID;
        string refName = refs[refId].RefName;
        int32_t refLength = refs[refId].RefLength;
        int32_t lStart = al.Position;        
        int32_t lEnd = lStart;
        int32_t rStart = lStart;
        int32_t rEnd = lStart;
        size_t opStart = 0;
        bool lastGood = false;
        bool ab = true;
        
        shared_ptr<BamAlignment> modifiedAlignment = make_shared<BamAlignment>(al);
        const char modOp = clipMode == ClipMode::HARD ? Constants::BAM_CIGAR_HARDCLIP_CHAR :
            clipMode == clipMode == ClipMode::SOFT ? Constants::BAM_CIGAR_SOFTCLIP_CHAR : 
                Constants::BAM_CIGAR_DEL_CHAR;

        for(size_t i = 0; i < al.CigarData.size(); i++) {

            CigarOp op = al.CigarData[i];
            
            if (op.Type == Constants::BAM_CIGAR_REFSKIP_CHAR) {
                rStart = lEnd + op.Length;
                
                // Create the intron
                shared_ptr<Intron> location(new Intron(RefSeq(refId, refName, refLength), lEnd, rStart, 
                    strandFromBool(modifiedAlignment->IsReverseStrand())));                        

                if (js.getJunction(*location) != nullptr) {
                    // Found a good junction, so region from start should be left as is, reset start to after junction
                    ab = false;
                    lastGood = true;
                }
                else {
                    if (lastGood) {
                        opStart = i;
                    }
                    
                    for(size_t j = opStart; j < i; j++) {
                        modifiedAlignment->CigarData[j] = CigarOp(modOp, al.CigarData[j].Length);
                    }
                    lastGood = false;
                }
                
                opStart = i+1;
            }
            else if (BamUtils::opFollowsReference(op.Type)) {
                lEnd += op.Length;                
            }
        }
        
        if (!lastGood) {
            for(size_t j = opStart; j < al.CigarData.size(); j++) {
                modifiedAlignment->CigarData[j] = CigarOp(modOp, al.CigarData[j].Length);
            }
        }
        
        allBad = ab;
        return modifiedAlignment;
    }
       

public:
    
    string getBamFile() const {
        return bamFile;
    }

    void setBamFile(string bamFile) {
        this->bamFile = bamFile;
    }

    string getJunctionFile() const {
        return junctionFile;
    }

    void setJunctionFile(string junctionFile) {
        this->junctionFile = junctionFile;
    }

    string getOutputBam() const {
        return outputBam;
    }

    void setOutputBam(string outputBam) {
        this->outputBam = outputBam;
    }

    StrandSpecific getStrandSpecific() const {
        return strandSpecific;
    }

    void setStrandSpecific(StrandSpecific strandSpecific) {
        this->strandSpecific = strandSpecific;
    }
    
    ClipMode getClipMode() const {
        return clipMode;
    }

    void setClipMode(ClipMode clipMode) {
        this->clipMode = clipMode;
    }
    
    bool isSaveMSRs() const {
        return saveMSRs;
    }

    void setSaveMSRs(bool saveMSRs) {
        this->saveMSRs = saveMSRs;
    }



    
    void filter() {
        
        
        
        cout << "Loading junctions from: " << junctionFile << endl;
        
        // Load junction system
        JunctionSystem js(junctionFile);
        
        cout << " - Found " << js.size() << " junctions" << endl << endl;
        
        // Sam header and refs info from the input bam    
        SamHeader header;    
        RefVector refs;
        
        BamReader reader;
        
        if (!reader.Open(bamFile)) {
            BOOST_THROW_EXCEPTION(BamFilterException() << BamFilterErrorInfo(string(
                    "Could not open BAM reader for input: ") + bamFile));
        }
        
        // Sam header and refs info from the input bam
        header = reader.GetHeader();
        refs = reader.GetReferenceData();

        js.setRefs(refs);
       
        cout << " - Processing alignments from: " << bamFile << endl;
        
        string indexFile = bamFile + BAM_INDEX_EXTENSION;
        
        // Opens the index for this BAM file
        if ( !reader.OpenIndex(indexFile) ) {            
            BOOST_THROW_EXCEPTION(BamFilterException() << BamFilterErrorInfo(string(
                    "Could not open index for BAM: ") + indexFile));             
        }
        
        cout << " - Using BAM index: " << indexFile << endl;
        
        BamWriter writer;
        
        if (!writer.Open(outputBam, header, refs)) {
            BOOST_THROW_EXCEPTION(BamFilterException() << BamFilterErrorInfo(string(
                    "Could not open BAM writer for output file: ") + outputBam));
        }
        
        cout << " - Saving filtered alignments to: " << outputBam << endl;
        
        BamWriter mod;
        BamWriter unmod;
        
        if (saveMSRs) {
        
            if (!mod.Open(outputBam + ".mod.bam", header, refs)) {
                BOOST_THROW_EXCEPTION(BamFilterException() << BamFilterErrorInfo(string(
                        "Could not open BAM writer for modified MSRs output file: ") + outputBam + ".mod.bam"));
            }
            
            if (!unmod.Open(outputBam + ".unmod.bam", header, refs)) {
                BOOST_THROW_EXCEPTION(BamFilterException() << BamFilterErrorInfo(string(
                        "Could not open BAM writer for unmodified MSRs output file: ") + outputBam + ".unmod.bam"));
            }

            cout << " - Saving modified MSRs to: " << outputBam << ".mod.bam" << endl;
        }
        
        BamAlignment al;
        uint64_t nbReadsIn = 0;
        uint64_t nbReadsOut = 0;
        uint64_t nbReadsModifiedOut = 0;
        uint32_t nbJunctionsFound = 0;
        
        while(reader.GetNextAlignment(al)) {
            
            nbReadsIn++;
            bool write = false;
                
            if (BamUtils::isSplicedRead(al)) {
                
                // If we are in complete clip mode, or this is a single spliced read, then keep the alignment
                // if its junction is found in the junctions system, otherwise discard it
                if (clipMode == COMPLETE || !BamUtils::isMultiplySplicedRead(al)) {
                    if (containsJunctionInSystem(al, refs, js)) {
                        writer.SaveAlignment(al);
                        nbReadsOut++;
                    }
                }
                // Else we are in HARD or SOFT clip mode and this is an MSR
                else {
                    bool allBad = false;
                    shared_ptr<BamAlignment> clipped = clipMSR(al, refs, js, allBad);
                    if (!allBad) {
                        writer.SaveAlignment(*clipped);
                        if (saveMSRs) {
                            mod.SaveAlignment(*clipped);
                            unmod.SaveAlignment(al);
                        }
                        nbReadsModifiedOut++;
                        nbReadsOut++;
                    }
                }
                
            }
            else {  // Unspliced read so add it to the output
                writer.SaveAlignment(al);
                nbReadsOut++;
            }            
        }
        
        reader.Close();
        writer.Close();
        
        if (saveMSRs) {
            mod.Close();
            unmod.Close();
        }
        
        cout << "done." << endl;
        
        uint32_t diff = nbReadsIn - nbReadsOut;
        
        cout << "Filtered out " << diff << " alignments.  In: " << nbReadsIn << "; Out: " << nbReadsOut << " (Modified: " << nbReadsModifiedOut << ");" << endl << endl;        
    }
  
    static string helpMessage() {
        return string("\nPortcullis BAM Filter Mode Help.\n\n") +
                      "Removes alignments associated with bad junctions from BAM file\n\n" + 
                      "Usage: portcullis bamfilt [options] <junction-file> <bam-file>\n\n" +
                      "Allowed options";
    }
    
    static int main(int argc, char *argv[]) {
        
        // Portcullis args
        string junctionFile;
        string bamFile;
        string output;
        string strandSpecific;
        string clipMode;
        bool saveMSRs; 
        bool verbose;
        bool help;
        
        // Declare the supported options.
        po::options_description generic_options(helpMessage());
        generic_options.add_options()
                ("output,o", po::value<string>(&output)->default_value("filtered.bam"), 
                    "Output BAM file generated by this program.")
                ("strand_specific,ss", po::value<string>(&strandSpecific)->default_value(SSToString(StrandSpecific::UNSTRANDED)), 
                    "Whether BAM alignments were generated using a strand specific RNAseq library: \"unstranded\" (Standard Illumina); \"firststrand\" (dUTP, NSR, NNSR); \"secondstrand\" (Ligation, Standard SOLiD, flux sim reads)  Default: \"unstranded\"")
                ("clip_mode,c", po::value<string>(&clipMode)->default_value(clipToString(ClipMode::HARD)), 
                    "How to clip reads associated with bad junctions: \"HARD\" (Hard clip reads at junction boundary - suitable for cufflinks); \"SOFT\" (Soft clip reads at junction boundaries); \"COMPLETE\" (Remove reads associated exclusively with bad junctions, MSRs covering both good and bad junctions are kept)  Default: \"HARD\"")
                ("save_msrs,m", po::bool_switch(&saveMSRs)->default_value(false), 
                    "Whether or not to output modified MSRs to a separate file.  If true will output to a file with name specified by output with \".msr.bam\" extension")
                ("verbose,v", po::bool_switch(&verbose)->default_value(false), 
                    "Print extra information")
                ("help", po::bool_switch(&help)->default_value(false), "Produce help message")
                ;

        // Hidden options, will be allowed both on command line and
        // in config file, but will not be shown to the user.
        po::options_description hidden_options("Hidden options");
        hidden_options.add_options()
                ("junction-file,g", po::value<string>(&junctionFile), "Path to the junction file containing good junctions.")
                ("bam-file,g", po::value<string>(&bamFile), "Path to the BAM file to filter.")
                ;

        // Positional option for the input bam file
        po::positional_options_description p;
        p.add("junction-file", 1);
        p.add("bam-file", 1);
        
        
        // Combine non-positional options
        po::options_description cmdline_options;
        cmdline_options.add(generic_options).add(hidden_options);

        // Parse command line
        po::variables_map vm;
        po::store(po::command_line_parser(argc, argv).options(cmdline_options).positional(p).run(), vm);
        po::notify(vm);

        // Output help information the exit if requested
        if (help || argc <= 1) {
            cout << generic_options << endl;
            return 1;
        }
        
        

        auto_cpu_timer timer(1, "\nPortcullis BAM filter completed.\nTotal runtime: %ws\n\n");        

        cout << "Running portcullis in BAM filter mode" << endl
             << "--------------------------------" << endl << endl;
        
        // Create the prepare class
        BamFilter filter(junctionFile, bamFile, output, verbose);
        filter.setStrandSpecific(SSFromString(strandSpecific));
        filter.setClipMode(clipFromString(clipMode));
        filter.setSaveMSRs(saveMSRs);
        filter.filter();
        
        return 0;
    }
};
}