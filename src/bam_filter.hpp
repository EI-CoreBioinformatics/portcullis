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


class BamFilter {

private:
    
    string junctionFile;
    string bamFile;
    string outputBam;
    StrandSpecific strandSpecific;
    bool verbose;
    
public:
    
    BamFilter(const string& _junctionFile, const string& _bamFile, const string& _outputBam, bool _verbose) {
        junctionFile = _junctionFile;
        bamFile = _bamFile;
        outputBam = _outputBam;
        verbose = _verbose;
        strandSpecific = StrandSpecific::UNSTRANDED;
        
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

    
    void filter() {
        
        
        
        cout << "Loading junctions from: " << junctionFile << endl;
        
        // Load junction system
        JunctionSystem js(junctionFile);
        
        cout << " - Found " << js.getJunctions().size() << " potential valid junctions" << endl << endl;
        
        // Sam header and refs info from the input bam    
        SamHeader header;    
        RefVector refs;
        
        BamReader reader;
        
        if (!reader.Open(bamFile)) {
            BOOST_THROW_EXCEPTION(JunctionBuilderException() << JunctionBuilderErrorInfo(string(
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
            BOOST_THROW_EXCEPTION(JunctionBuilderException() << JunctionBuilderErrorInfo(string(
                    "Could not open index for BAM: ") + indexFile));             
        }
        
        cout << " - Using BAM index: " << indexFile << endl;
        
        BamWriter writer;
        
        if (!writer.Open(outputBam, header, refs)) {
            BOOST_THROW_EXCEPTION(JunctionBuilderException() << JunctionBuilderErrorInfo(string(
                    "Could not open BAM writer for non-spliced file: ") + outputBam));
        }

        cout << " - Saving filtered alignments to: " << outputBam << endl;
        
        BamAlignment al;
        uint64_t nbReadsIn = 0;
        uint64_t nbReadsOut = 0;
        while(reader.GetNextAlignment(al)) {
            
            nbReadsIn++;
            bool write = true;
                
            if (BamUtils::isSplicedRead(al)) {
                
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
                            //strandSpecific ? strandFromBool(al.IsReverseStrand()) : UNKNOWN));
                                UNKNOWN));
                        
                        if (js.getJunction(*location) != nullptr) {
                            // Break out of the loop leaving write set to true
                            break;
                        }
                    }
                    else if (BamUtils::opFollowsReference(op.Type)) {
                        lEnd += op.Length;                
                    }
                }
                
                write = false;
            }            
            
            if (write) {
                writer.SaveAlignment(al);
                nbReadsOut++;
            }
        }
        
        reader.Close();
        writer.Close();
        cout << "done." << endl;
        
        uint32_t diff = nbReadsIn - nbReadsOut;
        
        cout << "Filtered out " << diff << " alignments.  In: " << nbReadsIn << "; Out: " << nbReadsOut << ";" << endl << endl;
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
        bool verbose;
        bool help;
        
        // Declare the supported options.
        po::options_description generic_options(helpMessage());
        generic_options.add_options()
                ("output,o", po::value<string>(&output)->default_value(DEFAULT_FILTER_OUTPUT_DIR), 
                    "Output BAM file generated by this program.")
                ("strand_specific,ss", po::value<string>(&strandSpecific)->default_value(SSToString(StrandSpecific::UNSTRANDED)), 
                    "Whether BAM alignments were generated using a strand specific RNAseq library: \"unstranded\" (Standard Illumina); \"firststrand\" (dUTP, NSR, NNSR); \"secondstrand\" (Ligation, Standard SOLiD, flux sim reads)  Default: \"unstranded\"")
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
        filter.filter();
        
        return 0;
    }
};
}