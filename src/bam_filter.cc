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

#include <fstream>
#include <string>
#include <iostream>
#include <vector>
using std::boolalpha;
using std::cout;
using std::endl;        
using std::ifstream;
using std::string;
using std::vector;

#include <boost/exception/all.hpp>
#include <boost/timer/timer.hpp>
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
using boost::timer::auto_cpu_timer;
using boost::lexical_cast;
namespace bfs = boost::filesystem;
namespace po = boost::program_options;

#include "bam/bam_reader.hpp"
#include "bam/bam_writer.hpp"
using namespace portcullis::bam;

#include "bam_filter.hpp"

    
portcullis::BamFilter::BamFilter(const path& _junctionFile, const path& _bamFile, const path& _outputBam, bool _verbose) {
    junctionFile = _junctionFile;
    bamFile = _bamFile;
    outputBam = _outputBam;
    verbose = _verbose;
    strandSpecific = Strandedness::UNKNOWN;
    clipMode = ClipMode::HARD;
    saveMSRs = false;

    // Test if provided genome exists
    if (!bfs::exists(junctionFile)) {
        BOOST_THROW_EXCEPTION(BamFilterException() << BamFilterErrorInfo(string(
                    "Could not find junction file at: ") + junctionFile.string()));
    }

    // Test if provided genome exists
    if (!bfs::exists(bamFile)) {
        BOOST_THROW_EXCEPTION(BamFilterException() << BamFilterErrorInfo(string(
                    "Could not find BAM file at: ") + bamFile.string()));
    }
}
    
/**
 * Checks a given alignment to see if it exists in the given junction system
 * @param al Alignment to check
 * @param refs References
 * @param js The junction system containing good junctions to keep
 * @return Whether or not the alignment contains a junction found in the junction system
 */
bool portcullis::BamFilter::containsJunctionInSystem(const BamAlignment& al, const RefSeqPtrList& refs, JunctionSystem& js) {

    int32_t refId = al.getReferenceId();
    string refName = refs[refId]->name;
    int32_t refLength = refs[refId]->length;
    int32_t lStart = al.getPosition();        
    int32_t lEnd = lStart;
    int32_t rStart = lStart;
    int32_t rEnd = lStart;

    for(size_t i = 0; i < al.getNbCigarOps(); i++) {

        CigarOp op = al.getCigarOpAt(i);
        if (op.type == BAM_CIGAR_REFSKIP_CHAR) {
            rStart = lEnd + op.length;

            // Create the intron
            shared_ptr<Intron> location = make_shared<Intron>(RefSeq(refId, refName, refLength), lEnd, rStart - 1);                        

            if (js.getJunction(*location) != nullptr) {
                // Break out of the loop leaving write set to true
                return true;
            }                    
        }
        else if (CigarOp::opConsumesReference(op.type)) {
            lEnd += op.length;                
        }
    }

    return false;
}
    
BamAlignmentPtr portcullis::BamFilter::clipMSR(const BamAlignment& al, const RefSeqPtrList& refs, JunctionSystem& js, bool& allBad) {

    int32_t refId = al.getReferenceId();
    string refName = refs[refId]->name;
    int32_t refLength = refs[refId]->length;
    int32_t lStart = al.getPosition();        
    int32_t lEnd = lStart;
    int32_t rStart = lStart;
    int32_t rEnd = lStart;
    size_t opStart = 0;
    bool lastGood = false;
    bool ab = true;

    BamAlignmentPtr modifiedAlignment = make_shared<BamAlignment>(al);
    const char modOp = clipMode == ClipMode::HARD ? BAM_CIGAR_HARDCLIP_CHAR :
        clipMode == clipMode == ClipMode::SOFT ? BAM_CIGAR_SOFTCLIP_CHAR : 
            BAM_CIGAR_DEL_CHAR;

    for(size_t i = 0; i < al.getNbCigarOps(); i++) {

        CigarOp op = al.getCigarOpAt(i);

        if (op.type == BAM_CIGAR_REFSKIP_CHAR) {
            rStart = lEnd + op.length;

            // Create the intron
            shared_ptr<Intron> location = make_shared<Intron>(RefSeq(refId, refName, refLength), lEnd, rStart - 1);     

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
                    modifiedAlignment->setCigarOpAt(j, CigarOp(modOp, al.getCigarOpAt(j).length));
                }
                lastGood = false;
            }

            opStart = i+1;
        }
        else if (CigarOp::opConsumesReference(op.type)) {
            lEnd += op.length;                
        }
    }

    if (!lastGood) {
        for(size_t j = opStart; j < al.getNbCigarOps(); j++) {
            modifiedAlignment->setCigarOpAt(j, CigarOp(modOp, al.getCigarOpAt(j).length));
        }
    }

    allBad = ab;
    return modifiedAlignment;
}
       
    
void portcullis::BamFilter::filter() {

    cout << "Loading junctions from: " << junctionFile << endl;

    // Load junction system
    JunctionSystem js(junctionFile);

    cout << " - Found " << js.size() << " junctions" << endl << endl;

    BamReader reader(bamFile);
    reader.open();

    shared_ptr<RefSeqPtrList> refs = reader.createRefList();
    js.setRefs(refs);

    if (!bfs::exists(outputBam.parent_path())) {
        if (!bfs::create_directory(outputBam.parent_path())) {
            BOOST_THROW_EXCEPTION(BamFilterException() << BamFilterErrorInfo(string(
                    "Could not create output directory: ") + outputBam.parent_path().string()));
        }
    }
    
    cout << " - Processing alignments from: " << bamFile << endl;

    BamWriter writer(outputBam);
    writer.open(reader.getHeader());

    cout << " - Saving filtered alignments to: " << outputBam << endl;

    BamWriter mod(outputBam.string() + ".mod.bam");
    BamWriter unmod(outputBam.string() + ".unmod.bam");

    if (saveMSRs) {

        mod.open(reader.getHeader());
        unmod.open(reader.getHeader());

        cout << " - Saving modified MSRs to: " << outputBam << ".mod.bam" << endl;
        cout << " - Saving unmodified MSRs to: " << outputBam << ".unmod.bam" << endl;
    }

    uint64_t nbReadsIn = 0;
    uint64_t nbReadsOut = 0;
    uint64_t nbReadsModifiedOut = 0;
    uint32_t nbJunctionsFound = 0;

    while(reader.next()) {

        const BamAlignment& al = reader.current();
    
        nbReadsIn++;
        bool write = false;

        if (al.isSplicedRead()) {

            // If we are in complete clip mode, or this is a single spliced read, then keep the alignment
            // if its junction is found in the junctions system, otherwise discard it
            if (clipMode == COMPLETE || !al.isMultiplySplicedRead()) {
                if (containsJunctionInSystem(al, *refs, js)) {
                    writer.write(al);
                    nbReadsOut++;
                }
            }
            // Else we are in HARD or SOFT clip mode and this is an MSR
            else {
                bool allBad = false;
                BamAlignmentPtr clipped = clipMSR(al, *refs, js, allBad);
                if (!allBad) {
                    writer.write(*clipped);
                    if (saveMSRs) {
                        mod.write(*clipped);
                        unmod.write(al);
                    }
                    nbReadsModifiedOut++;
                    nbReadsOut++;
                }
            }

        }
        else {  // Unspliced read so add it to the output
            writer.write(al);
            nbReadsOut++;
        }            
    }

    reader.close();
    writer.close();

    if (saveMSRs) {
        mod.close();
        unmod.close();
    }

    cout << "done." << endl;

    uint32_t diff = nbReadsIn - nbReadsOut;

    cout << "Filtered out " << diff << " alignments.  In: " << nbReadsIn << "; Out: " << nbReadsOut << " (Modified: " << nbReadsModifiedOut << ");" << endl << endl;        
    
    cout << "Indexing:" << endl;
    cout << " - filtered alignments ... ";
    cout.flush();

    // Create BAM index
    string indexCmd = BamHelper::createIndexBamCmd(outputBam, false);                

    int indexExitCode = system(indexCmd.c_str());                    

    cout << "done." << endl;
}

    
int portcullis::BamFilter::main(int argc, char *argv[]) {

    // Portcullis args
    path junctionFile;
    path bamFile;
    path outputBam;
    string strandSpecific;
    string clipMode;
    bool saveMSRs; 
    bool verbose;
    bool help;

    // Declare the supported options.
    po::options_description generic_options(helpMessage(), 100);
    generic_options.add_options()
            ("output,o", po::value<path>(&outputBam)->default_value("filtered.bam"), 
                "Output BAM file generated by this program.")
            ("strand_specific,ss", po::value<string>(&strandSpecific)->default_value(strandednessToString(Strandedness::UNKNOWN)), 
                "Whether BAM alignments were generated using a strand specific RNAseq library: \"unstranded\" (Standard Illumina); \"firststrand\" (dUTP, NSR, NNSR); \"secondstrand\" (Ligation, Standard SOLiD, flux sim reads)  Default: \"unstranded\".  By default we assume the user does not know the strand specific protocol used for this BAM file.  This has the affect that strand information is derived from splice site information alone, assuming junctions are either canonical or semi-canonical in form.  Default: \"unknown\"")
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
            ("junction-file,g", po::value<path>(&junctionFile), "Path to the junction file containing good junctions.")
            ("bam-file,g", po::value<path>(&bamFile), "Path to the BAM file to filter.")
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
         << "-------------------------------------" << endl << endl;

    // Create the prepare class
    BamFilter filter(junctionFile, bamFile, outputBam, verbose);
    filter.setStrandSpecific(strandednessFromString(strandSpecific));
    filter.setClipMode(clipFromString(clipMode));
    filter.setSaveMSRs(saveMSRs);
    filter.filter();

    return 0;
}
