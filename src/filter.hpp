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
    
typedef boost::error_info<struct FilterError,string> FilterErrorInfo;
struct FilterException: virtual boost::exception, virtual std::exception { };

const string DEFAULT_FILTER_OUTPUT_DIR = "portcullis_filter_out";
const string DEFAULT_FILTER_OUTPUT_PREFIX = "portcullis";
const double DEFAULT_MIN_ENTROPY = 1.0;
const uint32_t DEFAULT_MIN_DISTINCT_ALIGNMENTS = 3;
const uint32_t DEFAULT_MIN_RELIABLE_ALIGNMENTS = 0;
const uint32_t DEFAULT_MIN_MAXMMES = 5;
const uint16_t DEFAULT_MIN_HAMMING = 2;


class Filter {

private:
    
    string junctionFile;
    //string bamFile;
    string outputDir;
    string outputPrefix;
    uint32_t minDistinctAlignments;
    uint32_t minReliableAlignments;
    double minEntropy;
    uint16_t minHamming;
    uint32_t minMaxMMES;
    bool verbose;
    
public:
    
    Filter(const string& _junctionFile, const string& _outputDir, const string& _outputPrefix, bool _verbose) {
        junctionFile = _junctionFile;
        //bamFile = _bamFile;
        outputDir = _outputDir;
        outputPrefix = _outputPrefix;
        verbose = _verbose;
        
        // Set default thresholds
        minEntropy = DEFAULT_MIN_ENTROPY;
        minDistinctAlignments = DEFAULT_MIN_DISTINCT_ALIGNMENTS;
        minReliableAlignments = DEFAULT_MIN_RELIABLE_ALIGNMENTS;
        minMaxMMES = DEFAULT_MIN_MAXMMES;
        minHamming = DEFAULT_MIN_HAMMING;
        
        // Test if provided genome exists
        if (!exists(junctionFile)) {
            BOOST_THROW_EXCEPTION(FilterException() << FilterErrorInfo(string(
                        "Could not find junction file at: ") + junctionFile));
        }
        
        // Test if provided genome exists
        /*if (!exists(bamFile)) {
            BOOST_THROW_EXCEPTION(FilterException() << FilterErrorInfo(string(
                        "Could not find BAM file at: ") + bamFile));
        }*/
    }
    
    virtual ~Filter() {
    }
    
       

public:
    
    /*string getBamFile() const {
        return bamFile;
    }

    void setBamFile(string bamFile) {
        this->bamFile = bamFile;
    }*/

    string getJunctionFile() const {
        return junctionFile;
    }

    void setJunctionFile(string junctionFile) {
        this->junctionFile = junctionFile;
    }

    uint32_t getMinDistinctAlignments() const {
        return minDistinctAlignments;
    }

    void setMinDistinctAlignments(uint32_t minDistinctAlignments) {
        this->minDistinctAlignments = minDistinctAlignments;
    }
    
    uint16_t getMinHamming() const {
        return minHamming;
    }

    void setMinHamming(uint16_t minHamming) {
        this->minHamming = minHamming;
    }

    uint32_t getMinReliableAlignments() const {
        return minReliableAlignments;
    }

    void setMinReliableAlignments(uint32_t minReliableAlignments) {
        this->minReliableAlignments = minReliableAlignments;
    }


    double getMinEntropy() const {
        return minEntropy;
    }

    void setMinEntropy(double minEntropy) {
        this->minEntropy = minEntropy;
    }

    uint32_t getMinMaxMMES() const {
        return minMaxMMES;
    }

    void setMinMaxMMES(uint32_t minMaxMMES) {
        this->minMaxMMES = minMaxMMES;
    }
    
    string getOutputDir() const {
        return outputDir;
    }

    void setOutputDir(string outputDir) {
        this->outputDir = outputDir;
    }

    string getOutputPrefix() const {
        return outputPrefix;
    }

    void setOutputPrefix(string outputPrefix) {
        this->outputPrefix = outputPrefix;
    }


    
    
    void filter() {
        
        cout << "Filtering obvious false positive junctions or junctions without sufficient evidence" << endl
             << " - Minimum junction entropy required: " << minEntropy << endl
             << " - Minimum number of distinct alignments required: " << minDistinctAlignments << endl
             << " - Minimum number of reliable alignments required: " << minReliableAlignments << endl
             << " - Minimum required value for the maximum of the minimal match on either side of the exon junction (MaxMMES): " << minMaxMMES << endl
             << " - Minimum hamming distance between regions around 5' and 3' splice sites: " << minHamming << endl << endl;
        
        cout << "Loading junctions from: " << junctionFile << endl;
        
        // Load junction system
        JunctionSystem js(junctionFile);
        
        cout << " - Found " << js.getJunctions().size() << " potential junctions" << endl << endl;
        
        JunctionSystem filtered;
        
        uint32_t belowMinEntropyThreshold = 0;
        uint32_t belowDistinctAlignmentThreshold = 0;
        uint32_t belowReliableAlignmentThreshold = 0;
        uint32_t belowMaxMMESThreshold = 0;
        uint32_t below5pHammingThreshold = 0;
        uint32_t below3pHammingThreshold = 0;
        
        for(JunctionPtr j : js.getJunctions()) {
            
            bool filter = false;
            
            if (j->getEntropy() < minEntropy) {
                belowMinEntropyThreshold++;                
                filter = true;
            }
            
            if (j->getNbDistinctAlignments() < minDistinctAlignments) {
                belowDistinctAlignmentThreshold++;
                filter = true;
            }
            
            if (j->getNbReliableAlignments() < minReliableAlignments) {
                belowReliableAlignmentThreshold++;
                filter = true;
            }
                
            if (j->getMaxMMES() < minMaxMMES) {
                belowMaxMMESThreshold++;
                filter = true;
            }
            
            if (j->getHammingDistance5p() < minHamming) {
                below5pHammingThreshold++;
                filter = true;
            }
            
            if (j->getHammingDistance3p() < minHamming) {
                below3pHammingThreshold++;
                filter = true;
            }
            
            if (!filter) {
                filtered.addJunction(j);
            }
        }
        
        cout << "Filtering results: " << endl
             << " - Below entropy threshold: " << belowMinEntropyThreshold << endl
             << " - Below distinct alignment threshold: " << belowDistinctAlignmentThreshold << endl
             << " - Below reliable alignment threshold: " << belowReliableAlignmentThreshold << endl
             << " - Below maxMMES threshold: " << belowMaxMMESThreshold << endl 
             << " - Below hamming threshold at 5' splice site: " << below5pHammingThreshold << endl
             << " - Below hamming threshold at 3' splice site: " << below3pHammingThreshold << endl << endl;
        
        uint32_t diff = js.getJunctions().size() - filtered.getJunctions().size();
        
        cout << "Filtered out " << diff << " junctions.  " << filtered.getJunctions().size() << " junctions remaining" << endl << endl;
        
        if (!exists(outputDir)) {
            create_directory(outputDir);
        }
        
        filtered.saveAll(outputDir + "/" + outputPrefix);        
    }
  
    static string helpMessage() {
        return string("\nPortcullis Filter Mode Help.\n\n") +
                      "Usage: portcullis filter [options] <junction-file>\n\n" +
                      "Allowed options";
    }
    
    static int main(int argc, char *argv[]) {
        
        // Portcullis args
        string junctionFile;
        //string bamFile;
        uint32_t minDistinctAlignments;
        uint32_t minReliableAlignments;
        double minEntropy;
        uint32_t minMaxMMES;
        uint16_t minHamming;
        
        string outputDir;
        string outputPrefix;
        bool verbose;
        bool help;
        
        // Declare the supported options.
        po::options_description generic_options(helpMessage());
        generic_options.add_options()
                ("output_dir,o", po::value<string>(&outputDir)->default_value(DEFAULT_FILTER_OUTPUT_DIR), 
                    "Output directory for files generated by this program.")
                ("output_prefix,p", po::value<string>(&outputPrefix)->default_value(DEFAULT_FILTER_OUTPUT_PREFIX), 
                    "File name prefix for files generated by this program.")
                ("min_entropy,e", po::value<double>(&minEntropy)->default_value(DEFAULT_MIN_ENTROPY), 
                    "The minimum entropy required for a junction")
                ("min_distinct_alignments,d", po::value<uint32_t>(&minDistinctAlignments)->default_value(DEFAULT_MIN_DISTINCT_ALIGNMENTS), 
                    "The minimum number of distinct alignment required for a junction")
                ("min_reliable_alignments,r", po::value<uint32_t>(&minReliableAlignments)->default_value(DEFAULT_MIN_RELIABLE_ALIGNMENTS), 
                    "The minimum number of reliable alignments (alignments uniquely mapping to genome) required for a junction.  WARNING: putting this over 0 will probably produce false negatives.")
                ("min_maxmmes,m", po::value<uint32_t>(&minMaxMMES)->default_value(DEFAULT_MIN_MAXMMES), 
                    "The minimum value required for the maximum of minimal match on either side of exon junction")
                ("min_hamming,h", po::value<uint16_t>(&minHamming)->default_value(DEFAULT_MIN_HAMMING), 
                    "The minimum hamming distance allowed between regions around 5' and 3' splice sites.  (Low hamming score indicates potential repeat region).")
                ("verbose,v", po::bool_switch(&verbose)->default_value(false), 
                    "Print extra information")
                ("help", po::bool_switch(&help)->default_value(false), "Produce help message")
                ;

        // Hidden options, will be allowed both on command line and
        // in config file, but will not be shown to the user.
        po::options_description hidden_options("Hidden options");
        hidden_options.add_options()
                ("junction-file,g", po::value<string>(&junctionFile), "Path to the junction file to process.")
                //("bam-file,g", po::value<string>(&bamFile), "Path to the BAM file to filter.")
                ;

        // Positional option for the input bam file
        po::positional_options_description p;
        p.add("junction-file", 1);
        //p.add("bam-file", 1);
        
        
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
        
        

        auto_cpu_timer timer(1, "\nPortcullis filter completed.\nTotal runtime: %ws\n\n");        

        cout << "Running portcullis in filter mode" << endl
             << "--------------------------------" << endl << endl;
        
        // Create the prepare class
        Filter filter(junctionFile, outputDir, outputPrefix, verbose);
        filter.setMinDistinctAlignments(minDistinctAlignments);
        filter.setMinReliableAlignments(minReliableAlignments);
        filter.setMinEntropy(minEntropy);
        filter.setMinMaxMMES(minMaxMMES);
        filter.setMinHamming(minHamming);
        filter.filter();
        
        return 0;
    }
};
}