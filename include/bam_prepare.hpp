//  ********************************************************************
//  This file is part of Portculis.
//
//  Portculis is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  Portculis is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with Portculis.  If not, see <http://www.gnu.org/licenses/>.
//  *******************************************************************

#pragma once

#include <string>
#include <iostream>
#include <vector>

#include <boost/filesystem.hpp>

#include "bam_utils.hpp"

using std::string;
using std::vector;

using portculis::sortBam;
using portculis::mergeBams;
using portculis::indexBam;

namespace portculis {
    
std::pair<string, string> bamPrep(vector<string> bamFiles, string outputPrefix, bool forcePrep) {
    
    string indexedBamFile;
    string sortedBamFile;
    if (bamFiles.size() > 1) {
        string mergedBam = outputPrefix + string(".merged.bam");
        string sortedBam = outputPrefix + string(".merged.sorted.bam");
        string indexedBam = outputPrefix + string(".merged.sorted.bam.bti");

        bool mergedBamExists = boost::filesystem::exists(mergedBam);
        
        if (mergedBamExists) {
            cout << "Pre-merged BAM detected: " << mergedBam << endl;
            
            if (forcePrep) {
                cout << "Forcing merge anyway due to user request." << endl;
            }
        }
        
        if (!mergedBamExists || forcePrep) {
            cout << "Found " << bamFiles.size() << " BAM files." << endl 
                 << "Merging BAMS...";
            cout.flush();
            mergeBams(bamFiles, mergedBam);
            cout << "done." << endl;
        }

        bool mergedAndSortedBamExists = boost::filesystem::exists(sortedBam);
        
        if (mergedAndSortedBamExists) {
            cout << "Pre-merged and sorted BAM detected: " << sortedBam << endl;
            
            if (forcePrep) {
                cout << "Forcing sort anyway due to user request." << endl;
            }
        }
        
        if (!mergedAndSortedBamExists || forcePrep) {
            cout << "Sorting " << mergedBam << " ... ";
            cout.flush();
            sortBam(mergedBam, sortedBam, false);
            cout << "done." << endl;        
        }

        bool indexedBamExists = boost::filesystem::exists(indexedBam);
        
        if (indexedBamExists) {
            cout << "Pre-merged and sorted BAM index detected: " << indexedBam << endl;
            
            if (forcePrep) {
                cout << "Forcing index creation anyway due to user request." << endl;
            }
        }
        
        if (!indexedBamExists || forcePrep) {
            cout << "Indexing " << sortedBam << " ... ";
            cout.flush();
            indexBam(sortedBam, indexedBam);
            cout << "done." << endl;
        }

        indexedBamFile = indexedBam;
        sortedBamFile = sortedBam;
    }
    else {  // 1 bam File

        string bamFile = bamFiles[0];
        string bamFileToIndex = bamFile;

        // Sort the BAM file if necessary
        bool bamIsSorted = portculis::isSortedBam(bamFile);

        if (bamIsSorted) {
            cout << "Pre-sorted BAM detected: " << bamFile << endl;
            
            if (forcePrep) {
                cout << "Forcing sort anyway due to user request." << endl;
            }
        }
        
        if (!bamIsSorted || forcePrep) {
            string sortedBam = bamFile + ".sorted.bam";

            cout << "Sorting " << bamFile << " ... ";
            cout.flush();
            sortBam(bamFile, sortedBam, false);
            cout << "done" << endl;

            bamFileToIndex = sortedBam;
        }

        string indexedBam = (!bamIsSorted || forcePrep) ? bamFile + string(".sorted.bam.bti") : bamFile + string(".bti");

        bool bamIsIndexed = boost::filesystem::exists(indexedBam);

        if (bamIsIndexed) {
            cout << "Pre-indexed BAM detected: " << indexedBam << endl;
            
            if (forcePrep) {
                cout << "Forcing index creation anyway due to user request." << endl;
            }
        }

        // Index the BAM file if necessary.
        if (!bamIsIndexed || forcePrep) {
            cout << "Indexing " << bamFileToIndex << " ... ";
            cout.flush();
            indexBam(bamFileToIndex, indexedBam);
            cout << "done." << endl;
        }
        indexedBamFile = indexedBam;
        sortedBamFile = bamFileToIndex;
    }
    
    return std::pair<string, string>(sortedBamFile, indexedBamFile);
}
    
}
