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
#include <vector>

#include <api/BamMultiReader.h>
#include <api/BamWriter.h>

using std::string;
using std::vector;

using namespace BamTools;

namespace portculis {
    
    void mergeBams(vector<string>& bamFiles, string bamOut) {
        
        BamMultiReader reader;
        if (!reader.Open(bamFiles)) {
            throw "Could not open input BAM files";
        }

        const SamHeader header = reader.GetHeader();
        const RefVector refs = reader.GetReferenceData();
        
        BamWriter writer;
        if (!writer.Open(bamOut, header, refs)) {
            throw "Could not open output BAM file for merged results";
        }
        
        BamAlignment al;
        while (reader.GetNextAlignment(al)) {
            writer.SaveAlignment(al);
        }
        
        reader.Close();
        writer.Close();
    }
}