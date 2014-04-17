
//  ********************************************************************
//  This file is part of portculis.
//
//  Portculis is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  KAT is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with portculis.  If not, see <http://www.gnu.org/licenses/>.
//  *******************************************************************

#pragma once

#include <vector>

#include <seqan/bam_io.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>

using std::vector;

class ChunkProcessor
{
    private:
        uint16_t threads;

    protected:
    public:
        ChunkProcessor(uint16_t _threads) :
            threads(_threads)
        {}
        virtual ~ChunkProcessor() {}

        void process(   vector<seqan::BamAlignmentRecord*>& chunkIn,
                        vector<seqan::BamAlignmentRecord*>& unsplicedOut,
                        vector<seqan::BamAlignmentRecord*>& splicedOut,
                        vector<seqan::BamAlignmentRecord*>& weirdOut)
        {}
};
