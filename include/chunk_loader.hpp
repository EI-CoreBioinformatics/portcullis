
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

class ChunkLoader
{
    private:

        // Set in the constructor
        uint32_t chunkSize;
        seqan::BamStream bamIn;
        uint32_t gapSize;

        // Internal vars
        bool done;

    protected:


    public:
        ChunkLoader(seqan::BamStream _bamIn, uint32_t _chunkSize, uint32_t _gapSize) :
            bamIn(_bamIn), chunkSize(_chunkSize), gapSize(_gapSize)
        {}

        virtual ~ChunkLoader() {}

        void load(vector<seqan::BamAlignmentRecord*>* chunk)
        {
            bool suitableBreak = false;
            uint32_t lastRecordEnd;
            for(uint32_t i = 0; (i < chunkSize || !suitableBreak) && !atEnd(bamIn); i++)
            {
                // Create a new record... remember to clean this up later!!
                seqan::BamAlignmentRecord* record = new seqan::BamAlignmentRecord();

                readRecord(*record, bamIn);

                suitableBreak = hasFlagUnmapped(*record) || (record->beginPos > lastRecordEnd + gapSize);

                chunk->push_back(record);

                lastRecordEnd = record->beginPos + length(record->seq);
            }

            if (atEnd(bamIn))
                done = true;
        }

        bool isDone()
        {
            return done;
        }
};
