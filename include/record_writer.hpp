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

#include <boost/foreach.hpp>

using std::vector;

class RecordWriter
{
    private:
        seqan::BamStream outStream;

    protected:

    public:
        RecordWriter(seqan::BamStream _outStream) :
            outStream(_outStream)
        {}

        virtual ~RecordWriter() {}

        void write(vector<seqan::BamAlignmentRecord*>& records)
        {
            BOOST_FOREACH ( seqan::BamAlignmentRecord* record, records )
            {
                writeRecord(outStream, *record);
            }
        }


};
