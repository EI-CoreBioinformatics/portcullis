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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string.h>
#include <iostream>
#include <fstream>

#include <seqan/bam_io.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>

#include <boost/filesystem.hpp>
#include <boost/foreach.hpp>

#include "portculis_args.hpp"

using std::string;

using portculis::PortculisArgs;

/**
 * Start point for portculis.
 */
int main(int argc, char *argv[])
{
    // Parse args
    PortculisArgs args(argc, argv);

    // Test bam file exists
    if ( !boost::filesystem::exists( args.bam_file ) )
    {
        std::cerr << "ERROR: Specified BAM file " << args.bam_file << " does not exist!" << std::endl;
    }

    seqan::BamStream bamStreamIn(args.bam_file.c_str());

    // Open output stream, "-" means stdin on if reading, else stdout.
    seqan::BamStream bamStreamOut("-", seqan::BamStream::WRITE);

    // Copy header.  The header is automatically written out before
    // the first record.
    bamStreamOut.header = bamStreamIn.header;

    seqan::BamAlignmentRecord record;

    BOOST_FOREACH ( char ch, args.bam_file )
    {
        std::cout << ch << std::endl;

        readRecord(record, bamStreamIn);
        writeRecord(bamStreamOut, record);
    }

    return 0;
}
