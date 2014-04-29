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

#include <fstream>
#include <iostream>
#include <vector>

#include <bam.h>

#include <seqan/basic.h>
#include <seqan/bam_io.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <seqan/stream.h>

#include <boost/unordered_map.hpp>

#include "chunk_loader.hpp"
#include "chunk_processor.hpp"
#include "record_writer.hpp"

using std::string;
using std::cout;
using std::cerr;
using std::endl;
using std::vector;



namespace portculis {

const string DEFAULT_OUTPUT_PREFIX = "portculis_out";
const uint16_t DEFAULT_THREADS = 4;
const uint32_t DEFAULT_CHUNK_SIZE_PER_THREAD = 10000;
const uint32_t DEFAULT_GAP_SIZE = 100;

class Portculis {
private:

    // Can set these from the outside via the constructor
    string bamInputFile;
    string genomeFile;
    string outputPrefix;
    uint16_t threads;
    bool threadedIO;
    uint32_t chunkSizePerThread;
    uint32_t gapSize;
    bool verbose;

    // Created by the constructor
    seqan::BamStream* splicedOut = NULL;
    seqan::BamStream* weirdOut = NULL;
    uint32_t chunkSize;
    
    
    void init(  string _bamInputFile, string _genomeFile, string _outputPrefix, 
                uint16_t _threads, bool _threadedIO, uint32_t _chunkSizePerThread, 
                uint32_t _gapSize, bool _verbose) {
        
        bamInputFile = _bamInputFile;
        genomeFile = _genomeFile;
        outputPrefix = _outputPrefix;
        threads = _threads;
        threadedIO = _threadedIO;
        chunkSizePerThread = _chunkSizePerThread;
        gapSize = _gapSize;
        verbose = _verbose;
        
        chunkSize = (uint32_t) _threads * _chunkSizePerThread;
        
        if (verbose)
            cerr << "Initialised Portculis instance" << endl;
    }
    

protected:

    void findSeedsSeqan(vector<seqan::BamAlignmentRecord>& seeds) {
        // Open reading stream for bam file
        seqan::BamStream bamStreamIn(bamInputFile.c_str(), seqan::BamStream::READ);

        seqan::BamAlignmentRecord record;
        while (!atEnd(bamStreamIn)) {
            readRecord(record, bamStreamIn);

            /*if (!isDefinitelyUnspliced(record.cigar))
            {
                seeds.push_back(record);
            }*/
        }
    }

    void findSeedsSamtools(vector<bam1_t*>& seeds) {
        if (verbose)
            cerr << "Finding seeds" << endl;

        bamFile bamIn = bam_open(bamInputFile.c_str(), "r");

        if (bamIn == 0) {
            throw "Could not open bam file";
        }

        if (verbose)
            cerr << "Opened input file" << endl;

        bam_header_t *header = bam_header_read(bamIn);

        if (verbose)
            cerr << "Processed header" << endl;

        uint64_t alignmentCount = 0;


        // Loop until end of file
        for (;;) {
            bam1_t* b = bam_init1();

            // Exit loop if end of file
            if ((bam_read1(bamIn, b)) < 0) {
                if (verbose)
                    cerr << "EOF" << endl;

                free(b);
                break;
            }

            if (isDefinitelyUnsplicedSamtools(b)) {
                seeds.push_back(b);
            }

            alignmentCount++;
        }

        if (verbose)
            cerr << "Read a total of " << alignmentCount << " alignments" << endl
                << "Found " << seeds.size() << " seeds." << endl;

        bam_header_destroy(header);
        bam_close(bamIn);

    }

    bool isDefinitelyUnsplicedSamtools(bam1_t* bamRecord) {
        uint32_t* cigar = bam1_cigar(bamRecord);
        uint32_t nbCigarOps = bamRecord->core.n_cigar;
        uint32_t crefSkip = bam_cigar_op(BAM_CREF_SKIP);

        // Brute force pattern matching for every position
        for (uint32_t i = 0; i < nbCigarOps; i++) {
            if ((*(cigar + i)) == crefSkip)
                return false;
        }

        return true;
    }

    /*bool isDefinitelyUnsplicedSeqan(seqan::String<seqan::CigarElement<>>& cigar)
    {
        // Brute force pattern matching for every position
        for (unsigned i = 0; i < seqan::length(cigar); i++)
        {
            //if (cigar[i] == 'N')
            //   return false;
        }

        return true;
    }*/

    void linearExecution() {
    }

    void linearChunkedExecution() {
        /*ChunkLoader chunkLoader(*bamStreamIn, chunkSize, gapSize);
        ChunkProcessor chunkProcessor(threads);
        RecordWriter unsplicedWriter(*unsplicedOut);
        RecordWriter splicedWriter(*splicedOut);
        RecordWriter weirdWriter(*weirdOut);

        vector<seqan::BamAlignmentRecord*> chunk;

        while(!chunkLoader.isDone())
        {
            // Clear chunk
            chunk.clear();

            // Load chunk
            chunkLoader.load(&chunk);

            // Process chunk (might use threads)
            vector<seqan::BamAlignmentRecord*> unspliced;
            vector<seqan::BamAlignmentRecord*> spliced;
            vector<seqan::BamAlignmentRecord*> weird;

            chunkProcessor.process(chunk, unspliced, spliced, weird);

            // Append to outputs
            unsplicedWriter.write(unspliced);
            splicedWriter.write(spliced);
            weirdWriter.write(weird);
        }*/
    }

    void threadedChunkedExecution() {
        // Start thread to fetch data from bam
        // Start thread to process data chunks passed to it from fetcher thread
        // Start thread to write data

        // Have to carefully manage communication between threads to ensure we don't get into any nasty situations

    }

public:

    Portculis() {
        init(
                "", 
                "", 
                DEFAULT_OUTPUT_PREFIX, 
                DEFAULT_THREADS, 
                false, 
                DEFAULT_CHUNK_SIZE_PER_THREAD,
                DEFAULT_GAP_SIZE,
                false);
    }
    
    Portculis(  string _bamInputFile, string _genomeFile, string _outputPrefix, 
                uint16_t _threads, bool _threadedIO, uint32_t _chunkSizePerThread, 
                uint32_t _gapSize, bool _verbose) {
        init(   _bamInputFile, _genomeFile, _outputPrefix,
                _threads, _threadedIO, _chunkSizePerThread,
                _gapSize, _verbose);
    }
    
    virtual ~Portculis() {
        /*if (splicedOut != NULL)
            delete splicedOut;

        if (weirdOut != NULL)
            delete weirdOut;*/
    }

    void process() {
        vector<seqan::BamAlignmentRecord> seqanSeeds;
        vector<bam1_t*> samtoolSeeds;

        // Filter out all obvious unspliced reads
        findSeedsSamtools(samtoolSeeds);


        // Tidy up

        BOOST_FOREACH(bam1_t* bamRecord, samtoolSeeds) {
            free(bamRecord);
        }

        // Create a map of the genome
        /*map genomeMap;
        mapGenome(genomeMap);

        // Copy BAM header from input to all BAM output streams
        unsplicedOut->header = bamStreamIn->header;
        splicedOut->header = bamStreamIn->header;
        weirdOut->header = bamStreamIn->header;

        // Decide what mode to work in
        if (threadedIO)
            threadedChunkedExecution();
        else
            linearChunkedExecution();*/

    }
};
}
