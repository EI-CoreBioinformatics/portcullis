
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

#include <fstream>
#include <iostream>
#include <vector>

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
using std::cerr;
using std::endl;
using std::vector;

typedef boost::unordered_map<string, seqan::Dna5String> map;

class Portculis
{
    private:

        // Can set these from the outside via the constructor
        string bamFile;
        string genomeFile;
        string outputPrefix;
        uint16_t threads;
        bool threadedIO;
        uint32_t chunkSizePerThread;
        uint32_t gapSize;
        bool verbose;

        // Created by the constructor
        seqan::BamStream* bamStreamIn = NULL;
        seqan::BamStream* unsplicedOut = NULL;
        seqan::BamStream* splicedOut = NULL;
        seqan::BamStream* weirdOut = NULL;
        uint32_t chunkSize;


    protected:


        void mapGenome(map& genomeMap)
        {
            // Open file, create RecordReader and check all is well
            std::ifstream in(genomeFile.c_str());
            seqan::RecordReader<std::ifstream, seqan::SinglePass<> > reader(in);

            // Read file record-wise.
            seqan::CharString id;
            seqan::Dna5String seq;
            while (!atEnd(reader))
            {
                if (readRecord(id, seq, reader, seqan::Fasta()) != 0)
                    throw "ERROR: Could not read record from genome file";

                if (verbose)
                    cerr << "Adding to genome map: " << id << endl;

                string ssid;
                seqan::append(ssid, id);

                genomeMap[ssid] = seq;
            }
        }

        void linearExecution()
        {
        }

        void linearChunkedExecution()
        {
            ChunkLoader chunkLoader(*bamStreamIn, chunkSize, gapSize);
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
            }
        }

        void threadedChunkedExecution()
        {
            // Start thread to fetch data from bam
            // Start thread to process data chunks passed to it from fetcher thread
            // Start thread to write data

            // Have to carefully manage communication between threads to ensure we don't get into any nasty situations

        }

    public:
        Portculis(string _bamFile, string _genomeFile, string _outputPrefix, uint16_t _threads, bool _threadedIO, uint32_t _chunkSizePerThread, uint32_t _gapSize, bool _verbose) :
            bamFile(_bamFile), genomeFile(_genomeFile), outputPrefix(_outputPrefix), threads(_threads), threadedIO(_threadedIO), chunkSizePerThread(_chunkSizePerThread), gapSize(_gapSize), verbose(_verbose),
            chunkSize((uint32_t)_threads * _chunkSizePerThread)
        {
            // Open reading stream for bam file
            bamStreamIn = new seqan::BamStream(bamFile.c_str(), seqan::BamStream::READ);

            // Create filenames to output files
            std::stringstream sstmUnsplicedOut; sstmUnsplicedOut << outputPrefix << "_unspliced.bam"; string unsplicedOutFile = sstmUnsplicedOut.str();
            std::stringstream sstmSplicedOut; sstmSplicedOut << outputPrefix << "_spliced.bam"; string splicedOutFile = sstmSplicedOut.str();
            std::stringstream sstmWeirdOut; sstmWeirdOut << outputPrefix << "_weird.bam"; string weirdOutFile = sstmWeirdOut.str();

            // Create and open output streams
            unsplicedOut = new seqan::BamStream(unsplicedOutFile.c_str(), seqan::BamStream::WRITE);
            splicedOut = new seqan::BamStream(splicedOutFile.c_str(), seqan::BamStream::WRITE);
            weirdOut = new seqan::BamStream(weirdOutFile.c_str(), seqan::BamStream::WRITE);

            if (verbose)
                cerr << "Initialised Portculis instance" << endl;
        }

        virtual ~Portculis()
        {
            if (bamStreamIn != NULL)
                delete bamStreamIn;

            if (unsplicedOut != NULL)
                delete unsplicedOut;

            if (splicedOut != NULL)
                delete splicedOut;

            if (weirdOut != NULL)
                delete weirdOut;
        }

        void process()
        {
            // Create a map of the genome
            map genomeMap;
            mapGenome(genomeMap);

            // Copy BAM header from input to all BAM output streams
            unsplicedOut->header = bamStreamIn->header;
            splicedOut->header = bamStreamIn->header;
            weirdOut->header = bamStreamIn->header;

            // Decide what mode to work in
            if (threadedIO)
                threadedChunkedExecution();
            else
                linearChunkedExecution();

        }
};
