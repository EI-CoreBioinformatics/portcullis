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
using std::string;
using std::vector;

#include <api/BamReader.h>
using namespace BamTools;

#include <bam.h>

namespace portculis {
    
typedef struct {     // auxiliary data structure
	bamFile fp;      // the file handler
	bam_iter_t iter; // NULL if a region not specified
	int min_mapQ, min_len; // mapQ filter; length filter
} aux_t;

typedef struct {
    int ref;
    int32_t pos;
    uint32_t depth;    
} depth;

class DepthParser {
private:
 
    // Path to the original genome file in fasta format
    string bamFile;
    bool strandSpecific;
    
    bam_header_t *header;
    aux_t** data;
    bam_mplp_t mplp;
        
    depth last;
    bool start;
    int res;
    
protected:
    
    // This function reads a BAM alignment from one BAM file.
    static int read_bam(void *data, bam1_t *b) // read level filters better go here to avoid pileup
    {
        aux_t *aux = (aux_t*)data; // data in fact is a pointer to an auxiliary structure
        int ret = aux->iter? bam_iter_read(aux->fp, aux->iter, b) : bam_read1(aux->fp, b);
        if (!(b->core.flag&BAM_FUNMAP)) {
                if ((int)b->core.qual < aux->min_mapQ) b->core.flag |= BAM_FUNMAP;
                else if (aux->min_len && bam_cigar2qlen(&b->core, bam1_cigar(b)) < aux->min_len) b->core.flag |= BAM_FUNMAP;
        }
        return ret;
    }
    
    
public:
    
    DepthParser(string _bamFile, bool _strandSpecific) : bamFile(_bamFile), strandSpecific(_strandSpecific) {
    
        data = (aux_t**)calloc(1, sizeof(aux_t**));
        data[0] = (aux_t*)calloc(1, sizeof(aux_t));
        data[0]->fp = bam_open(bamFile.c_str(), "r");
        data[0]->min_mapQ = 0;                    
        data[0]->min_len  = 0;                 
             
        header = bam_header_read(data[0]->fp);
        mplp = bam_mplp_init(1, read_bam, (void**)data);
        
        res = 0;
        start = true;
    }
    
    virtual ~DepthParser() {
        
        bam_mplp_destroy(mplp);

        bam_header_destroy(header);
            
        bam_close(data[0]->fp);
        if (data[0]->iter) { 
            bam_iter_destroy(data[0]->iter);
        }
        free(data);
    }
    
     
    string getCurrentRefName() const {
        return string(header->target_name[last.ref]);
    }
    
    int32_t getCurrentRefIndex() const {
        return last.ref;
    }
    
    

    bool loadNextBatch(vector<uint32_t>& depths) {
        
        if (res == 0 && !start) {
            return false;
        }
        
        depths.clear();
        
        int pos = 0;
	int tid = -1;
        int n_plp = 0; // n_plp is the number of covering reads from the i-th BAM
        
        // the core multi-pileup loop
        const bam_pileup1_t** plp = (const bam_pileup1_t**)calloc(1, sizeof(void*)); // plp points to the array of covering reads (internal in mplp)
        
        if (start) {
            start = false;
            
            if ((res = bam_mplp_auto(mplp, &tid, &pos, &n_plp, plp)) > 0) {
                
                int m = 0;
                for (int j = 0; j < n_plp; ++j) {
                    const bam_pileup1_t *p = plp[0] + j;
                    if (p->is_del || p->is_refskip) ++m;
                }
                last.ref = tid;
                last.pos = pos+1;
                last.depth = n_plp - m; 
            }
            else {
                free(plp);
                return false;
            }
        }
        
        // Create the vector
        depths.resize(header->target_len[last.ref], 0);
        
        // Use the details from the last run
        depths[last.pos] = last.depth;
        
        while ((res = bam_mplp_auto(mplp, &tid, &pos, &n_plp, plp)) > 0) {
            
            int m = 0;
            for (int j = 0; j < n_plp; ++j) {
                const bam_pileup1_t *p = plp[0] + j;
                if (p->is_del || p->is_refskip) ++m;
            }

            int32_t rpos = pos+1;
            uint32_t cnt = n_plp - m;

            if (last.ref == tid) {
                // Set the depth
                depths[rpos] = cnt;
            }
            else {

                last.ref = tid;
                last.pos = rpos;
                last.depth = cnt;
                break;
            }
        }
        
        free(plp);
        
        return true;
    }
    
};
}