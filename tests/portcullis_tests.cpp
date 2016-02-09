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

#include <gtest/gtest.h>

#include <boost/filesystem.hpp>
#include <boost/filesystem/operations.hpp>
using boost::filesystem::remove_all;
using boost::filesystem::create_directory;

#include <prepare.hpp>
#include <portcullis_fs.hpp>
#include <junction_builder.hpp>
using portcullis::Prepare;
using portcullis::PortcullisFS;
using portcullis::JunctionBuilder;

TEST(interface, test1) {
    vector<path> bams;
    bams.push_back("resources/clipped3.bam");
    
    path outdir = "resources/temp";
    create_directory(outdir);
    
    path prepdir = "resources/temp/prep";
    create_directory(prepdir);
    
    path juncdir = "resources/temp/junc";
    create_directory(juncdir);
    
    Prepare prep(prepdir);
    prep.prepare(bams, path("resources/artha_chr4.fa"));
    prep.outputDetails(); // Outputs settings file
    
    //JunctionBuilder(prepdir, juncdir, "test", 1, false, true).process();
    remove_all(outdir);
    
    EXPECT_EQ(true, true);
}

TEST(interface, glob1) {
    vector<path> bams;
    bams.push_back("resources/*.bam");
    
    vector<path> globbedFiles = portcullis::Prepare::globFiles(bams);
    
    EXPECT_GE(globbedFiles.size(), 5);
    EXPECT_LT(globbedFiles.size(), 10);
}

TEST(interface, glob2) {
    vector<path> bams;
    bams.push_back("resources/*.bam");
    bams.push_back("resources/*.bam");
    
    vector<path> globbedFiles = portcullis::Prepare::globFiles(bams);
    
    EXPECT_GE(globbedFiles.size() >= 10);
}
