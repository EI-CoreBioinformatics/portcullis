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

#define BOOST_TEST_DYN_LINK
#ifdef STAND_ALONE
#define BOOST_TEST_MODULE PORTCULLIS
#endif
#include <boost/test/unit_test.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/operations.hpp>
using boost::filesystem::remove_all;
using boost::filesystem::create_directory;

#include <prepare.hpp>
#include <junction_builder.hpp>

BOOST_AUTO_TEST_SUITE(interface)

BOOST_AUTO_TEST_CASE(test1)
{
    vector<string> bams;
    bams.push_back("resources/clipped3.bam");
    
    string outdir = "resources/temp";
    create_directory(outdir);
    
    string prepdir = "resources/temp/prep";
    create_directory(prepdir);
    
    string juncdir = "resources/temp/junc";
    create_directory(juncdir);
    
    
    portcullis::Prepare prep(prepdir);
    prep.prepare(bams, "resources/artha_chr4.fa");
    prep.outputDetails(); // Outputs settings file
    
    portcullis::JunctionBuilder jb(prepdir, juncdir, "test", 1, false, true);
    jb.process();
    
    remove_all(outdir);
    
    BOOST_CHECK(true);
}

BOOST_AUTO_TEST_CASE(glob1)
{
    vector<string> bams;
    bams.push_back("resources/*.bam");
    
    vector<string> globbedFiles = portcullis::Prepare::globFiles(bams);
    
    BOOST_CHECK(globbedFiles.size() >= 5);
    BOOST_CHECK(globbedFiles.size() < 10);
}

BOOST_AUTO_TEST_CASE(glob2)
{
    vector<string> bams;
    bams.push_back("resources/*.bam");
    bams.push_back("resources/*.bam");
    
    vector<string> globbedFiles = portcullis::Prepare::globFiles(bams);
    
    BOOST_CHECK(globbedFiles.size() >= 10);
}

BOOST_AUTO_TEST_SUITE_END()