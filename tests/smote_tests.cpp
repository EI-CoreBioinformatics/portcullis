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

#include <iostream>
#include <fstream>
#include <string>
using std::cout;
using std::cerr;
using std::endl;
using std::stringstream;


#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/path.hpp>
namespace bfs = boost::filesystem;
using bfs::path;

#include <portcullis/ml/smote.hpp>
using portcullis::ml::Smote;

        
TEST(smote, simple) {
    
    double data[] = {
        0.2, 0.4, 1.5, 2.3, 0.0,
        0.3, 0.3, 2.6, 5.2, 0.1,
        0.3, 0.3, 2.4, 5.2, 0.1,
        0.1, 0.2, 0.5, 3.1, 0.3,
        0.3, 0.3, 2.6, 5.2, 0.9,
        0.3, 0.3, 2.6, 5.2, 0.1,
        0.3, 0.7, 2.6, 5.2, 0.1,
        0.2, 0.3, 2.6, 4.2, 0.1,
        0.3, 0.8, 2.6, 2.2, 0.1,
        1.3, 1.3, 2.6, 8.2, 0.8
    };
    Smote smote(2, 2, 1, data, 10, 5);
    smote.setVerbose(true);
    //smote.execute();
    //EXPECT_EQ(smote.getNbSynthRows(), 20);
}

