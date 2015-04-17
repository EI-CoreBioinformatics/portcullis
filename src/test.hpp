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

#pragma once

#include <iostream>
#include <vector>
#include <dlib/svm.h>
#include <dlib/rand.h>

#include "junction_system.hpp"

using namespace std;
using namespace dlib;

namespace portcullis {
    
typedef boost::error_info<struct TestError,string> TestErrorInfo;
struct TestException: virtual boost::exception, virtual std::exception { };

class Test {

private:
    
    string junctionFile;
    string modelFile;
    string resultsFile;
    bool verbose;
    
public:
    
    Test() : Test("", "", "", false) {}
    
    Test(const string& _junctionFile, const string& _modelFile, const string& _resultsFile, bool _verbose) {
        junctionFile = _junctionFile;
        modelFile = _modelFile;
        resultsFile = _resultsFile;
        verbose = _verbose;
        
        // Test if provided genome exists
        if (!exists(junctionFile)) {
            BOOST_THROW_EXCEPTION(TestException() << TestErrorInfo(string(
                        "Could not find junction file at: ") + junctionFile));
        }
        
        // Test if provided genome exists
        /*if (!exists(modelFile)) {
            BOOST_THROW_EXCEPTION(TestException() << TestErrorInfo(string(
                        "Could not find model file at: ") + modelFile));
        }*/
    }
    
    virtual ~Test() {
    }
    
       

public:
    
    
    void test() {
        
        // Load junction system
        JunctionSystem js(junctionFile);
        
        cout << " - Loaded " << js.getJunctions().size() << " samples to test" << endl;
        
        // Load one-class SVM decision function
        typedef matrix<double,0,1> sample_type;
        typedef radial_basis_kernel<sample_type> kernel_type;
        decision_function<linear_kernel<sample_type> > df;
        deserialize(modelFile + ".df") >> df;
        empirical_kernel_map<kernel_type> ekm;
        deserialize(modelFile + ".ekm") >> ekm;
        
        uint32_t no = 0;
        uint32_t yes = 0;
        sample_type bias(1);
        bias = -1;
        
        JunctionSystem good;
        JunctionSystem bad;
        
        for(JunctionPtr j : js.getJunctions()) {
            sample_type m(6);
        
            // Calculate the metrics and add them into the vector
            m(0) = j->getCoverage();
            m(1) = j->getEntropy();
            m(2) = j->getIntronSize();
            m(3) = j->getMaxMMES();
            m(4) = j->getHammingDistance3p();
            m(5) = j->getHammingDistance5p();
            
            if (df(join_cols(ekm.project(m), bias)) > -1) {
                good.addJunction(j);                
            }
            else {
                bad.addJunction(j);
            }
            
        }
        
        cout << "Found " << bad.getJunctions().size() << " bad junctions to discard" << endl;
        cout << "Found " << good.getJunctions().size() << " good junctions to keep" << endl;
        
        // Print junction stats to file
        string goodFile = resultsFile + ".good";
        string badFile = resultsFile + ".bad";
        
        good.saveAll(goodFile);
        bad.saveAll(badFile);        
    }
    
    
    static string helpMessage() {
        return string("\nPortcullis Testing Mode Help.\n\n") +
                      "This is mode is intended to test potential junctions against a pre-made one class SVM decision function.\n" +
                      "Usage: portcullis test [options] <junction-file> <model-file>\n\n" +
                      "Allowed options";
    }
    
    static int main(int argc, char *argv[]) {
        
        // Portcullis args
        string junctionFile;
        string modelFile;
        string resultsFile;
        bool verbose;
        bool help;
        
        // Declare the supported options.
        po::options_description generic_options(helpMessage());
        generic_options.add_options()
                ("output,o", po::value<string>(&resultsFile)->default_value(""), 
                    "File name for tabular results file generated by this program.")
                ("verbose,v", po::bool_switch(&verbose)->default_value(false), 
                    "Print extra information")
                ("help", po::bool_switch(&help)->default_value(false), "Produce help message")
                ;

        // Hidden options, will be allowed both on command line and
        // in config file, but will not be shown to the user.
        po::options_description hidden_options("Hidden options");
        hidden_options.add_options()
                ("junction-file,g", po::value<string>(&junctionFile), "Path to the junction file to process.")
                ("model-file,m", po::value<string>(&modelFile), "Path to the one-class SVM model file to process.")
                ;

        // Positional options
        po::positional_options_description p;
        p.add("junction-file", 1);
        p.add("model-file", 1);
        
        
        // Combine non-positional options
        po::options_description cmdline_options;
        cmdline_options.add(generic_options).add(hidden_options);

        // Parse command line
        po::variables_map vm;
        po::store(po::command_line_parser(argc, argv).options(cmdline_options).positional(p).run(), vm);
        po::notify(vm);

        // Output help information the exit if requested
        if (help || argc <= 1) {
            cout << generic_options << endl;
            return 1;
        }
        
        

        auto_cpu_timer timer(1, "\nPortcullis training completed.\nTotal runtime: %ws\n\n");        

        cout << "Running portcullis in training mode" << endl
             << "--------------------------------" << endl << endl;
        
        // Create the prepare class
        Test tester(junctionFile, modelFile, resultsFile, verbose);
        tester.test();
        
        return 0;
    }
};
}


