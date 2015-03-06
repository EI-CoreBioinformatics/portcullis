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

#include <fstream>
#include <string>
#include <iostream>
#include <vector>
using std::boolalpha;
using std::ifstream;
using std::string;
using std::vector;

#include <boost/exception/all.hpp>
#include <boost/timer/timer.hpp>
#include <boost/filesystem.hpp>
using boost::timer::auto_cpu_timer;
using boost::lexical_cast;
using boost::filesystem::absolute;
using boost::filesystem::copy_file;
using boost::filesystem::remove;
using boost::filesystem::exists;
using boost::filesystem::create_symlink;
using boost::filesystem::create_directory;
using boost::filesystem::symbolic_link_exists;
namespace po = boost::program_options;

#include "junction_system.hpp"

#include <dlib/clustering.h>
using namespace dlib;

namespace portcullis {
    
typedef boost::error_info<struct ClusterError,string> ClusterErrorInfo;
struct ClusterException: virtual boost::exception, virtual std::exception { };

class Cluster {

private:
    
    string junctionFile;
    string resultsFile;
    bool verbose;
    
public:
    
    Cluster() : Cluster("", "", false) {}
    
    Cluster(const string& _junctionFile, const string& _resultsFile, bool _verbose) {
        junctionFile = _junctionFile;
        resultsFile = _resultsFile;
        verbose = _verbose;
        
        // Test if provided genome exists
        if (!exists(junctionFile)) {
            BOOST_THROW_EXCEPTION(FilterException() << FilterErrorInfo(string(
                        "Could not find junction file at: ") + junctionFile));
        }
    }
    
    virtual ~Cluster() {
    }
    
       

public:
    
    
    void cluster() {
        
        // Load junction system
        JunctionSystem js(junctionFile);
        
        cout << " - Loaded junctions from: " << junctionFile << endl;
        
        // Here the matrix type which we will need to populate with junction data
        typedef matrix<double,3,1> sample_type;
        
        // Now we are making a typedef for the kind of kernel we want to use.  I picked the
        // radial basis kernel because it only has one parameter and generally gives good
        // results without much fiddling.
        typedef radial_basis_kernel<sample_type> kernel_type;

        // Here we declare an instance of the kcentroid object.  It is the object used to 
        // represent each of the centers used for clustering.  The kcentroid has 3 parameters 
        // you need to set.  The first argument to the constructor is the kernel we wish to 
        // use.  The second is a parameter that determines the numerical accuracy with which 
        // the object will perform part of the learning algorithm.  Generally, smaller values 
        // give better results but cause the algorithm to attempt to use more dictionary vectors 
        // (and thus run slower and use more memory).  The third argument, however, is the 
        // maximum number of dictionary vectors a kcentroid is allowed to use.  So you can use
        // it to control the runtime complexity.  
        kcentroid<kernel_type> kc(kernel_type(0.1), 0.01, 50);
        
        // Now we make an instance of the kkmeans object and tell it to use kcentroid objects
        // that are configured with the parameters from the kc object we defined above.
        kkmeans<kernel_type> test(kc);

        std::vector<sample_type> samples;
        std::vector<sample_type> initial_centers;
        
        // Load feature vector from junction file into matrix
        for(JunctionPtr j : js.getJunctions()) {
            sample_type m;
            m(0) = j->getNbJunctionAlignments();
            m(1) = j->getEntropy();
            m(2) = j->getIntronSize();
            samples.push_back(m);
        }
        
        cout << " - Populated matrix with junction data" << endl;

        // tell the kkmeans object we made that we want to run k-means with k set to 3. 
        // (i.e. we want 2 clusters)
        test.set_number_of_centers(2);

        // You need to pick some initial centers for the k-means algorithm.  So here
        // we will use the dlib::pick_initial_centers() function which tries to find
        // n points that are far apart (basically).  
        pick_initial_centers(2, initial_centers, samples, test.get_kernel());

        cout << " - Picked initial centers" << endl;
        
        // now run the k-means algorithm on our set of samples.  
        test.train(samples,initial_centers);

        cout << " - Trained model" << endl;
        
        unsigned long results[samples.size()];
        uint32_t nbClass0 = 0;
        uint32_t nbClass1 = 0;
        
        // now loop over all our samples and print out their predicted class.  In this example
        // all points are correctly identified.
        for (unsigned long i = 0; i < samples.size(); ++i) {
            
            results[i] = test(samples[i]);
            
            if (results[i] == 0) {
                nbClass0++;
            }
            else if (results[i] == 1) {
                nbClass1++;
            }
        }
        
        // Now print out how many dictionary vectors each center used.  Note that 
        // the maximum number of 8 was reached.  If you went back to the kcentroid 
        // constructor and changed the 8 to some bigger number you would see that these
        // numbers would go up.  However, 8 is all we need to correctly cluster this dataset.
        cout << "num dictionary vectors for center 0: " << test.get_kcentroid(0).dictionary_size() << endl;
        cout << "num dictionary vectors for center 1: " << test.get_kcentroid(1).dictionary_size() << endl;
        
        cout << "num class 0: " << nbClass0 << endl;
        cout << "num class 1: " << nbClass1 << endl;
        
        
        if (!resultsFile.empty()) {
            cout << "Saving results to: " << resultsFile << endl;
            std::ofstream out(resultsFile.c_str());
            for (unsigned long i = 0; i < samples.size(); ++i) {
                out << results[i] << "\t" << i << "\t" << samples[i](0) << "\t" << samples[i](1) << "\t" << samples[i](2) << endl;
            }
            out.close();
        }
    }
    
    void pca() {
        
    }
  
    static string helpMessage() {
        return string("\nPortcullis Cluster Mode Help.\n\n") +
                      "Usage: portcullis cluster [options] <junction-file>\n\n" +
                      "Allowed options";
    }
    
    static int main(int argc, char *argv[]) {
        
        // Portcullis args
        string junctionFile;
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
                ;

        // Positional option for the input bam file
        po::positional_options_description p;
        p.add("junction-file", 1);
        
        
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
        
        

        auto_cpu_timer timer(1, "\nPortcullis cluster completed.\nTotal runtime: %ws\n\n");        

        cout << "Running portcullis in cluster mode" << endl
             << "--------------------------------" << endl << endl;
        
        // Create the prepare class
        Cluster cluster(junctionFile, resultsFile, verbose);
        cluster.cluster();
        
        return 0;
    }
};
}