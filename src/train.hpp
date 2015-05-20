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

using namespace std;
using namespace dlib;

namespace portcullis {
    
typedef boost::error_info<struct TrainError,string> TrainErrorInfo;
struct TrainException: virtual boost::exception, virtual std::exception { };

class Train {

private:
    
    string junctionFile;
    string resultsFile;
    double width;
    bool verbose;
    
public:
    
    Train() : Train("", "", 4.0, false) {}
    
    Train(const string& _junctionFile, const string& _resultsFile, double _width, bool _verbose) {
        junctionFile = _junctionFile;
        resultsFile = _resultsFile;
        width = _width;
        verbose = _verbose;
        
        // Test if provided genome exists
        if (!exists(junctionFile)) {
            BOOST_THROW_EXCEPTION(TrainException() << TrainErrorInfo(string(
                        "Could not find junction file at: ") + junctionFile));
        }
    }
    
    virtual ~Train() {
    }
    
       

public:
    
    
    void train() {
        
        // Load junction system
        JunctionSystem js(junctionFile);
        
        // But putting the empirical_kernel_map aside, the most important step in turning a
        // linear SVM into a one-class SVM is the following.  We append a -1 value onto the end
        // of each feature vector and then tell the trainer to force the weight for this
        // feature to 1.  This means that if the linear SVM assigned all other weights a value
        // of 0 then the output from a learned decision function would always be -1.  The
        // second step is that we ask the SVM to label each training sample with +1.  This
        // causes the SVM to set the other feature weights such that the training samples have
        // positive outputs from the learned decision function.  But the starting bias for all
        // the points in the whole feature space is -1.  The result is that points outside our
        // training set will not be affected, so their outputs from the decision function will
        // remain close to -1.

        // We will use column vectors to store our points.  Here we make a convenient typedef
        // for the kind of vector we will use.
        typedef matrix<double,0,1> sample_type;

        // Use the radial basis kernel, which is normally quite effective.
        typedef radial_basis_kernel<sample_type> kernel_type;

        
        std::vector<sample_type> samples;
        sample_type m(6);
        
        // Load feature vector from junction file into matrix
        for(JunctionPtr j : js.getJunctions()) {
        
            // Calculate the metrics and add them into the vector
            m(0) = j->getCoverage();
            m(1) = j->getEntropy();
            m(2) = j->getIntronSize();
            m(3) = j->getMaxMMES();
            m(4) = j->getHammingDistance3p();
            m(5) = j->getHammingDistance5p();

            if (verbose)
                cout << m << endl;
            
            samples.push_back(m);
        }
        
        cout << " - Created feature vector from " << samples.size() << " samples" << endl;

        // Set the width of the radial basis kernel to 4 (may need adjusting later).
        // Larger values make the width smaller and give the radial basis kernel more 
        // resolution. 
        linearly_independent_subset_finder<kernel_type> lisf(kernel_type(width), 500);
        // populate lisf with samples.  We have configured it to allow at most 50 samples but this function 
        // may determine that fewer samples are necessary to form a good basis.  In this example program
        // it will select only 26.
        fill_lisf(lisf, samples);
    
        cout << " - Created linearly independent subset (LIS) from samples using radial kernel of width 4.0.  LIS contains: " << lisf.size() << " samples" << endl;
        
        // Use an empirical kernel map for this problem.  These are normally effective
        // when number of dimension are < 100.
        empirical_kernel_map<kernel_type> ekm;
        
        // Load LISF samples into the EKM.
        ekm.load(lisf);
        
        cout << " - Loaded LIS as basis for Empirical Kernel Map" << endl;
        
        // Set the label vector.  These will all be set to 1 as they should all be genuine.
        std::vector<double> labels;
        
        // make a vector with just 1 element in it equal to -1.
        sample_type bias(1);
        bias = -1;
        sample_type augmented;
        std::vector<sample_type> augmentedSamples;
        
        for(sample_type s : samples) {
            
            // Apply the empirical_kernel_map transformation on the existing sample data and 
            // then append the -1 value as the last column
            augmented = join_cols(ekm.project(s), bias);
            augmentedSamples.push_back(augmented);
            labels.push_back(+1);       // All these examples are real junctions!             
        }
        
        // Save memory
        samples.clear();
        
        cout << " - Created augmented and labelled feature vector suitable for one-class SVM" << endl;

        // The svm_c_linear_dcd_trainer is a very fast SVM solver which only works with the
        // linear_kernel.  Use "force_last_weight_to_1" to make this a one class SVM.  As mentioned before
        // the last weight will always be -1.  Implying that we need weighting applied from other weights
        // for a real junction to be predicted.
        svm_c_linear_dcd_trainer<linear_kernel<sample_type> > linear_trainer;
        linear_trainer.force_last_weight_to_1(true);

        // Train the one-class SVM
        decision_function<linear_kernel<sample_type> > df = linear_trainer.train(augmentedSamples, labels);

        cout << " - Trained one-class SVM" << endl;
                
        // Save the trained decision function to disk for use later
        serialize(resultsFile + ".df") << df;
        cout << " - Saved decision function to file: " << resultsFile << ".df" << endl;
        
        serialize(resultsFile + ".ekm") << ekm;
        cout << " - Saved Empirical Kernel Map to file: " << resultsFile << ".ekm" << endl;
        
        
        uint32_t no = 0;
        uint32_t yes = 0;
        
        for(JunctionPtr j : *(js.getJunctions()) {
            sample_type m(6);
        
            // Calculate the metrics and add them into the vector
            m(0) = j->getNbDistinctAlignments();
            m(1) = j->getEntropy();
            m(2) = j->getIntronSize();
            m(3) = j->getMaxMMES();
            m(4) = j->getHammingDistance3p();
            m(5) = j->getHammingDistance5p();
            
            if (df(join_cols(ekm.project(m), bias)) > -1) {
                yes++;                
            }
            else {
                no++;
            }
            
        }
        
        cout << "No: " << no << endl;
        cout << "Yes: " << yes << endl;
        
    }
        
    
    static string helpMessage() {
        return string("\nPortcullis Training Mode Help.\n\n") +
                      "This is mode is intended to train a one class SVM with genuine junctions.\n" +
                      "The trained model can then be used to filter potential junctions from real datasets\n\n" +
                      "Usage: portcullis train [options] <junction-file>\n\n" +
                      "Allowed options";
    }
    
    static int main(int argc, char *argv[]) {
        
        // Portcullis args
        string junctionFile;
        string resultsFile;
        double width;
        bool verbose;
        bool help;
        
        // Declare the supported options.
        po::options_description generic_options(helpMessage());
        generic_options.add_options()
                ("output,o", po::value<string>(&resultsFile)->default_value(""), 
                    "File name for tabular results file generated by this program.")
                ("width,w", po::value<double>(&width)->default_value(4.0), 
                    "The width to give the radial basis kernel.")
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
        
        

        auto_cpu_timer timer(1, "\nPortcullis training completed.\nTotal runtime: %ws\n\n");        

        cout << "Running portcullis in training mode" << endl
             << "--------------------------------" << endl << endl;
        
        // Create the prepare class
        Train trainer(junctionFile, resultsFile, width, verbose);
        trainer.train();
        
        return 0;
    }
};
}


