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


namespace portculis {
    
typedef boost::error_info<struct FilterError,string> FilterErrorInfo;
struct FilterException: virtual boost::exception, virtual std::exception { };

class Filter {

private:
    
    string junctionFile;
    string bamFile;
    bool verbose;
    
public:
    
    Filter() : Filter("", "", false) {}
    
    Filter(const string& _junctionFile, const string& _bamFile, bool _verbose) {
        junctionFile = _junctionFile;
        bamFile = _bamFile;
        verbose = _verbose;
        
        // Test if provided genome exists
        if (!exists(junctionFile)) {
            BOOST_THROW_EXCEPTION(FilterException() << FilterErrorInfo(string(
                        "Could not find junction file at: ") + junctionFile));
        }
        
        // Test if provided genome exists
        if (!exists(bamFile)) {
            BOOST_THROW_EXCEPTION(FilterException() << FilterErrorInfo(string(
                        "Could not find BAM file at: ") + bamFile));
        }
    }
    
    virtual ~Filter() {
    }
    
       

public:
    
    
    void filter() {
        
        // Load junction system
        JunctionSystem js(junctionFile);
        
        
    }
  
    static string helpMessage() {
        return string("\nPortculis Filter Mode Help.\n\n") +
                      "Usage: portculis filter [options] <junction-file> <bam-file> \n\n" +
                      "Allowed options";
    }
    
    static int main(int argc, char *argv[]) {
        
        // Portculis args
        string junctionFile;
        string bamFile;
        bool verbose;
        bool help;
        
        // Declare the supported options.
        po::options_description generic_options(helpMessage());
        generic_options.add_options()
                ("verbose,v", po::bool_switch(&verbose)->default_value(false), 
                    "Print extra information")
                ("help", po::bool_switch(&help)->default_value(false), "Produce help message")
                ;

        // Hidden options, will be allowed both on command line and
        // in config file, but will not be shown to the user.
        po::options_description hidden_options("Hidden options");
        hidden_options.add_options()
                ("junction-file,g", po::value<string>(&junctionFile), "Path to the junction file to process.")
                ("bam-file,g", po::value<string>(&bamFile), "Path to the BAM file to filter.")
                ;

        // Positional option for the input bam file
        po::positional_options_description p;
        p.add("junction-file", 1);
        p.add("bam-file", 1);
        
        
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
        
        

        auto_cpu_timer timer(1, "\nPortculis filter completed.\nTotal runtime: %ws\n\n");        

        cout << "Running portculis in filter mode" << endl
             << "--------------------------------" << endl << endl;
        
        // Create the prepare class
        Filter filter(junctionFile, bamFile, verbose);
        filter.filter();
        
        return 0;
    }
};
}