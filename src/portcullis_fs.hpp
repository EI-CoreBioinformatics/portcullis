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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif 

#include <iostream>
#include <string>
using std::cout;
using std::endl;
using std::string;

#include <boost/filesystem.hpp>
#include <boost/filesystem/path.hpp>
namespace bfs = boost::filesystem;
using bfs::exists;
using bfs::path;

namespace portcullis {
    
    typedef boost::error_info<struct FileSystemError,string> FileSystemErrorInfo;
    struct FileSystemException: virtual boost::exception, virtual std::exception { };

 
    class PortcullisFS {
    private:
        
        // Executables
        path portcullisExe;
        path samtoolsExe;
        path filterJuncsPy;
        
        // Directories
        path binDir;
        path scriptsDir;
        path etcDir;
        path rootDir;
        
        // Info
        string version;
        
        bool isAbsolute;
        bool isRelative;
        bool isOnPath;
        
        
    public:
       
        /**
         * Assume on PATH by default
         */
        PortcullisFS() : PortcullisFS("portcullis", PACKAGE_VERSION) {}
        
        /**
         * 
         * @param exe Full path to the exe, probably derived from argv0.
         */
        PortcullisFS(const char* argv, string _version) {
            
            isAbsolute = false;
            isRelative = false;
            isOnPath = false;
            
            path exe(argv);
            version = _version;
            
            if(exe.is_absolute()) {
                
                // Easy job... nothing special to do, resolve symlink then take two levels up
                portcullisExe = bfs::canonical(exe);
                isAbsolute=true;
            }
            else if (exe.string().find('/') != string::npos) {
                
                // Relative with some parent paths... get absolute path, resolving symlinks then take two levels up
                portcullisExe = bfs::canonical(bfs::system_complete(exe));
                isRelative=true;
            }
            else {

                portcullisExe = exe;
                isOnPath=true;
            }
            
            // We assume scripts are on the path if exe was on the path
            if (isAbsolute || isRelative) {
                // Check to see if scripts are adjacent to exe first
                path pb(portcullisExe.parent_path());
                pb /= "bed12.py";
                if (exists(pb)) {
                    scriptsDir = portcullisExe.parent_path();
                }
                else {
                
                    // Not 100% sure how far back we need to go (depends on whether using KAT exe or tests) 
                    // so try 2, 3 and 4 levels.
                    scriptsDir = portcullisExe.parent_path().parent_path();
                    scriptsDir /= "scripts";                 

                    if (!exists(scriptsDir)) {
                        scriptsDir = portcullisExe.parent_path().parent_path().parent_path();
                        scriptsDir /= "scripts";       
                        
                        if (!exists(scriptsDir)) {
                            scriptsDir = portcullisExe.parent_path().parent_path().parent_path().parent_path();
                            scriptsDir /= "scripts";        

                            if (!exists(scriptsDir)) {
                                BOOST_THROW_EXCEPTION(FileSystemException() << FileSystemErrorInfo(string(
                                    "Could not find suitable directory containing Portcullis scripts relative to provided exe: ") + portcullisExe.c_str()));
                            }
                        }
                    }

                    // Also double check the kat_distanalysis.py script exists
                    pb = scriptsDir;
                    pb /= "bed12.py";

                    if (!exists(pb)) {
                        BOOST_THROW_EXCEPTION(FileSystemException() << FileSystemErrorInfo(string(
                            "Found the scripts directory where expected") + scriptsDir.string() + 
                                ". However, could not find the \"bed12.py\" script inside."));
                    }
                }
            
            }
            
            binDir = path(rootDir);
            binDir /= "bin";
            
            etcDir = path(rootDir);
            etcDir /= "etc";
            
            scriptsDir = path(rootDir);
            scriptsDir /= "scripts";
            
            path srcDir = path(rootDir);
            srcDir /= "src";
            
            path testDir = path(rootDir);
            testDir /= "tests";
            
                
            if (portcullisExe.parent_path() == srcDir || portcullisExe.parent_path() == testDir) {
                samtoolsExe = path(rootDir);
                samtoolsExe /= "deps/samtools-1.2/portcullis_samtools"; 
                
                filterJuncsPy = path(scriptsDir);
                filterJuncsPy /= "filter_junctions.py";
            }
            else {
                samtoolsExe = path(binDir);
                samtoolsExe /= "samtools";
                
                filterJuncsPy = path(binDir);
                filterJuncsPy /= "filter_junctions.py";
            }
            
            if (!exists(samtoolsExe)) {
                BOOST_THROW_EXCEPTION(FileSystemException() << FileSystemErrorInfo(string(
                    "Could not find samtools executable at: ") + samtoolsExe.c_str()));
            }
            
            if (!exists(filterJuncsPy)) {
                BOOST_THROW_EXCEPTION(FileSystemException() << FileSystemErrorInfo(string(
                    "Could not find filter_junctions.py script at: ") + filterJuncsPy.c_str()));
            }
        }
        
        
        
        
        // **** Destructor ****
        virtual ~PortcullisFS() {

        }
        
        path getBinDir() const {
            return binDir;
        }

        path getEtcDir() const {
            return etcDir;
        }

        path getFilterJuncsPy() const {
            return filterJuncsPy;
        }

        path getPortcullisExe() const {
            return portcullisExe;
        }

        path getRootDir() const {
            return rootDir;
        }

        path getSamtoolsExe() const {
            return samtoolsExe;
        }

        path getScriptsDir() const {
            return scriptsDir;
        }
        
        string getVersion() const {
            return version;
        }


        
        friend std::ostream& operator<<(std::ostream &strm, const PortcullisFS& pfs) {
            
            return strm << "Directories: "<< endl
                        << " - Root: " << pfs.rootDir << endl
                        << " - Bin: " << pfs.binDir << endl
                        << " - Etc: " << pfs.etcDir << endl
                        << " - Scripts: " << pfs.scriptsDir << endl << endl
                        << "Executables: " << endl
                        << " - portcullis: " << pfs.portcullisExe << endl
                        << " - samtools: " << pfs.samtoolsExe << endl
                        << " - filter_junctions.py: " << pfs.filterJuncsPy << endl
                        << "Info: << endl"
                        << " - Version: " << pfs.version << endl;
        }     
    };
    
       
    
}



