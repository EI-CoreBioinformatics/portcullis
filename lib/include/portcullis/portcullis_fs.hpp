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
        
        // Directories
        path dataDir;
        path scriptsDir;
        
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
                dataDir = path(DATADIR "/portcullis/");
            }
            
            // We assume scripts are on the path if exe was on the path
            if (isAbsolute || isRelative) {
                // Check to see if scripts are adjacent to exe first
                path pb(portcullisExe.parent_path());
                pb /= "bed12.py";
                if (exists(pb)) {
                    scriptsDir = portcullisExe.parent_path();
                    dataDir = path(DATADIR "/portcullis/");
                }
                else {
                
                    // Not 100% sure how far back we need to go (depends on whether 
                    // using portcullis exe or tests) 
                    // so try 2, 3 and 4 levels.
                    scriptsDir = portcullisExe.parent_path().parent_path();
                    scriptsDir /= "scripts";
                    
                    dataDir = portcullisExe.parent_path().parent_path();
                    dataDir /= "data";                 


                    if (!exists(scriptsDir)) {
                        scriptsDir = portcullisExe.parent_path().parent_path().parent_path();
                        scriptsDir /= "scripts";       
                        
                        dataDir = portcullisExe.parent_path().parent_path().parent_path();
                        dataDir /= "data";       
                        
                        if (!exists(scriptsDir)) {
                            scriptsDir = portcullisExe.parent_path().parent_path().parent_path().parent_path();
                            scriptsDir /= "scripts";        

                            dataDir = portcullisExe.parent_path().parent_path().parent_path().parent_path();
                            dataDir /= "data";        

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
                    
                    path df = dataDir;
                    df /= "default_filter.json";
                    
                    if (!exists(df)) {
                        BOOST_THROW_EXCEPTION(FileSystemException() << FileSystemErrorInfo(string(
                            "Found the data directory where expected") + dataDir.string() + 
                                ". However, could not find the \"default_filter.json\" configuraton file inside."));
                    }
                    
                }
            
            }
        }
        
        
        
        
        // **** Destructor ****
        virtual ~PortcullisFS() {

        }
        
        path getPortcullisExe() const {
            return portcullisExe;
        }

        path getScriptsDir() const {
            return scriptsDir;
        }
        
        path getDataDir() const {
            return dataDir;
        }
        
        string getVersion() const {
            return version;
        }


        
        friend std::ostream& operator<<(std::ostream &strm, const PortcullisFS& pfs) {
            
            return strm << "Executables: " << endl
                        << " - portcullis: " << pfs.portcullisExe << endl
                        << "Directories: " << endl
                        << " - Data: " << pfs.dataDir << endl
                        << " - Scripts: " << pfs.scriptsDir << endl
                        << "Info:" << endl
                        << " - Version: " << pfs.version << endl;
        }     
    };
    
       
    
}



