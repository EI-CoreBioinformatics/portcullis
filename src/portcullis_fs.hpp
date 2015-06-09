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
        
        
        string exec(const char* cmd) {
            FILE* pipe = popen(cmd, "r");
            if (!pipe) return "ERROR";
            char buffer[512];
            string result = "";
            while(!feof(pipe)) {
                if(fgets(buffer, 512, pipe) != NULL)
                        result += buffer;
            }
            pclose(pipe);
            return result;
        }
    
    public:
       
        PortcullisFS() {}
        
        /**
         * 
         * @param exe Full path to the exe, probably derived from argv0.
         */
        PortcullisFS(const char* argv) {
            
            path exe(argv);
            
            cout << exe << endl;
            
            if(exe.is_absolute()) {
                
                //cout << "Absolute" << endl;
                
                // Easy job... nothing special to do, resolve symlink then take two levels up
                portcullisExe = bfs::canonical(exe);
                rootDir = portcullisExe.parent_path().parent_path();
            }
            else if (exe.string().find('/') != string::npos) {
                
                //cout << "Relative" << endl;
                
                // Relative with some parent paths... get absolute path, resolving symlinks then take two levels up
                portcullisExe = bfs::canonical(bfs::system_complete(exe));
                rootDir = portcullisExe.parent_path().parent_path();
            }
            else {

                //cout << "name only" << endl;
                
                // Tricky one... just exe name, no indication of where if comes from. Now we have to resort to using which.
                string cmd = string("which ") + exe.string();
                string res = exec(cmd.c_str());
                string fullpath = res.substr(0, res.length() - 1);

                //cout << "fullpath" << fullpath << endl;
                portcullisExe = bfs::canonical(path(fullpath));
                rootDir = portcullisExe.parent_path().parent_path();
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
                samtoolsExe /= "deps/samtools-1.2/samtools"; 
                
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

        
        friend std::ostream& operator<<(std::ostream &strm, const PortcullisFS& pfs) {
            
            return strm << "Directories: "<< endl
                        << " - Root: " << pfs.rootDir << endl
                        << " - Bin: " << pfs.binDir << endl
                        << " - Etc: " << pfs.etcDir << endl
                        << " - Scripts: " << pfs.scriptsDir << endl << endl
                        << "Executables: " << endl
                        << " - portcullis: " << pfs.portcullisExe << endl
                        << " - samtools: " << pfs.samtoolsExe << endl
                        << " - filter_junctions.py: " << pfs.filterJuncsPy << endl;
        }     
    };
    
       
    
}



