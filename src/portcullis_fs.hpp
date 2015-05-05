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

#include <boost/filesystem.hpp>
#include <boost/filesystem/path.hpp>
using boost::filesystem::exists;
using boost::filesystem::path;

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
        
    
    public:
       
        PortcullisFS() {}
        
        /**
         * 
         * @param exe Full path to the exe, probably derived from argv0.
         */
        PortcullisFS(path exe) {
            
            portcullisExe = exe;
            
            rootDir = exe.parent_path().parent_path();
            
            binDir = path(rootDir);
            binDir /= "bin";
            
            etcDir = path(rootDir);
            etcDir /= "etc";
            
            scriptsDir = path(rootDir);
            scriptsDir /= "scripts";
            
            path srcDir = path(rootDir);
            srcDir /= "src";
            
                
            if (exe.parent_path() == srcDir) {
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
        
        path GetBinDir() const {
            return binDir;
        }

        path GetEtcDir() const {
            return etcDir;
        }

        path GetFilterJuncsPy() const {
            return filterJuncsPy;
        }

        path GetPortcullisExe() const {
            return portcullisExe;
        }

        path GetRootDir() const {
            return rootDir;
        }

        path GetSamtoolsExe() const {
            return samtoolsExe;
        }

        path GetScriptsDir() const {
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



