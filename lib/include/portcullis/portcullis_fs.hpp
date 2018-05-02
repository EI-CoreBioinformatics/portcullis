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

#ifdef OS_MAC
#include <mach-o/dyld.h>
#endif

#include <boost/filesystem.hpp>
#include <boost/filesystem/path.hpp>
namespace bfs = boost::filesystem;
using bfs::exists;
using bfs::path;

namespace portcullis
{

typedef boost::error_info<struct FileSystemError, string> FileSystemErrorInfo;
struct FileSystemException: virtual boost::exception, virtual std::exception { };


class PortcullisFS
{
private:

    // Executables
    path portcullisExe;
    path canonicalExe;

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
    PortcullisFS() {}

    /**
     *
     * @param exe Full path to the exe, probably derived from argv0.
     */
    PortcullisFS(const char* argv)
    {
        version = "X.X.X";
        isAbsolute = false;
        isRelative = false;
        isOnPath = false;
        portcullisExe = argv;
        if (portcullisExe.is_absolute())
        {
            // Absolute path provided...  Easy job... nothing special to do,
            // resolve symlink then take two levels up
            canonicalExe = bfs::canonical(portcullisExe);
            isAbsolute = true;
        }
        else if (portcullisExe.string().find('/') != string::npos)
        {
            // Relative with some parent paths... get absolute path, resolving
            // symlinks then take two levels up
            canonicalExe = bfs::canonical(bfs::system_complete(portcullisExe));
            isRelative = true;
        }
        else
        {
            // Only name provided
            // In this case we just have to assume everything is properly installed
#ifdef OS_LINUX
            canonicalExe = do_readlink();
#elif OS_MAC
            canonicalExe = get_mac_exe();
#else
            canonicalExe = do_readlink(); // Assume linux
#endif
            isOnPath = true;
        }

#ifdef DATADIR
		path customDataDir(DATADIR);
#else
		path customDataDir("");
#endif


		// If python is installed we need to figure out where the scripts are located relative to the
		// running executable.  This can be in various different places depending on how everything is
		// setup: installed kat, running compiled binary from source directory, or running unit tests.

        // First get the executable directory
        path exe_dir(canonicalExe.parent_path());
			
		if (exe_dir.stem().string() == "bin") {

			// Ok, so we are in a installed location.  Figuring out the scripts directory isn't as straight
			// forward as it may seem because we might have installed to a alternate root.  So wind back the 
			// exec_prefix to get to root (or alternate root) directory.
			path ep(EXECPREFIX);
			path root = ep;
			path altroot = exe_dir.parent_path();
			while (root.has_parent_path()) {
				root = root.parent_path();
				altroot = altroot.parent_path();					
			}
			this->dataDir = altroot / customDataDir / "portcullis";
   	        this->scriptsDir = this->dataDir / "scripts";
		}
		else if (exe_dir.stem().string() == "src") {
			// Presumably if we are here then we are running the kat executable from the source directory
    	    if (exists(exe_dir / "portcullis.cc") || exists(exe_dir / "check_portcullis.cc")) {
        		this->dataDir = exe_dir.parent_path() / "data";
        		this->scriptsDir = exe_dir.parent_path() / "scripts" / "portcullis";
			}
		}
		else {
			this->dataDir = customDataDir;
			this->scriptsDir = this->dataDir / "scripts";
		}


		// Validate the existence of the scripts directory and scripts file.				
        if (!exists(this->dataDir)) {
        	BOOST_THROW_EXCEPTION(FileSystemException() << FileSystemErrorInfo(string(
                	"Could not find suitable directory containing Portcullis data files at the expected location: ") + this->dataDir.string()));
        }	
        path df = this->dataDir / "default_filter.json";
        if (!exists(df)) {
            BOOST_THROW_EXCEPTION(FileSystemException() << FileSystemErrorInfo(string(
                                  "Found the data directory where expected ") + dataDir.string() +
                                  ". However, could not find the \"default_filter.json\" at: " + df.string()));
        }

        if (!exists(this->scriptsDir)) {
        	BOOST_THROW_EXCEPTION(FileSystemException() << FileSystemErrorInfo(string(
                	"Could not find suitable directory containing Portcullis scripts at the expected location: ") + this->scriptsDir.string()));
        }	
        path prf = this->scriptsDir / "portcullis/rule_filter.py";
        if (!exists(prf)) {
            BOOST_THROW_EXCEPTION(FileSystemException() << FileSystemErrorInfo(string(
                                    "Found the scripts directory where expected: ") + this->scriptsDir.string() +
                                    ". However, could not find \"rule_filter.py\" script at: " + prf.string()));
        }

    }


    std::string do_readlink() {
        char buff[PATH_MAX];
        ssize_t len = ::readlink("/proc/self/exe", buff, sizeof(buff) - 1);
        if (len != -1)
        {
            buff[len] = '\0';
            return std::string(buff);
        }
        BOOST_THROW_EXCEPTION(FileSystemException() << FileSystemErrorInfo(string(
                                  "Could not find locations of executable from /proc/self/exe")));
    }

#ifdef OS_MAC
    std::string get_mac_exe() {
        char path[1024];
        uint32_t size = sizeof(path);
        _NSGetExecutablePath(path, &size);
        return path;
    }
#endif


    // **** Destructor ****
    virtual ~PortcullisFS() {}

    path getPortcullisExe() const
    {
        return portcullisExe;
    }

    path getDataDir() const
    {
        return dataDir;
    }

    path getScriptsDir() const
    {
        return scriptsDir;
    }

    void setVersion(string version)
    {
        this->version = version;
    }

    string getVersion() const
    {
        return version;
    }

    path GetCanonicalExe() const
    {
        return canonicalExe;
    }

    bool IsAbsolute() const
    {
        return isAbsolute;
    }

    bool IsOnPath() const
    {
        return isOnPath;
    }

    bool IsRelative() const
    {
        return isRelative;
    }

    friend std::ostream& operator<<(std::ostream &strm, const PortcullisFS& pfs)
    {
        return strm << "Executables: " << endl
               << " - portcullis: " << pfs.portcullisExe << " - " << pfs.canonicalExe << endl
               << "Directories: " << endl
               << " - Data: " << pfs.dataDir << endl
               << " - Scripts: " << pfs.scriptsDir << endl
               << "Info:" << endl
               << " - Version: " << pfs.version << endl;
    }
};


// Make available everywhere
extern PortcullisFS pfs;
}



