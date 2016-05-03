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

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
using std::cout;
using std::endl;
using std::string;
using std::stringstream;

#include <boost/filesystem.hpp>
using boost::filesystem::path;

#include <Python.h>

#include <portcullis/python_exe.hpp>

wchar_t* portcullis::PythonExe::convertCharToWideChar(const char* c) {
    const size_t cSize = strlen(c)+1;
    wchar_t* wc = new wchar_t[cSize];
    mbstowcs (wc, c, cSize);    
    return wc;
}       

void portcullis::PythonExe::executePythonScript(const path& scripts_dir, const string& script_name, const vector<string>& args) {
    
    const path full_script_path = path(scripts_dir.string() + "/" + script_name);
    
    stringstream ss;
    
    // Create wide char alternatives
    wchar_t* wsn = PythonExe::convertCharToWideChar(script_name.c_str());
    wchar_t* wsp = PythonExe::convertCharToWideChar(full_script_path.c_str());    
    wchar_t* wargv[10]; // Can't use variable length arrays!
    wargv[0] = wsp;
    ss << full_script_path.c_str();
    for (size_t i = 1; i <= args.size(); i++) {
        wargv[i] = PythonExe::convertCharToWideChar(args[i-1].c_str());
        ss << " " << args[i-1];
    }
        
    cout << endl << "Effective command line: " << ss.str() << endl << endl;
    
    std::ifstream script_in(full_script_path.c_str());
    string contents((std::istreambuf_iterator<char>(script_in)), std::istreambuf_iterator<char>());

    // Run python script
    Py_Initialize();
    Py_SetProgramName(wsp);
    PySys_SetArgv(4, wargv);
    PyRun_SimpleString(contents.c_str());
    Py_Finalize();

    // Cleanup
    delete wsn;
    // No need to free up "wsp" as it is element 0 in the array
    for(int i = 0; i < 4; i++) {
        delete wargv[i];
    }
}
