/*-------------------------------------------------------------------------------
This file is part of Ranger.
    
Ranger is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Ranger is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Ranger. If not, see <http://www.gnu.org/licenses/>.

Written by: 

Marvin N. Wright
Institut für Medizinische Biometrie und Statistik
Universität zu Lübeck
Ratzeburger Allee 160
23562 Lübeck 

http://www.imbs-luebeck.de
wright@imbs.uni-luebeck.de
#-------------------------------------------------------------------------------*/

#include <limits.h>
#include <math.h>
#include <iostream>
#include <vector>

#include <ranger/DataChar.h>

DataChar::DataChar() :
data(0) {
}

DataChar::DataChar(double* data_double, std::vector<std::string> variable_names, 
        size_t num_rows, size_t num_cols, bool& error) : 
        DataChar(variable_names, num_rows, num_cols) {

    // Save data and report errors
    for (size_t i = 0; i < num_cols; ++i) {
        for (size_t j = 0; j < num_rows; ++j) {
            double value = data_double[i * num_rows + j];
            if (value > CHAR_MAX || value < CHAR_MIN) {
                error = true;
            }
            if (floor(value) != ceil(value)) {
                error = true;
            }
            data[i * num_rows + j] = (char) value;
        }
    }
}

DataChar::~DataChar() {
    delete[] data;
}
