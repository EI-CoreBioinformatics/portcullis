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

#include <ranger/DataDouble.h>

DataDouble::DataDouble() :
data(0) {
}

DataDouble::DataDouble(double* data_double, std::vector<std::string> variable_names, size_t num_rows, size_t num_cols) : 
        DataDouble(variable_names, num_rows, num_cols) {
    
    for (size_t i = 0; i < num_cols; ++i) {
        for (size_t j = 0; j < num_rows; ++j) {
            data[i * num_rows + j] = data_double[i * num_rows + j];
        }
    }
}

DataDouble::~DataDouble() {
    delete[] data;
}
