/**************************************************************************
 *                                                                        *
 *   Author: Ivo Filot <i.a.w.filot@tue.nl>                               *
 *                                                                        *
 *   DEN2BIN is free software:                                            *
 *   you can redistribute it and/or modify it under the terms of the      *
 *   GNU General Public License as published by the Free Software         *
 *   Foundation, either version 3 of the License, or (at your option)     *
 *   any later version.                                                   *
 *                                                                        *
 *   DEN2BIN is distributed in the hope that it will be useful,           *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty          *
 *   of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.              *
 *   See the GNU General Public License for more details.                 *
 *                                                                        *
 *   You should have received a copy of the GNU General Public License    *
 *   along with this program.  If not, see http://www.gnu.org/licenses/.  *
 *                                                                        *
 **************************************************************************/

#ifndef _DENSITY_H
#define _DENSITY_H

#include <string>
#include <vector>
#include <iostream>
#include <ios>
#include <sstream>
#include <fstream>
#include <math.h>

#include <boost/regex.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <glm/glm.hpp>

#include "float_parser.h"

/**
 * @brief      class that can read density file and store it internally
 */
class Density{
private:
    std::string filename;               //!<input filename
    float scalar;                       //!<scalar
    float mat[3][3];                    //!<matrix dimensions
    float imat[3][3];                   //!<inverse of matrix
    glm::mat3 mat33;                    //!<matrix
    float volume;                       //!<unit cell volume

    unsigned int grid_dimensions[3];    //!<dimensions of the density grid
    std::vector<unsigned int> nrat;     //!<number of atoms in density file
    std::string gridline;
    std::vector<float> gridptr;         //!<grid to first pos of float array
    std::vector<float> gridptr2;        //!<grid to first pos of float array
    unsigned int gridsize;              //!<size of the grid
    bool vasp5_input;                   //!<whether file is vasp 5 input
    bool has_read;                      //!<whether file had been read
    bool header_read;                   //!<whether header has been read
    std::ifstream infile;               //!<infile stream
    unsigned int value_count;

public:

    /**
     * @brief      default constructor
     *
     * @param[in]  _filename  input filename
     */
    Density(const std::string &_filename);

    /**
     * @brief      get the unitcell matrix
     *
     * @return     The unitcell matrix.
     */
    glm::mat3 get_unitcell_matrix() {
        glm::mat3 out;
        for(unsigned int i=0; i<3; i++) {
            for(unsigned int j=0; j<3; j++) {
                out[i][j] = mat[i][j];
            }
        }

        return out;
    }

    /**
     * @brief      read number of lines
     *
     * @param[in]  lines  number of lines to read
     *
     * @return     percentage of file read
     */
    float read(unsigned int lines);

    /**
     * @brief      get minimum value in the density
     *
     * @return     minimum value
     */
    float minval() const;

    /**
     * @brief      get maximum value in the density
     *
     * @return     maximum value
     */
    float maxval() const;

private:

    /**
     * @brief      assert whether input is vasp 5 format
     */
    void test_vasp5();

    /**
     * @brief      read the scalar
     *
     *             Read the scalar value from the 2nd line of the CHGCAR file.
     *             Note that all read_* functions can be used seperately,
     *             although they may depend on each other and have to be used in
     *             some consecutive order as is done in the read() wrapper
     *             function.
     */
    void read_scalar();

    /**
     * @brief      read the unitcell matrix
     *
     *             Reads the matrix that defines the unit cell in the CHGCAR
     *             file. The inverse of that matrix is automatically
     *             constructed.
     *
     *             Note that all read_* functions can be used seperately,
     *             although they may depend on each other and have to be used in
     *             some consecutive order as is done in the read() wrapper
     *             function.
     */
    void read_matrix();

    /**
     * @brief      read the grid dimensions
     *
     *             Read the number of gridpoints in each direction.
     *
     *             Note that all read_* functions can be used seperately,
     *             although they may depend on each other and have to be used in
     *             some consecutive order as is done in the read() wrapper
     *             function.
     */
    void read_grid_dimensions();

    /**
     * @brief      read the atoms
     *
     *             Read the number of atoms of each element. These numbers are
     *             used to skip the required amount of lines.
     *
     *             Note that all read_* functions can be used seperately,
     *             although they may depend on each other and have to be used in
     *             some consecutive order as is done in the read() wrapper
     *             function.
     */
    void read_atoms();

    /**
     * @brief      read the grid
     *
     * @param[in]  lines  number of lines to read
     *
     * @return     percentage of file read
     */
    float read_grid(unsigned int lines);

public:

    /**
     * @brief      copy the grid dimensions
     *
     * @param      _grid_dimensions  array of grid dimensions
     */
    void copy_grid_dimensions(unsigned int _grid_dimensions[]) const;

    /**
     * @brief      Gets the matrix unitcell.
     *
     * @return     The matrix unitcell.
     */
    inline const glm::mat3& get_mat_unitcell() const {
        return this->mat33;
    }

    /**
     * @brief      Gets the grid pointer.
     *
     * @return     The grid pointer.
     */
    inline const float* get_grid_ptr() const {
        return &this->gridptr[0];
    }

    /**
     * @brief      Get a reference to the grid point vector.
     *
     * @return     The grid pointer.
     */
    inline const std::vector<float>& get_grid_vec() const {
        return this->gridptr;
    }

    /**
     * @brief      get the size of the grid array
     *
     * @return     size of the grid array
     */
    unsigned int get_size() const {
        return this->gridptr.size();
    }

    /**
     * @brief      Gets the filename.
     *
     * @return     The filename.
     */
    inline const std::string& get_filename() const {
        return this->filename;
    }

private:

    /**
     * @brief      Calculates the inverse.
     */
    void calculate_inverse();

    /**
     * @brief      calculate the volume
     *
     *             Calculates the inverse of a 3x3 matrix. This is a convenience
     *             function for the read_matrix() function.
     */
    void calculate_volume();
};

#endif //_DENSITY_H
