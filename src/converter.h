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

#ifndef _CONVERTER_H
#define _CONVERTER_H

#include <string>
#include <stdexcept>
#include <fstream>
#include <iostream>
#include <chrono>
#include <unistd.h>
#include <stdio.h>

#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/format.hpp>

#include "density.h"
#include "compressor.h"
#include "config.h"

/**
 * @brief      converter class that outputs density to binary file
 */
class Converter {
private:
    std::vector<float> coeff;
    std::vector<float> data;

    size_t blocksize;
    size_t coeffsize;

    unsigned int griddim[3];
    glm::mat3 mat;

    float minval;       //!< minimum value in dataset
    float maxval;       //!< maximum value in dataset

public:
    /**
     * @brief      default constructor for Converter object
     */
    Converter();

    /**
     * @brief      extract information of binary file and output to terminal
     *
     * @param[in]  filename  path to binary file
     */
    void get_info(const std::string& filename);

    /**
     * @brief      write density class to binary file
     *
     * @param[in]  comments    comments to store
     * @param[in]  density     density object
     * @param[in]  outputfile  output path
     *
     */
    void write_to_binary(const std::string& comments, const Density& sf, const std::string& outputfile);

    /**
     * @brief      write density class to binary file in a lossy format (using
     *             DCT)
     *
     * @param[in]  comments    comments to store
     * @param[in]  density     density object
     * @param[in]  outputfile  output path
     * @param[in]  blocksize   The blocksize
     * @param[in]  coeffsize   The coeffsize
     */
    void write_to_binary_dct(const std::string& comments, const Density& density, const std::string& outputfile, size_t blocksize, size_t coeffsize, bool check = false);

    /**
     * @brief      write density class to binary file (raw)
     *
     * @param[in]  density     density object
     * @param[in]  outputfile  output path
     *
     */
    void write_to_binary_raw(const Density& density, const std::string& outputfile);

    /**
     * @brief      rebuild CHGCAR / LOCPOT from data
     *
     * @param[in]  filename  path to CHGCAR / LOCPOT file
     */
    void build_density(const std::string& filename);
private:

    /**
     * @brief      normalize data to [-0.5,0.5]
     *
     * @param      data  The data
     * @param[in]  min   minimum value
     * @param[in]  max   maximum value
     */
    void normalize_data(std::vector<float>* data, float min, float max);

    /**
     * @brief      denormalize data to original range
     *
     * @param      data  The data
     * @param[in]  min   minimum value
     * @param[in]  max   maximum value
     */
    void denormalize_data(std::vector<float>* data, float min, float max);

    /**
     * @brief      calculate difference between the dct coefficient and the data
     *
     * @param[in]  orig       original data
     * @param[in]  coeff      dct coefficients
     * @param[in]  blocksize  The blocksize
     * @param[in]  coeffsize  The coeffsize
     *
     * @return     the difference
     */
    float calculate_difference(const std::vector<float>& orig, const std::vector<float>& coeff, size_t blocksize, size_t coeffsize);
};

#endif // _CONVERTER_H
