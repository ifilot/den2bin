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

#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/bzip2.hpp>

#include "density.h"

/**
 * @brief      converter class that outputs density to binary file
 */
class Converter {
private:

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
	 * @return     filesize
	 */
	void write_to_binary(const std::string& comments, const Density& sf, const std::string& outputfile);
private:
};

#endif // _CONVERTER_H