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

#include "converter.h"

/**
 * @brief      default constructor for Converter object
 */
Converter::Converter() {}

/**
 * @brief      extract information of binary file and output to terminal
 *
 * @param[in]  filename  path to binary file
 */
void Converter::get_info(const std::string& filename) {
	std::ifstream f(filename);
    if(f.is_open()) {

        std::stringstream compressed;
        compressed << f.rdbuf();
        f.close();

        std::stringstream decompressed;
        boost::iostreams::filtering_streambuf<boost::iostreams::input> out;
        out.push(boost::iostreams::bzip2_decompressor());
        out.push(compressed);
        boost::iostreams::copy(out, decompressed);

        // create read buffer
        char* buffer = new char[sizeof(float) * 3];

        // read message size
        decompressed.read(buffer, sizeof(uint32_t));
        const uint32_t msg_size = *(uint32_t *)buffer;

        // read and output header message
        char* msg = new char[msg_size+1];
        decompressed.read(msg, msg_size);
        msg[msg_size] = '\0';		// write end of string character
        std::string msgstr(msg);	// parse char array to std::string
        std::cout << "Header message: " << msgstr << std::endl;
        delete msg;					// we no longer need the memory

        // read the dimensions of the unit cell
        glm::mat3 mat;
        for(unsigned int i=0; i<3; i++) {
        	for(unsigned int j=0; j<3; j++) {
        		decompressed.read(buffer, sizeof(float));
        		mat[i][j] = *(float *)buffer;
        	}
        }

        // read the grid sizes
        unsigned int gridsize[3];
        for(unsigned int i=0; i<3; i++) {
    		decompressed.read(buffer, sizeof(uint32_t));
    		gridsize[i] = *(uint32_t *)buffer;
        }
        std::cout << "Grid dimension: " << gridsize[0] << "x" << gridsize[1] << "x" << gridsize[2] << std::endl;

        unsigned int gridsz = gridsize[0] * gridsize[1] * gridsize[2];
        float min = 1e50;
        float max = -1e50;
        std::cout << "Checking " << gridsz << " values." << std::endl;
        for(unsigned int i=0; i<gridsz; i++) {
            decompressed.read(buffer, sizeof(float));
            const float val = *(float *)buffer;
            min = std::min(min, val);
            max = std::max(max, val);
        }
        std::cout << "Minimum value: " << min << std::endl;
        std::cout << "Maximum value: " << max << std::endl;

        delete buffer;  // clean up the mess :-)

    } else {
        std::cerr << "Cannot read file " << filename << std::endl;
        throw std::runtime_error("Cannot open file");
    }
}

/**
 * @brief      write density class to binary file
 *
 * @param[in]  comments    comments to store
 * @param[in]  density     density object
 * @param[in]  outputfile  output path
 *
 * @return     filesize
 */
void Converter::write_to_binary(const std::string& comments, const Density& density, const std::string& outputfile) {
	std::fstream f(outputfile, std::ios_base::binary | std::ios::out);

    if(f.good()) {

        std::stringstream compressed;
        std::stringstream origin;

        // write a leading zero for comments
        const uint32_t header_length = comments.size();
        origin.write((char*)&header_length, sizeof(uint32_t));
        origin.write((char*)&comments.c_str()[0], sizeof(char) * header_length);

        // write the matrix; 9 floats
        for(unsigned int i=0; i<3; i++) {
        	for(unsigned int j=0; j<3; j++) {
        		const float val = density.get_mat_unitcell()[i][j];
        		origin.write((char*)&val, sizeof(float));
        	}
        }

        // write grid dimensions
        unsigned int griddim[3];
        density.copy_grid_dimensions(griddim);
        for(unsigned int i=0; i<3; i++) {
        	const uint32_t val = griddim[i];
        	origin.write((char*)&val, sizeof(uint32_t));
        }

        // write the grid
        uint32_t gridsize = density.get_size();
        const float* gridptr = density.get_grid_ptr();
        for(unsigned int i=0; i<gridsize; i++) {
        	origin.write((char*)&gridptr[i], sizeof(float));
        }

        boost::iostreams::filtering_streambuf<boost::iostreams::input> out;
        out.push(boost::iostreams::bzip2_compressor());
        out.push(origin);
        boost::iostreams::copy(out, compressed);

        f << compressed.str();
        f.close();

    } else {
        std::cerr << "Cannot write to file " << outputfile << std::endl;
        throw std::runtime_error("Could not write to file");
    }
}