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
        msg[msg_size] = '\0';       // write end of string character
        std::string msgstr(msg);    // parse char array to std::string
        std::cout << "File comments: " << msgstr << std::endl;
        delete[] msg;                 // we no longer need the memory

        // read header size
        decompressed.read(buffer, sizeof(uint32_t));
        const uint32_t header_size = *(uint32_t *)buffer;

        // read and output header message
        char* header = new char[header_size+10];
        decompressed.read(header, header_size);
        header[header_size] = '\0';       // write end of string character
        std::string headerstr(header);    // parse char array to std::string
        std::cout << "File header: " << headerstr << std::endl;
        delete[] header;               // we no longer need the memory

        // read the dimensions of the unit cell
        for(unsigned int i=0; i<3; i++) {
            for(unsigned int j=0; j<3; j++) {
                decompressed.read(buffer, sizeof(float));
                this->mat[i][j] = *(float *)buffer;
            }
        }

        // read the grid sizes
        for(unsigned int i=0; i<3; i++) {
            decompressed.read(buffer, sizeof(uint32_t));
            this->griddim[i] = *(uint32_t *)buffer;
        }
        std::cout << "Grid dimension: " << griddim[0] << "x" << griddim[1] << "x" << griddim[2] << std::endl;

        if(headerstr.substr(0,3).compare("DCT") == 0) {
            std::cout << "This is a DCT binary." << std::endl;

            this->blocksize = boost::lexical_cast<size_t>(headerstr.substr(4,4));
            this->coeffsize = boost::lexical_cast<size_t>(headerstr.substr(9,4));

            std::cout << "Blocksize: " << this->blocksize << std::endl;
            std::cout << "Coeffsize: " << this->coeffsize << std::endl;

            // obtain coeff size
            decompressed.read(buffer, sizeof(uint32_t));
            const uint32_t coeff_size = *(uint32_t *)buffer;

            // obtain min value
            decompressed.read(buffer, sizeof(float));
            this->minval = *(float *)buffer;

            // obtain max value
            decompressed.read(buffer, sizeof(float));
            this->maxval = *(float *)buffer;

            std::cout << "Data contains " << coeff_size << " DCT coefficients." << std::endl;
            std::cout << "Minimum value: " << this->minval << std::endl;
            std::cout << "Maximum value: " << this->maxval << std::endl;

            this->coeff.clear();
            this->data.clear();
            for(unsigned int i=0; i<coeff_size; i++) {
                decompressed.read(buffer, sizeof(float));
                const float val = *(float *)buffer;
                this->coeff.push_back(val);
            }

            // de-normalize data
            this->denormalize_data(&this->data, this->minval, this->maxval);
        }

        if(headerstr.substr(0,3).compare("BIN") == 0) {
            unsigned int gridsz = griddim[0] * griddim[1] * griddim[2];
            float min = 1e50;
            float max = -1e50;
            std::cout << "Checking " << gridsz << " values." << std::endl;
            this->coeff.clear();
            this->data.clear();
            for(unsigned int i=0; i<gridsz; i++) {
                decompressed.read(buffer, sizeof(float));
                const float val = *(float *)buffer;
                this->data.push_back(val);
                min = std::min(min, val);
                max = std::max(max, val);
            }
            std::cout << "Minimum value: " << min << std::endl;
            std::cout << "Maximum value: " << max << std::endl;
        }

        delete[] buffer;  // clean up the mess :-)

    } else {
        std::cerr << "Cannot read file " << filename << std::endl;
        throw std::runtime_error("Cannot open file");
    }
}

/**
 * @brief      rebuild CHGCAR / LOCPOT from data
 *
 * @param[in]  filename  path to CHGCAR / LOCPOT file
 */
void Converter::build_density(const std::string& filename) {
    // open file
    std::ofstream f(filename);
    std::stringstream out;

    if(this->coeff.size() > 0) {
        std::cout << "Decompressing DCT; this might take a while" << std::endl;
        Compressor compressor;
        auto start = std::chrono::system_clock::now();
        this->data = compressor.decompress_3d(this->coeff, this->griddim[0], this->griddim[1], this->griddim[2], this->blocksize, this->coeffsize);
        auto end = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = end-start;
        std::cout << "Decompressed DCT in " << elapsed_seconds.count() << " seconds." << std::endl;

        // de-normalize data
        this->denormalize_data(&this->data, this->minval, this->maxval);
    }

    if(this->data.size() > 0) {
        std::cout << "Writing data to " << filename << std::endl;
        auto start = std::chrono::system_clock::now();

        // write message
        out << "CHGCAR" << "\n";

        // write scalar
        out << boost::format("    %10.6f") % 1.0 << "\n";

        // write matrix
        out << boost::format("  %10.6f  %12.6f  %12.6f") % mat[0][0] % mat[0][1] % mat[0][2] << "\n";
        out << boost::format("  %10.6f  %12.6f  %12.6f") % mat[1][0] % mat[1][1] % mat[1][2] << "\n";
        out << boost::format("  %10.6f  %12.6f  %12.6f") % mat[2][0] % mat[2][1] % mat[2][2] << "\n";

        // write atoms
        out << boost::format("    %i") % 0 << "\n";
        out << "Direct" << "\n" << "\n";

        out << boost::format("  %i  %i  %i") % griddim[0] % griddim[1] % griddim[2] << "\n";

        size_t gridsz = griddim[0] * griddim[1] * griddim[2];
        float volume = glm::dot(glm::cross(this->mat[0], this->mat[1]), this->mat[2]);

        size_t nrthreads = omp_get_max_threads();
        std::stringstream local[nrthreads];

        #pragma omp parallel
        {
            size_t threadnum = omp_get_thread_num();

            size_t start = gridsz / nrthreads * threadnum;
            size_t stop = gridsz / nrthreads * (threadnum + 1);

            stop = std::min(gridsz, stop);

            char buffer[50];

            #pragma omp parallel for
            for(size_t i=start; i<stop; i++) {

                sprintf(buffer, "% 11.10E", this->data[i] * volume);
                local[threadnum] << buffer;

                if((i+1) % 5 == 0) {
                    local[threadnum] << "\n";
                } else {
                    local[threadnum] << " ";
                }
            }
        }

        // merge results
        for(unsigned int i=0; i<nrthreads; i++) {
            out << local[i].str();
        }

        auto end = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = end-start;
        std::cout << "Wrote file in " << elapsed_seconds.count() << " seconds." << "\n";
    }

    // write to file
    f << out.str();

    // close file
    f.close();
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
        const uint32_t comments_length = comments.size();
        origin.write((char*)&comments_length, sizeof(uint32_t));
        origin.write((char*)&comments.c_str()[0], sizeof(char) * comments_length);

        // construct header
        const std::string header = "BIN|BZIP2|" + PROGRAM_VERSION;

        // store header with file information
        const uint32_t header_length = header.size();
        origin.write((char*)&header_length, sizeof(uint32_t));
        origin.write((char*)&header.c_str()[0], sizeof(char) * header_length);

        // write the matrix; 9 floats
        for(unsigned int i=0; i<3; i++) {
            for(unsigned int j=0; j<3; j++) {
                const float val = density.get_mat_unitcell()[i][j];
                origin.write((char*)&val, sizeof(float));
            }
        }

        // write grid dimensions
        density.copy_grid_dimensions(this->griddim);
        for(unsigned int i=0; i<3; i++) {
            const uint32_t val = this->griddim[i];
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

/**
 * @brief      write density class to binary file in a lossy format (using DCT)
 *
 * @param[in]  comments    comments to store
 * @param[in]  density     density object
 * @param[in]  outputfile  output path
 *
 * @return     filesize
 */
void Converter::write_to_binary_dct(const std::string& comments, const Density& density, const std::string& outputfile, size_t blocksize, size_t coeffsize) {
    std::fstream f(outputfile, std::ios_base::binary | std::ios::out);

    if(f.good()) {

        std::stringstream compressed;
        std::stringstream origin;

        // write a leading zero for comments
        const uint32_t comments_length = comments.size();
        origin.write((char*)&comments_length, sizeof(uint32_t));
        origin.write((char*)&comments.c_str()[0], sizeof(char) * comments_length);

        // construct header
        const std::string header = (boost::format("DCT %04i %04i|BZIP2|%s") % blocksize % coeffsize % PROGRAM_VERSION).str();

        // store header with file information
        const uint32_t header_length = header.size();
        origin.write((char*)&header_length, sizeof(uint32_t));
        origin.write((char*)&header.c_str()[0], sizeof(char) * header_length);

        // write the matrix; 9 floats
        for(unsigned int i=0; i<3; i++) {
            for(unsigned int j=0; j<3; j++) {
                const float val = density.get_mat_unitcell()[i][j];
                origin.write((char*)&val, sizeof(float));
            }
        }

        // write grid dimensions
        density.copy_grid_dimensions(this->griddim);
        for(unsigned int i=0; i<3; i++) {
            const uint32_t val = this->griddim[i];
            origin.write((char*)&val, sizeof(uint32_t));
        }

        // normalize data
        this->minval = density.minval();
        this->maxval = density.maxval();
        auto data = density.get_grid_vec();
        this->normalize_data(&data, this->minval, this->maxval);

        // write the DCT coefficients
        Compressor compressor;
        auto coeff = compressor.compress_3d(data, this->griddim[0], this->griddim[1], this->griddim[2], blocksize, coeffsize);

        // write coefficients size
        const uint32_t coeff_length = coeff.size();
        origin.write((char*)&coeff_length, sizeof(uint32_t));

        // write min value
        origin.write((char*)&this->minval, sizeof(float));

        // write max value
        origin.write((char*)&this->maxval, sizeof(float));

        // write the individual coefficients
        for(unsigned int i=0; i<coeff.size(); i++) {
            origin.write((char*)&coeff[i], sizeof(float));
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

/**
 * @brief      write density class to binary file in raw format (mainly used for debugging and other purposes)
 *
 * @param[in]  comments    comments to store
 * @param[in]  density     density object
 * @param[in]  outputfile  output path
 *
 * @return     filesize
 */
void Converter::write_to_binary_raw(const Density& density, const std::string& outputfile) {
    std::fstream f(outputfile, std::ios_base::binary | std::ios::out);

    if(f.good()) {

        std::stringstream compressed;
        std::stringstream origin;

        // write the grid
        uint32_t gridsize = density.get_size();
        const float* gridptr = density.get_grid_ptr();
        for(unsigned int i=0; i<gridsize; i++) {
            origin.write((char*)&gridptr[i], sizeof(float));
        }

        f << origin.str();
        f.close();

    } else {
        std::cerr << "Cannot write to file " << outputfile << std::endl;
        throw std::runtime_error("Could not write to file");
    }
}

/**
 * @brief      normalize data to [-0.5,0.5]
 *
 * @param      data  The data
 * @param[in]  min   minimum value
 * @param[in]  max   maximum value
 */
void Converter::normalize_data(std::vector<float>* data, float min, float max) {
    auto start = std::chrono::system_clock::now();
    const float range = max - min;

    #pragma omp parallel for
    for(unsigned int i=0; i<data->size(); i++) {
        data->at(i) = (data->at(i) - min) / range - 0.5f;
    }

    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    std::cout << "Normalized data in " << elapsed_seconds.count() << " seconds." << std::endl;
}

/**
 * @brief      denormalize data to original range
 *
 * @param      data  The data
 * @param[in]  min   minimum value
 * @param[in]  max   maximum value
 */
void Converter::denormalize_data(std::vector<float>* data, float min, float max) {
    auto start = std::chrono::system_clock::now();
    const float range = max - min;

    #pragma omp parallel for
    for(unsigned int i=0; i<data->size(); i++) {
        data->at(i) = (data->at(i) + 0.5f) * range + min;
    }

    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    std::cout << "Normalized data in " << elapsed_seconds.count() << " seconds." << std::endl;
}
