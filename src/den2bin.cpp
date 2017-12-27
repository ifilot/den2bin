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

#include <chrono>
#include <tclap/CmdLine.h>
#include <boost/filesystem.hpp>

#include "density.h"
#include "converter.h"
#include "config.h"

int main(int argc, char* argv[]) {

    try {
        TCLAP::CmdLine cmd("Converts VASP density file (CHGCAR/LOCPOT/etc.) to a highly compressed binary", ' ', PROGRAM_VERSION);

        //**************************************
        // declare values to be parsed
        //**************************************

        // input file
        TCLAP::ValueArg<std::string> arg_input_filename("i","input","Input file (i.e. sphere.obj)",true,"__NONE__","filename");
        cmd.add(arg_input_filename);

        // output file
        TCLAP::ValueArg<std::string> arg_output_filename("o","output","Output file (i.e. sphere.mesh)",true,"__NONE__","filename");
        cmd.add(arg_output_filename);

        // message to store
        TCLAP::ValueArg<std::string> arg_message("m","message","Message",false,"","filename");
        cmd.add(arg_message);

        // whether to store raw
        TCLAP::SwitchArg arg_raw("r","raw","Raw",false);
        cmd.add(arg_raw);

        // whether to store in lossy format
        TCLAP::SwitchArg arg_lossy("l","lossy","Lossy",false);
        cmd.add(arg_lossy);

        // block size
        TCLAP::ValueArg<unsigned int> arg_block_size("b","blocksize","Block size",false,4,"blocksize");
        cmd.add(arg_block_size);

        // quality
        TCLAP::ValueArg<unsigned int> arg_quality("q","quality","Quality",false,2,"quality");
        cmd.add(arg_quality);

        // whether to check data loss
        TCLAP::SwitchArg arg_check("c","check","Check",false);
        cmd.add(arg_check);

        // whether to extract data in lossy format
        TCLAP::SwitchArg arg_extract("x","extract","Extract",false);
        cmd.add(arg_extract);

        cmd.parse(argc, argv);

        std::string input_filename = arg_input_filename.getValue();
        std::string output_filename = arg_output_filename.getValue();
        std::string message = arg_message.getValue();
        bool flag_write_raw = arg_raw.getValue();
        bool flag_write_lossy = arg_lossy.getValue();
        bool flag_extract = arg_extract.getValue();
        bool flag_check = arg_check.getValue();
        unsigned int blocksize = arg_block_size.getValue();
        unsigned int quality = arg_quality.getValue();

        std::cout << "-----------------------------------------" << std::endl;
        std::cout << "Executing DEN2BIN v." << PROGRAM_VERSION << std::endl;
        std::cout << "Author: Ivo Filot <i.a.w.filot@tue.nl>" << std::endl;
        std::cout << "-----------------------------------------" << std::endl;
        std::cout << std::endl;

        // build converter object
        Converter cv;

        if(flag_extract) {
            /*******************************************************************
             *
             * DECOMPRESSING / EXTRACTING
             *
             *******************************************************************/

            cv.get_info(input_filename);
            cv.build_density(output_filename);

        } else {

            /*******************************************************************
             *
             * COMPRESSING / PACKAGING
             *
             *******************************************************************/

            // read scalar field
            std::cout << "Building binary file" << std::endl;
            std::cout << "-----------------------------------------" << std::endl;
            auto start = std::chrono::system_clock::now();
            Density den(input_filename);
            while(den.read(1000) < 1.0) {}
            auto end = std::chrono::system_clock::now();
            std::chrono::duration<double> elapsed_seconds = end-start;
            std::cout << "Read " << input_filename << " in " << elapsed_seconds.count() << " seconds." << std::endl;

            if(flag_write_raw) {
                std::cout << "Creating raw binary file" << std::endl;
                cv.write_to_binary_raw(den, output_filename);
            } else if(flag_write_lossy) {
                std::cout << "Creating dct-compressed binary file" << std::endl;
                cv.write_to_binary_dct(message, den, output_filename, blocksize, quality, flag_check);
            }else {
                std::cout << "Creating lossless binary file" << std::endl;
                cv.write_to_binary(message, den, output_filename);
            }

            end = std::chrono::system_clock::now();
            elapsed_seconds = end-start;
            std::cout << "Constructed " << output_filename << " in " << elapsed_seconds.count() << " seconds." << std::endl;
            std::cout << "-----------------------------------------" << std::endl;
            std::cout << std::endl;

            if(!flag_write_raw) {
                // calculate compression statistics
                std::cout << "Compression statistics" << std::endl;
                std::cout << "-----------------------------------------" << std::endl;
                boost::filesystem::path ipath(input_filename);
                boost::filesystem::path opath(output_filename);
                size_t input_filesize = boost::filesystem::file_size(ipath);
                size_t output_filesize = boost::filesystem::file_size(opath);
                std::cout << "Filesize: " << (float)output_filesize / (1024 * 1024) << " mb." << std::endl;
                const float ratio = (float)input_filesize / (float)output_filesize;
                std::cout << "Compression ratio: " << ratio << std::endl;
                std::cout << "-----------------------------------------" << std::endl;
                std::cout << std::endl;

                std:: cout << "Verifying result..." << std::endl;
                std::cout << "-----------------------------------------" << std::endl;
                start = std::chrono::system_clock::now();
                cv.get_info(output_filename);
                end = std::chrono::system_clock::now();
                elapsed_seconds = end-start;
                std::cout << "Verification of results in " << elapsed_seconds.count() << " seconds." << std::endl;
            }
        }

        std::cout << "-----------------------------------------" << std::endl;
        std::cout << "Done" << std::endl << std::endl;

    } catch (TCLAP::ArgException &e) {
        std::cerr << "error: " << e.error() <<
                     " for arg " << e.argId() << std::endl;
        return -1;
    }
}
