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

#include "density.h"

/**
 * @brief      default constructor
 *
 * @param[in]  _filename  input filename
 */
Density::Density(const std::string &_filename) {
    this->filename = _filename;
    this->scalar = -1;
    this->vasp5_input = false;
    this->has_read = false;
    this->header_read = false;
}

/**
 * @brief      read number of lines
 *
 * @param[in]  lines  number of lines to read
 *
 * @return     percentage of file read
 */
float Density::read(unsigned int lines) {
    if(has_read) {
        return 1.0;
    }

    if(!this->header_read) {
        this->test_vasp5();
        this->read_scalar();
        this->read_matrix();
        this->read_atoms();
        this->read_grid_dimensions();
    }

    return this->read_grid(lines);
}

/**
 * @brief      get minimum value in the density
 *
 * @return     minimum value
 */
float Density::minval() const {
    return *std::min_element(this->gridptr.begin(), this->gridptr.end());
}

/**
 * @brief      get maximum value in the density
 *
 * @return     maximum value
 */
float Density::maxval() const {
    return *std::max_element(this->gridptr.begin(), this->gridptr.end());
}

/**
 * @brief      assert whether input is vasp 5 format
 */
void Density::test_vasp5() {
    std::ifstream infile(this->filename.c_str());
    std::string line;
    for(unsigned int i=0; i<5; i++) { // discard first two lines
        std::getline(infile, line);
    }
    std::getline(infile, line);

    // check if this line contains atomic information (i.e. alpha-characters)
    boost::regex regex_vasp_version("^(.*[A-Za-z]+.*)$");
    boost::smatch what;
    if(boost::regex_match(line, what, regex_vasp_version)) {
        this->vasp5_input = true;
    }
}

/**
  * @brief      read the scalar
  *
  *             Read the scalar value from the 2nd line of the density file. Note
  *             that all read_* functions can be used seperately, although they
  *             may depend on each other and have to be used in some consecutive
  *             order as is done in the read() wrapper function.
  */
void Density::read_scalar() {
    std::ifstream infile(this->filename.c_str());
    std::string line;
    std::getline(infile, line); // discard this line

    std::getline(infile, line);
    boost::regex regex_scalar("^\\s*([0-9.-]+)\\s*$");
    boost::smatch what;
    if(boost::regex_match(line, what, regex_scalar)) {
        this->scalar = boost::lexical_cast<float>(what[1]);
    } else {
        this->scalar = -1;
    }
}

/**
  * @brief      read the unitcell matrix
  *
  *             Reads the matrix that defines the unit cell in the density file.
  *             The inverse of that matrix is automatically constructed.
  *
  *             Note that all read_* functions can be used seperately, although
  *             they may depend on each other and have to be used in some
  *             consecutive order as is done in the read() wrapper function.
  */
void Density::read_matrix() {
    std::ifstream infile(this->filename.c_str());
    std::string line;
    for(unsigned int i=0; i<2; i++) { // discard first two lines
        std::getline(infile, line);
    }

    // setup match pattern
    boost::regex regex_vasp_matrix_line("^\\s*([0-9.-]+)\\s+([0-9.-]+)\\s+([0-9.-]+)\\s*$");
    for(unsigned int i=0; i<3; i++) {
        std::getline(infile, line);
        boost::smatch what;
        if(boost::regex_match(line, what, regex_vasp_matrix_line)) {
            for(unsigned int j=0; j<3; j++) {
                mat[i][j] = boost::lexical_cast<float>(what[j+1]);
            }
        }
    }

    for(unsigned int i=0; i<3; i++) {
        for(unsigned int j=0; j<3; j++) {
            this->mat[i][j] *= this->scalar;
        }
    }

    // also construct inverse matrix
    this->calculate_inverse();

    // calculate matrix volume
    this->calculate_volume();
}

/**
  * @brief      read the atoms
  *
  *             Read the number of atoms of each element. These numbers are used
  *             to skip the required amount of lines.
  *
  *             Note that all read_* functions can be used seperately, although
  *             they may depend on each other and have to be used in some
  *             consecutive order as is done in the read() wrapper function.
  */
void Density::read_atoms() {
    std::ifstream infile(this->filename.c_str());
    std::string line;
    for(unsigned int i=0; i < (this->vasp5_input ? 7 : 6); i++) { // discard first two lines
        std::getline(infile, line);
    }

    std::vector<std::string> pieces;
    boost::trim(line);
    boost::split(pieces, line, boost::is_any_of("\t "), boost::token_compress_on);
    for(unsigned int i=0; i<pieces.size(); i++) {
        boost::trim(pieces[i]);
        this->nrat.push_back(boost::lexical_cast<unsigned int>(pieces[i]));
    }
}

/**
  * @brief      read the grid dimensions
  *
  *             Read the number of gridpoints in each direction.
  *
  *             Note that all read_* functions can be used seperately, although
  *             they may depend on each other and have to be used in some
  *             consecutive order as is done in the read() wrapper function.
  */
void Density::read_grid_dimensions() {
    std::ifstream infile(this->filename.c_str());
    std::string line;
    // skip lines that contain atoms
    for(unsigned int i=0; i<(this->vasp5_input ? 10 : 9); i++) {
        std::getline(infile, line);
    }
    for(unsigned int i=0; i<this->nrat.size(); i++) {
        for(unsigned int j=0; j<this->nrat[i]; j++) {
                std::getline(infile, line);
        }
    }

    boost::trim(line);
    this->gridline = line;

    std::vector<std::string> pieces;
    boost::split(pieces, line, boost::is_any_of("\t "), boost::token_compress_on);
    for(unsigned int i=0; i<pieces.size(); i++) {
        this->grid_dimensions[i] = boost::lexical_cast<unsigned int>(pieces[i]);
    }
}

/*
 * void read_grid()
 *
 * Read all the grid points. This function depends
 * on the the gridsize being set via the
 * read_grid_dimensions() function.
 *
 * Note that all read_* functions can
 * be used seperately, although they may depend
 * on each other and have to be used in some
 * consecutive order as is done in the read()
 * wrapper function.
 *
 */
float Density::read_grid(unsigned int lines) {
    if(!this->header_read) {
        this->header_read = true;
        this->infile.open(this->filename.c_str());
        std::string line;
        // skip lines that contain atoms
        unsigned int skiplines = 0;
        for(unsigned int i=0; i<(this->vasp5_input ? 10 : 9); i++) {
            std::getline(this->infile, line);
            skiplines++;
        }
        for(unsigned int i=0; i<this->nrat.size(); i++) {
            for(unsigned int j=0; j<this->nrat[i]; j++) {
                    std::getline(this->infile, line);
                    skiplines++;
            }
        }

        this->gridsize = this->grid_dimensions[0] * this->grid_dimensions[1] * this->grid_dimensions[2];
    }

    float_parser p;

    /* read spin up */
    unsigned int linecounter=0; // for the counter
    static const boost::regex regex_augmentation("augmentation.*");
    std::string line;

    for(unsigned int i=0; i<lines; i++) {
        std::getline(this->infile, line);
        // stop looping when a second gridline appears (this
        // is where the spin down part starts)
        if(line.compare(this->gridline) == 0) {
            std::cout << "I am breaking the loop" << std::endl;
            break;
        }

        boost::smatch what;
        if(boost::regex_match(line, what, regex_augmentation)) {
            std::cout << "Augmentation break encountered" << std::endl;
            break;
        }

        // set iterators
        std::string::const_iterator b = line.begin();
        std::string::const_iterator e = line.end();

        // parse
        std::vector<float> floats;
        boost::spirit::qi::phrase_parse(b, e, p, boost::spirit::ascii::space, floats);

        // expand gridptr with the new size
        unsigned int cursize = this->gridptr.size();
        this->gridptr.resize(cursize + floats.size());

        // set the number of threads equal to the size of the pieces on the line
        for(unsigned int j=0; j<floats.size(); j++) {
            this->gridptr[cursize + j] = floats[j] / this->volume;
        }

        linecounter++;

        if(this->gridptr.size() >= this->gridsize) {
            this->has_read = true;
            return 1.0f;
        }
    }

    return (float)this->gridptr.size() / (float)this->gridsize;
}

/*
 * void calculate_inverse()
 *
 * Calculates the inverse of a 3x3 matrix. This is a convenience
 * function for the read_matrix() function.
 *
 */
void Density::calculate_inverse() {
    float det = 0;
    for(unsigned int i=0;i<3;i++) {
        det += (this->mat[0][i]*(this->mat[1][(i+1)%3]*this->mat[2][(i+2)%3] - this->mat[1][(i+2)%3]*this->mat[2][(i+1)%3]));
    }

    for(unsigned int i=0;i<3;i++){
            for(unsigned int j=0;j<3;j++) {
                     this->imat[i][j] = ((this->mat[(i+1)%3][(j+1)%3] * this->mat[(i+2)%3][(j+2)%3]) - (this->mat[(i+1)%3][(j+2)%3]*this->mat[(i+2)%3][(j+1)%3]))/ det;
            }
     }
}

/**
 * @brief      calculate the volume
 *
 *             Calculates the inverse of a 3x3 matrix. This is a convenience
 *             function for the read_matrix() function.
 */
void Density::calculate_volume() {
    for(unsigned int i=0; i<3; i++) {
        for(unsigned int j=0; j<3; j++) {
            this->mat33[i][j] = this->mat[i][j];
        }
    }

    this->volume = glm::dot(glm::cross(this->mat33[0], this->mat33[1]), this->mat33[2]);
}

/**
  * @brief      copy the grid dimensions
  *
  * @param      _grid_dimensions  array of grid dimensions
  */
void Density::copy_grid_dimensions(unsigned int _grid_dimensions[]) const {
    for(unsigned int i=0; i<3; i++) {
        _grid_dimensions[i] = this->grid_dimensions[i];
    }
}
