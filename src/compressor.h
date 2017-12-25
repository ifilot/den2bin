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

#ifndef _COMPRESSOR_H
#define _COMPRESSOR_H

#include <memory>
#include <fstream>

#include "dct.h"

class Compressor {
private:
    std::unique_ptr<DCT> dct;

public:
    Compressor();

    std::vector<float> compress_3d(const std::vector<float>& data_3d, size_t nx, size_t ny, size_t nz, size_t block_size, size_t coeff_size);
    std::vector<float> decompress_3d(const std::vector<float>& coeff, size_t nx, size_t ny, size_t nz, size_t block_size, size_t coeff_size);

private:
};

#endif // _COMPRESSOR_H
