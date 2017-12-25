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

#ifndef _DCT_H
#define _DCT_H

#include <iostream>
#include <vector>
#include <cmath>
#include <fftw3.h>
#include <complex>

class DCT {
private:
    float cc[4][4]; // preallocate cosine coefficients

public:

    DCT();

    std::vector<float> naive_dct_1d(const std::vector<float> &vec);
    std::vector<float> naive_dct_2d(const std::vector<float> &vec);
    std::vector<float> naive_dct_3d(const std::vector<float> &vec);

    std::vector<float> inv_dct_1d(const std::vector<float> &vec);
    std::vector<float> inv_dct_2d(const std::vector<float> &vec);
    std::vector<float> inv_dct_3d(const std::vector<float> &vec);

private:

};

#endif // _DCT_H
