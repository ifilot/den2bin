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

#include "compressor.h"

Compressor::Compressor() {
    this->dct = std::make_unique<DCT>();
}

/**
 * @brief      compress 2d data using dct
 *
 * @param[in]  data_3d     vector holding 3d data
 * @param[in]  nx          number of data points in x direction
 * @param[in]  ny          number of data points in y direction
 * @param[in]  nz          number of data points in y direction
 * @param[in]  block_size  block size (tile size)
 * @param[in]  coeff_size  number of dct coefficients
 *
 * @return     set of dct coefficients
 */
std::vector<float> Compressor::compress_3d(const std::vector<float>& data_3d, size_t nx, size_t ny, size_t nz, size_t block_size, size_t coeff_size) {
    // calculate number of tiles in each direction
    size_t nrtx = nx / block_size;
    size_t nrty = ny / block_size;
    size_t nrtz = nz / block_size;

    if(nx % block_size != 0) {
        nrtx++;
    }

    if(ny % block_size != 0) {
        nrty++;
    }

    if(nz % block_size != 0) {
        nrtz++;
    }

    // calculate compression
    size_t coeff_block_size = 0;
    for(size_t i = 0; i < block_size; i++) {
        for(size_t j = 0; j < block_size; j++) {
            for(size_t k = 0; k < block_size; k++) {
                if(i + j + k < coeff_size) {
                    coeff_block_size++;
                }
            }
        }
    }

    std::cout << "Building DCT; this might take a while. Target deflation factor: " << (float)coeff_block_size / (float)(block_size * block_size * block_size) << std::endl;

    auto npthreads = omp_get_max_threads();
    omp_set_num_threads(npthreads);
    std::cout << "Using " << npthreads << " threads for the DCT." << std::endl;

    // build coefficient vector
    std::vector<float> coeff(coeff_block_size * nrtx * nrty * nrtz, 0.0f);

    #pragma omp parallel for collapse(3) schedule(dynamic)
    for(size_t tz=0 ; tz < nrtz; tz++) {
        for(size_t ty=0 ; ty < nrty; ty++) {
            for(size_t tx=0 ; tx < nrtx; tx++) {

                std::vector<float> tile(block_size * block_size * block_size, 0.0f);
                std::vector<float> coeff_tile(block_size * block_size * block_size, 0.0f);

                // loop over values to store inside a tile
                for(size_t i = 0; i < block_size; i++) {
                    for(size_t j = 0; j < block_size; j++) {
                        for(size_t k = 0; k < block_size; k++) {

                            size_t x = tx * block_size + k;
                            size_t y = ty * block_size + j;
                            size_t z = tz * block_size + i;

                            // employ padding by using periodicity
                            x = x % nx;
                            y = y % ny;
                            z = z % nz;

                            // store data inside tile
                            tile[i * block_size * block_size + j * block_size + k] = data_3d[z * nx * ny + y * nx + x];
                        }
                    }
                }

                // calculate dct coefficients
                coeff_tile = this->dct->naive_dct_3d(tile);

                size_t idx = 0;
                size_t pos = (tz * nrtx * nrty + ty * nrtx + tx) * coeff_block_size;
                for(size_t i = 0; i < block_size; i++) {
                    for(size_t j = 0; j < block_size; j++) {
                        for(size_t k = 0; k < block_size; k++) {
                            if((i + j + k) < coeff_size) {
                                coeff[pos + idx] = coeff_tile[i * block_size * block_size + j * block_size + k];
                                idx++;
                            }
                        }
                    }
                }
            }
        }
    }

    return coeff;
}

/**
 * @brief      decompress 3d data using idct
 *
 * @param[in]  data_3d     vector holding 3d coefficients
 * @param[in]  nx          number of data points in x direction
 * @param[in]  ny          number of data points in y direction
 * @param[in]  nz          number of data points in y direction
 * @param[in]  block_size  block size (tile size)
 * @param[in]  coeff_size  number of dct coefficients
 *
 * @return     data
 */
std::vector<float> Compressor::decompress_3d(const std::vector<float>& coeff, size_t nx, size_t ny, size_t nz, size_t block_size, size_t coeff_size) {
    // calculate number of tiles in each direction
    size_t nrtx = nx / block_size;
    size_t nrty = ny / block_size;
    size_t nrtz = nz / block_size;

    if(nx % block_size != 0) {
        nrtx++;
    }

    if(ny % block_size != 0) {
        nrty++;
    }

    if(nz % block_size != 0) {
        nrtz++;
    }

    // calculate compression
    size_t coeff_block_size = 0;
    for(size_t i = 0; i < block_size; i++) {
        for(size_t j = 0; j < block_size; j++) {
            for(size_t k = 0; k < block_size; k++) {
                if(i + j + k < coeff_size) {
                    coeff_block_size++;
                }
            }
        }
    }

    std::vector<float> out(nx * ny * nz, 0.0f);

    #pragma omp parallel for collapse(3) schedule(dynamic)
    for(size_t tz=0 ; tz < nrtz; tz++) {
        for(size_t ty=0 ; ty < nrty; ty++) {
            for(size_t tx=0 ; tx < nrtx; tx++) {

                std::vector<float> tile(block_size * block_size * block_size, 0.0f);
                std::vector<float> coeff_tile(block_size * block_size * block_size, 0.0f);

                size_t idx = 0;
                size_t pos = (tz * nrtx * nrty + ty * nrtx + tx) * coeff_block_size;
                for(size_t i = 0; i < block_size; i++) {
                    for(size_t j = 0; j < block_size; j++) {
                        for(size_t k = 0; k < block_size; k++) {
                            if(i + j + k < coeff_size) {
                                coeff_tile[i * block_size * block_size + j * block_size + k] = coeff[pos + idx];
                                idx++;
                            }
                        }
                    }
                }

                // calculate inverse
                tile = this->dct->inv_dct_3d(coeff_tile);

                for(size_t i = 0; i < block_size; i++) {
                    for(size_t j = 0; j < block_size; j++) {
                        for(size_t k = 0; k < block_size; k++) {

                            size_t x = tx * block_size + k;
                            size_t y = ty * block_size + j;
                            size_t z = tz * block_size + i;

                            if(x >= nx || y >= ny || z >= nz) {
                                continue;
                            }

                            const size_t idx = z * nx * ny + y * nx + x;
                            out[idx] = tile[i * block_size * block_size + j * block_size + k];
                        }
                    }
                }
            }
        }
    }

    return out;
}
