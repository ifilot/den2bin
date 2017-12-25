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

#include "dct.h"

DCT::DCT() {
    size_t len = 4;
    const float factor = M_PI / (float)len;
    for (size_t i = 0; i < len; i++) {
        for (size_t j = 0; j < len; j++) {
            cc[i][j] = std::cos(((float)i + 0.5) * (float)j * factor);
        }
    }
}

/**
 * @brief      naive 1D DCT
 *
 * @param[in]  set of input values
 *
 * @return     DCT coefficients
 */
std::vector<float> DCT::naive_dct_1d(const std::vector<float> &vec) {
    // construct output vector
    std::vector<float> result;
    size_t len = vec.size();    // obtain length from vector
    result.reserve(len);        // pre-allocate memory

    float factor = M_PI / (float)len;

    // loop over dct coefficients
    for (size_t i = 0; i < len; i++) {
        float sum = 0.0;

        // loop over x positions
        for (size_t j = 0; j < len; j++) {
            sum += vec.at(j) * std::cos(((float)j + 0.5) * (float)i * factor);
        }

        if(i == 0) {
            sum *= 1.0 / (float)len;
        } else {
            sum *= 2.0 / (float)len;
        }

        result.push_back(sum);
    }

    return result;
}

/**
 * @brief      calculate 1d dct inverse
 *
 * @param[in]  vec   set of dct coefficients
 *
 * @return     output
 */
std::vector<float> DCT::inv_dct_1d(const std::vector<float> &vec) {
    // construct output vector
    std::vector<float> result;
    size_t len = vec.size();    // obtain length from vector
    result.reserve(len);        // pre-allocate memory

    float factor = M_PI / (float)len;

    // loop over x positions
    for (size_t i = 0; i < len; i++) {
        float sum = 0.0;

        // loop over dct coefficients
        for (size_t j = 0; j < len; j++) {
            sum += vec.at(j) * std::cos(((float)i + 0.5) * (float)j * factor);
        }

        result.push_back(sum);
    }

    return result;
}

/**
 * @brief      naive 2D DCT
 *
 * @param[in]  set of input values
 *
 * @return     DCT coefficients
 */
std::vector<float> DCT::naive_dct_2d(const std::vector<float> &vec) {
    // construct output vector
    std::vector<float> result;
    size_t len =  std::sqrt(vec.size());   // obtain length from vector
    result.reserve(vec.size());            // pre-allocate memory

    const float factor = M_PI / (float)len;

    // loop over dct periods
    for (size_t i = 0; i < len; i++) {
        for (size_t j = 0; j < len; j++) {

            float sum = 0.0;
            float ci, cj;

            // loop over xy coordinates
            for (size_t k = 0; k < len; k++) {
                for (size_t l = 0; l < len; l++) {

                    sum += vec.at(k * len + l) *
                           std::cos(((float)k + 0.5) * (float)i * factor) *
                           std::cos(((float)l + 0.5) * (float)j * factor);
                }
            }

            if(i == 0) {
                ci = 1.0 / (float)len;
            } else {
                ci = 2.0 / (float)len;
            }

            if(j == 0) {
                cj = 1.0 / (float)len;
            } else {
                cj = 2.0 / (float)len;
            }

            result.push_back(ci * cj * sum);
        }
    }

    return result;
}

/**
 * @brief      calculate 2d dct inverse
 *
 * @param[in]  vec   set of dct coefficients
 *
 * @return     output
 */
std::vector<float> DCT::inv_dct_2d(const std::vector<float> &vec) {
    // construct output vector
    std::vector<float> result;
    size_t len =  std::sqrt(vec.size());   // obtain length from vector
    result.reserve(vec.size());            // pre-allocate memory

    const float factor = M_PI / (float)len;

    // loop over xy coordinates
    for (size_t i = 0; i < len; i++) {
        for (size_t j = 0; j < len; j++) {

            float sum = 0.0;

            // loop over dct coefficients
            for (size_t k = 0; k < len; k++) {
                for (size_t l = 0; l < len; l++) {

                    if(vec.at(k * (float)len + l) == 0.0f) {
                        continue;
                    }

                    sum += vec.at(k * (float)len + l) *
                           std::cos(((float)i + 0.5) * (float)k * factor) *
                           std::cos(((float)j + 0.5) * (float)l * factor);
                }
            }

            result.push_back(sum);
        }
    }

    return result;
}

/**
 * @brief      naive 3D DCT
 *
 * @param[in]  set of input values
 *
 * @return     DCT coefficients
 */
std::vector<float> DCT::naive_dct_3d(const std::vector<float> &vec) {
    // construct output vector
    std::vector<float> result;
    size_t len =  std::cbrt(vec.size());   // obtain length from vector
    result.reserve(vec.size());            // pre-allocate memory

    const float factor = M_PI / (float)len;

    // loop over dct periods
    for (size_t i = 0; i < len; i++) {
        for (size_t j = 0; j < len; j++) {
            for (size_t k = 0; k < len; k++) {

                float sum = 0.0;
                float ci, cj, ck;

                // loop over xyz coordinates
                for (size_t m = 0; m < len; m++) {  // z
                    for (size_t n = 0; n < len; n++) {  // y
                        for (size_t o = 0; o < len; o++) {  // x

                            // sum += vec.at(m * (len * len) + n * len + o) *
                            //        std::cos(((float)m + 0.5) * (float)i * factor) *
                            //        std::cos(((float)n + 0.5) * (float)j * factor) *
                            //        std::cos(((float)o + 0.5) * (float)k * factor);

                            // use pre-generated cosines
                            sum += vec.at(m * (len * len) + n * len + o) * cc[m][i] * cc[n][j] * cc[o][k];
                        }
                    }
                }

                if(i == 0) {
                    ci = 1.0 / (float)len;
                } else {
                    ci = 2.0 / (float)len;
                }

                if(j == 0) {
                    cj = 1.0 / (float)len;
                } else {
                    cj = 2.0 / (float)len;
                }

                if(k == 0) {
                    ck = 1.0 / (float)len;
                } else {
                    ck = 2.0 / (float)len;
                }

                result.push_back(ci * cj * ck * sum);
            }
        }
    }

    return result;
}

/**
 * @brief      calculate 3d dct inverse
 *
 * @param[in]  vec   set of dct coefficients
 *
 * @return     output
 */
std::vector<float> DCT::inv_dct_3d(const std::vector<float> &vec) {
    // construct output vector
    std::vector<float> result;
    size_t len =  std::cbrt(vec.size());   // obtain length from vector
    result.reserve(vec.size());            // pre-allocate memory

    const float factor = M_PI / (float)len;

    // loop over xyz coordinates
    for (size_t i = 0; i < len; i++) {  // z
        for (size_t j = 0; j < len; j++) {  // y
            for (size_t k = 0; k < len; k++) {  // x

                float sum = 0.0;

                // loop over dct coefficients
                for (size_t m = 0; m < len; m++) {
                    for (size_t n = 0; n < len; n++) {
                        for (size_t o = 0; o < len; o++) {

                            if(vec.at(m * (len * len) + n * len + o) == 0.0f) {
                                continue;
                            }

                            // sum += vec.at(m * (len * len) + n * len + o) *
                            //        std::cos(((float)i + 0.5) * (float)m * factor) *
                            //        std::cos(((float)j + 0.5) * (float)n * factor) *
                            //        std::cos(((float)k + 0.5) * (float)o * factor);

                            sum += vec.at(m * (len * len) + n * len + o) * cc[i][m] * cc[j][n] * cc[k][o];
                        }
                    }
                }

                result.push_back(sum);
            }
        }
    }

    return result;
}
