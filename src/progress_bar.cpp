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

#include "progress_bar.h"

#define LENGTH_OF_PROGRESS_BAR 55
#define PERCENTAGE_BIN_SIZE (100.0/LENGTH_OF_PROGRESS_BAR)

std::string ProgressBar::generate_progress_bar(unsigned int percentage) {
    const int progress = static_cast<int>(percentage/PERCENTAGE_BIN_SIZE);
    std::ostringstream ss;
    ss << " " << std::setw(3) << std::right << percentage << "% ";
    std::string bar("[" + std::string(LENGTH_OF_PROGRESS_BAR-2, ' ') + "]");

    unsigned int numberOfSymbols = std::min(
            std::max(0, progress - 1),
            LENGTH_OF_PROGRESS_BAR - 2);

    bar.replace(1, numberOfSymbols, std::string(numberOfSymbols, '|'));

    ss << bar;
    return ss.str();
}

ProgressBar::ProgressBar(uint32_t expectedIterations, const std::string& initialMessage):
    total_iterations(expectedIterations),
    number_of_ticks(0),
    flag_ended(false) {

    // check if tty
    this->has_tty = isatty(fileno(stdout));

    if(this->has_tty) {
        std::cout << initialMessage << "\n";
        length_of_last_message_printed = initialMessage.size();
    }

    // set start point
    this->start = std::chrono::system_clock::now();

    if(this->has_tty) {
        std::cout << generate_progress_bar(0) << "\r" << std::flush;
    }
}

ProgressBar::~ProgressBar() {
    end_progress_bar();
}

void ProgressBar::operator++() {
    if (flag_ended) {
        throw std::runtime_error("Attempted to use progress bar after having terminated it");
    }

    if(!this->has_tty) {
        return;
    }

    number_of_ticks = std::min(total_iterations, number_of_ticks+1);
    const unsigned int percentage = static_cast<unsigned int>(number_of_ticks*100.0/total_iterations);

    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> diff = end - this->start;
    const double time_per_tick = diff.count() / (double)number_of_ticks;
    double expected_total_time = total_iterations * time_per_tick;
    this->update_last_printed_message((boost::format("%6.2f / %6.2f") % diff.count() % expected_total_time).str());

    std::cout << generate_progress_bar(percentage) << "\r" << std::flush;
}

void ProgressBar::print_new_message(const std::string& message) {
    if (flag_ended) {
        throw std::runtime_error("Attempted to use progress bar after having terminated it");
    }

    std::cout << "\r"
              << std::left
              << std::setw(LENGTH_OF_PROGRESS_BAR + 6)
              << message << "\n";

    length_of_last_message_printed = message.size();

    const unsigned int percentage = static_cast<unsigned int>(number_of_ticks*100.0/total_iterations);

    std::cout << generate_progress_bar(percentage) << "\r" << std::flush;

}

void ProgressBar::update_last_printed_message(const std::string& message) {
    if (flag_ended) {
        throw std::runtime_error("Attempted to use progress bar after having terminated it");
    }

    std::cout << "\r\033[F"
              << std::left
              << std::setw(length_of_last_message_printed)
              << message << "\n";

    length_of_last_message_printed = message.size();
}

void ProgressBar::end_progress_bar() {
    if (!flag_ended && this->has_tty) {
        std::cout << std::string(2, '\n');
    }

    flag_ended = true;
}
