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

#ifndef PROGRESS_BAR_H
#define PROGRESS_BAR_H

#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <chrono>

#include <boost/core/noncopyable.hpp>
#include <boost/format.hpp>

/**
 * RAII implementation of a progress bar.
 *
 * Implementation based on:
 * https://stackoverflow.com/questions/46717885/make-boostprogress-display-work-when-more-information-is-being-printed-to-the
 */
class ProgressBar : private boost::noncopyable {
private:
    unsigned int total_iterations;
    unsigned int number_of_ticks;
    bool flag_ended;
    size_t length_of_last_message_printed;

    std::chrono::time_point<std::chrono::system_clock> start;

public:
    /**
     * Constructor.
     * It takes two values: the expected number of iterations whose progress we
     * want to monitor and an initial message to be displayed on top of the bar
     * (which can be updated with update_last_printed_message()).
     */
    ProgressBar(uint32_t expectedIterations, const std::string& initialMessage="Time remaining: unknown");

    /**
     * Destructor to guarantee RAII.
     */
    ~ProgressBar();

    /**
     * Must be invoked when the progress bar is no longer needed to restore the
     * position of the cursor to the end of the output.
     * It is automatically invoked when the object is destroyed.
     */
    void end_progress_bar();

    /**
     * Prints a new message under the last printed message, without overwriting
     * it. This moves the progress bar down to be placed under the newly
     * written message.
     */
    void print_new_message(const std::string& message);

    /**
     * Prints a message while the progress bar is on the screen on top on the
     * last printed message. Since the cursor is right at the beginning of the
     * progress bar, it moves the cursor up by one line before printing, and
     * then returns it to its original position.
     */
    void update_last_printed_message(const std::string& message);

    /**
     * Overloaded prefix operator, used to indicate that the has been a new
     * iteration.
     */
    void operator++();

private:
    std::string generate_progress_bar(unsigned int percentage);
};

#endif /* PROGRESS_BAR_H */
