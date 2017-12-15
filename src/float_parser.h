/****************************************************************************
#                                                                           #
#    INORGANIC MATERIALS CHEMISTRY CONFIDENTIAL                             #
#                                                                           #
#    Copyright 2016 Inorganic Materials Chemistry                           #
#                   Eindhoven University of Technology (TU/e)               #
#                                                                           #
#    All Rights Reserved.                                                   #
#                                                                           #
#    NOTICE:  All information contained herein is, and remains              #
#    the property of Inorganic Materials Chemistry (TU/e) and its suppliers,#
#    if any.  The intellectual and technical concepts contained             #
#    herein are proprietary to Inorganic Materials Chemistry (TU/e)         #
#    and its suppliers and may be covered by U.S. and Foreign Patents,      #
#    patents in process, and are protected by trade secret or copyright law.#
#    Dissemination of this information or reproduction of this Materials    #
#    is strictly forbidden unless prior written permission is obtained      #
#    from Inorganic Materials Chemistry (TU/e).                             #
#                                                                           #
#    Authors: Ivo Filot       <i.a.w.filot@tue.nl>                          #
#             Emiel Hensen    <e.j.m.hensen@tue.nl>                         #
#                                                                           #
*****************************************************************************/

#ifndef _FLOAT_PARSER_H
#define _FLOAT_PARSER_H

#include <vector>
#include <boost/spirit/include/qi.hpp>

struct float_parser : boost::spirit::qi::grammar<std::string::const_iterator,
                                                  std::vector<float>(),
                                                  boost::spirit::ascii::space_type> {

  float_parser() : float_parser::base_type( vector ) {
    vector  %= +(boost::spirit::qi::float_);
  }

  boost::spirit::qi::rule<std::string::const_iterator, std::vector<float>(), boost::spirit::ascii::space_type> vector;
};

#endif //_FLOAT_PARSER_H
