/*
 * newton.h
 *
 * Newton interpolating polynomials for free energy estimates
 * Copyright (C) 2008   Conrad Shyu (conradshyu at hotmail.com)
 * Department of Physics, University of Idaho
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * Author's comments
 * -----------------
 * written by Conrad Shyu (conradshyu at hotmail.com)
 * Department of Physics
 * University of Idaho, Moscow, ID 83844
 *
 * first created on December 18, 2007
 * revised on September 3, 2008
 * revised on March 11, 2014
 */

#ifndef _NEWTON_H
#define _NEWTON_H

#include <list>
#include <cmath>
#include <bitset>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>

const unsigned int NEWTON_DEGREE = sizeof( unsigned int ) * 8;

typedef struct
{
    double x;   // positions on the x-axis
    double y;   // values on the y-axis, y=f(x)
} stNEWTON;

class Newton
{
public:
    Newton();
    Newton( const std::list<stNEWTON>& );
    Newton( const std::vector<double>&, const std::vector<double>& );
    ~Newton() {};

    bool GetEstimate( const std::string&, unsigned int ) const;

    double DoIntegral( bool = false ) const;
    double DoQuadrature( bool = false ) const;

    const std::vector<double>& GetPolynomial( bool = false ) const;
    const std::list<stNEWTON>& LoadData( const std::list<stNEWTON>& );
    const std::list<stNEWTON>& LoadData( const std::vector<double>&, const std::vector<double>& );
private:
    std::list<stNEWTON> sample;
    std::vector<double> factor;

    void ClearData();
    void DoPolynomial();        // construct the newton polynomial
    const std::vector<double>& GetFDD( std::vector<double>& );
    const std::vector<double>& GetPermute( std::vector<double>& );
};  // class definition for lagragne interpolating polynomial

#endif  // _NEWTON_H
