/*
 * newton.cpp
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

#include <newton.h>

/*
 * default class constructor
*/
Newton::Newton()
{
    ClearData();
}   // end of class constructor

/*
 * class constructor
*/
Newton::Newton(
    const std::list<stNEWTON>& _sample )
{
    LoadData( _sample );
}   // end of class constructor

/*
 * class constructor
*/
Newton::Newton(
    const std::vector<double>& _x,
    const std::vector<double>& _y )
{
    LoadData( _x, _y );
}   // end of class constructor

/*
 * reset and initialize essential variables
*/
const std::list<stNEWTON>& Newton::LoadData(
    const std::list<stNEWTON>& _sample )
{
    stNEWTON unit; ClearData();

    for ( std::list<stNEWTON>::const_iterator i = _sample.begin(); !( i == _sample.end() ); i++ )
    {
        unit.x = ( *i ).x; unit.y = ( *i ).y; sample.push_back( unit );
    }   // save a local copy of the data

    // perform interpolation with newton polynomials
    DoPolynomial(); return( sample );
}   // end of LoadData()

/*
 * reset and initialize essential variables
*/
const std::list<stNEWTON>& Newton::LoadData(
    const std::vector<double>& _x,
    const std::vector<double>& _y )
{
    stNEWTON unit; ClearData();

    for ( unsigned int i = 0; i < _x.size(); ++i )
    {
        unit.x = _x[ i ]; unit.y = _y[ i ]; sample.push_back( unit );
    }   // save a local copy of the data

    // perform interpolation with newton polynomials
    DoPolynomial(); return( sample );
}   // end of LoadData()

/*
 * clear all contents
*/
void Newton::ClearData()
{
    sample.clear(); factor.clear();
}   // end of ClearData()

/*
 * calculate the divided differences
 * note: the calculation of the forward divided differences has been verified correcly
 * on december 22, 2007
*/
const std::vector<double>& Newton::GetFDD(
    std::vector<double>& _fdd )
{
    std::vector<double> next;
    std::vector<double> x;
    std::vector<double> y;
    _fdd.clear(); y.clear(); x.clear();

    for ( std::list<stNEWTON>::iterator i = sample.begin(); !( i == sample.end() ); i++ )
    {
        x.push_back( ( *i ).x ); y.push_back( ( *i ).y );
    }   // make a local copy of the x and y elements

    _fdd.push_back( y[ 0 ] );       // only save the first element in the list

    for ( unsigned int s = 1; s < sample.size(); ++s )
    {
        next.clear();

        for ( unsigned int t = 0; t < ( y.size() - 1 ); ++t )
        {
            next.push_back( ( y[ t + 1 ] - y[ t ] ) / ( x[ ( t + s ) ] - x[ t ] ) );
        }   // calculate the forward divided differences

        y = next; _fdd.push_back( y[ 0 ] );     // update forward divided difference
    }   // iterate through all elements in the vector

    return( _fdd );
}   // end of GetTerm()

/*
 * permute the coefficients for the construction of the polynomial
*/
const std::vector<double>& Newton::GetPermute(
    std::vector<double>& _x )
{
    unsigned int bit_value = ( 0x1 << _x.size() );
    std::vector<double> term( ( _x.size() + 1 ), 0.0 );
    std::bitset<NEWTON_DEGREE> permute;
    double unit;

    for ( unsigned int i = 0; i < bit_value; ++i )
    {
        unit = 1.0; permute = i;

        for ( unsigned int j = 0; j < _x.size(); ++j )
        {
            unit *= ( permute[ j ] ) ? ( -1.0 * _x[ j ] ) : 1.0;
        }   // calculate combinatoric terms; -a_1 * -a_2 * ...

        term[ permute.count() ] += unit;
    }   // iterate through all possible combinations of factors

    _x = term; return( _x );    // invoke copy constructor and overwrite contents
}   // end of GetPermute()

/*
 * return the newton interpolating polynomial
*/
const std::vector<double>& Newton::GetPolynomial(
    bool _print ) const
{
    if ( _print )
    {
        printf( "Degree, Coefficients\n" );

        for ( unsigned int i = 0; i < factor.size(); ++i )
        {
            printf( "%6d, %.8f\n", i, factor[ i ] );
        }   // print out the newton polynomial
    }   // print out the coefficients of the polynomial, if necessary

    return( factor );
}   // end of GetPolynomial()

/*
 * construct the newton polynomial
 * note: the construction of newton polynomial has been verified to work correctly on
 * december 22, 2007
*/
void Newton::DoPolynomial()
{
    std::vector<double> fdd;
    std::vector<double> x;
    std::vector<double> term;
    x.clear(); term.clear(); factor.clear();    // empty the lists
    std::list<stNEWTON>::iterator i = sample.begin();
    factor.resize( sample.size(), 0.0 );
    term.push_back( 1.0 );
    GetFDD( fdd );      // calculate the newton forward divided differences

    for ( unsigned int s = 0; s < fdd.size(); ++s )
    {
        for ( unsigned int t = 0; t < term.size(); ++t )
        {
            factor[ t ] += ( fdd[ s ] * term[ term.size() - t - 1 ] );
        }   // accumulate the coefficients with the same power

        x.push_back( ( *i ).x ); i++;   // iteratively add coefficient to the list
        term = x; GetPermute( term );   // calcualte the coefficients for polynomial
    }   // construct the newton polynomial
}   // end of DoPolynomial()

/*
 * perform integration on the newton polynomial
*/
double Newton::DoIntegral(
    bool _print ) const
{
    double lower = ( sample.front() ).x;
    double upper = ( sample.back() ).x;
    double area = 0.0;
    double power;

    for ( unsigned int i = 0; i < factor.size(); ++i )
    {
        power = static_cast<double>( i + 1.0 );
        area += ( ( pow( upper, power ) / power ) * factor[ i ] -
            ( pow( lower, power ) / power ) * factor[ i ] );
    }   // perform the integration with given upper and lower bounds

    if ( _print )
    {
        printf( "area under the curve: %.8f\n", area );
    }   // print out the integration result

    return( area );
}   // end of DoIntegral()

/*
 * calculate the area under the curve using quadrature
*/
double Newton::DoQuadrature(
    bool _print ) const
{
    std::list<stNEWTON>::const_iterator a = sample.begin();
    std::list<stNEWTON>::const_iterator b = sample.begin(); b++;
    double area = 0.0;

    while ( !( b == sample.end() ) )
    {
        area += ( ( *b ).y + ( *a ).y ) * 0.5 * ( ( *b ).x - ( *a ).x );
        a = b; b++;
    }   // iterate through the entire list

    if ( _print )
    {
        printf( "area under the curve: %.8f\n", area );
    }   // print out the integration result

    return( area );
}   // end of DoQuadrature()

/*
 * get the estimate of f(x) with a given value
*/
bool Newton::GetEstimate(
    const std::string& _file, unsigned int _step ) const
{
    std::ofstream ofs( _file.c_str(), std::ios::trunc );

    if ( ofs.bad() )
    {
        std::cout << "file " << _file << "cannot be opened" << std::endl; return( false );
    }   // make sure the file stream has been opened successfully

    char buffer[ 80 ]; double y = 0.0; double x = 0.0;
    double step = 1.0 / static_cast<double>( _step );

    for ( unsigned int s = 0; !( s > _step ); ++s )
    {
        for ( unsigned int i = 0; i < factor.size(); ++i )
        {
            y += ( pow( x, static_cast<double>( i ) ) * factor[ i ] );
        }   // substitute the given value and calculate the estimate

        sprintf( buffer, "%.4f, %.8f", x, y ); ofs << buffer << std::endl;
        y = 0.0; x += step;
    }   // iterate through the entire interval

    ofs.close(); return( true );
}   // end of GetEstimate()

/*
int main( void )
{
    list<stNEWTON> sample;
    stNEWTON unit;

    sample.clear();

    unit.x = 0.0; unit.y =  51.49866347; sample.push_back( unit );
    unit.x = 0.1; unit.y =  23.92508775; sample.push_back( unit );
    unit.x = 0.2; unit.y =  10.35390700; sample.push_back( unit );
    unit.x = 0.3; unit.y =   2.58426990; sample.push_back( unit );
    unit.x = 0.4; unit.y =  -2.18351656; sample.push_back( unit );
    unit.x = 0.5; unit.y =  -5.41745387; sample.push_back( unit );
    unit.x = 0.6; unit.y =  -7.62452181; sample.push_back( unit );
    unit.x = 0.7; unit.y =  -9.25455804; sample.push_back( unit );
    unit.x = 0.8; unit.y = -10.45592989; sample.push_back( unit );
    unit.x = 0.9; unit.y = -11.39244138; sample.push_back( unit );
    unit.x = 1.0; unit.y = -12.12433704; sample.push_back( unit );

    Newton a( sample );
    a.DoIntegral( true ); a.GetPolynomial( true );

    return( 1 );
}   // end of main()
*/
