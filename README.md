# Newton interpolating polynomial for free energy estimates

Written by [Conrad Shyu](mailto:conradshyu@hotmail.com)

*release 2, September 8, 2008*<br>
*release 3, March 11, 2014*<br>
*revised on July 4, 2020*

## Synopsis
This archive contains the implementation of Newton interpolating polynomial for estimates of free energy
differences using thermodynamic integration. It also contains an example data file to jumpstart the use of
software. The implementation makes heavy use of STL (standard template library) which requires a C++ compiler.
This program is self-contained and should not require any additional libraries.

To compile the program, simply type:

`make`

The `make` program will invoke the Makefile and compiles the source code. Alternatively, to compile the files
manually, issue the command:

```
g++ -I. driver.cpp newton.cpp -o newton -lm
```

The archive lists three implementation files and an example file.

| Filename | Description |
| --- | --- |
| `driver.cpp` | the user interface program |
| `newton.cpp` | the main implementation of Newton interpolation polynomial |
| `newton.h` | the header file for Lagrange interpolation polynomial |
| `example.csv` | example file of simulations data from thermodynamic integration |
| `Makefile` | the makefile to compile the source code |
| `README.md` | this file |

The purpose of this program is to construct a Newton interpolating polynomial that fits the thermodynamic
integration data. The implementation also contains an estimate function that allows the interpolating polynomial
be used to extrapolate points. This feature may be useful to determine the accuracy of the interpolating
polynomial.

The Newton interpolation program takes three command-line parameters: `input_file`, `plot_file` (optional), and
`data_points` (optional). The followings detail the format and meaning of each command-line parameter.

### `input_file`
The input file must be saved in the plain text format. The first column contains the lambda values and second the
ensemble averages from the simulations. Each column of data must be separated by either a comma or a white space.
For example, the followings list the contents of the example input file:

```
0.0,  51.49866347
0.1,  23.92508775
0.2,  10.35390700
0.3,   2.58426990
0.4,  -2.18351656
0.5,  -5.41745387
0.6,  -7.62452181
0.7,  -9.25455804
0.8, -10.45592989
0.9, -11.39244138
1.0, -12.12433704
```

To run the analysis program with only the required parameter, type:

`newton example.csv`

The program will output the interpolating polynomial and the estimates of free energy using the Newton
interpolating polynomial and trapezoidal rule. For example, the example file should generate the following output:

```
Degree, Coefficients
     0,      51.49866347
     1,    -440.06771988
     2,    2724.49944200
     3,  -16355.52325204
     4,   74074.45611540
     5, -226406.51705359
     6,  454849.19406749
     7, -590543.95494381
     8,  476225.01322752
     9, -216633.75369268
    10,   42443.03080908

Free energy difference
 Lagrange: 0.65644834
Trapezoid: 1.02220063
```

### `plot_file`
This parameter specifies the name of the file that contains the estimates obtained from the interpolating
polynomial. This parameter is optional. If this parameter is not given, no plot data will be generated by the
program. For example, to generate a file of plot data with the default number of data points extrapolated from the
interpolating polynomial, issue the command:

`newton example.csv plot.csv`

The file, `example.csv`, is the input file that contains the thermodynamic integration data. The file, `plot.csv`,
is the output file that contains the estimates from the interpolating polynomial.

### `data_points`
This parameter specifies the number of data points extrapolated from the interpolating polynomial. This parameter
is optional. If this parameter is not given, the number of data points will be the same as the input file. For
example, to generate a file of plot data with 100 data points extrapolated from the interpolating polynomial,
issue the command:

`newton example.csv plot.csv 100`

## Author's Comments
It is strongly recommended to limit the number of lambda values to less than a dozen. High degree of Newton
polynomials suffers greatly from numerical instability. The analysis from the program can fluctuate widely and
becomes completely unreliable (see the reference for more information). Please report any problems or send
comments to [me](mailto:conradshyu@hotmail.com).

> **Note**: Lagrange and Newton polynomials are mathematically equivalent. The numerical differences are only
notable after the eight digit. The Newton polynomial is generally considered more computationally efficient than
Lagrange.

## Reference
C. Shyu and F.M. Ytreberg (2011). Accurate estimation of solvation free energy using polynomial fitting
techniques. *Journal of Computational Chemistry*, 32(1): 134-141, doi: 10.1002/jcc.21609.

---

Copyright (C) 2014 [Conrad Shyu](mailto:conradshyu@hotmail.com)<br>
Department of Physics<br>
University of Idaho, Moscow, ID 83844

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public
License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any
later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
details.

You should have received a copy of the GNU General Public License along with this program. If not, see
<http://www.gnu.org/licenses/>.