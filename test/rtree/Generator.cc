/******************************************************************************
 * Project:  libspatialindex - A C++ library for spatial indexing
 * Author:   Marios Hadjieleftheriou, mhadji@gmail.com
 ******************************************************************************
 * Copyright (c) 2002, Marios Hadjieleftheriou
 *
 * All rights reserved.
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
 * THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
******************************************************************************/


#include <spatialindex/tools/Tools.h>
#include <cmath>
#include <limits>

#include <set>

#define INSERT 1
#define DELETE 0
#define QUERY 2

class Region
{
public:
	double m_xmin, m_ymin, m_xmax, m_ymax;

	Region(double x1, double y1, double x2, double y2)
	{
		m_xmin = (x1 < x2) ? x1 : x2;
		m_ymin = (y1 < y2) ? y1 : y2;
		m_xmax = (x1 > x2) ? x1 : x2;
		m_ymax = (y1 > y2) ? y1 : y2;
	}
};

int main(int argc, char** argv)
{
	if (argc != 3)
	{
		std::cerr << "Usage: " << argv[0] << " number_of_data skewness." << std::endl;
		return -1;
	}

	size_t skewness = atol(argv[2]);
	size_t numberOfObjects = atol(argv[1]);
	// std::map<size_t, Region> data;
	Tools::Random rnd;

	for (size_t i = 0; i < numberOfObjects; i++)
	{
		double x = rnd.nextUniformDouble();
		double y = rnd.nextUniformDouble();
		double y_ = std::pow(y, skewness); 



		// double dx = rnd.nextUniformDouble(0.0001, 0.1);
		// double dy = rnd.nextUniformDouble(0.0001, 0.1);
		// Region r = Region(x, y, x + dx, y + dy);

		// Region r = Region(x, y, x, y );


		// data.insert(std::pair<size_t, Region>(i, r));

		std::cout << INSERT << " " << i << " " << x << " " << y_ << " "
			<< x << " " << y_ << std::endl;
	}


	for (size_t i = 0; i < 1000; i++)
	{
		double stx = rnd.nextUniformDouble();
		double sty = rnd.nextUniformDouble();
		std::cout << QUERY << " 9999999 " << stx << " " << sty << " " << (stx + 0.01) << " " << (sty + 0.01) << std::endl;
	}

	return 0;
}
