// Copyright 2023 Rik Essenius
// 
// Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file
// except in compliance with the License. You may obtain a copy of the License at
// 
//     http://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software distributed under the License
// is distributed on an "AS IS" BASIS WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and limitations under the License.

// Fits a number of points to an ellipse. Uses the Eigen library which is quite large.
// Also it is quite a bit of work, so we optimize where we can to minimize calculations during fits, which need to complete within the time of a sample interval.  

#ifndef ELLIPSEFIT_H
#define ELLIPSEFIT_H

#include "SolverMatrix.h"

class EllipseFit {

public:
	explicit EllipseFit();
	bool addMeasurement(const double x, const double y);
	void begin();
	std::vector<double> fit();
	bool bufferIsFull() const { return _size >= BUFFER_SIZE; }
	unsigned int pointCount() const { return _size; }
    static unsigned int size() { return BUFFER_SIZE; }

private:
	static constexpr unsigned int BUFFER_SIZE = 32;
	Matrix _c1Inverse;
	Matrix _design1;
	Matrix _design2;

	unsigned int _size = 0;
};
#endif
