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

#include "EllipseFit.h"
#include <iostream>

// see https://autotrace.sourceforge.net/WSCG98.pdf for an explanation of the algorithm
// optimized version of https://github.com/mericdurukan/ellipse-fitting

// fixed buffer size as we don't want to fragment the heap
constexpr unsigned int EllipseFit::BUFFER_SIZE;  // NOLINT(readability-redundant-declaration) -- using C++ 11 on device

EllipseFit::EllipseFit() : _c1Inverse(3, 3), _design1(BUFFER_SIZE, 3), _design2(BUFFER_SIZE, 3) {
	// Calculating the constant vectors/matrices once to save a bit of time when fitting
	Matrix c1({{0, 0, 2},{ 0, -1, 0},{ 2, 0, 0}});
	_c1Inverse = c1.inverse();
	_design2.setColumn(2, 1); // set the last column to all ones
}

bool EllipseFit::addMeasurement(double x, double y) {
	if (_size >= BUFFER_SIZE) return false;
	// Design matrix (D in the article). Building up incrementally to minimize compute at fit
	_design1( _size, 0) = x * x;
	_design1(_size, 1) = x * y;
	_design1(_size, 2) = y * y;
	_design2( _size, 0) = x;
	_design2( _size, 1) = y;
	_size++;
	return true;
}

void EllipseFit::begin() {
	_size = 0;
}

std::vector<double> EllipseFit::fit() {

    // Scatter matrix (S in the article)
 	const Matrix scatter1 = _design1.transpose() * _design1;
	const Matrix scatter2 = _design1.transpose() * _design2;
	const Matrix scatter3 = _design2.transpose() * _design2;

	// reduced scatter matrix (M in the article)
	const Matrix reducedScatter = _c1Inverse * (scatter1 - scatter2 * scatter3.inverse() * scatter2.transpose());

	const auto eigenvectors = reducedScatter.getEigenvectors();

	// 4ac - b^2 must be positive for an ellipse
	Eigen::VectorXd condition = 4 * (eigenvector.row(0).array() * eigenvector.row(2).array()) - eigenvector.row(1).array().pow(2);

	// a1 and a2 are the coefficients: a1 = (a,b,c), a2 = (d,f,g).
	Eigen::VectorXd a1;

	// there should be one where the condition is positive, that's the one we need

	for (int i = 0; i < 3; i++) {
		if (condition(i) > 0) {
			a1 = eigenvector.col(i);
			break;
		}
		// we didn't find it (should not happen)
		if (i == 2) {
			/* TODO: remove this print statement */
			printf("No positive eigenvector found");
			return {0, 0, 0, 0, 0, 0};
		}
	}

    const Eigen::VectorXd a2 = -1 * scatter3.inverse() * scatter2.transpose() * a1;

	return {a1(0), a1(1), a1(2), a2(0), a2(1), a2(2)};
}
