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

#include "Array.h"
#include <cassert>
#include <cmath>

namespace RixMatrix {
	Array::Array(const Dimension rows, const Dimension columns) :
		_data(rows* columns),
		_rows(rows),
		_columns(columns),
		_arraySize(rows* columns) {
	}

	Array::Array(const std::initializer_list<std::initializer_list<double>> list) :
		_rows(static_cast<Dimension>(list.size())),
		_columns(static_cast<Dimension>(list.begin()->size())),
		_arraySize(_rows* _columns) {
		_data = std::vector<double>(_rows * _columns);
		Dimension row = 0;
		for (auto rowList : list) {
			Dimension column = 0;
			for (const auto value : rowList) {
				_data[row * _columns + column] = value;
				column++;
			}
			row++;
		}
	}

	Array::Array(Array&& other) noexcept {
		_data = std::move(other._data);
		_rows = other._rows;
		_columns = other._columns;
		_arraySize = other._arraySize;
		other._rows = 0;
		other._columns = 0;
		other._arraySize = 0;
		other._data.clear();
	}

	Array::Array(const Array& other) {
		_data = other._data;
		_rows = other._rows;
		_columns = other._columns;
		_arraySize = other._arraySize;
	}

	double& Array::operator[](const Dimension cell) {
		assert(cell < _arraySize);
		return _data[cell];
	}

	const double& Array::operator[](const Dimension cell) const {
		assert(cell < _arraySize);
		return _data[cell];
	}

	double& Array::operator()(const Dimension row, const Dimension column) {
		assert(row < _rows && column < _columns);
		return _data[row * _columns + column];
	}

	const double& Array::operator()(const Dimension row, const Dimension column) const {
		assert(row < _rows && column < _columns);
		return _data[row * _columns + column];
	}

	void Array::operator+=(const Array& other) {
		assert(other.sizeIsEqual(*this));
		for (Dimension cell = 0; cell < _arraySize; cell++) {
			_data[cell] += other[cell];
		}
	}

	void Array::operator-=(const Array& other) {
		assert(other.sizeIsEqual(*this));
		for (Dimension cell = 0; cell < _arraySize; cell++) {
			_data[cell] -= other[cell];
		}
	}

	void Array::operator*=(const Array& other) {
		assert(other.sizeIsEqual(*this));
		for (Dimension cell = 0; cell < _arraySize; cell++) {
			_data[cell] *= other[cell];
		}
	}

	void Array::operator/=(const double other) {
		for (Dimension cell = 0; cell < _arraySize; cell++) {
			_data[cell] /= other;
		}
	}

	void Array::operator+=(const double other) {
		for (Dimension cell = 0; cell < _arraySize; cell++) {
			_data[cell] += other;
		}
	}

	void Array::operator-=(const double other) {
		for (Dimension cell = 0; cell < _arraySize; cell++) {
			_data[cell] -= other;
		}
	}

	void Array::operator*=(const double other) {
		for (Dimension cell = 0; cell < _arraySize; cell++) {
			_data[cell] *= other;
		}
	}

	bool Array::operator==(const Array& other) const {
		if (!sizeIsEqual(other)) {
			return false;
		}
		for (Dimension cell = 0; cell < _arraySize; cell++) {
			if (abs(_data[cell] - other[cell]) > Epsilon) {
				return false;
			}
		}
		return true;
	}

	Dimension Array::columnCount() const {
		return _columns;
	}

	Array Array::getColumn(const Dimension column) const {
		assert(column < _columns);
		Array result(_rows, 1);
		for (Dimension row = 0; row < _rows; row++) {
			result(row, 0) = me(row, column);
		}
		return result;
	}

	Array Array::getRow(const Dimension row) const {
		assert(row < _rows);
		Array result(1, _columns);
		for (Dimension column = 0; column < _columns; column++) {
			result(0, column) = me(row, column);
		}
		return result;
	}

	bool Array::isSquare() const {
		return _rows == _columns;
	}

	double Array::me(const Dimension row, const Dimension column) const {
		assert(row < _rows && column < _columns);
		return _data[row * _columns + column];
	}

	Array Array::pow2() const {
		Array result(*this);
		for (Dimension cell = 0; cell < _arraySize; cell++) {
			result[cell] *= _data[cell];
		}
		return result;
	}

	Dimension Array::rowCount() const {
		return _rows;
	}

	void Array::setColumn(const Dimension column, const Array& input) {
		assert(column < _columns && input.rowCount() == _rows);
		for (Dimension row = 0; row < _rows; row++) {
			(*this)(row, column) = input(row, 0);
		}
	}

	void Array::setColumn(const Dimension column, const double value) {
		assert(column < _columns);
		for (Dimension row = 0; row < _rows; row++) {
			(*this)(row, column) = value;
		}
	}

	void Array::setColumnCount(const Dimension columns) {
		if (columns == _columns) {
			return;
		}
		Array result(_rows, columns);
		const Dimension maxColumns = std::min(_columns, columns);
		for (Dimension row = 0; row < _rows; row++) {
			for (Dimension column = 0; column < maxColumns; column++) {
				result(row, column) = me(row, column);
			}
		}
		*this = result;
	}

	void Array::setRow(const Dimension row, const Array& input) {
		assert(row < _rows && input.columnCount() == _columns);
		for (Dimension column = 0; column < _columns; column++) {
			(*this)(row, column) = input(0, column);
		}
	}

	void Array::setRow(const Dimension row, const double value) {
		assert(row < _rows);
		for (Dimension column = 0; column < _columns; column++) {
			(*this)(row, column) = value;
		}
	}

	void Array::setRowCount(const Dimension rows) {
		if (rows == _rows) {
			return;
		}
		Array result(rows, _columns);
		const Dimension maxRows = std::min(_rows, rows);
		for (Dimension row = 0; row < maxRows; row++) {
			for (Dimension column = 0; column < _columns; column++) {
				result(row, column) = me(row, column);
			}
		}
		*this = result;
	}

	Dimension Array::size() const {
		return _arraySize;
	}

	bool Array::sizeIsEqual(const Array& other) const {
		return _rows == other._rows && _columns == other._columns;
	}

	void Array::swapColumns(const Dimension column1, const Dimension column2) {
		assert(column1 < _columns && column2 < _columns);
		if (column1 == column2) return;
		for (Dimension row = 0; row < _rows; row++) {
			const double temp = me(row, column1);
			(*this)(row, column1) = me(row, column2);
			(*this)(row, column2) = temp;
		}
	}

	void Array::swapRows(const Dimension row1, const Dimension row2) {
		assert(row1 < _rows && row2 < _rows);
		if (row1 == row2) return;
		for (Dimension column = 0; column < _columns; column++) {
			const double temp = me(row1, column);
			(*this)(row1, column) = me(row2, column);
			(*this)(row2, column) = temp;
		}
	}

	Array Array::transposed() const {
		Array result(columnCount(), rowCount());
		for (Dimension row = 0; row < rowCount(); row++) {
			for (Dimension column = 0; column < columnCount(); column++) {
				result(column, row) = me(row, column);
			}
		}
		return result;
	}

	Array operator+(Array left, const Array& right) {
		left += right;
		return left;
	}

	Array operator-(Array left, const Array& right) {
		left -= right;
		return left;
	}

	Array operator*(Array left, const Array& right) {
		left *= right;
		return left;
	}

	Array operator*(Array left, const double right) {
		left *= right;
		return left;
	}

	Array operator/(Array left, const double right) {
		left /= right;
		return left;
	}

	Array operator*(const double left, Array right) {
		right *= left;
		return right;
	}
}