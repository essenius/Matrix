#include "Array.h"
#include <cassert>

Array::Array(const Dimension rows, const Dimension columns) : 
    _data(rows * columns), 
    _rows(rows),
    _columns(columns), 
    _arraySize(rows * columns) {}

Array::Array(const std::initializer_list<std::initializer_list<double>> list):
    _rows(static_cast<Dimension>(list.size())), 
    _columns(static_cast<Dimension>(list.begin()->size())), 
    _arraySize(_rows * _columns) {
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

double& Array::operator[](Dimension cell) {
    assert(cell < _arraySize);
    return _data[cell];
}

const double &Array::operator[](Dimension cell) const {
    assert(cell < _arraySize);
    return _data[cell];
}

double& Array::operator()(Dimension row, Dimension column) {
    assert(row < _rows && column < _columns);
    return _data[row * _columns + column];
}

const double& Array::operator()(Dimension row, Dimension column) const {
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

void Array::operator+=(double other) {
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
        if (abs(_data[cell]- other[cell]) > EPSILON) {
            return false;
        }
    }
    return true;
}


Array Array::getColumn(const Dimension column) const {
    assert(column < _columns);
    Array result(_rows, 1);
    for (Dimension row = 0; row < _rows; row++) {
        result(row, 0) = me(row, column);
    }
    return result;
}

Dimension Array::columnCount() const {
    return _columns;
}

Array Array::getRow(const Dimension row) const {
    assert(row < _rows);
    Array result(1, _columns);
    for (Dimension column = 0; column < _columns; column++) {
        result(0, column) = me(row, column);
    }
    return result;
}

Dimension Array::rowCount() const {
    return _rows;
}

bool Array::isSquare() const {
    return _rows == _columns;
}

double Array::me(const Dimension row, const Dimension column) const {
    assert(row < _rows && column < _columns);
    return _data[row * _columns + column];
}

Array Array::pow2() {
    Array result(*this);
    for (Dimension cell = 0; cell < _arraySize; cell++) {
        result[cell] *= _data[cell];
    }
    return result;
}

void Array::setColumn(Dimension column, const Array& input) {
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

void Array::swapColumns(Dimension column1, Dimension column2) {
    assert(column1 < _columns && column2 < _columns);
    if (column1 == column2) return;
    for (Dimension row = 0; row < _rows; row++) {
        double temp = me(row, column1);
        (*this)(row, column1) = me(row, column2);
        (*this)(row, column2) = temp;
    }
}

void Array::swapRows(Dimension row1, Dimension row2) {
    assert(row1 < _rows && row2 < _rows);
    if (row1 == row2) return;
    for (Dimension column = 0; column < _columns; column++) {
        double temp = me(row1, column);
        (*this)(row1, column) = me(row2, column);
        (*this)(row2, column) = temp;
    }
}

Dimension Array::size() const {
    return _arraySize;
}

bool Array::sizeIsEqual(const Array& other) const {
    return _rows == other._rows && _columns == other._columns;
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

Array operator*(double left, Array right) {
    right *= left;
    return right;
}
