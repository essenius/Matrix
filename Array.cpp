#include "Array.h"

const double Array::_epsilon = 1e-12; //std::numeric_limits<double>::epsilon();

Array::Array(const Dimension rows, const Dimension columns) : _data(rows * columns) {
    _rows = rows;
    _columns = columns;
    _arraySize = rows * columns;
}

Array::Array(const std::initializer_list<std::initializer_list<double>> list) {
    _rows = static_cast<Dimension>(list.size());
    _columns = static_cast<Dimension>(list.begin()->size());
    _arraySize = _rows * _columns;
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

Array::Array(const Array& other) {
    _rows = other._rows;
    _columns = other._columns;
    _arraySize = _rows * _columns;
    _data = other._data;
}

double& Array::operator()(Dimension row, Dimension column) {
    return _data[row * _columns + column];
}

const double& Array::operator()(Dimension row, Dimension column) const {
    return _data[row * _columns + column];
}

void Array::operator+=(const Array& other) {
    assertMatchingSize(other);
    for (Dimension cell = 0; cell < _arraySize; cell++) {
        _data[cell] += other._data[cell];
    }
}

void Array::operator-=(const Array& other) {
    assertMatchingSize(other);
    for (Dimension cell = 0; cell < _arraySize; cell++) {
        _data[cell] -= other._data[cell];
    }
}

void Array::operator*=(const Array& other) {
    assertMatchingSize(other);
    for (Dimension cell = 0; cell < _arraySize; cell++) {
        _data[cell] *= other._data[cell];
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
    if (!equalSize(other)) {
        return false;
    }
    for (Dimension cell = 0; cell < _arraySize; cell++) {
        if (abs(_data[cell]- other._data[cell]) > _epsilon) {
            return false;
        }
    }
    return true;
}

void Array::assertMatchingSize(const Array &other) const {
    if (!equalSize(other)) {
        throw std::invalid_argument("Matrix dimensions do not match");
    }
}

Dimension Array::columns() const {
    return _columns;
}

bool Array::equalSize(const Array& other) const {
    return _rows == other._rows && _columns == other._columns;
}

Array Array::getColumn(const Dimension column) const {
    Array result(_rows, 1);
    for (Dimension row = 0; row < _rows; row++) {
        result(row, 0) = me(row, column);
    }
    return result;
}

Array Array::getRow(const Dimension row) const {
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
    return _data[row * _columns + column];
}

Array Array::pow2() const {
    Array result(_rows, _columns);
    for (Dimension cell = 0; cell < _arraySize; cell++) {
        result._data[cell] = _data[cell] * _data[cell];
    }
    return result;
}

Dimension Array::rows() const {
    return _rows;
}

void Array::setColumn(Dimension column, const Array& input) {
    for (Dimension row = 0; row < _rows; row++) {
        (*this)(row, column) = input(row, 0);
    }
}

void Array::setColumn(const Dimension column, const double value) {
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

Array operator/(Array left, const double right) {
    left /= right;
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

Array operator*(double left, Array right) {
    right *= left;
    return right;
}

Array operator+(Array left, const Array& right) {
    left += right;
    return left;
}

Array operator-(Array left, const Array& right) {
    left -= right;
    return left;
}
