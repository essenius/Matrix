#include "BasicMatrix.h"

#include <stdexcept>
#include <cmath>
#include <iostream>

BasicMatrix::BasicMatrix(const Dimension rows, const Dimension columns) : Array(rows, columns) {} 

BasicMatrix::BasicMatrix(const BasicMatrix& other) : Array(other){}

BasicMatrix::BasicMatrix(const Array other) : Array(other){}

BasicMatrix::BasicMatrix(const std::initializer_list<std::initializer_list<double>> list) : Array(list) {}


void BasicMatrix::operator*=(const BasicMatrix& other) {
    if (_columns != other._rows) {
        throw std::invalid_argument("Matrix dimensions do not match");
    }
    BasicMatrix result(_rows, other._columns);
    for (Dimension row = 0; row < _rows; row++) {
        for (Dimension otherColumn = 0; otherColumn < other._columns; otherColumn++) {
            double sum = 0;
            for (Dimension column = 0; column < _columns; column++) {
                sum += me(row, column) * other(column, otherColumn);
            }
            result(row, otherColumn) = sum;
        }
    }
    *this = result;
}

BasicMatrix BasicMatrix::adjoint() const {
    BasicMatrix result(_rows, _columns);
    for (Dimension row = 0; row < _rows; row++) {
        for (Dimension column = 0; column < _columns; column++) {
            result(row, column) = cofactor(row, column);
        }
    }
    return result;
}

BasicMatrix BasicMatrix::adjugate() const {
    return adjoint().transpose();
}


double BasicMatrix::cofactor(Dimension row, Dimension column) const {
    BasicMatrix minor(_rows - 1, _columns - 1);
    for (Dimension subRow = 0; subRow < _rows; subRow++) {
        for (Dimension subColumn = 0; subColumn < _columns; subColumn++) {
            if (subRow < row && subColumn < column) {
                minor(subRow, subColumn) = me(subRow, subColumn);
            }
            if (subRow < row && subColumn > column) {
                minor(subRow, subColumn - 1) = me(subRow, subColumn);
            }
            if (subRow > row && subColumn < column) {
                minor(subRow - 1, subColumn) = me(subRow, subColumn);
            }
            if (subRow > row && subColumn > column) {
                minor(subRow - 1, subColumn - 1) = me(subRow, subColumn);
            }
        }
    }
    return minor.getDeterminant() * pow(-1, row + column);
}

BasicMatrix BasicMatrix::cubic() const {
    auto result = squared();
    result *= *this;
    return result;
}

double BasicMatrix::getDeterminant() const {
    if (!isSquare()) {
        throw std::invalid_argument("Matrix is not square");
    }
    if (_rows == 1) {
        return _data[0];
    }
    if (_rows == 2) {
        return _data[0] * _data[3] - _data[1] * _data[2];
    }
    double result = 0;
    for (Dimension row = 0; row < _rows; row++) {
        BasicMatrix minor(_rows - 1, _columns - 1);
        for (Dimension subRow = 1; subRow < _rows; subRow++) {
            for (Dimension column = 0; column < _columns; column++) {
                if (column < row) {
                    minor(subRow - 1, column) = me(subRow, column);  //_data[subRow * _columns + column];
                }
                else if (column > row) {
                    minor(subRow - 1, column - 1) = me(subRow, column);
                }
            }
        }
        result += _data[row] * minor.getDeterminant() * (row % 2 == 0 ? 1 : -1);
    }
    return result;
}

BasicMatrix BasicMatrix::diagonal(const std::initializer_list<double> list) {
    const auto size = static_cast<Dimension>(list.size());
    BasicMatrix result(size, size);
    Dimension diagonalCell = 0;
    for (const auto value : list) {
        result(diagonalCell, diagonalCell) = value;
        diagonalCell++;
    }
    return result;
}

/*BasicMatrix BasicMatrix::getColumn(const Dimension column) const {
    BasicMatrix result(_rows, 1);
    for (Dimension row = 0; row < _rows; row++) {
        result(row, 0) = me(row, column);
    }
    return result;
}

BasicMatrix BasicMatrix::getRow(const Dimension row) const {
    BasicMatrix result(1, _columns);
    for (Dimension column = 0; column < _columns; column++) {
        result(0, column) = me(row, column);
    }
    return result;
}*/

double BasicMatrix::getTrace() const {
    if (!isSquare()) {
        throw std::invalid_argument("Matrix is not square");
    }
    double result = 0;
    for (Dimension row = 0; row < _rows; row++) {
        result += me(row, row);
    }
    return result;
}

BasicMatrix BasicMatrix::identity(Dimension size) {
    BasicMatrix result(size, size);
    for (Dimension diagonalCell = 0; diagonalCell < size; diagonalCell++) {
        result(diagonalCell, diagonalCell) = 1;
    }
    return result;
}

BasicMatrix BasicMatrix::inverse() const {
    if (!isInvertible()) {
        throw std::invalid_argument("Matrix is not invertible");
    }
    return adjugate() / getDeterminant();
}

bool BasicMatrix::isInvertible() const {
    return isSquare() && abs(getDeterminant()) > _epsilon;
}



BasicMatrix BasicMatrix::normalize() const {
    double norm = 0;
    for (const auto entry : _data) {
        norm += entry * entry;
    }
    norm = sqrt(norm);
    if (abs(norm) > _epsilon) return BasicMatrix(*this / norm);
    return *this;
}

/*void BasicMatrix::setAt(const BasicMatrix& input, const Dimension startRow, const Dimension startColumn) {
    for (Dimension row = 0; row < input._rows; row++) {
        for (Dimension column = 0; column < input._columns; column++) {
            (*this)(startRow + row, startColumn + column) = input(row, column);
        }
    }
}

void BasicMatrix::setAll(const double value) {
    for (double& entry : _data) {
        entry = value;
    }
}

void BasicMatrix::setRow(const Dimension row, const BasicMatrix& input) {
    for (Dimension column = 0; column < _columns; column++) {
        (*this)(row, column) = input(0, column);
    }
}*/

BasicMatrix BasicMatrix::squared() const {
    return *this * *this;
}

BasicMatrix BasicMatrix::transpose() const {
    BasicMatrix result(_columns, _rows);
    for (Dimension row = 0; row < _rows; row++) {
        for (Dimension column = 0; column < _columns; column++) {
            result(column, row) = me(row, column);
        }
    }
    return result;
}

BasicMatrix operator*(BasicMatrix left, const BasicMatrix& right) {
    left *= right;
    return left;
}
