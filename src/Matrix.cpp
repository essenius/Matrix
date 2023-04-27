#include "Matrix.h"
#include <stdexcept>
#include <cmath>

Matrix::Matrix(const Dimension rows, const Dimension columns) : Array(rows, columns) {} 

Matrix::Matrix(const Matrix& other) : Array(other){}

Matrix::Matrix(const Array& other) : Array(other){}

Matrix::Matrix(const std::initializer_list<std::initializer_list<double>> list) : Array(list) {}

void Matrix::operator*=(const Matrix& other) {
    assert(_columns == other._rows);
    Matrix result(_rows, other._columns);
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

void Matrix::operator*=(double other){
    Array::operator*=(other);
}

Matrix Matrix::adjoint() const {
    Matrix result(_rows, _columns);
    for (Dimension row = 0; row < _rows; row++) {
        for (Dimension column = 0; column < _columns; column++) {
            result(row, column) = cofactor(row, column);
        }
    }
    return result;
}

Matrix Matrix::adjugate() const {
    return adjoint().transpose();
}

double Matrix::cofactor(Dimension row, Dimension column) const {
    assert(row < _rows && isSquare());
    Matrix minor(_rows - 1, _columns - 1);
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

double Matrix::getDeterminant() const {
    assert(isSquare());
    if (_rows == 1) {
        return _data[0];
    }
    if (_rows == 2) {
        return _data[0] * _data[3] - _data[1] * _data[2];
    }
    double result = 0;
    for (Dimension row = 0; row < _rows; row++) {
        Matrix minor(_rows - 1, _columns - 1);
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

double Matrix::getTrace() const {
    assert(isSquare());
    double result = 0;
    for (Dimension row = 0; row < _rows; row++) {
        result += me(row, row);
    }
    return result;
}

Matrix Matrix::identity(Dimension size) {
    Matrix result(size, size);
    for (Dimension diagonalCell = 0; diagonalCell < size; diagonalCell++) {
        result(diagonalCell, diagonalCell) = 1;
    }
    return result;
}

Matrix Matrix::inverse() const {
    assert(isInvertible());
    return adjugate() / getDeterminant();
}

bool Matrix::isInvertible() const {
    return isSquare() && abs(getDeterminant()) > EPSILON;
}

Matrix Matrix::normalize() const {
    double norm = 0;
    for (const auto entry : _data) {
        norm += entry * entry;
    }
    norm = sqrt(norm);
    if (norm < EPSILON) return *this;
    if (_data[0] < 0) norm = -norm;
    return Matrix(*this / norm);
}

Matrix Matrix::squared() const {
    return *this * *this;
}

Array Matrix::toArray() const {
    Array result(_rows, _columns);
    
    for (Dimension row = 0; row < _rows; row++) {
        for (Dimension column = 0; column < _columns; column++) {
            result(row, column) = me(row, column);
        }
    }
    return result; 
}

Matrix Matrix::transpose() const {
    Matrix result(_columns, _rows);
    for (Dimension row = 0; row < _rows; row++) {
        for (Dimension column = 0; column < _columns; column++) {
            result(column, row) = me(row, column);
        }
    }
    return result;
}

Matrix operator+(Matrix left, const Matrix &right) {
    left += right;
    return left;
}

Matrix operator-(Matrix left, const Matrix &right) {
    left -= right;
    return left;
}

Matrix operator*(Matrix left, const Matrix &right) {
    left *= right;
    return left;
}

Matrix operator*(Matrix left, double right) {
    left *= right;
    return left;
}

Matrix operator*(double left, Matrix right) {
    right *= left;
    return right;
}
