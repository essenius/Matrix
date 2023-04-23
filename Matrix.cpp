#include "Matrix.h"
#include <stdexcept>
#include <cmath>

Matrix::Matrix(const Dimension rows, const Dimension columns) : Array(rows, columns) {} 

Matrix::Matrix(const Matrix& other) : Array(other){}

Matrix::Matrix(const Array& other) : Array(other){}

Matrix::Matrix(const std::initializer_list<std::initializer_list<double>> list) : Array(list) {}


void Matrix::operator*=(const Matrix& other) {
    if (_columns != other._rows) {
        throw std::invalid_argument("Matrix dimensions do not match");
    }
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

Matrix Matrix::cubic() const {
    auto result = squared();
    result *= *this;
    return result;
}

double Matrix::getDeterminant() const {
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

Matrix Matrix::diagonal(const std::initializer_list<double> list) {
    const auto size = static_cast<Dimension>(list.size());
    Matrix result(size, size);
    Dimension diagonalCell = 0;
    for (const auto value : list) {
        result(diagonalCell, diagonalCell) = value;
        diagonalCell++;
    }
    return result;
}

/*Matrix Matrix::getColumn(const Dimension column) const {
    Matrix result(_rows, 1);
    for (Dimension row = 0; row < _rows; row++) {
        result(row, 0) = me(row, column);
    }
    return result;
}

Matrix Matrix::getRow(const Dimension row) const {
    Matrix result(1, _columns);
    for (Dimension column = 0; column < _columns; column++) {
        result(0, column) = me(row, column);
    }
    return result;
}*/

double Matrix::getTrace() const {
    if (!isSquare()) {
        throw std::invalid_argument("Matrix is not square");
    }
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
    if (!isInvertible()) {
        throw std::invalid_argument("Matrix is not invertible");
    }
    return adjugate() / getDeterminant();
}

bool Matrix::isInvertible() const {
    return isSquare() && abs(getDeterminant()) > _epsilon;
}

Matrix Matrix::normalize() const {
    double norm = 0;
    for (const auto entry : _data) {
        norm += entry * entry;
    }
    norm = sqrt(norm);
    if (abs(norm) > _epsilon) return Matrix(*this / norm);
    return *this;
}

/*void Matrix::setAt(const Matrix& input, const Dimension startRow, const Dimension startColumn) {
    for (Dimension row = 0; row < input._rows; row++) {
        for (Dimension column = 0; column < input._columns; column++) {
            (*this)(startRow + row, startColumn + column) = input(row, column);
        }
    }
}

void Matrix::setAll(const double value) {
    for (double& entry : _data) {
        entry = value;
    }
}

void Matrix::setRow(const Dimension row, const Matrix& input) {
    for (Dimension column = 0; column < _columns; column++) {
        (*this)(row, column) = input(0, column);
    }
}*/

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
