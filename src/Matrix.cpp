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

#include "Matrix.h"
#include <cassert>
#include <stdexcept>
#include <cmath>

namespace RixMatrix {

    Matrix::Matrix(const Dimension rows, const Dimension columns) : Array(rows, columns) {}

    Matrix::Matrix(const Array& other) : Array(other) {}

    Matrix::Matrix(const std::initializer_list<std::initializer_list<double>> list) : Array(list) {}

    void Matrix::operator*=(const Matrix& other) {
        assert(columnCount() == other.rowCount());
        Matrix result(rowCount(), other.columnCount());
        for (Dimension row = 0; row < rowCount(); row++) {
            for (Dimension otherColumn = 0; otherColumn < other.columnCount(); otherColumn++) {
                double sum = 0;
                for (Dimension column = 0; column < columnCount(); column++) {
                    sum += me(row, column) * other(column, otherColumn);
                }
                result(row, otherColumn) = sum;
            }
        }
        *this = result;
    }

    void Matrix::operator*=(const double other) {
        Array::operator*=(other);
    }

    Matrix Matrix::getAdjoint() const {
        Matrix result(rowCount(), columnCount());
        for (Dimension row = 0; row < rowCount(); row++) {
            for (Dimension column = 0; column < columnCount(); column++) {
                result(row, column) = getCofactor(row, column);
            }
        }
        return result;
    }

    Matrix Matrix::getAdjugate() const {
        return getAdjoint().transposed();
    }

    double Matrix::getCofactor(const Dimension row, const Dimension column) const {
        const auto determinant = getMinor(row, column).getDeterminant();
        return determinant * pow(-1, row + column);
    }

    double Matrix::getDeterminant() const {
        assert(isSquare());
        if (rowCount() == 1) {
            return me(0, 0);
        }
        if (rowCount() == 2) {
            return me(0, 0) * me(1, 1) - me(0, 1) * me(1, 0);
        }
        double result = 0;
        for (Dimension column = 0; column < columnCount(); column++) {
            const auto cofactor = getCofactor(0, column);
            result += me(0, column) * cofactor;
        }
        return result;
    }

    Matrix Matrix::getMinor(const Dimension row, const Dimension column) const {
        assert(row < rowCount() && column < columnCount() &&
            rowCount() > 1 && columnCount() > 1);
        Matrix result(rowCount() - 1, columnCount() - 1);

        for (Dimension subRow = 0; subRow < rowCount(); subRow++) {
            if (subRow == row) continue;
            for (Dimension subColumn = 0; subColumn < columnCount(); subColumn++) {
                if (subColumn == column) continue;
                const auto targetRow = subRow < row ? subRow : subRow - 1;
                const auto targetColumn = subColumn < column ? subColumn : subColumn - 1;
                result(targetRow, targetColumn) = me(subRow, subColumn);
            }
        }
        return result;
    }

    double Matrix::getTrace() const {
        assert(isSquare());
        double result = 0;
        for (Dimension row = 0; row < rowCount(); row++) {
            result += me(row, row);
        }
        return result;
    }

    Matrix Matrix::getIdentity(const Dimension size) {
        Matrix result(size, size);
        for (Dimension diagonalCell = 0; diagonalCell < size; diagonalCell++) {
            result(diagonalCell, diagonalCell) = 1;
        }
        return result;
    }

    Matrix Matrix::inverted() const {
        assert(isInvertible());
        return Matrix(getAdjugate() / getDeterminant());
    }

    bool Matrix::isInvertible() const {
        return isSquare() && abs(getDeterminant()) > Epsilon;
    }

    Matrix Matrix::normalized() const {
        double norm = 0;
        for (Dimension cell = 0; cell < size(); cell++) {
            norm += (*this)[cell] * (*this)[cell];
        }
        norm = sqrt(norm);
        if (norm < Epsilon) return *this;
        if ((*this)[0] < 0) norm = -norm;
        return Matrix(*this / norm);
    }

    Matrix Matrix::squared() const {
        return *this * *this;
    }

    Array Matrix::toArray() const {
        Array result(rowCount(), columnCount());

        for (Dimension row = 0; row < rowCount(); row++) {
            for (Dimension column = 0; column < columnCount(); column++) {
                result(row, column) = me(row, column);
            }
        }
        return result;
    }

    Matrix Matrix::transposed() const {
        Matrix result(columnCount(), rowCount());
        for (Dimension row = 0; row < rowCount(); row++) {
            for (Dimension column = 0; column < columnCount(); column++) {
                result(column, row) = me(row, column);
            }
        }
        return result;
    }

    Matrix operator+(Matrix left, const Matrix& right) {
        left += right;
        return left;
    }

    Matrix operator-(Matrix left, const Matrix& right) {
        left -= right;
        return left;
    }

    Matrix operator*(Matrix left, const Matrix& right) {
        left *= right;
        return left;
    }

    Matrix operator*(Matrix left, const double right) {
        left *= right;
        return left;
    }

    Matrix operator*(const double left, Matrix right) {
        right *= left;
        return right;
    }
}