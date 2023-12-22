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

#ifndef ARRAY_H
#define ARRAY_H

#include <vector>

namespace RixMatrix {
    using Dimension = unsigned int;

    /// Class for array manipulations (coefficient wise)
    class Array {
    public:
        Array(Dimension rows, Dimension columns);
        explicit Array(std::initializer_list<std::initializer_list<double>> list);

        Array(Array&& other) noexcept;
        Array(const Array& other);
        Array& operator=(const Array& other) = default;

        double& operator[](Dimension cell);
        const double& operator[](Dimension cell) const;

        double& operator()(Dimension row, Dimension column);
        const double& operator()(Dimension row, Dimension column) const;

        void operator+=(const Array& other);
        void operator+=(double other);

        void operator-=(const Array& other);
        void operator-=(double other);

        void operator*=(const Array& other);
        void operator*=(double other);
        void operator/=(double other);
        bool operator==(const Array& other) const;

        Dimension columnCount() const;
        Array getColumn(Dimension column) const;
        Array getRow(Dimension row) const;

        bool isSquare() const;
        double me(Dimension row, Dimension column) const;
        Array pow2() const;
        Dimension rowCount() const;

        void setColumn(Dimension column, const Array& input);
        void setColumn(Dimension column, double value);
        void setColumnCount(Dimension columns);

        void setRow(Dimension row, const Array& input);
        void setRow(Dimension row, double value);
        void setRowCount(Dimension rows);

        void swapRows(Dimension row1, Dimension row2);
        void swapColumns(Dimension column1, Dimension column2);

        Dimension size() const;
        bool sizeIsEqual(const Array& other) const;

        friend Array operator+(Array left, const Array& right);
        friend Array operator-(Array left, const Array& right);
        friend Array operator*(Array left, const Array& right);
        friend Array operator*(Array left, double right);
        friend Array operator*(double left, Array right);
        friend Array operator/(Array left, double right);

        // not using std::numeric_limits<double>::epsilon() because it is too small
        static constexpr double Epsilon = 1e-12;

    private:
        std::vector<double> _data;

        Dimension _rows;
        Dimension _columns;
        Dimension _arraySize;
    };
}
#endif