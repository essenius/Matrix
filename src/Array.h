#ifndef ARRAY_H
#define ARRAY_H

#include <vector>

using Dimension = unsigned int;

/// Class for array manipulations (coefficient wise)
class Array {
public:
    Array(Dimension rows, Dimension columns);
    explicit Array(const std::initializer_list<std::initializer_list<double>> list);

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
    bool operator==(const Array &other) const;

    Dimension columnCount() const;
    Array getColumn(const Dimension column) const;
    Array getRow(const Dimension row) const;

    bool isSquare() const;
    double me(Dimension row, Dimension column) const;
    Array pow2();
    Dimension rowCount() const;

    void setColumn(Dimension column, const Array &input);
    void setColumn(const Dimension column, const double value);
    void setColumnCount(const Dimension columns);

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
    static constexpr double EPSILON = 1e-12; 

private:
    std::vector<double> _data;

    Dimension _rows;
    Dimension _columns;
    Dimension _arraySize;
};
#endif