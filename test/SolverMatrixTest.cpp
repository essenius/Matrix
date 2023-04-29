#include <gtest/gtest.h>
#include "MatrixTest.h"
#include "../src/SolverMatrix.h"

class SolverMatrixTest : public MatrixTest {
protected:
    bool normalizedContains(const Matrix& expected, const Matrix& actual) const {
        auto expectedNormalized = expected.normalize();
        for (Dimension column = 0; column < actual.columns(); column++) {
            // actual should already be normalized
            Array vector = actual.getColumn(column);
            if (isEqual(expectedNormalized, vector, SolverMatrix::EIGEN_EPSILON)) return true;
        }
        return false;
    }
};

// *** Row echelon format ***

TEST_F(SolverMatrixTest, toRowEchelonForm10x) {
    SolverMatrix m({ {10, 20, 30, 0}, {30, 40, 50, 0}, {40, 50, 60, 0} });
    const Matrix permutation = m.toReducedRowEchelonFormWithPivot();
    const Array expected({ {1, 0, 0.5, 0}, {0, 1, 0.5, 0}, {0, 0, 0, 0} });
    expectEqual(expected, m);
    expectEqual(Array({{0,1,0,0}, {0,0,1,0}, {1,0,0,0}, {0,0,0,1}}), permutation);
}


TEST_F(SolverMatrixTest, toRowEchelonForm1) {
    SolverMatrix m({ {1, 2, 3, 0}, {3, 4, 5, 0}, {4, 5, 6, 0} });
    const Matrix permutation = m.toReducedRowEchelonFormWithPivot();
    const Array expected({ {1, 0, 0.5, 0}, {0, 1, 0.5, 0}, {0, 0, 0, 0} });
    expectEqual(expected, m);
    expectEqual(Array({{0,1,0,0}, {0,0,1,0}, {1,0,0,0}, {0,0,0,1}}), permutation);
    m.toReducedRowEchelonFormWithPivot();
    expectEqual(expected, m);
}

TEST_F(SolverMatrixTest, toRowEchelonForm2) {
    SolverMatrix m({ {-5, -4, 2, 0}, {-2, -2, 2, 0}, {4, 2, 2, 0} });
    const Matrix permutation = m.toReducedRowEchelonFormWithPivot();
    const Array expected({ {1, 0, 2.0/3.0, 0}, {0, 1, -1.0/3.0, 0}, {0, 0, 0, 0} });
    expectEqual(expected, m);
    expectEqual(Array({{1,0,0,0}, {0,0,1,0}, {0,1,0,0}, {0,0,0,1}}), permutation);
}

TEST_F(SolverMatrixTest, toRowEchelonForm3) {
    SolverMatrix m({ {0, 0, 0}, {1, 0, 1}, {0, 1, -1} });
    const Matrix permutation = m.toReducedRowEchelonFormWithPivot();
    const Array expected({ {1, 0, 1}, {0, 1, -1}, {0, 0, 0} });
    expectEqual(expected, m);
    expectEqual(Array({{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}), permutation);
}

TEST_F(SolverMatrixTest, toRowEchelonForm4) {
    SolverMatrix m({ {1, 0, 0}, {0, 2, 1}, {0, 0, 2} });
    const Matrix permutation = m.toReducedRowEchelonFormWithPivot();
    const Array expected({ {1, 0, 0}, {0, 1, 0}, {0, 0, 1} });
    expectEqual(expected, m);
    expectEqual(Array({ {0, 0, 1}, {1, 0, 0}, {0, 1, 0} }), permutation);
}

TEST_F(SolverMatrixTest, toRowEchelonFormSpecial) {
    SolverMatrix m({ {8192.10277, -0.0341611, 8192.03422}, {0.0686703, -6.3222615e-7, 0.0686707}, {8192.03422, -0.03416085, 8192.10277} });
    const Matrix permutation = m.toReducedRowEchelonFormWithPivot();
    const Array expected({ {1, 0, -3.9e-6}, {0, 1, -2.62e-7}, {0, 9.74e-7, -3.46e-7} });
    expectEqual(expected, m, "toRowEchelonFormSpecial", 1e-7);
    expectEqual(Array({ {1, 0, 0}, {0, 0, 1}, {0, 1, 0} }), permutation);
}

// *** Null spaces and free variables ***

TEST_F(SolverMatrixTest, oneFreeVariableColumn) {
    SolverMatrix m({{0,1,1}, {0,0,1}, {0,0,0}});
    const auto actual = m.getFreeVariables();
    EXPECT_EQ(1, actual.size());
    EXPECT_EQ(0, actual[0]);
}

TEST_F(SolverMatrixTest, twoFreeVariableColumns) {
    SolverMatrix m({{0,1,0}, {0,0,0}, {0,0,0}});
    const auto actual = m.getFreeVariables();
    EXPECT_EQ(2, actual.size());
    EXPECT_EQ(0, actual[0]);
    EXPECT_EQ(2, actual[1]);
}

TEST_F(SolverMatrixTest, noFreeVariableColumns) {
    SolverMatrix m({{1,0,0}, {0,1,0}, {0,0,1}});
    const auto actual = m.getFreeVariables();
    EXPECT_EQ(0, actual.size());
}

TEST_F(SolverMatrixTest, nullSpaceNoFreeVariables) {
    SolverMatrix m({ {1, 0, 0}, {0, 1, 0}, {0, 0, 1} });
    const auto actual = m.getNullSpace();
    EXPECT_EQ(0, actual.columns());
    EXPECT_EQ(3, actual.rows());
}

TEST_F(SolverMatrixTest, nullSpaceForOneFreeVariable) {
    SolverMatrix m({ {1, 0, 0}, {0, 1, 0}, {0, 0, 0} });
    const auto actual = m.getNullSpace();
    EXPECT_EQ(1, actual.columns());
    expectNormalizedEqual(Matrix({ { 0 }, { 0 }, { 1 } }), actual);
}

TEST_F(SolverMatrixTest, nullSpaceForTwoFreeVariable) {
    SolverMatrix m({ { 0, 1, 0}, {0, 0, 0}, {0, 0, 0} });
    const auto actual = m.getNullSpace();
    
    EXPECT_EQ(2, actual.columns());
    const auto expected1 = Matrix({ { 1 }, { 0 }, { 0 } });
    const auto expected2 = Matrix({ { 0 }, { 0 }, { 1 } });
    expectNormalizedEqual(expected1, Matrix(actual.getColumn(0)), "expected1");
    expectNormalizedEqual(expected2, Matrix(actual.getColumn(1)), "expected2"); 
}

TEST_F(SolverMatrixTest, nullSpaceSpecial) {
    SolverMatrix m({{1, -3.9e-6, 0}, {0, -2.62e-7, 1}, {0, -3.46e-7, 9.74e-7}});
    const auto actual = m.getNullSpace();
    
    EXPECT_EQ(1, actual.columns());
    auto outcome = m * actual;
    const auto expected = Matrix({ { 3.9e-6 }, { 1 }, { 3.46e-7 } });
    expectNormalizedEqual(expected, actual, "null space", SolverMatrix::EIGEN_EPSILON);
    expectEqual(Array({{0},{0},{0}}), outcome, "matrix * null space", SolverMatrix::EIGEN_EPSILON);
} 

// *** Eigenvalues, eigen vectors ***

TEST_F(SolverMatrixTest, getEigenvalues1d) {
    SolverMatrix m({ {6} });
    const auto actual = m.getEigenvalues();
    expectEqual(m, actual);
}

TEST_F(SolverMatrixTest, getEigenvalues2d) {
    SolverMatrix m({ {6, -1}, {2, 3} });
    const Matrix expected({ {5}, {4} });
    const auto actual = m.getEigenvalues();
    expectEqual(expected, actual);
}

TEST_F(SolverMatrixTest, getEigenvalues2dNoSolution) {
    SolverMatrix m({ {0, 1}, {-1, 0} });
    const Matrix expected(0, 0);
    const auto actual = m.getEigenvalues();
    expectEqual(expected, actual);
}

TEST_F(SolverMatrixTest, getEigenvalues3dDup) {
    SolverMatrix m({ {2, 0, 0}, {1, 2, 1}, {-1, 0, 1} });
    const auto actual = m.getEigenvalues();
    EXPECT_EQ(2, actual.rows());
    EXPECT_EQ(1, actual.columns());
    EXPECT_TRUE(contains(actual, 2));
    EXPECT_TRUE(contains(actual, 1));
}

TEST_F(SolverMatrixTest, getEigenvalues3dOne) {
    // trace(M) = 0, trace(M^2) = 0, determinant = 0
    SolverMatrix m({ {0, 0, 1}, {0, 0, -1}, {1, 1, 0} });
    const auto actual = m.getEigenvalues();
    EXPECT_EQ(1, actual.rows());
    EXPECT_EQ(1, actual.columns());
    EXPECT_TRUE(contains(actual, 0));
}

TEST_F(SolverMatrixTest, getEigenvalues3dComplex) {
    SolverMatrix m({ {33, -23, 9}, {22, 33, -23}, {19, 14, 50} });
    m /= 29;
    const Array expected({ {2} });
    const auto actual = m.getEigenvalues();
    expectEqual(expected, actual);
}

TEST_F(SolverMatrixTest, getEigenvectorsFor3dReal) {
    SolverMatrix m({ {-2, -4, 2}, {-2, 1, 2}, {4, 2, 5} });
    const auto actual = m.getEigenvalues();
    EXPECT_EQ(3, actual.rows());
    EXPECT_EQ(1, actual.columns());
    EXPECT_TRUE(contains(actual, 3.0));
    EXPECT_TRUE(contains(actual, -5.0));
    EXPECT_TRUE(contains(actual, 6.0));

    const auto actualVector1 = m.getEigenvectorFor(3.0);
    expectNormalizedEqual(Matrix({ { -2 }, { 3 }, { 1 } }), actualVector1, "eigenvector 1");
    const auto actualVector2 = m.getEigenvectorFor(-5.0);
    expectNormalizedEqual(Matrix({ { -2 }, { -1 }, { 1 } }), actualVector2, "eigenvector 2");
    const auto actualVector3 = m.getEigenvectorFor(6.0);
    expectNormalizedEqual(Matrix({ { 1 }, { 6 }, { 16 } }), actualVector3, "eigenvector 3");
}

TEST_F(SolverMatrixTest, eigenvectorsForTwoFreeVariables) {
    SolverMatrix m({ {1, 0, 0}, {0, 0, 0}, {0, 0, 1} });
    const auto actual = m.getEigenvalues();
    EXPECT_EQ(2, actual.rows());
    EXPECT_EQ(1, actual.columns());
    EXPECT_TRUE(contains(actual, 1.0));
    EXPECT_TRUE(contains(actual, 0.0));

    const auto actualVectors1 = m.getEigenvectorFor(1.0);
    EXPECT_EQ(2, actualVectors1.columns());
    const auto expected1 = Matrix({ { 1 }, { 0 }, { 0 } });
    const auto expected2 = Matrix({ { 0 }, { 0 }, { 1 } });
    expectNormalizedEqual(expected1, Matrix(actualVectors1.getColumn(0)), "actualVector1-1");
    expectNormalizedEqual(expected2, Matrix(actualVectors1.getColumn(1)), "actualVector1-2"); 

    const auto actualVector2 = m.getEigenvectorFor(0.0);
    expectNormalizedEqual(Matrix({ { 0 }, { 1 }, { 0 } }), actualVector2, "actualVector2");
}

TEST_F(SolverMatrixTest, getEigenVectorsTwoFreeVariables) {
    SolverMatrix m({ {1, 0, 0}, {0, 0, 0}, {0, 0, 1} });
    const auto actual = m.getEigenvectors();
    EXPECT_EQ(3, actual.rows());
    EXPECT_EQ(3, actual.columns());

    const auto expected = Matrix({ {0, 1, 0}, {1, 0, 0}, {0, 0, 1} }).transpose();
    for (int i = 0; i < 3; ++i) {
        expectNormalizedEqual(Matrix(expected.getColumn(i)), Matrix(actual.getColumn(i)), "column " + std::to_string(i));
    }
}

TEST_F(SolverMatrixTest, getEigenvectors3dReal) {
    SolverMatrix m({ {-2, -4, 2}, {-2, 1, 2}, {4, 2, 5} });
    const auto actual = m.getEigenvectors();
    EXPECT_EQ(3, actual.rows());
    EXPECT_EQ(3, actual.columns());

    EXPECT_TRUE(normalizedContains(Matrix({ { -2 }, { 3 }, { 1 } }), actual)) << "eigenvector 1";
    EXPECT_TRUE(normalizedContains(Matrix({ { -2 }, { -1 }, { 1 } }), actual)) << "eigenvector 2";
    EXPECT_TRUE(normalizedContains(Matrix({ { 1 }, { 6 }, { 16 } }), actual)) << "eigenvector 3";
}