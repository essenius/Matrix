#include <gtest/gtest.h>
#include "BasicMatrixTest.h"

#include <gtest/gtest.h>
#include "Matrix.h"

class MatrixTest : public BasicMatrixTest {};

// *** Row echelon format ***

TEST_F(MatrixTest, toRowEchelonForm1) {
    Matrix m({ {1, 2, 3, 0}, {3, 4, 5, 0}, {4, 5, 6, 0} });
    m.toRowEchelonForm();
    const Matrix expected({ {1, 0, -1, 0}, {0, 1, 2, 0}, {0, 0, 0, 0} });
    expectEqual(expected, m);
}

TEST_F(MatrixTest, toRowEchelonForm2) {
    Matrix m({ {-5, -4, 2, 0}, {-2, -2, 2, 0}, {4, 2, 2, 0} });
    m.toRowEchelonForm();
    const Matrix expected({ {1, 0, 2, 0}, {0, 1, -3, 0}, {0, 0, 0, 0} });
    expectEqual(expected, m);
    m.toRowEchelonForm();
    expectEqual(expected, m);
}

TEST_F(MatrixTest, toRowEchelonForm3) {
    Matrix m({ {0, 0, 0}, {1, 0, 1}, {0, 1, -1} });
    m.toRowEchelonForm();
    const Matrix expected({ {1, 0, 1}, {0, 1, -1}, {0, 0, 0} });
    expectEqual(expected, m);
}

TEST_F(MatrixTest, toRowEchelonForm4) {
    Matrix m({ {1, 0, 0}, {0, 2, 1}, {0, 0, 2} });
    m.toRowEchelonForm();
    const Matrix expected({ {1, 0, 0}, {0, 1, 0}, {0, 0, 1} });
    expectEqual(expected, m);
}

// *** Null spaces and free variables ***

TEST_F(MatrixTest, oneFreeVariableColumn) {
    Matrix m({{0,1,1}, {0,0,1}, {0,0,0}});
    const auto actual = m.getFreeVariables();
    EXPECT_EQ(1, actual.size());
    EXPECT_EQ(0, actual[0]);
}

TEST_F(MatrixTest, twoFreeVariableColumns) {
    Matrix m({{0,1,0}, {0,0,0}, {0,0,0}});
    const auto actual = m.getFreeVariables();
    EXPECT_EQ(2, actual.size());
    EXPECT_EQ(0, actual[0]);
    EXPECT_EQ(2, actual[1]);
}

TEST_F(MatrixTest, noFreeVariableColumns) {
    Matrix m({{1,0,0}, {0,1,0}, {0,0,1}});
    const auto actual = m.getFreeVariables();
    EXPECT_EQ(0, actual.size());
}

TEST_F(MatrixTest, nullSpaceForOneFreeVariable) {
    Matrix m({ {1, 0, 0}, {0, 1, 0}, {0, 0, 0} });
    const auto actual = m.getNullSpace();
    EXPECT_EQ(1, actual.columns());
    expectNormalizedEqual(Matrix({ { 0 }, { 0 }, { 1 } }), actual);
}

TEST_F(MatrixTest, nullSpaceForTwoFreeVariable) {
    Matrix m({ { 0, 1, 0}, {0, 0, 0}, {0, 0, 0} });
    const auto actual = m.getNullSpace();
    
    EXPECT_EQ(2, actual.columns());
    const auto expected1 = Matrix({ { 1 }, { 0 }, { 0 } });
    const auto expected2 = Matrix({ { 0 }, { 0 }, { 1 } });
    expectNormalizedEqual(expected1, actual.getColumn(0), "expected1");
    expectNormalizedEqual(expected2, actual.getColumn(1), "expected2"); 
}

// *** Eigenvalues, eigen vectors ***

TEST_F(MatrixTest, getEigenvalues2d) {
    Matrix m({ {6, -1}, {2, 3} });
    const Matrix expected({ {5}, {4} });
    const auto actual = m.getEigenvalues();
    expectEqual(expected, actual);
}

TEST_F(MatrixTest, getEigenvalues3dDup) {
    Matrix m({ {2, 0, 0}, {1, 2, 1}, {-1, 0, 1} });
    const Matrix expected({ {1}, {2}, {2} });
    const auto actual = m.getEigenvalues();
    EXPECT_EQ(2, actual.rows());
    EXPECT_EQ(1, actual.columns());
    EXPECT_TRUE(contains(actual, 2));
    EXPECT_TRUE(contains(actual, 1));
}

TEST_F(MatrixTest, getEigenvalues3dComplex) {
    Matrix m({ {33, -23, 9}, {22, 33, -23}, {19, 14, 50} });
    m /= 29;
    const Matrix expected({ {2} });
    const auto actual = m.getEigenvalues();
    expectEqual(expected, actual);
}

TEST_F(MatrixTest, getEigenvectorsFor3dReal) {
    Matrix m({ {-2, -4, 2}, {-2, 1, 2}, {4, 2, 5} });
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

TEST_F(MatrixTest, eigenvectorsForTwoFreeVariables) {
    Matrix m({ {1, 0, 0}, {0, 0, 0}, {0, 0, 1} });
    const auto actual = m.getEigenvalues();
    EXPECT_EQ(2, actual.rows());
    EXPECT_EQ(1, actual.columns());
    EXPECT_TRUE(contains(actual, 1.0));
    EXPECT_TRUE(contains(actual, 0.0));

    const auto actualVectors1 = m.getEigenvectorFor(1.0);
    EXPECT_EQ(2, actualVectors1.columns());
    const auto expected1 = Matrix({ { 1 }, { 0 }, { 0 } });
    const auto expected2 = Matrix({ { 0 }, { 0 }, { 1 } });
    expectNormalizedEqual(expected1, actualVectors1.getColumn(0), "actualVector1-1");
    expectNormalizedEqual(expected2, actualVectors1.getColumn(1), "actualVector1-2"); 

    const auto actualVector2 = m.getEigenvectorFor(0.0);
    expectNormalizedEqual(Matrix({ { 0 }, { 1 }, { 0 } }), actualVector2, "actualVector2");
}

TEST_F(MatrixTest, getEigenVectorsTwoFreeVariables) {
    Matrix m({ {1, 0, 0}, {0, 0, 0}, {0, 0, 1} });
    const auto actual = m.getEigenvectors();
    EXPECT_EQ(3, actual.rows());
    EXPECT_EQ(3, actual.columns());

}