#include <SolverMatrix.h>

void printArray(const char* caption, const Array& array) {
    Serial.printf("%s:\n", caption);
    for (int row = 0; row < array.rowCount(); row++) {
        for (int column = 0; column < array.columnCount(); column++) {
            Serial.printf("%.3f ", array(row, column));
        }
        Serial.printf("\n");
    }
}

void setup() {
  Serial.begin(115200);
  delay(100);
  Serial.println("Array operations (cell by cell)");
  Array a({ {1, 2}, {3, 4} } );
  printArray("Array a Squared", a.pow2());
  // 1.000 4.000
  // 9.000 16.000
  printArray("Array a * 3", a * 3);
  // 3.000 6.000
  // 9.000 12.000
  printArray("Array a / 2", a / 2);
  // 0.500 1.000
  // 1.500 2.000
  Array b{{0, 2}, {3, 0}};
  printArray("Array a + Array b", a + b);
  // 1.000 4.000
  // 6.000 4.000
  printArray("Array a - Array b", a - b);
  // 1.000 0.000
  // 0.000 4.000  
  printArray("Array a * Array b", a * b);
  // 0.000 4.000
  // 9.000 0.000 
  Matrix c(b.getColumn(1));
  printArray("Matrix c", c);
  // 2.000
  // 0.000 
  b.setColumn(1, {{7}, {11}});
  printArray("New Array b", b);
  // 0.000 7.000
  // 2.000 11.000  
  b.setRow(1, c.transposed());
  printArray("New Array b", b);
  // 0.000 7.000
  // 2.000 0.000  
  
  Matrix m({ {1, 2}, {3, 4} });
  printArray("Matrix m Squared Transposed", m.squared().transposed());  
  // Should return:
  // Squared Transposed:
  // 7.000 15.000 
  // 10.000 22.000 

  const SolverMatrix solver({ {-2, -4, 2}, {-2, 1, 2}, {4, 2, 5} });
  printArray("Eigenvalues",solver.getEigenvalues());
  // Should return:
  // Eigenvalues:
  // 6.000
  // -5.000
  // 3.000
  
  printArray("Eigenvectors", solver.getEigenvectors());
  // should return:
  // Eigenvectors:
  // 0.058 0.816 0.535 
  // 0.351 0.408 -0.802 
  // 0.935 -0.408 -0.267 
}

void loop() {
}