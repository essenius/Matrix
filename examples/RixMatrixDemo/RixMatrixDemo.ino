#include <Matrix.h>

void printArray(const Array& array) {
    for (int row = 0; row < array.rowCount(); row++) {
        for (int column = 0; column < array.columnCount(); column++) {
            Serial.printf("%.0f ", array(row, column));
        }
        Serial.printf("\n");
    }
}


void setup() {
  Serial.begin(115200);
  delay(10);
  Matrix m({ {1, 2}, {3, 4} });
  printArray(m.squared().transposed());  
  // should return:
  // 7.000000 15.000000 
  // 10.000000 22.000000 
}

void loop() {
}