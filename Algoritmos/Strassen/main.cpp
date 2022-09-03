/**
 * @author Stephan Guingor
 * @date September 2 2022
*/

#include <iostream>

#include "MatrixSolver.h"

int main() {
    // Load matrices 16x16
    MatrixSolver solver = MatrixSolver::FromFiles("mat/01. Matrix_A_16_2_4.txt", "mat/02. Matrix_B_16_2_4.txt");
    solver.solve("Matrix_16_2_4_Test");

    solver.load("mat/03. Matrix_A_128_2_7.txt", "mat/04. Matrix_B_128_2_7.txt");
    solver.solve("Matrix_128_2_7_Test");

    solver.load("mat/05. Matrix_A_4096_2_12.txt", "mat/06. Matrix_B_4096_2_12.txt");
    solver.solve("Matrix_4096_2_12_Test");
}