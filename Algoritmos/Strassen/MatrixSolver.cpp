#include "MatrixSolver.h"
#include <fstream>
#include <chrono>
#include <cmath>
using namespace std::chrono;
 
// After function call
auto stop = high_resolution_clock::now();

// Constructors

// Default Initializer
MatrixSolver::MatrixSolver(): m_matrixA(), m_matrixB(), m_elementary(0) {

}

MatrixSolver MatrixSolver::FromFiles(const std::string& pathA, const std::string& pathB) {
    MatrixSolver solver = MatrixSolver();
    solver.load(pathA, pathB);
    
 
    return solver;
}

std::vector<std::vector<int> > MatrixSolver::slice(const std::vector<std::vector<int> > mat, int qx, int qy) {
    int qSize = std::floor(mat.size() / 2);
    std::vector<std::vector<int> > matC;
    matC.reserve(qSize);

    for (int itX = qSize*qx; itX <  qSize*qx + qSize; itX++) {  
        matC.emplace_back(mat[itX].begin() + qSize*qy, mat[itX].begin() + qSize*qy + qSize);
    }

    return matC;
}

std::vector<std::vector<int> > MatrixSolver::sum(const std::vector<std::vector<int> > matA,const std::vector<std::vector<int> > matB) {
    int qSize = matA.size();
    std::vector<std::vector<int> > matC;
    matC.reserve(qSize);

    for (int i = 0; i < qSize; i++) {
        std::vector<int> row;
        row.reserve(qSize);
        for (int j = 0; j < matA.size(); j++) {
            row.push_back(matA[i][j] + matB[i][j]);
        }

        matC.push_back(row);
    }

    return matC;
}

std::vector<std::vector<int> > MatrixSolver::sub(const std::vector<std::vector<int> > matA,const std::vector<std::vector<int> > matB) {
    int qSize = matA.size();
    std::vector<std::vector<int> > matC;
    matC.reserve(qSize);


    for (int i = 0; i < qSize; i++) {
        std::vector<int> row;
        row.reserve(qSize);
        for (int j = 0; j < qSize; j++) {
            row.push_back(matA[i][j] - matB[i][j]);
        }

        matC.push_back(row);
    }

    return matC;
}

std::vector<std::vector<int> > MatrixSolver::merge(std::vector<std::vector<int> > mA,std::vector<std::vector<int> > mB,std::vector<std::vector<int> > mC,  std::vector<std::vector<int> > mD) {
   int qSize = mA.size();

    for (int i = 0; i < qSize; i++ ) {
        for (int j = 0; j < qSize; j++) {
           mA[i].push_back(mB[i][j]);
           mC[i].push_back(mD[i][j]);
        }
        mA.push_back(mC[i]);
    }

    return mA;
}

// recursive strassen
std::vector<std::vector<int> > MatrixSolver::recursiveStrassen(const std::vector<std::vector<int> > matA,const std::vector<std::vector<int> > matB) {
    // base case
      if (matA.size() == 1){
        std::vector<std::vector<int> > c{{ matA[0][0] * matB[0][0] }};
        return c;
      }
    // recursive steps
   auto A00 = slice(matA, 0, 0);
   auto A10 = slice(matA, 1, 0);
   auto A01 = slice(matA, 0, 1);
   auto A11 = slice(matA, 1, 1);

   auto B00 = slice(matB, 0, 0);
   auto B10 = slice(matB, 1, 0);
   auto B01 = slice(matB, 0, 1);
   auto B11 = slice(matB, 1, 1);

    auto P1  = recursiveStrassen(sum(A00, A11), sum(B00, B11));
    auto P2  = recursiveStrassen(sum(A10,A11),B00);
    auto P3  = recursiveStrassen(A00, sub(B01, B11));
    auto P4  = recursiveStrassen(A11, sub(B10, B00));
    auto P5  = recursiveStrassen(sum(A00, A01), B11);
    auto P6  = recursiveStrassen(sub(A10, A00), sum(B00, B01));
    auto P7  = recursiveStrassen(sub(A01, A11), sum(B10, B11));

    auto C11 = sum(sub(sum(P1, P4), P5), P7);
    auto C12 = sum(P3,P5);
    auto C21 = sum(P2,P4);
    auto C22 = sum(sub(P1, P2), sum(P3, P6));
    return merge(C11, C12, C21,C22);
}

 std::vector<std::vector<int> > MatrixSolver::textBook(const std::vector<std::vector<int> > matA,const std::vector<std::vector<int> > matB) {
    std::vector<std::vector<int> > cMat;
    for (int r1 = 0; r1 < matA.size(); r1++){
        std::vector<int> row;
        for (int c2 = 0; c2 < matA.size(); c2++){
            int res=0;
            // dot product
            for (int c1 = 0; c1 < matA.size(); c1++){
                res += matA[r1][c1]*matB[c1][c2];
                m_elementary += 1;
            }
            row.push_back(res);
        }
        cMat.push_back(row);
    }
    return cMat;
 }

// helpers
void MatrixSolver::loadMatrix(std::ifstream& stream, std::vector<std::vector<int> >& outMat) {
         if (!stream.is_open()) return;

        int lineIdx = 0;
        std::string line;
        while( std::getline(stream, line, '\n') ) {
            outMat.emplace_back(std::vector<int>{});
            std::stringstream ss(line);
            std::string stringNumber;

            while (std::getline( ss, stringNumber, ',') ) {
                int number;

                if (!(std::stringstream(stringNumber) >> number)) {
                    std::cout << "loadMatrix::error " + stringNumber + " is not a int" << std::endl;
                    number = -1;
                }

                outMat[lineIdx].push_back(number);
            }

            lineIdx++;
        }
}

void MatrixSolver::printMatrix(const std::vector<std::vector<int> >& mat) {
    std::cout << "Matrix" << std::endl;
    for (auto vec : mat) {
        for (int num : vec) {
            std::cout << "\t" << num << " ";
        }
        std::cout << std::endl;
    }
}

void MatrixSolver::load(const std::string& pathA, const std::string& pathB) {

    m_matrixA = std::vector<std::vector<int> >();

    std::ifstream matrixA(pathA);
    loadMatrix(matrixA, m_matrixA);
    matrixA.close();

    m_matrixB = std::vector<std::vector<int> >();
    std::ifstream matrixB(pathB);
    loadMatrix(matrixB, m_matrixB);
    matrixB.close();
}

void MatrixSolver::matrixToCsv(const std::vector<std::vector<int> >& mat, const std::string& path) {
    std::ofstream f(path);

    for (int i = 0; i < mat.size(); i++){
        for (int j = 0; j < mat[i].size(); j++){
            f << mat[i][j] << ",";
        }
        f << "\n";
    }
    f.close();
    std::cout<<"File created at " << path << std::endl;
}

void MatrixSolver::solve(std::string testName) {
    // reset steps
    m_elementary = 0;
     // Get starting timepoint
    auto start = high_resolution_clock::now();
    auto res = recursiveStrassen(m_matrixA, m_matrixB);
    // Get ending timepoint
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
 
    std::cout << "--Strassen Ran--" << "\n" << "Number of multiplications: " << m_elementary << std::endl;
    std::cout << "Time taken by function: " << duration.count() << " microseconds" << std::endl;

    matrixToCsv(res, testName + "_strassen" + ".csv");

    // reset steps
    m_elementary = 0;
     // Get starting timepoint
    start = high_resolution_clock::now();
    res = textBook(m_matrixA, m_matrixB);
    // Get ending timepoint
    stop = high_resolution_clock::now();
    duration = duration_cast<microseconds>(stop - start);
 
   
    std::cout << "--Textbook Ran--" << "\n" << "Number of multiplications: " << m_elementary << std::endl;
    std::cout << "Time taken by function: " << duration.count() << " microseconds" << std::endl;

    matrixToCsv(res, testName + "_textbook" + ".csv");
}
