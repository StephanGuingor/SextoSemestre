#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include <string>

class MatrixSolver {
    private:
    std::vector<std::vector<int> > m_matrixA;
    std::vector<std::vector<int> > m_matrixB;
    int m_elementary;

    void loadMatrix(std::ifstream& stream, std::vector<std::vector<int> >& mat);
    void printMatrix(const std::vector<std::vector<int> >&);
    std::vector<std::vector<int> > slice(const std::vector<std::vector<int> >, int qx, int qy);
    std::vector<std::vector<int> > recursiveStrassen(const std::vector<std::vector<int> > matA,const std::vector<std::vector<int> > matB);
    std::vector<std::vector<int> > textBook(const std::vector<std::vector<int> > matA,const std::vector<std::vector<int> > matB);

    std::vector<std::vector<int> > sum(const std::vector<std::vector<int> > matA,const  std::vector<std::vector<int> > matB);
    std::vector<std::vector<int> > sub(const std::vector<std::vector<int> > matA, const std::vector<std::vector<int> > matB);
    std::vector<std::vector<int> > merge(std::vector<std::vector<int> > matA, std::vector<std::vector<int> > matB, std::vector<std::vector<int> > matC, std::vector<std::vector<int> > matD);
    void matrixToCsv(const std::vector<std::vector<int> >& mat,const std::string& path);

    public:
    MatrixSolver();
    static MatrixSolver FromFiles(const std::string&, const std::string&);

    void load(const std::string&, const std::string&);
    void solve(std::string testName);
    
    // void load(std::string&, std::string&);
};