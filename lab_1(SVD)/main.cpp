#include <iostream>
#include "Matrix.h"
#include "eigen.h"
using namespace std;

void main()
{
	//Input block
	int inputRows, inputColumns; 
	cout << "Enter count of rows and columns in your system: ";
	cin >> inputRows >> inputColumns;
	float** input = new float*[inputRows];
	cout << "Enter matrix" << endl;
	for (int i = 0; i < inputRows; i++) {
		input[i] = new float[inputColumns];
		for (int j = 0; j < inputColumns; j++)
			cin >> input[i][j];
	} 
	cout  << "Enter vector of right parts: " << endl;
	float* right = new float[inputRows];
	for (int i = 0; i < inputRows; i++)
		cin >> right[i];
	/********************************************/

	Matrix<float> B(right, inputRows, 1);
	Matrix <float> M((float**)input, inputRows, inputColumns);
	Matrix<float> D = M;
	Matrix<float> U, V;
	SVD(D, U, V);

	cout << endl << "Singular value decomposition matrices:" << endl;
	cout << "U = " << endl; U.print();
	cout << "V = " << endl; V.print();
	cout << "Sigma (D) matrix = " << endl; D.print();

	//Matrix for finding a solution of initial system
	Matrix<float> diagD(M.getRows(), D.getColumns());
	for (int i = 0; i < D.getRows(); i++)
		for (int j = 0; j < D.getColumns(); j++) {
			if (i == j) {
				diagD[i][j] = 1 / D[i][j];
			}
			else diagD[i][j] = 0;
		}

	//solution
	if (inputRows == inputColumns) {
		cout << "A solution of initial system: " << endl;
		Matrix<float> X = V * diagD * U.transpose() * B;
		X.print();
	}
	system("pause");
}   