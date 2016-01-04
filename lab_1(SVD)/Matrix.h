//Implementation of matrix on custom class List

#pragma once
#include<iostream>

template <typename U>
class Matrix
{
public:
	Matrix(int rows_, int columns_);
	Matrix();
	//Use this constructor with explicit type casting
	//e. g. ((int*)a, 3, 3);
	Matrix(U* source, int rows_, int columns_);
	Matrix(U** source, int rows_, int columns_);
	Matrix(const Matrix&);
	~Matrix();
	void print();
	Matrix transpose();
	int getRows() const;
	int getColumns() const;
	Matrix operator+(Matrix);
	Matrix operator*(Matrix);
	Matrix operator*(U);
	U* &operator[] (int i);
	Matrix& operator =(Matrix);
	operator U** ();
	
private:
	U** matrix;
	int rows;
	int columns;
};

template <typename U>
void Matrix<U>::print()
{
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < columns; j++)
			std::cout << matrix[i][j] << " ";
		std:: cout << std::endl;
	}
	std::cout << std::endl;
}

template <typename U>
Matrix<U> Matrix<U>::transpose()
{
	Matrix res(columns, rows);
	for (int i = 0; i < rows; i++)
		for (int j = 0; j < columns; j++) 
			res[j][i] = matrix[i][j];
	return res;
}

template<typename U>
inline int Matrix<U>::getRows() const
{
	return rows;
}

template<typename U>
int Matrix<U>::getColumns() const
{
	return columns;
}

template <typename U>
Matrix<U>::Matrix(int rows_, int columns_) :
	rows(rows_), columns(columns_)
{
	matrix = new U*[rows];
	for (int i = 0; i < rows; i++)
		matrix[i] = new U[columns];
}

template<typename U>
Matrix<U>::Matrix() :rows(0), columns(0){}

template <typename U>
Matrix<U>::Matrix(U* source, int rows_, int columns_)
{
	this->rows = rows_;
	this->columns = columns_;
	matrix = new U*[rows];
	for (int i = 0; i < rows; i++) {
		matrix[i] = new U[columns];
		for (int j = 0; j < columns; j++)
			matrix[i][j] = (*(source + i*columns + j));
	}
}

template<typename U>
Matrix<U>::Matrix(U ** source, int rows_, int columns_)
{
	this->rows = rows_;
	this->columns = columns_;
	matrix = new U*[rows];
	for (int i = 0; i < rows; i++) {
		matrix[i] = new U[columns];
		for (int j = 0; j < columns; j++)
			matrix[i][j] = source[i][j];
	}
}

template<typename U>
Matrix<U>::Matrix(const Matrix &another)
{
	rows = another.rows;
	columns = another.columns;
	matrix = new U*[rows];
	for (int i = 0; i < rows; i++) {
		matrix[i] = new U[columns];
		for (int j = 0; j < columns; j++)
			matrix[i][j] = another.matrix[i][j];
	}
}

template<typename U>
Matrix<U>::~Matrix()
{
	for (int i = 0; i < rows; i++)
		delete[] matrix[i];
	delete[] matrix;
}

template <typename U>
Matrix<U> Matrix<U>::operator + (Matrix<U> A)
{
	if (A.rows == rows && A.columns == columns)
	{
		Matrix<U> C = A;
		for (int i = 0; i < A.rows; i++)
			for (int j = 0; j < A.columns; j++)
				C[i][j] = C[i][j] + matrix[i][j];
		return C;
	}
	else std::cerr << "Error: adding matrices with different dimensions";
}

template <typename U>
Matrix<U> Matrix<U>::operator * (Matrix<U> A)
{
	if (columns == A.rows)
	{
		Matrix<U> C(rows, A.columns);
		U data = 0.0;
		for (int i = 0; i < rows; i++)
		{
			for (int j = 0; j < A.columns; j++)
			{
				for (int r = 0; r < columns; r++)
					data += matrix[i][r] * A[r][j];
				C[i][j] = data;
				data = 0;
			}
		}
		return C;
	}
	else std::cerr << "Error: multipling matrices whith inappropriate dimensions";
}

template <typename U>
Matrix<U> Matrix<U>::operator * (U alpha)
{
	Matrix<U> C = *this;
	for (int i = 0; i < rows; i++)
		for (int j = 0; j < columns; j++)
			C[i][j] *= alpha;
	return C;
}

template <typename U>
U* & Matrix<U>::operator[] (int i) {
	return matrix[i];
}

template<typename U>
Matrix<U> & Matrix<U>::operator=(Matrix<U> another)
{
	rows = another.rows;
	columns = another.columns;
	matrix = new U*[rows];
	for (int i = 0; i < rows; i++) {
		matrix[i] = new U[columns];
		for (int j = 0; j < columns; j++)
			matrix[i][j] = another[i][j];
	}
	return *this;
}

template<typename U>
inline Matrix<U>::operator U**()
{
	return matrix;
}
