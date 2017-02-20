//implement matrix class LU

#ifndef MatrixLU_H
#define MatrixLU_H

#include<Eigen/Dense>
#include<tuple>
#include<iostream>
#include<math.h>
#include<algorithm>

//switch row i and j
using Matrix = Eigen::MatrixXd;

void SwitchRow(Matrix& m, int i, int j) {
	if (i >= m.rows() || j >= m.rows()) { throw "Row overflow"; }
	Eigen::RowVectorXd rtemp = m.row(i);
	m.row(i) = m.row(j);
	m.row(j) = rtemp;
}


//conduct LU decomposition without pivoting
std::tuple<Matrix, Matrix> lu_no_pivoting(Matrix m) {
		Matrix L(m.rows(),m.cols()), U(m.rows(),m.cols());
		L.setZero(); U.setZero();
		//check if it is square matrix
		if (m.rows() != m.cols()) { throw "non-square matrix has no LU decomposition"; }

		for (int i = 0; i < m.rows();++i) {
			//get entries for U's ith row and L's ith column
			for (int j = i; j < m.cols(); ++j) {
				U(i, j) = m(i, j);
				L(j, i) = m(j, i) / m(i, i);
			}
			//update matrix m(i+1:Row,i+1:Col)
			for (int k = i + 1; k < m.rows(); ++k) {
				for (int q = i + 1; q < m.cols(); ++q) {
					m(k, q) = m(k, q) - L(k, i)*U(i, q);
				}
			}
		}
		//construct tuple with L and U
		std::tuple<Matrix, Matrix> result(L, U);
		return result;
	}

//conduct LU with pivoting
std::tuple<Matrix, Matrix, Matrix> lu_row_pivoting(Matrix m) {
		//check if it is square matrix
		if (m.rows() != m.cols()) { throw "non-square matrix has no such decomposition"; }
		Matrix L(m.rows(),m.cols()), U(m.rows(),m.cols()); //construct L and U
		L.setZero(); U.setZero();
		Matrix P(m.rows(),m.cols()); //construct P
		//set P to be identity matrix
		P.setIdentity();
		double max; int index;

		for (int i = 0; i < m.rows(); ++i) {
			//find entry with max absolute value
			max =m.cwiseAbs().col(i).tail(m.cols()-i).maxCoeff(&index);
			
			SwitchRow(P, i, index+i);
			SwitchRow(m, i, index+i);
			SwitchRow(L, i, index+i);
			
			//get entries for U's ith row and L's ith column
			for (int j = i; j < m.cols(); ++j) {
				U(i, j) = m(i, j);
				L(j, i) = m(j, i) / m(i, i);
			}
			//update matrix A
			for(int k=i+1;k<m.rows();++k)
				for (int q = i + 1; q <m.cols(); ++q) {
					m(k, q) -= L(k, i)*U(i, q);
				}
			
		}
		std::tuple<Matrix, Matrix, Matrix> result(P, L, U);
		return result;
	}





#endif