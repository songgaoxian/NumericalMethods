//implement LU for bandmatrix

#ifndef BandMatrix_HPP
#define BandMatrix_HPP
#include"MatrixLU.hpp"
#include <Eigen\src\Core\BandMatrix.h>

using BandMatrix = Eigen::internal::BandMatrix<double>;

//access B(i,j)
double& item(BandMatrix& B,unsigned int i,unsigned int j) {
	if (i >= B.rows()) { throw "row index overload"; }
	if (i > j + B.subs()) { throw "column index underflow"; }
	return B.coeffs()(B.supers() + i - j, j);
}
//transpose
BandMatrix Transpose(BandMatrix B) {
	BandMatrix result(B.cols(), B.rows(), B.subs(), B.supers());
	for (int i = 0; i < B.rows(); ++i) {
		for (int j = i; j < std::min(B.cols(), i + B.supers() + 1); ++j) {
			item(result, j, i) = item(B, i, j);
		}
		for (int j = std::max(0, i - B.subs()); j < i; ++j) {
			item(result, j, i) = item(B, i, j);
		}
	}
	return result;
}



//conduct LU decomposition without pivoting
std::tuple<BandMatrix, BandMatrix> lu_no_pivoting(BandMatrix m) {
	BandMatrix L(m.rows(), m.cols(),0,m.subs()), U(m.rows(), m.cols(),m.supers(),0);

	for (int i = 0; i < m.rows(); ++i) {
		//get entries for U's ith row and L's ith column
		for (int j = i; j < std::min(m.cols(),i+m.supers()+1); ++j) {
			item(U,i,j)= item(m,i, j);
		}
		for (int j = i; j < std::min(m.rows(), i + m.subs()+1); ++j) {
			item(L,j, i) = item(m, j, i) / item(m, i, i);
		}
		//update matrix m(i+1:Row,i+1:Col)
		for (int k = i + 1; k < std::min(m.rows(),i+m.subs()+1); ++k) {
			for (int q = i + 1; q < std::min(m.cols(),i+m.supers()+1); ++q) {
				item(m,k, q) = item(m,k, q) - item(L,k, i)*item(U,i, q);
			}
		}
	}
	//construct tuple with L and U
	std::tuple<BandMatrix, BandMatrix> result(L, U);
	return result;
}
//conduct LU decomposition without pivoting
std::tuple<BandMatrix, BandMatrix> lu_no_pivoting_tri(BandMatrix m) {
	BandMatrix L(m.rows(), m.cols(),0,1), U(m.rows(), m.cols(),1,0);

	for (int i = 0; i < m.rows(); ++i) {
		//get entries for U's ith row and L's ith column
		for (int j = i; j < std::min(m.cols(), i  + 2); ++j) {
			item(U,i, j) = item(m, i, j);
		}
		for (int j = i; j < std::min(m.rows(), i  + 2); ++j) {
			item(L,j, i) = item(m, j, i) / item(m, i, i);
		}
		//update matrix m(i+1:Row,i+1:Col)
		for (int k = i + 1; k < std::min(m.rows(), i + 2); ++k) {
			for (int q = i + 1; q < std::min(m.cols(), i + 2); ++q) {
				item(m, k, q) = item(m, k, q) - item(L,k, i)*item(U,i, q);
			}
		}
	}
	//construct tuple with L and U
	std::tuple<BandMatrix, BandMatrix> result(L, U);
	return result;
}

//conduct LU with pivoting
std::tuple<Matrix, Matrix, Matrix> lu_row_pivoting(BandMatrix A) {
	int sub = A.subs(), sup = A.supers();
	Matrix m = A.toDenseMatrix();
	Matrix L(m.rows(), m.cols()), U(m.rows(), m.cols()); //construct L and U
	L.setZero(); U.setZero();
	Matrix P(m.rows(), m.cols()); //construct P								  
	//set P to be identity matrix
	P.setIdentity();
	double max; int index;

	for (int i = 0; i < m.rows(); ++i) {
		//find the entry with max absolute value and exchange rows
		max = m.cwiseAbs().col(i).segment(i,std::min(m.cols()-i,sub+1)).maxCoeff(&index);
		SwitchRow(P, i, index+i);
		SwitchRow(m, i, index+i);
		SwitchRow(L, i, index+i);

		//get entries for U's ith row and L's ith column
		for (int j = i; j < std::min(m.cols(),i+sub+sup+1); ++j) {
			U(i, j) = m(i, j);		
		}
		for (int j = i; j < std::min(m.rows(), i + sub + 1); ++j) {
			L(j, i) = m(j, i) / m(i, i);
		}
		//update matrix A
		for (int k = i + 1; k<std::min(m.rows(),i+sub+1); ++k)
			for (int q = i + 1; q <std::min(m.cols(),i+sub+sup+1); ++q) {
				m(k, q) -= L(k, i)*U(i, q);
			}

	}
	std::tuple<Matrix, Matrix, Matrix> result(P, L, U);
	return result;
}


#endif