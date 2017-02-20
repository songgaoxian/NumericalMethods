#include<cmath>
#include<iostream>
#include"BinomialBSR.hpp"
#include"Binomial_BlackScholes.hpp"
#include"AverageBinomial.hpp"
#include"BinomialTree.hpp"
#include<sstream>
#include<string>
#include<fstream>
#include<iomanip>
#include<boost\algorithm\string.hpp>
#include"BinomialBarrierOption.hpp"
#include"BinomialDiscreteDiv.hpp"
#include"SecantMethod.hpp"
#include<utility>

void IniMatrixCsv(std::ifstream& file, Matrix& A) {
	std::string line, field;
	std::vector<std::string> dataline;
	while (std::getline(file, line)) {
		dataline.push_back(line);
	}

	std::stringstream is;
	for (int i = 0; i < dataline.size(); ++i) {
		std::vector<std::string> data;
		boost::split(data, dataline[i], boost::is_any_of(" ,\t"));
		for (int j = 0; j < data.size(); ++j) {
			is = std::stringstream(data[j]);
			is >> A(i, j);
		}
	}
}

void IniVectorCsv(std::ifstream& file, Vector& A) {
	std::string line, field;
	std::vector<std::string> dataline;
	std::stringstream is;
	int i = 0;
	while (std::getline(file, line)) {
		is = std::stringstream(line);
		is >> A(i);
		i++;
	}
}

void WriteToCsv(std::string filename, Matrix A) {
	std::ofstream file(filename);
	for (int i = 0; i < A.rows(); ++i) {
		for (int j = 0; j < A.cols(); ++j) {
			file << std::setprecision(18) << A(i, j) << ",";
		}
		file << "\n";
	}
}
void WriteToCsv(std::string filename, Vector b, double r) {
	std::ofstream file(filename);
	for (int i = 0; i < b.rows(); ++i) {
		file << std::setprecision(18) << b(i) << "\n";
	}
	file << "residualORcount\n";
	file << std::setprecision(18) << r << "\n";
}
template<typename T>
double TupleElement(T mytuple,int index) {
	if (index == 0) return std::get<0>(mytuple);
	if (index == 1) return std::get<1>(mytuple);
	if (index == 2) return std::get<2>(mytuple);
	if (index == 3) return std::get<3>(mytuple);
	if (index == 4) return std::get<4>(mytuple);
	if (index == 5) return std::get<5>(mytuple);
	if(index == 6) return std::get<6>(mytuple);
	if (index == 7) return std::get<7>(mytuple);
	return 0;
}

template<std::size_t I=0, typename... Tp>
inline typename std::enable_if<I==sizeof...(Tp),void>::type
  TupleToVector(Vector& v, std::tuple<Tp...>&t) {}

template<std::size_t I=0, typename... Tp>
inline typename std::enable_if<I<sizeof...(Tp),void>::type
	TupleToVector(Vector& v, std::tuple<Tp...>&t) {
	v(I) = std::get<I>(t);
	TupleToVector<I + 1, Tp...>(v,t);
}

int main() {
	//s,sigma,q,r,k,t
	double s = 42.0, sigma = 0.1, q = 0.02, r = 0.04, k = 45, t = 9.0/12.0;
	bool isCall = false;
	bool isEuro =false;
	std::tuple<double, double, double, double, double, double> inputData(s,sigma,q,r,k,t);
	BinomialTree bt(inputData, 2500, isCall, isEuro);
	double tol = 1E-4;
	SecantVolatility s_vol(bt, 0.1, 0.5, 4.09, tol);
	Vector v(1);
	v(0) = s_vol.ImpliedVolConsecutive();
	WriteToCsv("implied1.csv", v);
	/*
	double bs_price = BSPrice(s, sigma, q, r, k, t, isCall);
	SecantVolatility s_vol(bt, 0.4, 0.36, bs_price, 0.00001);
	std::cout << s_vol.ImpliedVolAbsolute() << "absolute:\n";
	std::cout << s_vol.ImpliedVolConsecutive() << "consecutive:\n";*/
	return 0;
}