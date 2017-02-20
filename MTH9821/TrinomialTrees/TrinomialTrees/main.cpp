#include<cmath>
#include<iostream>
#include"TrinomialBSR.hpp"
#include"Trinomial_BlackScholes.hpp"
#include"TrinomialTree.hpp"
#include<sstream>
#include<string>
#include<fstream>
#include<iomanip>
#include<boost\algorithm\string.hpp>
#include"TrinomialBarrierOption.hpp"
#include<utility>
#include"SecantMethod.hpp"

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
template<std::size_t I = 0, typename... Tp>
inline typename std::enable_if<I == sizeof...(Tp), void>::type
TupleToVector(Vector& v, std::tuple<Tp...>&t) {}

template<std::size_t I = 0, typename... Tp>
inline typename std::enable_if<I<sizeof...(Tp), void>::type
	TupleToVector(Vector& v, std::tuple<Tp...>&t) {
	v(I) = std::get<I>(t);
	TupleToVector<I + 1, Tp...>(v, t);
}

int main() {
	//s,sigma,q,r,k,t
	double s = 50.0, sigma = 0.3, q = 0.005, r = 0.02, k = 48.0, t = 3.0/12.0;
	bool iscall = true, iseuro = true;
	double Bup = 60, Blow = 40, pay = 2.0;
	std::tuple<double, double, double, double, double, double> inputData(s, sigma, q, r, k, t);
	int n = 10;
	Matrix Val2(8, 4);
	for (int i = 0; i < 8; ++i) {
		TrinomialDownOut tO(Bup, Blow, pay, inputData, n, iscall, iseuro);
		auto tup = tO.Pricing();
		Vector v(4);
		TupleToVector(v, tup);
		Val2.row(i) = v.transpose();
		std::cout << n << "\n";
		n *= 2;
	}
	WriteToCsv("Val2.csv", Val2);
	n = 20000;
	TrinomialDownOut tDO(Bup, Blow, pay, inputData, n, iscall, iseuro);
	Vector exact(1);
	exact(0) = std::get<0>(tDO.Pricing());
	WriteToCsv("ExactP1.csv", exact);

	return 0;
}