#include "HeatEquationSolver.hpp"
#include "AmericanOptionSolver.hpp"
#include "EuroBarrier.hpp"
#include<utility>
#include<vector>
#include<boost\algorithm\string.hpp>
#include"BarrierUpLow.hpp"
#include"EuroOptionSolver2.hpp"
#include"AmericanOptionSolver2.hpp"

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
	std::cout << std::setprecision(8);

	double tol = 1E-8;
	double w = 1.2;
	double sigma = 0.3, s0 = 50, q = 0.005, K = 48, T = 3.0/12.0, r = 0.02;
	double qdiv = 0.02, t_div = 5.0 / 12.0;
	double alpha = 0.5, alpha2 = 4.0;
	bool isCall = true;
	int M = 4, N = 12;
	double Blow = 40, Bup = 60;
	double payUp = 2.0, payLow = 2.0;
	Matrix Domain1(4, 6), Domain2(4, 6);
	for (int i = 0; i < 4; ++i) {
		BarrierUpLow bUL(payUp, payLow, Bup, Blow, isCall, K, sigma, s0, q, T, r, M, alpha, tol, w);
		BarrierUpLow bUL2(payUp, payLow, Bup, Blow, isCall, K, sigma, s0, q, T, r, M, alpha2, tol, w);
		auto tup = bUL.Domain();
		Vector v(6);
		TupleToVector(v, tup);
		Domain1.row(i) = v.transpose();
		tup = bUL2.Domain();
		TupleToVector(v, tup);
		Domain2.row(i) = v.transpose();
		M *= 4;
	}
	WriteToCsv("Domain1.csv", Domain1);
	WriteToCsv("Domain2.csv", Domain2);
	Matrix FE3(4, 6);
	M = 4;
	for (int i = 0; i < 4; ++i) {
		BarrierUpLow bUL(payUp, payLow, Bup, Blow, isCall, K, sigma, s0, q, T, r, M, alpha, tol, w);
		bUL.ForwardEuler();
		auto tup = bUL.GetUs();
		FE3(i, 0) = std::get<0>(tup);
		FE3(i, 1) = std::get<1>(tup);
		FE3(i, 2) = std::get<0>(bUL.Approximate0());
		auto tup2 = bUL.DeltaGammaTheta();
		FE3(i, 3) = std::get<0>(tup2);
		FE3(i, 4) = std::get<1>(tup2);
		FE3(i, 5) = std::get<2>(tup2);
		M *= 4;
	}
	WriteToCsv("FE3.csv", FE3);
	M = 4;
	for (int i = 0; i < 4; ++i) {
		BarrierUpLow bUL(payUp, payLow, Bup, Blow, isCall, K, sigma, s0, q, T, r, M, alpha, tol, w);
		bUL.CrankNicolson(true);
		auto tup = bUL.GetUs();
		FE3(i, 0) = std::get<0>(tup);
		FE3(i, 1) = std::get<1>(tup);
		FE3(i, 2) = std::get<0>(bUL.Approximate0());
		auto tup2 = bUL.DeltaGammaTheta();
		FE3(i, 3) = std::get<0>(tup2);
		FE3(i, 4) = std::get<1>(tup2);
		FE3(i, 5) = std::get<2>(tup2);
		M *= 4;
	}
	WriteToCsv("CN4.csv", FE3);
	M = 4;
	for (int i = 0; i < 4; ++i) {
		BarrierUpLow bUL(payUp, payLow, Bup, Blow, isCall, K, sigma, s0, q, T, r, M, alpha2, tol, w);
		bUL.CrankNicolson(true);
		auto tup = bUL.GetUs();
		FE3(i, 0) = std::get<0>(tup);
		FE3(i, 1) = std::get<1>(tup);
		FE3(i, 2) = std::get<0>(bUL.Approximate0());
		auto tup2 = bUL.DeltaGammaTheta();
		FE3(i, 3) = std::get<0>(tup2);
		FE3(i, 4) = std::get<1>(tup2);
		FE3(i, 5) = std::get<2>(tup2);
		M *= 4;
	}
	WriteToCsv("CN5.csv", FE3);
	M = 4;
	BarrierUpLow bUL(payUp, payLow, Bup, Blow, isCall, K, sigma, s0, q, T, r, M, alpha, tol, w);
	auto Rt = bUL.ForwardEuler();
	auto UFE = bUL.GetU();
	WriteToCsv("FEU6.csv", UFE);
	WriteToCsv("FEOpt7.csv", Rt);
	Rt = bUL.CrankNicolson(true);
	auto UCN = bUL.GetU();
	WriteToCsv("CNU8.csv", UCN);
	WriteToCsv("CNOpt9.csv", Rt);
	return 0;
}