#include"MCVarReduction.hpp"
#include"MCBasketOption.hpp"
#include"HestonModel.hpp"
#include<Eigen\Dense>
#include<sstream>
#include<string>
#include<fstream>
#include<iomanip>
#include<boost\algorithm\string.hpp>
#include"MonteCarloDiscreteDividend.hpp"
#include<utility>

using Matrix = Eigen::MatrixXd;
using Vector = Eigen::VectorXd;

void IniMatrix(std::string filename, Matrix& A) {
	std::ifstream file(filename);
	std::string line;
	std::vector<std::string> dataline;
	while (std::getline(file, line)) { dataline.push_back(line); }
	std::stringstream is;
	for (int i = 0; i < dataline.size(); ++i) {
		std::vector<std::string> data;
		boost::split(data, dataline[i], boost::is_any_of(",\t "));
		for (int j = 0; j < data.size(); ++j) {
			is = std::stringstream(data[j]);
			is >> A(i, j);
		}
	}
}
void IniVector(std::string filename, Vector& A) {
	std::ifstream file(filename);
	std::string line;
	int i = 0;
	std::stringstream is;
	while (std::getline(file, line)) {
		is = std::stringstream(line);
		is >> A(i++);
	}
}
void WriteToCsv(std::string filename, Matrix A) {
	std::ofstream file(filename);
	for (int i = 0; i < A.rows(); ++i) {
		for (int j = 0; j < A.cols(); ++j) {
			file <<std::setprecision(18)<<A(i, j) << ",";
		}
		file << "\n";
	}
}
void WriteToCsv(std::string filename, Vector b, double r) {
	std::ofstream file(filename);
	for (int i = 0; i < b.rows(); ++i) {
		file << std::setprecision(18) << b(i) << "\n";
	}
	file << "residualOrcount\n";
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
	std::cout << std::setprecision(18);
	long N =std::pow(2,10)*10000 ;
	double s0 =50;
	double k =48;
	double t =3.0/12.0;
	double sigma =0.3;
	double q = 0.005;
	double r = 0.02;
	std::tuple<double, double, double, double, double, double> stock(s0, k, t, sigma, q, r);
	double Bup = 60.0, Blow = 40.0, pay = 2.0;
	MCBarrierOption mc_b(N, Bup, Blow, pay, stock, 0);
	Matrix MCR(10, 3);
	long Nk = 10000;
	for (int i = 0; i < 10; ++i) {
		Vector v(3);
		auto tup = mc_b.DependentBarriertCall(Nk);
		TupleToVector(v, tup);
		MCR.row(i) = v.transpose();
		std::cout << Nk << "\n";
		Nk *= 2;
	}
	WriteToCsv("MCResult.csv", MCR);





	/*std::vector<double> times{ 2.0 / 12.0,4.0 / 12.0,6.0 / 12.0 };
	std::vector<DividendFunct> functs;
	DividendFunct dvf = [](double s) {return s - 0.5; };
	functs.push_back(dvf);
	dvf = [](double s) {return s*(1 - 0.01); };
	functs.push_back(dvf);
	dvf = [](double s) {return s - 0.75; };
	functs.push_back(dvf);
	MonteCarloDiscreteDividend mcd(N, stock, 4, times, functs,false);//need to adjust delta
	Matrix m(8, 2);
	for (int k = 0; k <= 7; ++k) {
		long n = 10000 * std::pow(2, k);
		auto result = mcd.PriceDelta(n);
		m(k, 0) = std::get<0>(result);
		m(k, 1) = std::get<1>(result);
		std::cout << k << "\n";
	}
	*/
	return 0;
}