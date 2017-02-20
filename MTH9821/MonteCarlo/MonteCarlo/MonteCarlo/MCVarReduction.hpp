//implement variance reduction for monteCarlo
#ifndef MCVarReduction_HPP
#define MCVarReduction_HPP

#include"MonteCarlo.hpp"

class MCVarReduction : public MonteCarlo
{
public:
	MCVarReduction(long size_, std::tuple<double, double, double, double, double, double>& stck, int type): MonteCarlo(size_,stck,type){}
	
	std::tuple<double, double> ControlVariate(long N) {
		if (N > normal.size()) throw "normal variables are not enough";
		double S_bar = 0,V_bar = 0, b_nomi = 0, b_denomi = 0, b, V_cv = 0;
		std::vector<double> Vi, Si;
		double s0 = std::get<0>(stock);
		double k = std::get<1>(stock);
		double t = std::get<2>(stock);
		double sigma = std::get<3>(stock);
		double q = std::get<4>(stock);
		double r = std::get<5>(stock);
		double S, V;
		for (int i = 0; i < N; ++i) {
			S = s0*std::exp((r - q - 0.5*sigma*sigma)*t + sigma*std::sqrt(t)*normal[i]);
			S_bar += S / double(N);
			Si.push_back(S);
			V = std::exp(-r*t)*std::max(k - S, 0.0);
			V_bar += V / double(N);
			Vi.push_back(V);
		}
		for (int i = 0; i < N; ++i) {
			b_nomi += (Si[i] - S_bar)*(Vi[i] - V_bar);
			b_denomi += (Si[i] - S_bar)*(Si[i] - S_bar);
		}
		b = b_nomi / b_denomi;
		for (int i = 0; i < N; ++i) {
			double W = Vi[i] - b*(Si[i] - std::exp(r*t)*s0);
			V_cv += W / double(N);
		}
		auto bs_info = BSPricing();
		double V_bs = std::get<3>(bs_info);
		double error = std::abs(V_cv - V_bs);
		std::tuple<double, double> result(V_cv, error);
		return result;
	}

	std::tuple<double, double> AntitheticVariable(long N) {
		if (N > normal.size()) throw "uniform variables are not enough";
		double s0 = std::get<0>(stock);
		double k = std::get<1>(stock);
		double t = std::get<2>(stock);
		double sigma = std::get<3>(stock);
		double q = std::get<4>(stock);
		double r = std::get<5>(stock);
		double S1, S2, V1, V2, V_av = 0, V_bs, Z1, Z2;
		auto bs_info = BSPricing();
		V_bs = std::get<3>(bs_info);
		for (int i = 0; i < N; ++i) {
			Z1 = normal[i];
			Z2 = -Z1;
			S1 = s0*std::exp((r - q - sigma*sigma*0.5)*t + sigma*std::sqrt(t)*Z1);
			S2 = s0*std::exp((r - q - 0.5*sigma*sigma)*t + sigma*std::sqrt(t)*Z2);
			V1 = std::exp(-r*t)*std::max(k - S1, 0.0);
			V2 = std::exp(-r*t)*std::max(k - S2, 0.0);
			V_av += 0.5*(V1 + V2) / double(N);
		}
		double error = std::abs(V_av - V_bs);
		std::tuple<double, double> result(V_av, error);
		return result;
	}

	std::tuple<double, double> MomentMatch(long N) {
		if (N > normal.size()) throw "normal variables are not enough";
		double S_bar = 0, Vi, Vmm = 0, S;
		std::vector<double> Si;
		double s0 = std::get<0>(stock);
		double k = std::get<1>(stock);
		double t = std::get<2>(stock);
		double sigma = std::get<3>(stock);
		double q = std::get<4>(stock);
		double r = std::get<5>(stock);
		for (int i = 0; i < N; ++i) {
			S = s0*std::exp((r - q - 0.5*sigma*sigma)*t + sigma*std::sqrt(t)*normal[i]);
			Si.push_back(S);
			S_bar += S / double(N);
		}
		for (int i = 0; i < N; ++i) {
			S = Si[i] * std::exp(r*t)*s0 / S_bar;
			Vi = std::exp(-r*t)*std::max(k - S, 0.0);
			Vmm += Vi / double(N);
		}
		auto bs_info = BSPricing();
		double V_bs = std::get<3>(bs_info);
		double error = std::abs(Vmm - V_bs);
		std::tuple<double, double> result(Vmm, error);
		return result;
	}

	std::tuple<double, double> MomentMatchControlVariate(long N) {
		if (N > normal.size()) throw "normal variables are not enough";
		double s0 = std::get<0>(stock);
		double k = std::get<1>(stock);
		double t = std::get<2>(stock);
		double sigma = std::get<3>(stock);
		double q = std::get<4>(stock);
		double r = std::get<5>(stock);
		double S, S_bar = 0, S_blur, V_blur, V_bar = 0, Wi, b, b_nomi = 0, b_denomi = 0, Vcvmm = 0;
		std::vector<double> Si, Si_blur_v, V_blur_v;
		for (int i = 0; i < N; ++i) {
			S = s0*std::exp((r - 0.5*sigma*sigma - q)*t + sigma*std::sqrt(t)*normal[i]);
			Si.push_back(S);
			S_bar += S / double(N);
		}
		for (int i = 0; i < N; ++i) {
			S_blur = Si[i] * std::exp(r*t)*s0 / S_bar;
			Si_blur_v.push_back(S_blur);
			V_blur = std::exp(-r*t)*std::max(k - S_blur, 0.0);
			V_bar += V_blur / double(N);
			V_blur_v.push_back(V_blur);
		}
		for (int i = 0; i < N; ++i) {
			b_nomi += (Si_blur_v[i] - std::exp(r*t)*s0)*(V_blur_v[i] - V_bar);
			b_denomi += (Si_blur_v[i] - std::exp(r*t)*s0)*(Si_blur_v[i] - std::exp(r*t)*s0);
		}
		b = b_nomi / b_denomi;
		for (int i = 0; i < N; ++i) {
			Wi = V_blur_v[i] - b*(Si_blur_v[i] - std::exp(r*t)*s0);
			Vcvmm += Wi / double(N);
		}
		auto bs_info = BSPricing();
		double V_bs = std::get<3>(bs_info);
		double error = std::abs(Vcvmm - V_bs);
		std::tuple<double, double> result(Vcvmm, error);
		return result;
	}

};


#endif