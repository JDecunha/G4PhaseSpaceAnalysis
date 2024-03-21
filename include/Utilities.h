#pragma once

class TH1F;

namespace Utilities
{
	//Utilities
	std::vector<double> Linspace(float min, float max, int nVals);

	//Histogram utilities
	std::pair<float,float> GetHistogramAverage(const TH1F& histogram, int start, int end);
	void PrintHistogram(const TH1& toPrint);

	void VerifyNormalization(const TH1& h);
	void Prepare_for_Semilog(TH1* h);

	void PMF_to_FrequencyFunction(TH1* h);
	void PMF_to_DoseFunction(TH1* h);

	//File parsing utilities
	std::string GetFileEnergy(std::string filename);
	std::string GetFileTargetSize(std::string filename);

	TH1D GetNy(std::string path, double Energy, std::string TargetSize);
	long long GetEffectiveNumberOfTracks(std::string path, double Energy, std::string TargetSize);
	TH1D GetFy(std::string path, double Energy, std::string TargetSize);
	TH1D GetDy(std::string path, double Energy, std::string TargetSize);

	double GetYD(TH1* h);
};
