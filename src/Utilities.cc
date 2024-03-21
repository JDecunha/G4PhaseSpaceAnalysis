//std
#include <iostream>
#include <filesystem>
#include <vector>
#include <utility>
#include <string>

//ROOT
#include "TH1F.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TH1D.h"
#include "TString.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TSystem.h"

//This project
#include "Utilities.h"

std::vector<double> Utilities::Linspace(float min, float max, int nVals)
{
	std::vector<double> outputVector;
	outputVector.reserve(nVals); //reserve memory for efficiency
	float stepSize = (max - min)/((float)nVals-1.);

	for (int i = 0; i < nVals; ++i)
	{
		outputVector.push_back(min+(stepSize*(float)i));
	}

	return outputVector;
}

std::pair<float,float> Utilities::GetHistogramAverage(const TH1F& histogram, int start, int end)
{
	float accumulator = 0;
	float binCenterAccumulator = 0;
	float nVals = end - start + 1;

	for (int i = start; i <= end; i++)
	{
		accumulator += histogram.GetBinContent(i);
		binCenterAccumulator += histogram.GetBinCenter(i);
	}

	return std::pair<float,float>(accumulator/nVals,binCenterAccumulator/nVals);
}

std::string Utilities::GetFileEnergy(std::string filename)
{
	std::stringstream pathsplitter(filename);
	std::string output;

	//This is a dumb workaround because of my naming convention on the new lineal energy library
	std::getline(pathsplitter,output,'_');
	std::getline(pathsplitter,output,'_');
	std::getline(pathsplitter,output,'_');

	while (std::getline(pathsplitter,output,'_'))
	{
		//Find and erase MeV
		auto n = output.find("MeV");
		if (n != -1)
		{
			output.erase(n,3);
			return output;
		}
	}

	return "failure :(";
}

std::string Utilities::GetFileTargetSize(std::string filename)
{
	std::stringstream pathsplitter(filename);
	std::string output;

	while (std::getline(pathsplitter,output,'_'))
	{
		//Find and erase nm
		auto n = output.find("nm");
		if (n != -1)
		{
			output.erase(n,2);
			return output;
		}
	}

	return "failure :(";
}

double Utilities::GetYD(TH1* h)
{
	//Input: takes a normalized d(y) function
	//In future: to get uncertainty out of these for the positions on the jig will require some calculations
	//To figure out how to propagate all the uncertainties.

	double normalization = 0;
	double yd = 0;
	int length = h->GetNbinsX(); //get length

	//Step 1.) Should be normalized, but we'll double check
	for (int i = 1; i <= length; i++) 
	{
		auto value = h->GetBinContent(i);
		double low_edge = h->GetBinLowEdge(i);
		double high_edge = h->GetBinLowEdge(i+1);
		double width = high_edge-low_edge;
		normalization += value*width;
	}
	for (int i = 1; i <= length; i++) 
	{
		auto value = h->GetBinContent(i);
		auto normalized = value/normalization;
		h->SetBinContent(i,normalized);
	}

	//Step 2.) Calculate y_d
	for (int i = 1; i <= length; i++) 
	{
		auto value = h->GetBinContent(i);
		double low_edge = h->GetBinLowEdge(i);
		double high_edge = h->GetBinLowEdge(i+1);
		double y = (high_edge+low_edge)/2;
		double width = high_edge-low_edge;

		yd += value*width*y;
	}	

	return yd;
}

void Utilities::PMF_to_FrequencyFunction(TH1* h)
{
	//This converts a probability mass function to a dose-weighted probability density function
	int length = h->GetNbinsX();

	double normalization = 0;

	for (int i = 1; i <= length; i++) //bins per bin width
	{
		double value = h->GetBinContent(i);
		double low_edge = h->GetBinLowEdge(i);
		double high_edge = h->GetBinLowEdge(i+1);
		double bin_middle = (high_edge+low_edge)/double(2);
		double width = h->GetBinWidth(i);
		double per_width = value/(width);

		//This works because each bin value = value/width now
		//So integration would be given by value*width/width
		normalization += value;
		h->SetBinContent(i,per_width);
	}
	
	for (int i = 1; i <= length; i++) //normalize
	{
		double value = h->GetBinContent(i);
		double normalized = value/normalization;

		h->SetBinContent(i,normalized);
	}
}

void Utilities::PMF_to_DoseFunction(TH1* h)
{
	//This converts a probability mass function to a dose-weighted probability density function
	int length = h->GetNbinsX();

	double normalization = 0;

	for (int i = 1; i <= length; i++) //bins per bin width
	{
		double value = h->GetBinContent(i);
		double low_edge = h->GetBinLowEdge(i);
		double high_edge = h->GetBinLowEdge(i+1);
		double bin_middle = (high_edge+low_edge)/double(2);
		double width = h->GetBinWidth(i);
		double per_width = bin_middle*value/(width);

		//This works because each bin value = bin_middle*value/width now
		//So integration would be given by bin_middle*value*width/width
		normalization += bin_middle*value;
		h->SetBinContent(i,per_width);
	}
	
	for (int i = 1; i <= length; i++) //normalize
	{
		double value = h->GetBinContent(i);
		double normalized = value/normalization;

		h->SetBinContent(i,normalized);
	}
}

TH1D Utilities::GetDy(std::string path, double Energy, std::string TargetSize)
{
	TH1::AddDirectory(false); //So that we own the TH1 and gROOT won't delete it on us
	TH1D output; //Output histogram
	bool fileFound = false; //Will throw an exception if we don't find the file

	//We have separated the files into different folders by target size. Specify that path.
	path = path + "/" + TargetSize + "nm"; 

	for (const auto &entry : std::filesystem::directory_iterator(path)) //Loop over all the files in the folder
	{
		double energy = std::stod(GetFileEnergy(entry.path().filename())); //Get the energy of the file

		//Check if the energy is what we requested
		if (std::fabs(Energy - energy) < 0.000001) //I wish this wasn't necessary. But subtraction of the energies never returns exactly zero
		{
			TFile f = TFile((TString)entry.path());
			output = std::move(*(TH1D*)f.Get("Lineal energy histogram")); //Move the value pointed to by the TH1D pointer on to the stack
			PMF_to_DoseFunction(&output); //Transform the N(y) to a d(y)
			fileFound = true;
			break;
		}
	}

	if (fileFound) { return output; } else {throw std::runtime_error("From Utilities::GetDy. Lineal energy histogram not found.");}
}

TH1D Utilities::GetFy(std::string path, double Energy, std::string TargetSize)
{
	TH1::AddDirectory(false); //So that we own the TH1 and gROOT won't delete it on us
	TH1D output; //Output histogram
	bool fileFound = false; //Will throw an exception if we don't find the file

	//We have separated the files into different folders by target size. Specify that path.
	path = path + "/" + TargetSize + "nm"; 

	for (const auto &entry : std::filesystem::directory_iterator(path)) //Loop over all the files in the folder
	{
		double energy = std::stod(GetFileEnergy(entry.path().filename())); //Get the energy of the file

		//Check if the energy is what we requested
		if (std::fabs(Energy - energy) < 0.000001) //I wish this wasn't necessary. But subtraction of the energies never returns exactly zero
		{
			TFile f = TFile((TString)entry.path());
			output = std::move(*(TH1D*)f.Get("Lineal energy histogram")); //Move the value pointed to by the TH1D pointer on to the stack
			PMF_to_FrequencyFunction(&output); //Transform the N(y) to a f(y)
			fileFound = true;
			break;
		}
	}

	if (fileFound) { return output; } else {throw std::runtime_error("From Utilities::GetDy. Lineal energy histogram not found.");}
}

TH1D Utilities::GetNy(std::string path, double Energy, std::string TargetSize)
{
	TH1::AddDirectory(false); //So that we own the TH1 and gROOT won't delete it on us
	TH1D output; //Output histogram
	bool fileFound = false; //Will throw an exception if we don't find the file

	//We have separated the files into different folders by target size. Specify that path.
	path = path + "/" + TargetSize + "nm"; 

	for (const auto &entry : std::filesystem::directory_iterator(path)) //Loop over all the files in the folder
	{
		double energy = std::stod(GetFileEnergy(entry.path().filename())); //Get the energy of the file

		//Check if the energy is what we requested
		if (std::fabs(Energy - energy) < 0.000001) //I wish this wasn't necessary. But subtraction of the energies never returns exactly zero
		{
			TFile f = TFile((TString)entry.path());
			output = std::move(*(TH1D*)f.Get("Lineal energy histogram")); //Move the value pointed to by the TH1D pointer on to the stack
			fileFound = true;
			break;
		}
	}

	if (fileFound) { return output; } else {throw std::runtime_error("From Utilities::GetNy. Lineal energy histogram not found.");}
}

long long Utilities::GetEffectiveNumberOfTracks(std::string path, double Energy, std::string TargetSize)
{
	bool fileFound = false; //Will throw an exception if we don't find the file
	long long effectiveNumTracks = 0;

	//We have separated the files into different folders by target size. Specify that path.
	path = path + "/" + TargetSize + "nm"; 

	for (const auto &entry : std::filesystem::directory_iterator(path)) //Loop over all the files in the folder
	{
		double energy = std::stod(GetFileEnergy(entry.path().filename())); //Get the energy of the file

		//Check if the energy is what we requested
		if (std::fabs(Energy - energy) < 1e-6) //I wish this wasn't necessary. But subtraction of the energies never returns exactly zero
		{
			TFile f = TFile((TString)entry.path());
			TNamed* tnNumTracks = (TNamed*)f.Get("Effective Number of Tracks"); 
			if (tnNumTracks)
			{
				effectiveNumTracks = std::strtoll(tnNumTracks->GetTitle(), nullptr, 10);
				fileFound = true;
			}
			break;
		}
	}

	if (fileFound) { return effectiveNumTracks; } else {throw std::runtime_error("From Utilities::GetEffectiveNumberOfTracks. Lineal energy histogram not found.");}
}

void Utilities::PrintHistogram(const TH1& toPrint)
{
	for (int i = 0; i <= toPrint.GetNbinsX(); ++i)
	{
		std::cout << "bin: " << i << " bin value: " << toPrint.GetBinLowEdge(i) << " value: " << toPrint.GetBinContent(i) << std::endl; 
	}
}

void Utilities::Prepare_for_Semilog(TH1* h)
{
	//This takes a PDF and normalizes it by an extra factor of y (i.e. to give y*f(y))to preserve the graphical properties
	//of a PDF on a semilog axis
	int length = h->GetNbinsX();
	for (int i = 1; i <= length; i++) //bins per bin width
	{
		auto value = h->GetBinContent(i);
		auto low_edge = h->GetBinLowEdge(i);
		auto high_edge = h->GetBinLowEdge(i+1);
		double mid_value = (high_edge+low_edge)/double(2);
		h->SetBinContent(i,value*mid_value);
	}
}

void Utilities::VerifyNormalization(const TH1& h)
{
	double totalVal = 0;

	for (int i = 1; i <= h.GetNbinsX(); i++) //Start at 1 to skip underflow bin
	{
		auto value = h.GetBinContent(i);
		auto width = h.GetBinWidth(i);
		//auto edge = h.GetBinLowEdge(i);
		//std::cout << "Edge: " << edge << " Value: " << value << " Width: " << width << std::endl;
		totalVal += value*width;
	}

	double differenceFromOne = std::fabs(1.-totalVal);
	std::cout << "Normalization deviation from 1: " << differenceFromOne << std::endl;
	assert(differenceFromOne < 1e-5); //close enough to 1.
}