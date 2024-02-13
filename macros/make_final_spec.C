/**
 * \file TBS.cxx
 * \brief Compute total beta spectra
 * \author L. Perisse
 * \date July 2018
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <cstdlib>
#include <dirent.h>

using namespace std;

#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TFile.h>
#include <TH1.h>
#include <TKey.h>
#include <TDirectory.h>
#include <TLatex.h>







TFile *file_save;
TFile *file_spectrum;
TFile *file_factor;
TFile *file_oscillation_up;
TFile *file_oscillation_down;

TH1D *nspec;
TH1D *nspec_OU;
TH1D *nspec_OD;
TH1D *oscillation_prefactor;

TH2F *ncov_Tot;
TH2F *ncor_Tot;

int Nbins   = 0;
double Emin = 0.;
double Emax = 0.;





////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Constant

//Physics constant
const double PI = M_PI;
const double MEC2 = 0.51099895000;  	//MeV.C-2, electron mass
const double MNC2 = 939.56542052;  	 	//MeV.C-2, neutron mass
const double MPC2 = 938.27208816;  	 	//MeV.C-2, proton mass
const double mu_p = 2.792847351; 	 	//Proton magnetic moment
const double mu_n = -1.9130427;		 	//Neutron magnetic moment
const double Ca   = -1.2756; 			//Axial-vector constant (from V-A theory), Ca=gA/gV
const double ALPHA = 0.0072973525693;	//QED fine structure constant, CODATA 2018
const double hBarC = 1.973269804e-13; 	//In MeV.m
const double Gf    = 1.1663788e-11; 	//Fermi cste, in MeV^(-2)
const double tau_n= 878.4; 				//Neutron average lifetime in seconds
const double fR   = 1.7152;				//Phase-space factor for beta-decay of the free neutron

//Conversion factor
const double sec2MeV   = 1.5192674e21;		//1 sec = 1.52e21 MeV-1
const double cm2ToMeV2 = 2.56819e21; 		//1cm² = 2.56819e21 MeV-²
const double uToG     = 1.66053906660e-24;	//CODATA 20.05.2019, uma to gram
const double uToMeV   = 931.49410242;		//CODATA 20.05.2019, uma to MeV



//Names of elements Z=0 to Z=118 of the Mendeleiev table.
const char element[][119] = {
	"n" ,"H" ,"He","Li","Be","B" ,"C" ,"N" ,"O" ,"F" ,
	"Ne","Na","Mg","Al","Si","P" ,"S" ,"Cl","Ar","K" ,
	"Ca","Sc","Ti","V" ,"Cr","Mn","Fe","Co","Ni","Cu",
	"Zn","Ga","Ge","As","Se","Br","Kr","Rb","Sr","Y" ,
	"Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In",
	"Sn","Sb","Te","I" ,"Xe","Cs","Ba","La","Ce","Pr",
	"Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm",
	"Yb","Lu","Hf","Ta","W" ,"Re","Os","Ir","Pt","Au",
	"Hg","Tl","Pb","Bi","Po","At","Rn","Fr","Ra","Ac",
	"Th","Pa","U" ,"Np","Pu","Am","Cm","Bk","Cf","Es", 
	"Fm","Md","No","Lr","Rf","Db","Sg","Bh","Hs","Mt",
	"Ds","Rg","Cn","Nh","Fl","Mc","Lv","Ts","Og"};
	

	
///////////////////////////////////////////////////////////////////////////////////////////////////

/**
 * \fn int n_files_in_dir(string path)
 * \brief Count the number of files located in the file pointed by "path"
 * \brief a folder will be counted as a file
 */	
 
int n_files_in_dir(string path)
{
    struct dirent *de;
    DIR *dir = opendir(path.data());
    if(!dir)
    {
        cout <<  path.data() << " not found ! Does it exist?\n";
        return 0;
    }

    int count=0;
    while(de = readdir(dir))
    {
		count++;
    }

    closedir(dir);

    return count-2;
}
	
	
///////////////////////////////////////////////////////////////////////////////////////////////////

/**
 * \fn string ListItem(string filecase)
 * \brief Display available items (file and case) from the given path
 * \param filecase Path located one step above from the TBS repertory
 */	

string ListItem(string folder)
{
	string list_case = "ls "; list_case += folder;
	string file = "" ;

	int cmd ;

	cout << n_files_in_dir(folder) <<" different databases have been found in the database repository...\n" << endl;
	cmd = system(list_case.c_str()) ;
	cout << "\n... which one should be used ?" << endl ;
	cin >> file ;
	cout << "\n => TBS will use " << file.data() << " as input" ;

	return file ;
}


///////////////////////////////////////////////////////////////////////////////////////////////////

void MakeCorrelationMatrix(TH2F *mat_cov, TH2F *mat_cor)
{
	vector<double> vec_error;

	for(int i=0 ; i<Nbins ; i++)
	{
		vec_error.push_back( sqrt(mat_cov->GetBinContent(i+1, i+1)) );
	}
	
	for(int i=0 ; i<Nbins ; i++)
	{
		for(int j=0 ; j<Nbins ; j++)
		{
			if(i==j)									 {mat_cor->SetBinContent(i+1, j+1, 1.);}
			else if(vec_error[i]!=0. && vec_error[j]!=0.){mat_cor->SetBinContent(i+1, j+1, mat_cov->GetBinContent(i+1, j+1) / vec_error[i] / vec_error[j]);}
			else							 			 {mat_cor->SetBinContent(i+1, j+1, 0.);}
		}
	}
}

void MakeCorrelationMatrix(TMatrixT<double> *mat_cov, TMatrixT<double> *mat_cor)
{
	vector<double> vec_error;

	for(int i=0 ; i<Nbins ; i++)
	{
		vec_error.push_back( sqrt((*mat_cov)(i,i)) );
	}
	
	for(int i=0 ; i<Nbins ; i++)
	{
		for(int j=0 ; j<Nbins ; j++)
		{
			if(i==j)									 {(*mat_cor)(i,j) = 1.;}
			else if(vec_error[i]!=0. && vec_error[j]!=0.){(*mat_cor)(i,j) = (*mat_cov)(i,j) / vec_error[i] / vec_error[j];}
			else										 {(*mat_cor)(i,j) = 0.;}			
		}
	}
}


///////////////////////////////////////////////////////////////////////////////////////////////////

double SumMatrix(TMatrixT<double> mat)
{
	int nbins = mat.GetNcols();
	double tmp_res = 0.;
	
	for(int i=0 ; i<nbins ; i++)
	{
		for(int j=0 ; j<nbins ; j++)
		{
			tmp_res += mat(i,j);
		}	
	}
	
	return tmp_res;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//xmin start of the intervalle where to look at best contributors
//xmax end of the intervalle where to look at best contributors

void Routine_cov(int label)
{
	string sufix = "";
	if     (label == 0){sufix = "FC";}
	else if(label == 1){sufix = "FF_RS";}
	else if(label == 2){sufix = "OS";}

	cout << endl << sufix << endl;


	/*****************************************/
	//Compute covariance
	//For FF + RS, we apply the oscillation prefactor to the covariance matrix of a single reactor spectrum
	
	TH2F *ncov;
	if(label == 1)
	{
		ncov = (TH2F *)file_spectrum->Get("ncov_Tot_spectrum");    ncov->SetName(Form("ncov_%s", sufix.c_str()));
		for(int i=0 ; i<Nbins ; i++)
		{
			int bin_i = oscillation_prefactor->FindBin(nspec->GetBinCenter(i+1));
			for(int j=0 ; j<i+1 ; j++)
			{
				int bin_j = oscillation_prefactor->FindBin(nspec->GetBinCenter(j+1));
				ncov->SetBinContent(i+1, j+1,   ncov->GetBinContent(i+1, j+1)  *  oscillation_prefactor->GetBinContent(bin_i)  *  oscillation_prefactor->GetBinContent(bin_j) );
				ncov->SetBinContent(j+1, i+1,   ncov->GetBinContent(i+1, j+1)  );
			}
		}		
	}
	else	   	       
	{
		ncov = new TH2F(Form("ncov_%s", sufix.c_str()), Form("ncov_%s", sufix.c_str()),   Nbins, Emin, Emax,   Nbins, Emin, Emax);
		for(int i=0 ; i<Nbins ; i++)
		{
			for(int j=0 ; j<i+1 ; j++)
			{
				if     (label == 0){ncov->SetBinContent(i+1, j+1,   0.005*0.005 * nspec->GetBinContent(i+1) * nspec->GetBinContent(j+1)  );}
				else if(label == 2){ncov->SetBinContent(i+1, j+1,   0.5 * ( (nspec_OU->GetBinContent(i+1) - nspec->GetBinContent(i+1)) * (nspec_OU->GetBinContent(j+1) - nspec->GetBinContent(j+1))   +  
																		    (nspec_OD->GetBinContent(i+1) - nspec->GetBinContent(i+1)) * (nspec_OD->GetBinContent(j+1) - nspec->GetBinContent(j+1)) )   );}

				ncov->SetBinContent(j+1, i+1,   ncov->GetBinContent(i+1, j+1)  );
			}
		}
	}

	ncov_Tot->Add(ncov, 1.);
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Type converter

void TH2toMatrix(TH2D *hist, TMatrixT<double> *mat)
{
    int size = hist->GetNbinsX();
    mat->ResizeTo(size, size);
    for(int i=0 ; i<size ; i++)
    {
        for(int j=0 ; j<size ; j++)
        {
            mat[0](i,j) = hist->GetBinContent(i+1, j+1);
        }
    }
}

void TH2toMatrix(TH2F *hist, TMatrixT<double> *mat)
{
    int size = hist->GetNbinsX();
    mat->ResizeTo(size, size);
    for(int i=0 ; i<size ; i++)
    {
        for(int j=0 ; j<size ; j++)
        {
            mat[0](i,j) = hist->GetBinContent(i+1, j+1);
        }
    }
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void WriteFile(string filename, TH1D *hspec, TH2F *hcov)
{
	ofstream file;
	file.open(filename.data(), ios_base::app);

	cout << "\nWritting to file: " << filename << endl;

	file << setprecision(3);
	file << fixed;
	file << "Reactor neutrino spectrum at SK from 0 to 16 MeV in 1 keV bins" << endl;
	file << "Bin center [MeV]          Spectrum [nu/fission/MeV]     Uncertainty [nu/fission/MeV]" << endl;
	file << setprecision(6);
	file << scientific;

	//Spectrum and uncertainty
	for(int i=0 ; i<hspec->GetNbinsX() ; i++)
	{
		file << hspec->GetBinCenter(i+1) << "                 " << hspec->GetBinContent(i+1) << "                 " << sqrt(hcov->GetBinContent(i+1, i+1)) << endl; 
	}

	file.close();
	cout << "End of file writting" << endl;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void macro() 
{
	// int clear = system("clear");
	int rbin = 10;
	string path_SKReact= "/mnt/c/Users/lolos/Documents/Soft/skreact/SKReact";
	string path_write  = path_SKReact + "/data/Prediction_SKReact_v2p0.dat";


	//Choose which ROOT to use
	string path_save               = path_SKReact + "/data/Prediction_SKReact_v2p0.root";
	string source_spectrum         = path_SKReact + "/data/ReactorModel_SKReact_v2p0.root";
	string source_oscillation_up   = path_SKReact + "/data/ReactorModel_SKReact_v2p0_oscillationUp.root";
	string source_oscillation_down = path_SKReact + "/data/ReactorModel_SKReact_v2p0_oscillationDown.root";

	file_save     = new TFile(path_save.data(), "recreate");
	file_spectrum = new TFile(source_spectrum.data());
	file_oscillation_up   = new TFile(source_oscillation_up.data());
	file_oscillation_down = new TFile(source_oscillation_down.data());


/********************************************************/
//Define histograms

	// file_factor = new TFile(source_spectrum.data());
	oscillation_prefactor = (TH1D *)file_spectrum->Get("oscillation_prefactor_SK");
		
	nspec    = (TH1D *)file_spectrum->Get("nspec_SK");          nspec->SetName("nspec");
	nspec_OU = (TH1D *)file_oscillation_up->Get("nspec_SK");    nspec_OU->SetName("nspec_OU");
	nspec_OD = (TH1D *)file_oscillation_down->Get("nspec_SK");  nspec_OD->SetName("nspec_OD");

	Nbins  = nspec->GetNbinsX(); 						
	Emin= nspec->GetBinLowEdge(1); 			//in MeV
	Emax= nspec->GetBinLowEdge(Nbins+1); 	//in MeV	

	ncov_Tot = new TH2F("ncov_Tot", "ncov_Tot",   Nbins, Emin, Emax,   Nbins, Emin, Emax);
	ncor_Tot = new TH2F("ncor_Tot", "ncor_Tot",   Nbins, Emin, Emax,   Nbins, Emin, Emax);

	cout << endl << endl;
	// cout << "nu flux   = " << nspec->Integral("width") << "  " << nspec_OU->Integral("width") << "  " << nspec_OD->Integral("width") << endl;
	

	//Make covariance
	Routine_cov(0);  //FC
	Routine_cov(1);  //FF + RS
	Routine_cov(2);  //OS


	//Make correlation
	nspec->Rebin(rbin);						nspec->Scale(1./rbin);
	ncov_Tot->Rebin2D(rbin, rbin);			ncov_Tot->Scale(1./rbin/rbin);
	ncor_Tot->Rebin2D(rbin, rbin);			ncor_Tot->Scale(0.);
	MakeCorrelationMatrix(ncov_Tot, ncor_Tot);


/********************************************************/
//Check positive definess
//Cholesky decomposition, computes Chol s.t. Chol * CholT = Cor

    TMatrixT<double> spec_cor_skreact_m;
    TH2toMatrix(ncor_Tot, &spec_cor_skreact_m);

	for(int i=0 ; i<nspec->GetNbinsX() ; i++)
	{
		nspec->SetBinError(i+1,  sqrt(ncov_Tot->GetBinContent(i+1, i+1))  );
		
		for(int j=0 ; j<nspec->GetNbinsX() ; j++)
		{
			if(i != j){spec_cor_skreact_m(i,j) *= 0.99;}  //Regularization
		}
	}

	TDecompChol finder(spec_cor_skreact_m);               //Cholesky decomposition
	bool IsSemiPositive = finder.Decompose();
	if(!IsSemiPositive)
	{
		cout << "Correlation matrix not semi-positive, you must recomputing it..." << endl;
		exit(1);
	}
	

/********************************************************/
//Write and save to file

	// WriteFile(path_write, nspec, ncov_Tot);


	file_save->cd();
	nspec->Write();
	ncov_Tot->Write();
	ncor_Tot->Write();

	file_save->Close();
	delete file_save;
}

