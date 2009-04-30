//CalcPack.h
//The following is a container class
//holds needed values to calculate various scores throughout the program


#ifndef CalcPack_H
#define CalcPack_H

#include <iostream>
#include <stdlib.h>
#include <string>
#define _USE_MATH_DEFINES
#include <math.h>
#include <fstream>
#include <map>
#include <sstream>
#include "Translate.h"
#include "GO.h"
#include "FastaRead.h"

using std::cout;
using std::cerr;
using std::string;
using std::ostream;
using std::ofstream;
using std::ifstream;
using std::map;
using std::stringstream;


//Values obtained from Glimmer 3.02
/*#define  DEFAULT_POS_ENTROPY_PROF  {0.08468,0.01606,0.05739,0.05752,0.04328,\
  0.07042,0.02942,0.05624,0.04442,0.05620,0.03029,0.03975,0.05116,0.04098,\
  0.05989,0.08224,0.05660,0.06991,0.02044,0.03310}
#define  DEFAULT_NEG_ENTROPY_PROF  {0.07434,0.03035,0.05936,0.04729,0.05662,\
  0.07704,0.05777,0.05328,0.03360,0.05581,0.01457,0.03718,0.04594,0.05977,\
  0.08489,0.05990,0.04978,0.07227,0.01050,0.01974}*/

class CalcPack{//open prototype
public:
	double Lambda;//constant in MaxBit calculation
	double K;//constant in MaxBit calculation
	map <string,int> ConsValue;//map for storing the max value of a conserved amino acid
	map<string,string> Genomes; //for storing the genome
	map<string,string>::iterator CurrentGenome;//iterator that points to the genome currently being used
        map<string, long> GenomeInfo;//stores genome information using the offset
        map<string,double> FStartCodons;//stores the start codons and their weight
        map<string,double> RStartCodons;
        set<string> FStopCodons;
        set<string> RStopCodons;
	string CurrentGenomeID;
	string TransFile;//the command for running entropy-calc
	string Matrix;
	string GenomeFile;
	int GenomeSize;
	int OStartCount;//Count for find starts that keeps track of the number of times original start has been processed
	double DefaultCProfile[20];
	double DefaultNCProfile[20];//for non coding genes
	double CProfile[20];
	double NCProfile[20];
	double SmallNCProfile[20];
	double SmallCProfile[20];
	Translate Translator;//for translating nucl. to AA
	bool UseSmallProf;
	int NumSmallC;//the number of small orfs used to calculate new EDP
	int NumLargeC;//the number of orfs >300bp used to calculate new EDP
	int NumLargeNC;
	int NumSmallNC;
	double* UseCProfile;//the coding profile in use
	double* UseNCProfile;//the non coding profile in use
	GO* GOAccess;
	FastaRead Reader; //For reading fasta files

	CalcPack();//default constructor
	

	CalcPack(string Matx, string GF, string TF, int TN);//parameterized constructor
		
	//Copy Constructor
	//Should not be USED!! pack contains GENOME
	 CalcPack(const CalcPack &Source);// open defintion
	

	//Assignment operator
	CalcPack& operator =(const CalcPack &Source);// open defintion
		
	//Function for setting ontological access
	int SetGOAccess(GO* Access);
	
	//Function to initialize parameters based on blast matrix used
	int InitCodes();//open definition
		


	int SetDefaultEDP();
	

	//this function returns the coordinate of the amino acid in the array
	//expects that any non-AA coding seq. enountered will be translated to '*'
	int MapAA(char AA);
		
	//This function calculates the entropy distance for the section of the genome specified
	double GetEntropy(int AACount[]);//open definition

	//create EDP from AACounts
	int CountToEDP(int AACount[], double ResultEDP[]);
	
        //Create new EDP for coding and non coding genes specific to this organism
	int CreateOrgEDPs();
		
	//function for finding setting the frequencies of the amino acids in a sequence
	int GetAACount(int AACount[],const long& LB, const long& HB, const bool& Reverse, const string& gid);
		
	
	//function for retrieving the Gene's sequence from the genome given the LowBase HighBase and Reverse status
	string GeneSequence(const string& GenomeSeq, const long& LB, const long& HB, const bool& Reverse);
        string GenomeSubseq(const bool& Reverse, const long& LB, const long& HB, const string& gid);
        string GetTrans(const bool& Reverse, const long& LB, const long& HB, const string& gid);
                
        int SetStarts(const string& filename);
        int SetStops();
        int SetupTrans(const int& TN, const string& TF);
        
        bool CheckStarts(); //checks to see if the reverse and forward starts have been initialized
        bool CheckStops();

	//Get the Entropy Distance Ratio
	double GetEDR(double EntropyProfile[]);
		
	//total up the EDPs from all the orfs
	//to come up with a new EDP for coding and noncoding
	int TotalEDP(int Counts[], const long& Length, const bool& Defeated);



//Add the EDP to the profile submitted so that it can be averaged
	int AddEDP(double TempEDP[], double Profile[]);
		

	//function to lower the case of all characters in a string
	int LowerTheCase(string & Seq);
		

	
	//This function is designed to check if submitted string is reverse codon
	bool ReverseStart(const string& Codon);//open definition
		
	//This function is designed to check if submitted string is forward start
	bool ForwardStart(const string& Codon);//open definition

	//This function is designed to check if submitted string a forward stop
	bool ForwardStop(const string& Codon);//open definition
		
	//This function is designed to check if submitted string a forward stop
	bool ReverseStop(const string& Codon);//open definition
		
	

	//This function is intended to give a likelihood of correctness relative to other starts
	double CalcSS(const string& Codon);//open definition
		



	//This function returns the complement of a nucleotide passed as a parameter
	char Complement(const char& Base);//open definition
		

	//This function gets the genome based on the genome file name
	int ReadGenome();
		

	//This function returns the reverse complement of a given sequence
	string ReverseComp(const string& Forward);
	


	//Calculate the RawBit score from sequence coordinates
	double CalcRawBit(const long& LowB, const long& HighB, const bool& Rev, const string& gid);
		
	//Function to select which genomic sequence to use based on its ID
	int SelectGenome(const string& gid);
        
        //Returns a unique, short id for a replicon
        string SmallGID(const string& gid);
	
	//Clear Genome free memory associated with genomes
	int ClearGenome();

	//This function continually adjusts the start site until its back at the original
	//assumes there is a query alignment offset to start at and an original start site to come back to
	bool FindStarts(long& St, const long& OSt, const long& Sp, const long& QAS, const bool& Reverse, double& StartScore, const string& gid);//open definition
        
        //Write out the genomes
        int WriteGenome(std::ostream& Out);
		


/*	//Function to Adjust the start site based on where the alignment begins
int AdjustStart1(long& St, const long& Sp, const long& QAS, const bool& Reverse, const string& Genome){//open definition
	long Start=St;
	long Stop=Sp;
	long QAlignStart=QAS;
	int StartSearch;
	double StartScore=0;
	double MaxStartScore=0;
	double Travel=0;
	int NuclDist=(3*QAlignStart);


	if(QAlignStart<2){
		return 1;
	}
	else if(Reverse){//if the pGene is reversed
		StartSearch=(Start-1)-NuclDist;// start search position

		for (int s=0; s<NuclDist; s+=3){//search codons in the upstream direction
			if(s%3==0){//if its the next codon
				if(ReverseStart(Genome.substr(StartSearch+s-2, 3))){//if its a start
					St=StartSearch+s+1;
					return 1;
					//}//close max start score
				}//close is start
			}//close next codon
		}//close search codons
	}//close reverse pGene

	else{//else its not reversed
		StartSearch=(Start-1)+NuclDist;

		for (int s=0; s<NuclDist; s+=3){//search 3 codons in the upstream direction
			if(s%3==0){//if its the next codon
				if(ForwardStart(Genome.substr(StartSearch-s, 3))){//if its a start
					St=StartSearch-s+1;
					return 1;
				}//close if start
			}//close if next codon
		}//close search next codons
	}//close not reversed
	return 1;
}//close defintion*/

};//close prototype

#endif
