//CalcPack.h
//The following is a container class
//holds needed values to calculate various scores throughout the program


#ifndef CalcPack_H
#define CalcPack_H

#include <iostream>
#include <string>
#include <fstream>
#include <map>
#include <sstream>

using std::cout;
using std::cerr;
using std::string;
using std::ostream;
using std::ofstream;
using std::ifstream;
using std::map;
using std::stringstream;

class CalcPack{//open prototype
public:
	double Lambda;//constant in MaxBit calculation
	double K;//constant in MaxBit calculation
	map <string,int> ConsValue;//map for storing the max value of a conserved amino acid
	string Genome; //for storing the genome
	string ECommand;//the command for running entropy-calc
	string Matrix;
	string GenomeFile;


	CalcPack(){//default constructor
		Lambda=0;
		K=0;
		Genome="";
		ECommand="";
		GenomeFile="";
	}

	CalcPack(string Matx, string GF, string ECom){//parameterized constructor
		ECommand=ECom;
		Matrix=Matx;
		int Status=InitCodes();//read in the values.
		if(Status!=0){
			cout<<"\ngrc_overlap unable to initialize BLAST parameters:exiting...\n";
			throw 20;
		}
		GenomeFile=GF;
		GetGenome();
	}

	//Copy Constructor
	//Should not be USED!! pack contains GENOME
	 CalcPack(const CalcPack &Source){// open defintion
		cerr<<"Error Trying to copy CalcPack object";
		throw 15;//throw exception
	}

	//Assignment operator
	CalcPack& operator =(const CalcPack &Source){// open defintion
		if (this!= &Source){// open non-self assignment consq.
			cerr<<"Error trying to copy calcpack object";
			throw 15;//throw exception
		}
		return *this;
	}


	//Function to initialize parameters based on blast matrix used
	int InitCodes(){//open definition
	//Get the Genetic Codes for translation
		ifstream InCode;
		InCode.open(Matrix.c_str());//open matrix specified
		if(!InCode.is_open()){
			cout<<"Error opening MaxMatrix file";
			return 1;
		}
		string Delim;
		string OpCase="";
		int AA; //Max value for conservation of the amino acid
		//InCode>>Delim;

		while(InCode){
			InCode>>Delim;
			if(isupper(Delim[0])){//create lower case in map
				OpCase=tolower(Delim[0]);
				for(int t=1; t<Delim.size(); t++){
					OpCase+=tolower(Delim[t]);
				}
			}
			else {
				OpCase=Delim;
			}

			InCode>>AA; //Read in the AA conservation value
			//CodeArray[TableNum].insert(map<string,char>::value_type(Delim,AA));
			ConsValue.insert(map<string,int>::value_type(OpCase,AA));
			
		}//close while loop
		InCode.close();

		//Set constants
		if(string::npos!=Matrix.find("MaxB62",0)){//if the matrix is BLOSUM62
			Lambda=0.267;
			K=0.041;
		}

		else if(string::npos!=Matrix.find("MaxB80",0)){//if the matrix is BLOSUM80
			Lambda=0.299;
			K=0.071;
		}

		else if(string::npos!=Matrix.find("MaxB45",0)){//if the matrix is BLOSUM45
			Lambda=0.195;
			K=0.032;
		}

		else if(string::npos!=Matrix.find("MaxP30",0)){//if the matrix is PAM30
			Lambda=0.309;
			K=0.15;
		}

		else if(string::npos!=Matrix.find("MaxP70",0)){//if the matrix is PAM30
			Lambda=0.270;
			K=0.060;
		}

		else{//default B62
			Lambda=0.267;
			K=0.041;
		}


		return 0;
	}//close definition


	//This function calculates the entropy for the section of the genome specified
	double GetEntropy(const long& LB, const long& HB, const bool& Reverse){//open definition
		long Length=HB-LB+1;
		long Begin=LB-1;
		FILE* TempF;
		char TempC[sizeof(double)];
		double Entropy=-1;
		string Command=ECommand;
		
		Command=Command+" ";//add space before sequence
		if(Reverse){//if the pGene is reversed
			Command=Command+ReverseComp(Genome.substr(Begin, Length));	
		}
		else{//not reverse complement
			Command=Command+Genome.substr(Begin, Length);
		}
		
	
		if ( ( TempF = popen ( Command.c_str(), "r") ) != NULL ){
			fgets (TempC, sizeof(double), TempF);
			pclose(TempF); 
			stringstream StoD;
			StoD<<TempC;
			StoD>>Entropy;
			return Entropy;
		}
		else{
			cerr<<"Could not run command with popen in GetEntropy.";
			return Entropy;
		}
	}//close function

	//function to lower the case of all characters in a string
	int LowerTheCase(string & Seq){
		for(int i=0; i<Seq.length(); i++){
			Seq[i]=tolower(Seq[i]);
		}
		return 1;
	}//close definition
		

	
	//This function is designed to check if submitted string is reverse codon
	bool ReverseStart(const string& Codon){//open definition
		if(Codon=="cat"||Codon=="cac"||Codon=="caa"){
			return true;
		}
		else return false;
	}
	


	//This function is designed to check if submitted string is reverse codon
	bool ForwardStart(const string& Codon){//open definition
		if(Codon=="atg"||Codon=="gtg"||Codon=="ttg"){
			return true;
		}
		else return false;
	}
	


	//This function is intended to give a likelihood of correctness relative to other starts
	double CalcSS(const string& Codon){//open definition
		if (Codon=="atg"||Codon=="cat"){
			return (.77);
		}
		else if (Codon=="gtg"||Codon=="cac"){
			return (.14);
		}
		else if (Codon=="ttg"||Codon=="caa"){
			return (.08);
		}
		return 0;
	}



	//This function returns the complement of a nucleotide passed as a parameter
	char Complement(const char& Base){//open definition
		switch(Base){
			case 'a':
				return 't';
				
			case 't':
				return 'a';
				
			case 'c':
				return 'g';
				
			case 'g':
				return 'c';
			default: return '*';
		}
	}//close definition

	//This function gets the genome based on the genome file name
	int GetGenome(){
		ifstream In2;//ofstream operator for reading in the genomic sequence
		In2.open(GenomeFile.c_str());//open up the translated file
		
		string GenomeID;//for reading in the id
		getline(In2,GenomeID); //get line for description
		string Seq;//for reading in the sequence
		Genome="";//Initialize to empty string
	
		while(In2){//read in the genome file
	
			In2>>Seq; //read in the genomic sequence
			LowerTheCase(Seq);
			Genome+=Seq;//concat. each line to Genome
	
			/*FindIt=HitList.find(ID.substr(1,(ID.length())-1));
	
			if(FindIt!=HitList.end()){//if its found then its a hit
				FindIt->second->Sequence=Seq;//assign the sequence
			}*/
		}
		In2.close();//close the input stream
	}

	//This function returns the reverse complement of a given sequence
	string ReverseComp(const string& Forward){
	
		string Comp="";
		for(int s= int(Forward.length())-1; s>=0; s--){
			Comp+=Complement(Forward[s]);
		}
		return Comp;
	}//close definition

};//close prototype

#endif