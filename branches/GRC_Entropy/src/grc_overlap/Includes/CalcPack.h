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
#include "Translate.h"

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
	string Genome; //for storing the genome
	string TransFile;//the command for running entropy-calc
	string Matrix;
	string GenomeFile;
	int OStartCount;//Count for find starts that keeps track of the number of times original start has been processed
	double DefaultCProfile[20];
	double DefaultNCProfile[20];//for non coding genes
	Translate Translator;//for translating nucl. to AA


	CalcPack(){//default constructor
		Lambda=0;
		K=0;
		Genome="";
		TransFile="";
		GenomeFile="";
		OStartCount=0;
	}

	CalcPack(string Matx, string GF, string TF){//parameterized constructor
		TransFile=TF;
		Matrix=Matx;
		OStartCount=0;
		int Status=InitCodes();//read in the values.
		if(Status!=0){
			cout<<"\ngrc_overlap unable to initialize BLAST parameters:exiting...\n";
			throw 20;
		}
		SetDefaultEDP();
		Translator.SetTransFile(TransFile);
		Translator.InitCodes();
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

		//Set constants Values obtained from FSABLAST
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


	int SetDefaultEDP(){
	//Annoying
	//no list initialization of non-static data member....
		DefaultCProfile[0]=0.08468;
		DefaultCProfile[1]=0.01606;
		DefaultCProfile[2]=0.05739;
		DefaultCProfile[3]=0.05752;
		DefaultCProfile[4]=0.04328;
		DefaultCProfile[5]=0.07042;
		DefaultCProfile[6]=0.02942;
		DefaultCProfile[7]=0.05624;
		DefaultCProfile[8]=0.04442;
		DefaultCProfile[9]=0.05620;
		DefaultCProfile[10]=0.03029;
		DefaultCProfile[11]=0.03975;
		DefaultCProfile[12]=0.05116;
		DefaultCProfile[13]=0.04098;
		DefaultCProfile[14]=0.05989;
		DefaultCProfile[15]=0.08224;
		DefaultCProfile[16]=0.05660;
		DefaultCProfile[17]=0.06991;
		DefaultCProfile[18]=0.02044;
		DefaultCProfile[19]=0.03310;

		//setup Default noncoding profile
		DefaultNCProfile[0]=0.07434;
		DefaultNCProfile[1]=0.03035;
		DefaultNCProfile[2]=0.05936;
		DefaultNCProfile[3]=0.04729;
		DefaultNCProfile[4]=0.05662;
		DefaultNCProfile[5]=0.07704;
		DefaultNCProfile[6]=0.05777;
		DefaultNCProfile[7]=0.05328;
		DefaultNCProfile[8]=0.03360;
		DefaultNCProfile[9]=0.05581;
		DefaultNCProfile[10]=0.01457;
		DefaultNCProfile[11]=0.03718;
		DefaultNCProfile[12]=0.04594;
		DefaultNCProfile[13]=0.05977;
		DefaultNCProfile[14]=0.08489;
		DefaultNCProfile[15]=0.05990;
		DefaultNCProfile[16]=0.04978;
		DefaultNCProfile[17]=0.07227;
		DefaultNCProfile[18]=0.01050;
		DefaultNCProfile[19]=0.01974;
		return 0;
	}


	//this function returns the coordinate of the amino acid in the array
	//expects that any non-AA coding seq. enountered will be translated to '*'
	int MapAA(char AA){
		int Value=21;
		if(AA=='*'){
			return Value;
		}
		else{
			Value=AA-'A';
		}
		//26 letters in the alphabet but only 20 AA convert to correct indecies
		if(Value>23){
			Value-=5;
		}
		else if(Value>20){
			Value-=4;
		}
		else if(Value>14){
			Value-=3;
		}
		else if(Value>9){
			Value-=2;
		}
		else if(Value>1){
			Value-=1;
		}
		
		return Value;
	}



	//This function calculates the entropy for the section of the genome specified
	double GetEntropy(int AACount[]){//open definition

		double Entropy=-1;
		int t=0;
		double TotalEntropy=0;
		int TotalCount=0;

		double TempFreq[21];

		for(int x=0; x<21; x++){
			TempFreq[x]=0;
		}

		for(t=0; t<20; t++){
			TotalCount+=AACount[t];
		}
		//convert counts to frequencies
		for(t=0; t<20; t++){
			TempFreq[t]=double(AACount[t])/(double(TotalCount));
		}
		//convert frequencies to entropies
		for(t=0; t<20; t++){
			if(TempFreq[t]!=0){
				TempFreq[t]=-1*TempFreq[t]*log10(TempFreq[t]);
				TotalEntropy+=TempFreq[t];
			}
		}
		for(t=0;t<20;t++){
			TempFreq[t]=TempFreq[t]/TotalEntropy;
		}

		Entropy=GetEDR(TempFreq);
	
		/*if ( ( TempF = popen ( Command.c_str(), "r") ) != NULL ){
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
		}*/
		return Entropy;
	}//close function


	//function for finding setting the frequencies of the amino acids in a sequence
	int GetAACount(int AACount[],const long& LB, const long& HB, const bool& Reverse){
		long Length=HB-LB+1;
		long Begin=LB-1;
		string Translation="";


		if(Reverse){//if the pGene is reversed
			Translation=Translator.TranslateSeq(ReverseComp(Genome.substr(Begin,Length)));
			//Command=Command+ReverseComp(Genome.substr(Begin, Length));	
		}
		else{//not reverse complement
			Translation=Translator.TranslateSeq(Genome.substr(Begin, Length));
			//Command=Command+Genome.substr(Begin, Length);
		}
		
		//for each amino acid in the sequence
		for (int t=0; t<Translation.size(); t++){
			AACount[MapAA(Translation[t])]++;
		}
		return 0;
	}


	//Get the Entropy Distance Ratio
	double GetEDR(double EntropyProfile[]){
		int x=0;
		double CodingDist=0;
		double NonCodingDist=0;
		for (x=0; x<20; x++){
			double CodingDif=EntropyProfile[x]-DefaultCProfile[x];
			CodingDist+=(CodingDif*CodingDif);
			double NonCodingDif=EntropyProfile[x]-DefaultNCProfile[x];
			NonCodingDist+=(NonCodingDif*NonCodingDif);
		}
		CodingDist=sqrt(CodingDist);
		NonCodingDist=sqrt(NonCodingDist);
		return CodingDist/NonCodingDist;
	}


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


	//Calculate the RawBit score from sequence coordinates
	double CalcRawBit(const long& LowB, const long& HighB, const bool& Rev){
		long LB=LowB;
		long HB=HighB;
		long StartSearch=LB-1;//Subtract one to convert to string coordinates
		long Length=HB-LB+1;
		string Codon;
		double RawBit;
		map <string,int>::iterator FindIt;		
	//even though its in reverse go forward through the sequence and get reverse complement
		if(Rev){
			for (long s=0; s<Length; s=s+3){//calc max possible score
				if(s%3==0){//if its the next codon
					Codon=ReverseComp(Genome.substr(StartSearch+s, 3));//get reverse complement
					FindIt=ConsValue.find(Codon);//find the score of this codon
					if(FindIt!=ConsValue.end()){
						RawBit=RawBit+FindIt->second;//add up score
					}
				}
			}//close max loop
		}//close if Reverse

		else {
			for (long s=0; s<Length; s=s+3){//calc max possible score
				if(s%3==0){//if its the next codon
					Codon=Genome.substr(StartSearch+s, 3);//codon
					FindIt=ConsValue.find(Codon);//find the score of this codon
					if(FindIt!=ConsValue.end()){
						RawBit=RawBit+(FindIt->second);//add up score
					}
				}
			}//close max loop
		}//close not Reverse
		return RawBit;
	}//close definition



	//This function continually adjusts the start site until its back at the original
	bool FindStarts(long& St, const long& OSt, const long& Sp, const long& QAS, const bool& Reverse, double& StartScore){//open definition
		long Start=St;
		long Stop=Sp;
		long OrigStart=OSt;
		long QAlignStart=QAS;
		int StartSearch=0;
		double MaxStartScore=0;
		double Travel=0;
		double NuclDist=(3*QAlignStart);
		int Halt=0;//defines when the search has gone back to original

		//if there is no room to search for a start between the aligned region and current start
		//increase orig start count and return
		//if have returned to Original start (OStartCount>0) then return false
		if(Start==OrigStart){
			if(OStartCount>0){	
				OStartCount=0;//reset for next use
				return false;
			}
			else if(QAlignStart<2){//if there isn't any room for finding another start return true this time and false next time
				OStartCount++;
				if(Reverse){
					StartScore=CalcSS(Genome.substr(Start-3, 3));//get the probability of the codon being start site
				}
				else{
					StartScore=CalcSS(Genome.substr(Start-1, 3));//get the probability of the codon being start site
				}
				return true;
			} 
		}
		if(Reverse){//if the pGene is reversed LOOKING FOR reverse starts in THE FORWARD DIRECTION!!!
			if(Start==OrigStart){
				OStartCount++;//increase count that originalstart has been processed
				StartSearch=(Start-3-1)-(3*QAlignStart);// start search position //minus 3 to account for subtraction difference of(+1*3) and -1 to convert to string coordinates
				Halt=(QAlignStart*3)+3;
			}
			else{//start from the position of the last start found
				StartSearch=Start+2;//next codon 
				Halt=OrigStart-Start-3;//room left between orig and current search position
			}
	
			for (int s=0; s<=Halt; s=s+3){//search codons in the upstream direction
				if(s%3==0){//if its the next codon
					if(ReverseStart(Genome.substr(StartSearch+s-2, 3))){//if its a start -2 because looking in forward direction
						StartScore=CalcSS(Genome.substr(StartSearch+s-2, 3));//get the probability of the codon being start site
						St=StartSearch+s+1;
						return true;//if back to original start, stop searching
						//}//close max start score
					}//close is start
				}//close next codon
			}//close search codons
		}//close reverse pGene
	
		else{//else its not reversed
			if(Start==OrigStart){
				OStartCount++;//increase the count that original start has been processed
				StartSearch=(Start-1)+(3*QAlignStart);
				Halt=(QAlignStart*3);
			}
			else{
				StartSearch=Start-4;//start at the next codon
				Halt=Start-OrigStart-3;//room left between orig and current search position
			}
	
			for (int s=0; s<=Halt; s=s+3){//search 3 codons in the upstream direction
				if(s%3==0){//if its the next codon
					if(ForwardStart(Genome.substr(StartSearch-s, 3))){//if its a start
						StartScore=CalcSS(Genome.substr(StartSearch-s, 3));//get the probability of the codon being start site
						St=StartSearch-s+1;
						return true;
						//}//close if max score
					}//close if start
				}//close if next codon
			}//close search next codons
		}//close not reversed
		return false;
	}//close defintion


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
