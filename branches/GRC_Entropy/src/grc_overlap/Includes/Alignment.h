// Alignment.h


//Programmer: Andrew Warren
//email: anwarren@vt.edu
//Date of last modification: 02/xx/07




#ifndef Alignment_H
#define Alignment_H

#include <iostream>
#include <string>
#include <fstream>
#include <math.h>
#include <cmath>
#include <map>
#include <list>
#include <sstream>

using std::cout;
using std::cerr;
using std::string;
using std::ostream;
using std::ofstream;
using std::ifstream;
using std::multimap;
using std::map;
using std::list;
using std::stringstream;

class Alignment {//open prototype
public:


	long Start; // the start position for the orf
	long Stop; // the stop position for the orf
	long HighBase; //the highest base number for the orf
	long LowBase; //the lowest base number for the orf
	double Bit; // the bit score for the hit
	string EScore; //the E-score for the hit
	double EValue;
	bool Reverse; //Is it in a - frame
	long ALength; //the length of the alignment
	bool Defeated; //marks whether this record has been knocked out(tombstone method)
	long QAlignStart;//offset to start alignment from Start
	long QAlignStop;//offset to stop alignment from Start
	double MaxBit;//the maximum possible bit score
	double RelBit;//bit score adjusted by bit/maxbit
	double BitFrac;//Bit/MaxBit
	double StartScore;//frequency of occurence of start score
	double AlignScore;//Conditional Probability of(BitFrac, StartSite Codon Prob, MaxFrac)



	
	Alignment(){//default constructor
		Start=0;
		Stop=0;
		HighBase=0;
		LowBase=0;
		Bit=0;
		EScore="Not Assigned";
		EValue=100000;
		Reverse =false;
		Defeated=false;
		ALength=0;
		QAlignStart=0;
		QAlignStop=0;
		MaxBit=0;
		BitFrac=0;
		RelBit=0;
		AlignScore=0;
		StartScore=0;
	}

	//parameterized constructor
	Alignment(long St=0, long Sp=0, double B=0, string ES="none", long AL=0, long QASt=0, long QASp=0, double MxBit=0, double SScore=0){ // parameterized constructor1
		Start=St;
		Stop=Sp;
		Bit=B;
		EScore=ES;
		Defeated=false;
		ALength=AL;
		QAlignStart=QASt;//query align start offset
		QAlignStop=QASp;
		MaxBit=MxBit;
		RelBit=0;
		AlignScore=0;
		StartScore=SScore;



		if (Start>Stop){ 
			Reverse=true;//see if the orf is reversed
			HighBase=Start;
			LowBase=Stop;
			Stop=Stop-3;//adjust for stop codon
		}
		else{
			Reverse=false;
			HighBase=Stop;
			LowBase=Start;
			Stop=Stop+3;
		}

		if(Bit==0){EValue=100000;}
		else {
			stringstream ss;
			if(EScore[0]=='e'){
				EScore="1"+EScore;
			}
			ss.str(EScore);
			ss>>EValue;
			//EValue=(double(HLength))*(double(QLength))*(1.0/pow(2.0,Bit));
			BitFrac=Bit/MaxBit;
			RelBit=Bit*BitFrac;
		}//assign the E-value 
	}


		//Copy Constructor
	 Alignment(const Alignment &Source){// open defintion

		 Start=Source.Start; // the start position for the orf
		 Stop=Source.Stop; // the stop position for the orf
		 HighBase=Source.HighBase; 
		 LowBase=Source.LowBase;
		 Defeated=Source.Defeated;
		 RelBit=Source.RelBit;
		 MaxBit=Source.MaxBit;
		BitFrac=Source.BitFrac;
		 Bit=Source.Bit; // the bit score for the hit
		 EScore=Source.EScore; //the E-score for the hit
		 Reverse=Source.Reverse; //Is it in a - frame
		 EValue=Source.EValue;//copy the evalue for the hit
		ALength=Source.ALength;//alignment length
		QAlignStart=Source.QAlignStart;
		QAlignStop=Source.QAlignStop;
		AlignScore=Source.AlignScore;
		StartScore=Source.StartScore;
	}// close definition



	 //Assignment Operator
	 Alignment& operator =(const Alignment &Source){// open defintion
		if (this!= &Source){// open non-self assignment consq.

			Start=Source.Start; // the start position for the orf
			Stop=Source.Stop; // the stop position for the orf
			HighBase=Source.HighBase;
			RelBit=Source.RelBit;
			Defeated=Source.Defeated;
			LowBase=Source.LowBase;
			MaxBit=Source.MaxBit;
			BitFrac=Source.BitFrac;
			Bit=Source.Bit; // the bit score for the hit
			EScore=Source.EScore; //the E-score for the hit
			Reverse=Source.Reverse; //Is it in a - frame;	
			EValue=Source.EValue;
			ALength=Source.ALength;//alignment length
			QAlignStart=Source.QAlignStart;
			QAlignStop=Source.QAlignStop;
			AlignScore=Source.AlignScore;
			StartScore=Source.StartScore;
		}// close self assignment
		return *this;
	}// close definition


	//> OPERATOR overload
	//Returns wheter Alignment LHS>RHS
	//Based on whether the (start site, alignment) pair is the best
	//See AlignScore Calculation for more detail
	//bool operator >(const Alignment& RHS)const{
	//	AlignScore>RHS.AlignScore;
	//}

	//> OPERATOR overload
	//Returns wheter Alignment LHS<RHS
	//Based on whether the (start site, alignment) pair is the best
	//See AlignScore Calculation for more detail
	//bool operator <(const Alignment& RHS)const{
	//	AlignScore<RHS.AlignScore;
	//}

	//report the alignment score
	double ReportScore(){
		return AlignScore;
	}

	//function for retrieving information from alignment
	int GetInfo(long& ALen, double& BScore, string& EScr, double& EVal, long& HB, long& LB, double& MxBit, long& QAStart, long& QAStop, double& BFrac, double& RelB, long& Strt, long& Stp){
		ALen=ALength;
		BScore=Bit;
		EVal=EValue;
		HB=HighBase;
		LB=LowBase;
		MxBit=MaxBit;
		QAStart=QAlignStart;
		QAStop=QAlignStop;
		RelB=RelBit;
		BFrac=BitFrac;
		Strt=Start;
		Stp=Stop;
		return 0;
	}

	//this function calculates the score for this alignment/start pair
	//this score is conditional probability (Bit/MaxBit=proximity of start to alignment=AX, Bit/HighScore=RelativePerformance=RP, StartScore=StartFrequency=SF)
	//Bayes: AlignScore=Prob(Correct|SF,AX,RP)=alpha
	//alpha/(1-alpha)=(SF/(1-SF))*(AX/(1-AX))*(RP/(1-RP))
	int SetScore(const double& HighScore){
		double RP=Bit/HighScore;//relative performance to all other alignments
		double RPDiv=1-RP;
		double BFDiv=1-BitFrac;
		double SSDiv=1-StartScore;
		//NOTE: Criteria for most important Score factor is determined by setting the minimum for the divisor
		//Currently from most to least important (1)Relative Performance (2)alignment proximity (3) start site frequency
		if(RPDiv<=0.01){
			RPDiv=.01;
		}
		if(BFDiv<=.1){
			BFDiv=.1;
		}
		double Beta=(BitFrac/BFDiv)*(StartScore/(1-StartScore))*(RP/RPDiv);
		AlignScore=(Beta/(1+Beta));//set align score to conditional probability
		return 0;
	}
		
		
};//close prototype

struct OrderAlign {
	bool operator()(Alignment* S1, Alignment*S2){
		if(S1==NULL){
			return true;
		}
		else if(S2==NULL){
			return false;
		}
		else return ((S1->AlignScore)<(S2->AlignScore));
	}//close def.
	
};//close prototype

#endif
