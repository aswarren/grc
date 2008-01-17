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
	long Length;
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
	double RP;//relative performance of this alignment to all other alignments
	double EDR;//entropy distance ratio



	
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
		RP=0;
		EDR=9999;
	}

	//parameterized constructor
	Alignment(long St=0, long Sp=0, double B=0, string ES="none", long AL=0, long QASt=0, long QASp=0, double MxBit=0, double SScore=0, double EDRatio=9999){ // parameterized constructor1
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
		RP=0;
		EDR=EDRatio;
		if (Start>Stop){ 
			Reverse=true;//see if the orf is reversed
			HighBase=Start;
			LowBase=Stop+3;
			//Stop=Stop-3;//adjust for stop codon
		}
		else{
			Reverse=false;
			HighBase=Stop-3;
			LowBase=Start;
			//Stop=Stop+3;
		}
		Length=HighBase-LowBase+1;

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
		RP=Source.RP;
		Length=Source.Length;
		EDR=Source.EDR;
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
			RP=Source.RP;
			Length=Source.Length;
			EDR=Source.EDR;
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

	//report adjusted bit value
	double ReportRelBit(){
		return RelBit;
	}

	//report the alignment score
	double ReportScore(){
		return AlignScore;
	}

	//function for retrieving information from alignment
	int GetInfo(long& ALen, double& BScore, string& EScr, double& EVal, long& HB, long& LB, double& MxBit, long& QAStart, long& QAStop, double& BFrac, double& RelB, long& Strt, long& Stp){
		ALen=ALength;
		BScore=Bit;
		EVal=EValue;
		EScr=EScore;
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
		RP=Bit/HighScore;//relative performance to all other alignments
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


//Given two alignments this function determines the allowable overlap threshold
	double OverlapThreshold(Alignment& RHS, const double& MxOLap){
		double OLapThreshold=12;//default overlap threshold is 12 nucl.
		double BitFrac1, BitFrac2;
		BitFrac1=BitFrac;
		BitFrac2=RHS.BitFrac;
		double MaxOverlap=MxOLap;

		if(BitFrac1>1){
			BitFrac1=1;
		}
		if(BitFrac2>1){
			BitFrac2=1;
		}

		//double OLapT1=-1*log(EValue)*OEFactor;
		//double OLapT2=-1*log(RHS.EValue)*OEFactor;
		//double ScoreFactor=0;
		//double Diff=0;
		//BitFrac1=BitFrac1*Bit;
		//BitFrac2=BitFrac2*RHS.Bit;
		/*if(Length<300 || RHS.Length<300){ //if the length of one of the orfs is less than 300 decrease maxoverlap
			if(Length<RHS.Length){
				MaxOverlap=MaxOverlap*(Length/300);
				//MaxOverlap=0.4*QLength;
			}
			else {
				MaxOverlap=MaxOverlap*(RHS.Length/300);
				//MaxOverlap=0.4*RHS.QLength;
			}
		}*/
			
		if(BitFrac1<BitFrac2){
			//Diff=(Bit/RHS.Bit)*2;
			//ScoreFactor=(Bit/1000)+1;
			OLapThreshold=((BitFrac1+(1-EDR))*MaxOverlap);
		}
		else{
			//Diff=(RHS.Bit/Bit)*2;
			//ScoreFactor=(RHS.Bit/1000)+1;
			OLapThreshold=((BitFrac2+(1-RHS.EDR))*MaxOverlap);
		}

		return OLapThreshold;
		//OLapThreshold=BitFrac*90;//for every .10 they have similar score give a 1% overlap
		//if(OLapT2<OLapT1){//set overlap threshold to be a function of the e-value score
		//	OLapThreshold=OLapT2;
		//}
		//else {OLapThreshold=OLapT1;}//take the min threshold
	}


	//Function that reports the bit ratio of LHS.RelBit to RHS.RelBit
	//returns a negative bit ratio if LHS.RelBit < RHS.RelBit (meaning RHS is better)
	double BitRatio(Alignment& RHS){//open def
		double BRatio=0;
		if (RelBit==0 || RHS.RelBit==0){
			BRatio=0;
		}
		else if(RelBit<RHS.RelBit){
			BRatio=((RHS.RelBit-RelBit)/RelBit)*-1;
		}
		else {BRatio=(RelBit-RHS.RelBit)/RHS.RelBit;}
		return BRatio;
	}

		//Function that reports the ratio of the two ratios
	//returns a negative  ratio if LHS.EDR > RHS.EDR (meaning RHS is better)
	double EDRRatio(Alignment& RHS){//open def
		double ERatio=0;
		if (EDR==0 || RHS.RelBit==0){
			ERatio=0;
		}
		else if(EDR<RHS.EDR){
			ERatio=((RHS.EDR-EDR)/EDR);
		}
		else {ERatio=(EDR-RHS.EDR)/RHS.EDR*-1;}
		return ERatio;
	}


	int DisplayInfo(std::ostream& Out){
		Out<<Bit<<"\t";
		Out<<EScore<<"\t";
		return 0;
	}

	//function for updating the EDR of the alignment
	int UpdateEDR(const double& UpdatedValue){
		EDR=UpdatedValue;
		return 0;
	}

	//This function returns the value of the lower part of the orf coordinates
	long ReportLowBase(){
		return LowBase;
	}

	//Reports the Bit/MaxBit value for this alignment
	double ReportBitFrac(){
		return BitFrac;
	}

	//Gives value of Alignment Length
	double GetALength(){
		return ALength;
	}

	//Gives value of ORF length
	double GetLength(){
		return Length;
	}

	//This function returns the value of the higher part of the orf coordinates
	long ReportHighBase(){
		return HighBase;
	}
		
		
};//close prototype


//return S1<S2
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
