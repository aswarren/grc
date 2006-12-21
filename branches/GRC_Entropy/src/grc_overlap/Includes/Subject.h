// Subject.h


//Programmer: Andrew Warren
//email: anwarren@vt.edu
//Date of last modification: 10/xx/06




#ifndef Subject_H
#define Subject_H

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

class Subject {//open prototype
public:
	int ID;
	string Hit;
	long Start; // the start position for the orf
	long Stop; // the stop position for the orf
	long HighBase; //the highest base number for the orf
	long LowBase; //the lowest base number for the orf
	double Bit; // the bit score for the hit
	string EScore; //the E-score for the hit
	double EValue;
	long HLength; //the length of the hit sequence
	bool Reverse; //Is it in a - frame
	bool Blank; //indicates whether the was a hit to this orf in the DB
	long QLength; //the length of the orf
	long ALength; //the length of the alignment
	bool Defeated; //marks whether this record has been knocked out(tombstone method)
	long QAlignStart;//offset to start alignment from Start
	long QAlignStop;//offset to stop alignment from Start
	double MaxBit;//the maximum possible bit score
	double RelBit;//bit score adjusted by bit/maxbit
	bool Hypot;//bool to tell whether description contains hypothetical
	string HitID;//id of the hit in db
	string HitOrg;//name of the organism in the db hit


	
	Subject(){//default constructor
		ID=-1;
		Hit="unassigned";
		Start=0;
		Stop=0;
		HighBase=0;
		LowBase=0;
		Bit=0;
		EScore="Not Assigned";
		EValue=100000;
		HLength=0;
		Reverse =false;
		Blank =true;
		Defeated=false;
		QLength=0;
		ALength=0;
		QAlignStart=0;
		QAlignStop=0;
		MaxBit=0;
		RelBit=0;
		HitID="none";
		HitOrg="none";
	}

	//parameterized constructor
	Subject(int I=1, long St=0, long Sp=0, string H="None", double B=0, string ES="none", long L=0, long A=0, long QASt=0, long QASp=0, double MxBit=0, string HID="none", string HOrg="none"){ // parameterized constructor1
		ID=I;
		Start=St;
		Stop=Sp;
		Hit=H;
		Bit=B;
		EScore=ES;
		HLength=L;
		Defeated=false;
		ALength=A;
		QAlignStart=QASt;//query align start offset
		QAlignStop=QASp;
		MaxBit=MxBit;
		RelBit=0;
		HitID=HID;
		HitOrg=HOrg;

		Blank=(B==0); //if the bit score is 0 then blank is true
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
		Hypot=(Hit.npos!=Hit.find("hypothetical"));
		QLength=labs(Start-Stop)+1;
		if(Bit==0){EValue=100000;}
		else {
			stringstream ss;
			if(EScore[0]=='e'){
				EScore="1"+EScore;
			}
			ss.str(EScore);
			ss>>EValue;
			//EValue=(double(HLength))*(double(QLength))*(1.0/pow(2.0,Bit));
			RelBit=(Bit*Bit)/MaxBit;
		}//assign the E-value 
	}


		//Copy Constructor
	 Subject(const Subject &Source){// open defintion
		 ID=Source.ID;
		 Hit=Source.Hit;
		 Start=Source.Start; // the start position for the orf
		 Stop=Source.Stop; // the stop position for the orf
		 HighBase=Source.HighBase; 
		 LowBase=Source.LowBase;
		 Defeated=Source.Defeated;
		 RelBit=Source.RelBit;
		 Hypot=Source.Hypot;
		 MaxBit=Source.MaxBit;
		 Bit=Source.Bit; // the bit score for the hit
		 EScore=Source.EScore; //the E-score for the hit
		 HLength=Source.HLength; //the length of the hit sequence
		 Reverse=Source.Reverse; //Is it in a - frame
		 Blank=Source.Blank;
		 EValue=Source.EValue;//copy the evalue for the hit
		ALength=Source.ALength;//alignment length
		QLength=Source.QLength;//Orf length
		QAlignStart=Source.QAlignStart;
		QAlignStop=Source.QAlignStop;
		HitOrg=Source.HitOrg;
		HitID=Source.HitID;
	}// close definition

	 	//> OPERATOR overload
		//adding e-score evaluation and inverting signs 04/04/06
	//bool operator>(const Subject& RHS)const{
	//	return(RelBit>RHS.RelBit);
		//return(EValue<RHS.EValue);
	//}

		//< OPERATOR overload
	//bool operator<(const Subject& RHS)const{

	//	return (RelBit<RHS.RelBit);
		//return(EValue>RHS.EValue);
	//}


	 //Assignment Operator
	 Subject& operator =(const Subject &Source){// open defintion
		if (this!= &Source){// open non-self assignment consq.
			ID=Source.ID;
			Hit=Source.Hit;
			Start=Source.Start; // the start position for the orf
			Stop=Source.Stop; // the stop position for the orf
			HighBase=Source.HighBase;
			RelBit=Source.RelBit;
			Hypot=Source.Hypot;
			Defeated=Source.Defeated;
			LowBase=Source.LowBase;
			MaxBit=Source.MaxBit;
			Bit=Source.Bit; // the bit score for the hit
			EScore=Source.EScore; //the E-score for the hit
			HLength=Source.HLength; //the length of the hit sequence
			Reverse=Source.Reverse; //Is it in a - frame
			Blank=Source.Blank;	
			EValue=Source.EValue;
			QLength=Source.QLength;//Orf length
			ALength=Source.ALength;//alignment length
			QAlignStart=Source.QAlignStart;
			QAlignStop=Source.QAlignStop;
			HitOrg=Source.HitOrg;
			HitID=Source.HitID;
		}// close self assignment
		return *this;
	}// close definition




	//bool Overlap
	//returns the length of there being an overlap of ORFs
	 int Overlap(Subject& RHS){//open def
		int OverLen=0;
		 
		if (RHS.LowBase>=LowBase && RHS.LowBase <=HighBase){
			if(HighBase>=RHS.HighBase){OverLen=RHS.QLength;}//if one frame encompasses the other
			else OverLen=HighBase-RHS.LowBase+1;
		}
		else if(RHS.LowBase <= LowBase && RHS.HighBase >= LowBase){
			 if(RHS.HighBase>=HighBase){OverLen=QLength;}//if one frame encompasses the other
			 else OverLen=RHS.HighBase-LowBase+1;
		}
		return OverLen;
	 }// close defintion





		
		


	bool operator ^(const Subject& RHS)const{//open definition
		return (QLength>RHS.QLength);
	}// close definition

	bool ReportHit() const{//open defintion
		return !Blank;
	}//close definiton

	//Reporter function that tells whether a record has been knocked out
	bool Dead(){return Defeated;}

	//Tombstones this record
	int KnockOut(){
		Defeated=true;
		return 0;
	}

	//string ReportAcc ()const{//open definition
	//	return Accession;
	//}// close definition

//	string ReportSeq ()const{//open definition
//		return Sequence;
//	}// close definition

};//close prototype

struct MoreBit {
	bool operator()(Subject* S1, Subject*S2){
		if(S1==NULL){
			return true;
		}
		else if(S2==NULL){
			return false;
		}
		else return ((S1->RelBit)<(S2->RelBit));
	}//close def.
	
};//close prototype

#endif