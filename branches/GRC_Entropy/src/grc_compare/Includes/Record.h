// Record.h


//Programmer: Andrew Warren
//email: anwarren@vt.edu
//Date of last modification: 04/xx/06      




#ifndef Record_H
#define Record_H

#include <iostream>
#include <string>
#include <fstream>
#include <math.h>
#include <cmath>
#include <list>
#include <iomanip>
#include <sstream>
#include <set>
#include <map>
#include "DirectHash.h"
#include "GO.h"
		
using std::setw;
using std::list;
using std::cerr;
using std::string;
using std::ostream;
using std::ofstream;
using std::ifstream;
using std::cout;
using std::stringstream;
using std::set;
using std::map;

	enum compete{WIN, LOSE, TIE, NONE};//keeps track of who won what
	enum result{TP, FP, TN, FN, NA, NRP, NRN};//keeps track status with respect to reference
	typedef set<string> StringSet;
	typedef map<int,StringSet> FunctionMap;

	//class Match;//foward declaration

class Record {//open prototype
	friend std::ostream& operator<<(std::ostream& ACOut, const Record& AC);


//private:
public:
	string Sequence;
	string ID; //unique for each record
	char Strand;//indicate the direction of the orf
	list<string> Description;
	string Hit;
	long Start; // the start position for the orf
	long Stop; // the stop position for the orf
	double Bit; // the bit score for the hit
	string EScore; //the E-score for the hit
	long HLength; //the length of the hit sequence
	bool Reverse; //Is it in a - frame
	bool Blank;//keeps track of whether a pGene has a hit in the DB 
	long OLength; //the length of the orf
	//long ALength; //the length of the alignment
	double QPercent;//from blast the % of the query that is aligned
	double HPercent;//from blast the % of the hsp that is aligned
	double MaxOLap; //the maximum sequence overlap % of this annotated gene from one set of annotations with the other set of annotations
	double MaxTerms;
	double EValue;
	//static list<string> FalseHit;
	bool Ref; //boolean variable to keep track of whether it is a reference record
	string DBID;
	string DBOrg;
	bool RefMatched;//boolean to keep track if the reference has been used in a match aka it keeps track of whether this reference gene has been predicted by any of the orfs generated
	FunctionMap GOTerms;//GOTerms,EvidenceCodes assigned to this prediction
	bool HasGO;
	double Entropy;//entropy distance ratio
	
//public:
	result Evaluation;//is this orf a TP,FP,FN,TN??
	Record* RefEval;//pointer so that Reference can be present in the FNAnalysis
	//this pointer can and should be removed in favor of a better designed Match class
	
	Record(){//default constructor
		Sequence="blah";
		ID="unassigned";
		//Hit="unassigned";
		Start=0;
		Stop=0;
		Bit=0;
		EScore="Not Assigned";
		EValue=999999;
		HLength=0;
		
		Reverse =false;
		HasGO=false;
		Blank =true;
		OLength=0;
		QPercent=0;
		Evaluation=NA;
		HPercent=0;
		MaxOLap=0;
		MaxTerms=0;
		Ref=false;
		Strand='#';
		RefEval=NULL;
		DBID="none";
		DBOrg="none";
		RefMatched=false;
		Entropy=9999;
	}


	//parameterized constructor
	Record(string TID="unassigned", long St=0, long Sp=0, string H="No_hits", double Ent=9999, bool R=false, double B=0, string ES="none", long L=0, double QP=0, double HP=0, string DID= "none", string DOrg="none"){ // parameterized constructor1
		ID=TID;
		Start=St;
		Stop=Sp;
		Hit=H;
		Bit=B;
		EScore=ES;
		HLength=L;
		Evaluation=NA;
		QPercent=QP;
		HPercent=HP;
		DBID=DID;
		DBOrg=DOrg;
		OLength=labs(Start-Stop)+1;
		RefMatched=false;
		HasGO=false;
		Entropy=Ent;//set entropy
		if(Bit==0){
			EValue=999999;
		}
		else{
			if(EScore[0]=='e'){
				EScore="1"+EScore;
			}
			stringstream TempSS(EScore);
			TempSS>>EValue;
		}
		ParseHit();//Parse the Hit description
		Ref=R;
		MaxOLap=0;
		RefEval=NULL;
		MaxTerms=0;
		Blank=(Bit<=0); //if the bit score is 0 then blank is true
		if (Start>Stop){
			Reverse=true;//see if the orf is reversed
			Strand='-';
		}
		else{
			Reverse=false;
			Strand='+';
		}
		//Taken out since this code was added to grc_overlap
		/*if(!Ref){
			if(Reverse){
				Stop=Stop-3;
			}
			else Stop=Stop+3;
		}*/
		
		
	}


		//Copy Constructor
	 Record(const Record &Source){// open defintion
		 Sequence=Source.Sequence;
		 Strand=Source.Strand;
		 ID=Source.ID; //unique for each record
		 Hit=Source.Hit;
		 Description=Source.Description;
		 Start=Source.Start; // the start position for the orf
		 Stop=Source.Stop; // the stop position for the orf
		 Bit=Source.Bit; // the bit score for the hit
		 EScore=Source.EScore; //the E-score for the hit
		 HLength=Source.HLength; //the length of the hit sequence
		 Reverse=Source.Reverse; //Is it in a - frame
		 Blank=Source.Blank;
		 EValue=Source.EValue;
		 //OLapOrfs=Source.OLapOrfs;//get overlapping orfs
		QPercent=Source.QPercent;
		Evaluation=Source.Evaluation;
		MaxOLap=Source.MaxOLap;//the % of max overlap for the other set of annotations
		MaxTerms=Source.MaxTerms;//The number of annotation terms associated with max overlap
		HPercent=Source.HPercent; //percent the blast 
		OLength=Source.OLength;//Orf length
		Ref=Source.Ref; //reference record indicator
		RefEval=Source.RefEval;
		DBID=Source.DBID;
		DBOrg=Source.DBOrg;
		RefMatched=Source.RefMatched;
		GOTerms=Source.GOTerms;
		HasGO=Source.HasGO;
		Entropy=Source.Entropy;
		
	}// close definition

	 	//> OPERATOR overload
	bool operator>(const Record& RHS)const{
		return(Bit>RHS.Bit);
	}

		//< OPERATOR overload
	bool operator<(const Record& RHS)const{
		return(Bit<RHS.Bit);
	}
	
	 //Assignment Operator
	 Record& operator =(const Record &Source){// open defintion
		if (this!= &Source){// open non-self assignment consq.
		 Sequence=Source.Sequence;
		 ID=Source.ID; //unique for each record
		 Hit=Source.Hit;
		 Description=Source.Description;
		 Strand=Source.Strand;
		 Start=Source.Start; // the start position for the orf
		 Stop=Source.Stop; // the stop position for the orf
		 EValue=Source.EValue;
		 Bit=Source.Bit; // the bit score for the hit
		 EScore=Source.EScore; //the E-score for the hit
		 HLength=Source.HLength; //the length of the hit sequence
		 Evaluation=Source.Evaluation;
		 Reverse=Source.Reverse; //Is it in a - frame
		 Blank=Source.Blank;
		 //OLapOrfs=Source.OLapOrfs;
		MaxOLap=Source.MaxOLap; //maximum overlap
		MaxTerms=Source.MaxTerms;//The number of annotation terms associated with max overlap
		OLength=Source.OLength;//Orf length
		QPercent=Source.QPercent;
		HPercent=Source.HPercent;
		Ref=Source.Ref;
		RefEval=Source.RefEval;
		DBID=Source.DBID;
		DBOrg=Source.DBOrg;
		RefMatched=Source.RefMatched;
		GOTerms=Source.GOTerms;
		HasGO=Source.HasGO;
		Entropy=Source.Entropy;
		}// close self assignment
		return *this;
	}// close definition

	//bool Overlap
	//returns the bool status of there being an overlap of ORFs
	//RHS is Setubal/reference annotation
	 bool Overlap(Record& RHS, double& OverLen){//open def
		bool OLap=false;
		//double OverLen=0;
		long LStart=0;
		long LStop=0;
		long RStart=0;
		long RStop=0;
		if (Reverse){LStart=Stop; LStop=Start;} //if the LHS is reverse frame
		else {LStart=Start;  LStop=Stop;} //if LHS is not
		if (RHS.Reverse){RStart= RHS.Stop; RStop =RHS.Start;} //if RHS is reverse
		else {RStart=RHS.Start; RStop=RHS.Stop;}

		 
		if (RStart>=LStart && RStart <=LStop){
			 if(LStop>RStop){OverLen=RHS.OLength;}//if one frame encompasses the other
			 else OverLen=LStop-RStart+1;
			OLap= true;}
		else if(RStart <= LStart && RStop >= LStart){
			 if(RStop>LStop){OverLen=OLength;}//if one frame encompasses the other
			 else OverLen=RStop-LStart+1;
			 OLap= true;}
		else {
			OLap= false;}
		 

		/*if (OLap){ //if the orfs overlap
			
			RefOLapPercent=OverLen/RHS.OLength; //setubal annotation overlap
			GRCOLapPercent=OverLen/OLength; //grc annotation overlap
			//adjust maxOverlap
			double CombinedOLap=RefOLapPercent*GRCOLapPercent; //keep track of combined overlap
			if(RHS.MaxOLap<CombinedOLap){//if the maximum overlap for the reference orf is less than the current update it
				RHS.MaxOLap=CombinedOLap;
			}
			if(MaxOLap<CombinedOLap){//if the overlap with respect to the ref. set orf is the greatest so far update it
				MaxOLap=CombinedOLap;
				
			}
		}
		else { RefOLapPercent=0; GRCOLapPercent=0;}*/
		return OLap;
	 }// close defintion


	//bool JustOverlap
	//returns the bool status of there being an overlap of ORFs
	//RHS is Setubal/reference annotation
	 bool JustOverlap(Record& RHS, bool& SameDirec, double& OverLen){//open def
		bool OLap=false;
		OverLen=0;
		long LStart=0;
		long LStop=0;
		long RStart=0;
		long RStop=0;
		SameDirec=(Reverse==RHS.Reverse);//set the same direction
		if (Reverse){LStart=Stop; LStop=Start;} //if the LHS is reverse frame
		else {LStart=Start;  LStop=Stop;} //if LHS is not
		if (RHS.Reverse){RStart= RHS.Stop; RStop =RHS.Start;} //if RHS is reverse
		else {RStart=RHS.Start; RStop=RHS.Stop;}

		 
		if (RStart>=LStart && RStart <=LStop){
			 if(LStop>RStop){OverLen=RHS.OLength;}//if one frame encompasses the other
			 else OverLen=LStop-RStart+1;
			OLap= true;}
		else if(RStart <= LStart && RStop >= LStart){
			 if(RStop>LStop){OverLen=OLength;}//if one frame encompasses the other
			 else OverLen=RStop-LStart+1;
			 OLap= true;}
		else {
			OLap= false;}
		 

		return OLap;
	 }// close defintion
//Function to determine if the predicted ORF is true with respect to the reference ORF it is overlapping with
	 bool SameFrame(Record& RHS){//open def
		bool Pos=false;
		//modified 10/19/06 due to observed inconsistency of including stop codons in chromosome tables
		//Pos=(Reverse==RHS.Reverse)&&(Stop==RHS.Stop);//true positive or false negative
		Pos=(Reverse==RHS.Reverse)&&(Stop==RHS.Stop||Stop==RHS.Stop+3||Stop==RHS.Stop-3);//true positive	
		return Pos;
	 }// close defintion

//Default destructor
	 ~Record(){
		 RefEval=NULL;
	 }

//Function for retrograde analysis
//This is just to evaluate whether this orf is a better match
//than the one that was selected by grc
//This program is starting to get unweildy
	 bool OverMax(Record& RHS){//open def
		bool OLap=false;
		double OverLen=0;
		long LStart=0;
		long LStop=0;
		long RStart=0;
		long RStop=0;

		if (Reverse){LStart=Stop; LStop=Start;} //if the LHS is reverse frame
		else {LStart=Start;  LStop=Stop;} //if LHS is not
		if (RHS.Reverse){RStart= RHS.Stop; RStop =RHS.Start;} //if RHS is reverse
		else {RStart=RHS.Start; RStop=RHS.Stop;}
		//if(RStart>6338000){
		//	cout<<"AHAH!\n";
		//}
		 
		if (RStart>=LStart && RStart <=LStop){
			 if(LStop>RStop){OverLen=RHS.OLength;}//if one frame encompasses the other
			 else OverLen=LStop-RStart+1;
			OLap= true;}
		else if(RStart <= LStart && RStop >= LStart){
			 if(RStop>LStop){OverLen=OLength;}//if one frame encompasses the other
			 else OverLen=RStop-LStart+1;
			 OLap= true;}
		else {
			OLap= false;}
		 

		if (OLap){ //if the orfs overlap
			
			double RefOLapPercent=OverLen/RHS.OLength; //retro OLength
			double GRCOLapPercent=OverLen/OLength; //reference annotation overlap
			//adjust maxOverlap
			double CombinedOLap=RefOLapPercent*GRCOLapPercent; //keep track of combined overlap
			//if(RHS.MaxOLap<CombinedOLap){//if the maximum overlap for the reference orf is less than the current update it
			//	RHS.MaxOLap=CombinedOLap;
			//}
			if(MaxOLap<CombinedOLap){//if the overlap with respect to the ref. set orf is the greatest so far update it
				RHS.MaxOLap=CombinedOLap;//if the overlap is better then store it in the retro record and return true
				return true;
			}
		}
		return false;//default for non overlapping and non impressive CombinedOLap
	 }// close defintion


	bool operator ^(const Record& RHS)const{//open definition
		return (OLength>RHS.OLength);
	}// close definition

	bool ReportHit() const{//open defintion
		return Blank;
	}//close definiton


	//Function for reporting record variables
	int RecordOut(ostream & Out){
		if(Ref){
			Out<<"Ref"<<'\t'<<ID<<'\t'<<Start<<'\t'<<Stop<<'\t'<<OLength<<'\t'<<Strand<<'\t';
			PrintDescription(Out);
			Out<<'\n';
		}
		else{
			Out<<"GRC"<<'\t'<<ID<<'\t'<<Start<<'\t'<<Stop<<'\t'<<OLength<<'\t'<<Strand<<'\t'<<Entropy<<'\t';
			PrintDescription(Out);
			Out<<'\t'<<Bit<<"\t"<<EScore<<"\t"<<HLength<<'\t'<<DBID<<'\t'<<DBOrg<<'\n';
		}
		return 0;
	}//close definition



	//This function prints the words in the descriptiion list
	int PrintDescription(ostream & Out){//open def
		for (list<string>:: iterator It1 =Description.begin(); It1!=Description.end(); It1++ ){
			Out<<*It1<<" ";
		}
		
		return 0;
	}//close def


	//This function breaks up the hit line and initializes the various data structures appropriately
	int ParseHit(){//open definition
		stringstream Breakup(Hit);//term for reading parts of the function
		string Term;
		FunctionMap::iterator FindIt;//iterator for finding go id in goterms map
		StringSet EmptySet;//empty set for initializing in HitParsing

		while (Breakup>>Term){ //read in terms
			
			if (GO::StringIsGO(Term)){//if the term is a go term
				int TempID=GO::StringToID(Term);//get integer value
				GOTerms.insert(FunctionMap::value_type(TempID,EmptySet));//insert ID (only unique)
				FindIt=GOTerms.find(TempID);//find the ID key for evidence code insertion
			}//close go term
			else if (GO::IsECode(Term)){//if its an evidence code
				if(FindIt!=GOTerms.end()){
					FindIt->second.insert(Term);//insert evidence code
				}
				else {
					cerr<<"WARNING: Possible parsing error in grc_compare at "<<Hit<<'\n';
					Description.push_back(Term);
				}
				FindIt=GOTerms.end();
			}//close if evidence code
			else{//else its a description Term
				Description.push_back(Term);
			}
		}//close while loop

		if (GOTerms.size()>0){
			HasGO=true;
		}
		return 0;
	}
		
	//string ReportAcc ()const{//open definition
	//	return Accession;
	//}// close definition

//	string ReportSeq ()const{//open definition
//		return Sequence;
//	}// close definition

}; // close prototype

#endif

