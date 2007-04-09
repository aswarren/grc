// Subject.h
//The subject class models the subjects from the db
//that the user provides
//each subject can have multiple alignments 

//Programmer: Andrew Warren
//email: anwarren@vt.edu
//Date of last modification: 10/xx/06




#ifndef Subject_H
#define Subject_H

#include "Alignment.h"
#include <iostream>
#include <string>
#include <fstream>
#include <math.h>
#include <cmath>
#include <map>
#include <list>
#include <sstream>
#include <queue>
#include <vector>
#include "GO.h"

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
using std::deque;
using std::priority_queue;
using std::vector;
	typedef set<string> StringSet;
	typedef map<int,StringSet> FunctionMap;

class Subject {//open prototype
public:
	int SubjectID;
	string Function;
	long HLength; //the length of the hit sequence
	bool Defeated; //marks whether this record has been knocked out(tombstone method)
	bool Hypot;//bool to tell whether description contains hypothetical
	string HitID;//id of the hit in db
	string HitOrg;//name of the organism in the db hit
	list<Alignment> AlignList;//the list of alignments for this subject
	priority_queue<Alignment*,vector<Alignment*>,OrderAlign> AlignQ;
	list<string> Description;
	FunctionMap GOTerms;//GOTerms,EvidenceCodes assigned to this prediction
	bool HasGO;
	double BestBitFrac;

	
	Subject(){//default constructor
		SubjectID=-1;
		HitID="unassigned";
		HLength=0;
		Defeated=false;
		Function="none";
		HitOrg="none";
		HasGO=false;
		BestBitFrac=0;
	}

	//parameterized constructor
	Subject(int I=-1, string HID="None", long HL=0, string Func="none", string HOrg="none"){ // parameterized constructor1
		SubjectID=I;
		Function=Func;
		HLength=HL;
		Defeated=false;
		HitID=HID;
		BestBitFrac=0;
		HitOrg=HOrg;
		HasGO=false;
		Hypot=(Function.npos!=Function.find("hypothetical"));
		ParseHit();

	}


		//Copy Constructor
	 Subject(const Subject &Source){// open defintion
		SubjectID=Source.SubjectID;
		Function=Source.Function;
		Defeated=Source.Defeated;
		Hypot=Source.Hypot;
		HLength=Source.HLength; //the length of the hit sequence
		HitOrg=Source.HitOrg;
		HitID=Source.HitID;
		AlignList=Source.AlignList;
		Description=Source.Description;
		GOTerms=Source.GOTerms;
		HasGO=Source.HasGO;
		BestBitFrac=Source.BestBitFrac;
		if (Source.AlignQ.size()>0){//if there is something to copy
			for(list<Alignment>::iterator It=AlignList.begin(); It!=AlignList.end(); It++){
				AlignQ.push(&(*It));
			}
		}
		
	}// close definition

	 //Assignment Operator
	 Subject& operator =(const Subject &Source){// open defintion
		if (this!= &Source){// open non-self assignment consq.
			SubjectID=Source.SubjectID;
			Function=Source.Function;
			Defeated=Source.Defeated;
			Hypot=Source.Hypot;
			HLength=Source.HLength; //the length of the hit sequence
			HitOrg=Source.HitOrg;
			HitID=Source.HitID;
			AlignList=Source.AlignList;
			Description=Source.Description;
			GOTerms=Source.GOTerms;
			HasGO=Source.HasGO;
			BestBitFrac=Source.BestBitFrac;
			if (Source.AlignQ.size()>0){//if there is something to copy
				for(list<Alignment>::iterator It=AlignList.begin(); It!=AlignList.end(); It++){
					AlignQ.push(&(*It));
				}
			}
		}// close self assignment
		return *this;
	}// close definition


	//>operator to tell whether one subject is better than another
	//this is judged on the current top alignment
	bool operator>(const Subject& RHS)const{
		double BitRatio=AlignQ.top()->BitRatio(*RHS.AlignQ.top());
		double EDRR=AlignQ.top()->EDRRatio(*RHS.AlignQ.top());
		double AbsBRatio=BitRatio;
		double AbsEDRR=EDRR;
		if(BitRatio<0){
			AbsBRatio*=-1;
		}
		if(EDRR<0){
			AbsEDRR*=-1;
		}
		//being a hypothetical costs you 20% versus another Hit
		//this isnt really fair to hypothetical
		//since others may be called putative, hypot. etc.
		//if the values are fairly close
		if(AbsBRatio<.20 && EDRR<.20){
			if(Hypot && !RHS.Hypot){
				return false;
			}
			else if(!Hypot && RHS.Hypot){
				return true;
			}
		}
		//else check which score is a stronger discriminator
		if(AbsBRatio>AbsEDRR){
			return(BitRatio>0);//BitRatio is positive if LHS>RHS
		}
		else{ //EDR ratio is positive if LHS>RHS (if LHS.EDR<RHS.EDR)
			return (EDRR>0);
		}
	}


	//<operator to tell whether one subject is better than another
	//this is judged on the current top alignment
	bool operator<(const Subject& RHS)const{
		double BitRatio=AlignQ.top()->BitRatio(*RHS.AlignQ.top());
		double EDRR=AlignQ.top()->EDRRatio(*RHS.AlignQ.top());
		double AbsBRatio=BitRatio;
		double AbsEDRR=EDRR;
		if(BitRatio<0){
			AbsBRatio*=-1;
		}
		if(EDRR<0){
			AbsEDRR*=-1;
		}
		//if the values are fairly close
		if(AbsBRatio<.20 && EDRR<.20){
			if(Hypot && !RHS.Hypot){
				return true;
			}
			else if(!Hypot && RHS.Hypot){
				return false;
			}
		}
		//else check which if Bit score is a stronger discriminator
		if(AbsBRatio>AbsEDRR){
			return (BitRatio<0);//Bit Ratio is negative if LHS<RHS
		}
		else{
			return (EDRR<0);
		}
	}


	//Function for adding alignment to the alignment list
	Alignment* AddAlign(long St, long Sp, double B, string ES, long AL, long QASt, long QASp, double MxBit, double SScore, double EDRatio){
		AlignList.push_back(Alignment(St,Sp,B,ES,AL,QASt,QASp,MxBit,SScore,EDRatio));//add Alignment
		if(AlignList.back().ReportBitFrac() > BestBitFrac){//record the highest bit fraction
			BestBitFrac=AlignList.back().ReportBitFrac();
		}
		return &(AlignList.back());//return address to Alignment
	}


	//This function updates the score for each (alignment, start site pair) aka each alignment object
	int UpdateScores(const double& HighScore){
		while(!AlignQ.empty()){//clear the alignQ
			AlignQ.pop();
		}
		for(list<Alignment>::iterator It=AlignList.begin(); It!=AlignList.end(); It++){
			It->SetScore(HighScore);//update the score
			AlignQ.push(&(*It));//add to the alignment Q
		}
		return 0;
	}
	

	
	//Function for returning the SubjectID of this subject
	string GetID(){
		return HitID;
	}

	//Report BestBitFrac
	double ReportTopBitFrac(){
		return BestBitFrac;
	}

	//Report GOTerms
	bool ReportGO(){
		return HasGO;
	}


	//This function returns the value of the lower part of the orf coordinates
	long ReportLowBase(){
		return AlignQ.top()->ReportLowBase();
	}


	//This function returns the value of the higher part of the orf coordinates
	long ReportHighBase(){
		return AlignQ.top()->ReportHighBase();
	}


	//Reporter function that tells whether a record has been knocked out
	bool Dead(){return Defeated;}

	//Tombstones this record
	int KnockOut(){
		Defeated=true;
		return 0;
	}
	
	//function for retrieving information from the top alignment
	int GetInfo(int& SID, string& Func, string& HID, string& HOrg, long& HLen, bool& Hyp, long& ALen, double& BScore, string& EScr, double& EVal, long& HB, long& LB, double& MxBit, long& QAStart, long& QAStop, double& BFrac, double& RelB, long& Strt, long& Stp){
		Alignment* TopA=AlignQ.top();
		SID=SubjectID;
		Func=Function;
		HID=HitID;
		HOrg=HitOrg;
		HLen=HLength;
		Hyp=Hypot;
		TopA->GetInfo(ALen, BScore, EScr, EVal, HB, LB, MxBit, QAStart, QAStop, BFrac, RelB, Strt, Stp);
		TopA=NULL;
		return 0;
	}

	//reports the size of the alignment queue
	int ReportAlignQSize(){
		return AlignQ.size();
	}

	//pops the top Alignment off the queue
	int PopTopAlign(){
		if(AlignQ.size()==0){
			cerr<<"Error trying to pop off empty queue\n";
		}
		else {
			AlignQ.pop();//remove top element
		}
		return AlignQ.size();
	}

	int RefreshAlignQ(){
		while(!AlignQ.empty()){
			AlignQ.pop();
		}
		//put each alignment back into the queue
		for(list<Alignment>::iterator It=AlignList.begin(); It!=AlignList.end(); It++){
			AlignQ.push(&(*It));
		}
		return 0;
	}


	//Return the overlap threshold for the subjects involved
	double OverlapThreshold(Subject& RHS, const double& MxOLap){
		return AlignQ.top()->OverlapThreshold(*RHS.AlignQ.top(), MxOLap);
	}


	//Reports the Bit/MaxBit value for this alignment
	double ReportBitFrac(){
		return AlignQ.top()->ReportBitFrac();
	}



	//This function breaks up the function description and initializes the various data structures appropriately
	int ParseHit(){//open definition
		stringstream Breakup(Function);//term for reading parts of the function
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
					cerr<<"WARNING: Possible parsing error in grc_overlap at "<<Function<<'\n';
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

	//Return GO Terms
	int GOContent(vector<int>& TermIDs){
		for(FunctionMap::iterator It=GOTerms.begin(); It!=GOTerms.end(); It++){
			TermIDs.push_back(It->first);
		}
		return 0;
	}

	//Display the Information about this Subject and its representative
	int DisplayInfo(std::ostream& Out){

			double ALength=AlignQ.top()->GetALength();
			double OLength=AlignQ.top()->GetLength();

			//Out<<Function<<"\t"; //print hit description
			for(FunctionMap::iterator It= GOTerms.begin(); It!=GOTerms.end(); It++){//Print GO terms if there are any
				Out<<GO::IDToString(It->first)<<" ";
				//print evidence codes
				for(StringSet::iterator It2=It->second.begin(); It2!=It->second.end(); It2++){
					Out<<*It2<<" ";
				}
			}
			//print out the description "fasta line function"
			for(list<string>::iterator DIt=Description.begin(); DIt!=Description.end(); DIt++){
				Out<<*DIt<<" ";
			}
			Out<<"\t";
			Out<<HitID<<"\t";
			Out<<HitOrg<<"\t";
			AlignQ.top()->DisplayInfo(Out);
			Out<<HLength<<"\t";
			Out<<(ALength/(OLength/3))*100<<"\t";
			Out<<(ALength/HLength)*100<<"\n";
	}


};//close prototype

//Struct for ordering priority of Subjects
//the priority of subjects is determined by their alignments
//Returns wheter Subject2 is better than Subject1
struct OrderSubject {
	bool operator()(Subject* S1, Subject*S2){
		if(S1==NULL){
			return true;
		}
		else if(S2==NULL){
			return false;
		}
		else if(S1->ReportAlignQSize()==0){
			return true;
		}
		else if(S2->ReportAlignQSize()==0){
			return false;
		}
		else return ((S1->AlignQ.top()->ReportScore()) < (S2->AlignQ.top()->ReportScore()));
	}//close def.
	
};//close prototype


string ltos(long i)	// convert long to string
	{
		stringstream s;
		s << i;
		return s.str();
	}
string btos(bool i)	// convert bool to string
	{
		stringstream s;
		s << i;
		return s.str();
	}


#endif
