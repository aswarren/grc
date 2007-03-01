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
#include<vector>


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


	
	Subject(){//default constructor
		SubjectID=-1;
		HitID="unassigned";
		HLength=0;
		Defeated=false;
		Function="none";
		HitOrg="none";
	
	}

	//parameterized constructor
	Subject(int I=-1, long St=0, long Sp=0, string Func="None", double B=0, string ES="none", long HL=0, long AL=0, long QASt=0, long QASp=0, double MxBit=0, string HID="none", string HOrg="none"){ // parameterized constructor1
		SubjectID=I;
		Function=Func;
		HLength=HL;
		Defeated=false;
		HitID=HID;
		HitOrg=HOrg;
		Hypot=(Function.npos!=Function.find("hypothetical"));
		AlignList.push_back(Alignment(St,Sp,B,ES,AL,QASt,QASp,MxBit));//add Alignment
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
		for(list<Alignment>::iterator It=AlignList.begin(); It!=AlignList.end(); It++){
			AlignQ.push(&(*It));
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
			for(list<Alignment>::iterator It=AlignList.begin(); It!=AlignList.end(); It++){
				AlignQ.push(&(*It));
			}
		}// close self assignment
		return *this;
	}// close definition




	//Reporter function that tells whether a record has been knocked out
	bool Dead(){return Defeated;}

	//Tombstones this record
	int KnockOut(){
		Defeated=true;
		return 0;
	}
	
	//function for retrieving information from the top alignment
	int GetInfo(int& SID, string& Func, string& HID, string& HOrg, long& HLen, bool& Hyp, long& ALen, double& BScore, string& EScr, double& EVal, long& HB, long& LB, double& MxBit, long& QAStart, long& QAStop, double& RelB, long& Strt, long& Stp){
		Alignment* TopA=AlignQ.top();
		SID=SubjectID;
		Func=Function;
		HID=HitID;
		HOrg=HitOrg;
		HLen=HLength;
		Hyp=Hypot;
		TopA->GetInfo(ALen, BScore, EScr, EVal, HB, LB, MxBit, QAStart, QAStop, RelB, Strt, Stp);
		TopA=NULL;
		return 0;
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
		else return ((*(S1->AlignQ.top()))<(*(S2->AlignQ.top())));
	}//close def.
	
};//close prototype

#endif
