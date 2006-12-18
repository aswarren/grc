// AARecord.h


//Programmer: Andrew Warren
//email: anwarren@vt.edu
//Date of last modification: 10/xx/06




#ifndef AARecord_H
#define AARecord_H


#include "Subject.h"
#include <queue>
#include <functional>
#include<vector>
#include <set>
using std::vector;
using std::deque;
using std::priority_queue;
using std::set;





class AARecord {//open prototype
	friend std::ostream& operator<<(std::ostream& ACOut, const AARecord& AC);
	friend std::ostream& operator<<(std::ostream& ACOut, AARecord* AC);

//private:
public:
	string Sequence;
	string ID; //unique for each record
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
	int SID;//identifier of subject that is current rep
	list<Subject> PrimaryHits;//other blast hits for this query orf
	list<Subject> SecondaryHits;//other blast hits for this query orf
	priority_queue<Subject*,vector<Subject*>,MoreBit> PrimaryQ;
	priority_queue<Subject*,vector<Subject*>,MoreBit> SecondaryQ;
	double LowComplexity;
//public:
	
	AARecord(){//default constructor
		Sequence="blah";
		ID="unassigned";
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
		SID=0;
		LowComplexity=0;
	}

	//parameterized constructor
	AARecord(string TID="unassigned", long St=0, long Sp=0, string H="None", double LC=0, double B=0, string ES="none", long L=0, long A=0, long QASt=0, long QASp=0, double MxBit=0, string HID="none", string HOrg="none"){ // parameterized constructor1
		ID=TID;
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
		SID=0;
		LowComplexity=LC;

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
		//This is a lazy addition. ToDo: Modify AArecord to use the values at the top of the BitQueue
		AddPrimary(St,Sp,H,B,ES,L,A,QASt,QASp,MxBit,HID,HOrg);//add Subject
	}


		//Copy Constructor
	 AARecord(const AARecord &Source){// open defintion
		 Sequence=Source.Sequence;
		 ID=Source.ID; //unique for each record
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
		PrimaryHits=Source.PrimaryHits;
		SecondaryHits=Source.SecondaryHits;
		SID=Source.SID;
		LowComplexity=Source.LowComplexity;
		for(list<Subject>::iterator It=PrimaryHits.begin(); It!=PrimaryHits.end(); It++){
			PrimaryQ.push(&(*It));
		}
		for(list<Subject>::iterator It=SecondaryHits.begin(); It!=SecondaryHits.end(); It++){
			SecondaryQ.push(&(*It));
		}
	}// close definition

	 	//> OPERATOR overload
		//adding e-score evaluation and inverting signs 04/04/06
	bool operator>(const AARecord& RHS)const{
		double BitRatio=0;
		if (RelBit==0 || RHS.RelBit==0){
			BitRatio=0;
		}
		else if(RelBit<RHS.RelBit){
			BitRatio=RelBit/RHS.RelBit;
		}
		else {BitRatio=RHS.RelBit/RelBit;}

		if(Hypot && !RHS.Hypot && BitRatio>.80){
			return false;
		}
		else if(!Hypot && RHS.Hypot && BitRatio>.80){
			return true;
		}
		else return(RelBit>RHS.RelBit);
		//return(EValue<RHS.EValue);
	}

		//< OPERATOR overload
	bool operator<(const AARecord& RHS)const{
		double BitRatio=0;
		if (RelBit==0 || RHS.RelBit==0){
			BitRatio=0;
		}
		else if(RelBit<RHS.RelBit){
			BitRatio=RelBit/RHS.RelBit;
		}
		else {BitRatio=RHS.RelBit/RelBit;}

		if(Hypot && !RHS.Hypot && BitRatio>.80){
			return true;
		}
		else if(!Hypot && RHS.Hypot && BitRatio>.80){
			return false;
		}
		return (RelBit<RHS.RelBit);
		//return(EValue>RHS.EValue);
	}
	
	 //Assignment Operator
	 AARecord& operator =(const AARecord &Source){// open defintion
		if (this!= &Source){// open non-self assignment consq.
		 Sequence=Source.Sequence;
		 ID=Source.ID; //unique for each record
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
		PrimaryHits=Source.PrimaryHits;
		SecondaryHits=Source.SecondaryHits;
		SID=Source.SID;
		LowComplexity=Source.LowComplexity;
		for(list<Subject>::iterator It=PrimaryHits.begin(); It!=PrimaryHits.end(); It++){
			PrimaryQ.push(&(*It));
		}
		for(list<Subject>::iterator It=SecondaryHits.begin(); It!=SecondaryHits.end(); It++){
			SecondaryQ.push(&(*It));
		}
		}// close self assignment
		return *this;
	}// close definition
	
	 //destructor
	 ~AARecord(){
		 PrimaryQ.empty();
		 SecondaryQ.empty();
	 }




	//bool Overlap
	//returns the length of there being an overlap of ORFs
	//OR if two no_hit orfs are being compared it returns the distance between them if no overlap
	 int Overlap(AARecord& RHS){//open def
		int OverLen=0;
		 
		if (RHS.LowBase>=LowBase && RHS.LowBase <=HighBase){
			if(HighBase>=RHS.HighBase){OverLen=RHS.QLength;}//if one frame encompasses the other
			else OverLen=HighBase-RHS.LowBase+1;
		}
		else if(RHS.LowBase <= LowBase && RHS.HighBase >= LowBase){
			 if(RHS.HighBase>=HighBase){OverLen=QLength;}//if one frame encompasses the other
			 else OverLen=RHS.HighBase-LowBase+1;
		}
		else if(RHS.Bit==0 && Bit==0){//neither have hits return distance between two orfs in possible intergenic region
			if(LowBase<RHS.LowBase){
				OverLen=HighBase-RHS.LowBase;//Negative overlap is distance between
			}
			else {
				OverLen=RHS.HighBase-LowBase;
			}
		}
		return OverLen;
	 }// close defintion




	//bool Knockout
	//returns whether two ORFs break overlapping threshold so that one should be knocked out
	//modified for different threshold overlaps based on blast score
	bool KnockOut(AARecord& RHS, int OverLen){//open def
		if(RHS.ID==ID){return true;}//if the same id then one has to go

		bool NoScore =(RHS.Bit==0 || Bit==0);//if either ORF does not have a score
		double BitFrac1;
		double BitFrac2;
		//double BitFrac=0;
		double MaxOverlap=120;
		float OEFactor=.05/3;//E-score to overlap threshold converstion factor 
		double OLapThreshold=12;//default overlap threshold is 12 nucl.

		if(!NoScore){
			BitFrac1=Bit/MaxBit;
			BitFrac2=RHS.Bit/RHS.MaxBit;
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
			if(QLength<300 || RHS.QLength<300){ //if the length of one of the orfs is less than 300 decrease maxoverlap
				if(QLength<RHS.QLength){
					MaxOverlap=MaxOverlap*(QLength/300);
					//MaxOverlap=0.4*QLength;
				}
				else {
					MaxOverlap=MaxOverlap*(RHS.QLength/300);
					//MaxOverlap=0.4*RHS.QLength;
				}
			}
				
			if(BitFrac1<BitFrac2){
				//Diff=(Bit/RHS.Bit)*2;
				//ScoreFactor=(Bit/1000)+1;
				OLapThreshold=(BitFrac1*MaxOverlap);
			}
			else{
				//Diff=(RHS.Bit/Bit)*2;
				//ScoreFactor=(RHS.Bit/1000)+1;
				OLapThreshold=(BitFrac2*MaxOverlap);
			}
			//OLapThreshold=BitFrac*90;//for every .10 they have similar score give a 1% overlap
			//if(OLapT2<OLapT1){//set overlap threshold to be a function of the e-value score
			//	OLapThreshold=OLapT2;
			//}
			//else {OLapThreshold=OLapT1;}//take the min threshold
		}//close if both have score
		else if(OverLen>12){return true;}//LAST ADDED 080306  If two putative ORFs overlap and one has no hits it must overlap very little
		else if(RHS.Bit==0 && Bit==0){
			if(OverLen>0){return true;}
			else{//else the overlap is negative which is distance between two no_hit orfs Possible intergenic region
				int Distance=OverLen*-1;
				return(LowComplexity>0.90||RHS.LowComplexity>0.90||QLength<300||RHS.QLength<300);
			}
		}

		//double RPercentOLap=double(OverLen)/RHS.QLength;
		//double LPercentOLap=double(OverLen)/QLength;

		if(OverLen>OLapThreshold && OverLen>12){
			return true;
		}
		else return false;
	}//close definition
		
		

	//If there is no hit information for either ORF evaluate the fraction of LowComplexity X's/Length
	bool operator ^(const AARecord& RHS)const{//open definition
		return (QLength>RHS.QLength);
		//return (LowComplexity<RHS.LowComplexity);
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

	//Adds Subjects to OtherHits and makes top RelBit value Subject to be Record Rep.
	int AddPrimary(long& St, long& Sp, string& H, double& B, string& ES, long& L, long& A, long& QASt, long& QASp, double& MxBit, string& HID, string& HOrg){
		int TempID=PrimaryHits.size();
		PrimaryHits.push_back(Subject(TempID,St,Sp,H,B,ES,L,A,QASt,QASp,MxBit,HID,HOrg));//add Subject
		PrimaryQ.push(&PrimaryHits.back());
		
		//if(BitQueue.top()->ID==OtherHits.back().ID && OtherHits.size()>0){//if the newest hit is the best
			//switch representative hit
		//}
		return 0;
	}//close def.

	//Adds Subjects to OtherHits and makes top RelBit value Subject to be Record Rep.
	int AddSecondary(long& St, long& Sp, string& H, double& B, string& ES, long& L, long& A, long& QASt, long& QASp, double& MxBit, string& HID, string& HOrg){
		int TempID=SecondaryHits.size();
		SecondaryHits.push_back(Subject(TempID,St,Sp,H,B,ES,L,A,QASt,QASp,MxBit,HID,HOrg));//add Subject
		SecondaryQ.push(&SecondaryHits.back());
		
		//if(BitQueue.top()->ID==OtherHits.back().ID && OtherHits.size()>0){//if the newest hit is the best
			//switch representative hit
		//}
		return 0;
	}//close def.

	//Switch the representative for the query orf
	//Assumes that the first hit is at the top of the OtherHit list
	int SwitchRep(){

		if(PrimaryHits.size()>1 && SID!=PrimaryQ.top()->ID){
			Subject* TempS=PrimaryQ.top();
			SID=TempS->ID;
			Hit=TempS->Hit;
			HitID=TempS->HitID;
			HitOrg=TempS->HitOrg;
			HLength=TempS->HLength;
			Hypot=TempS->Hypot;
			ALength=TempS->ALength;
			Bit=TempS->Bit;
			EScore=TempS->EScore;
			EValue=TempS->EValue;
			HighBase=TempS->HighBase;
			LowBase=TempS->LowBase;
			MaxBit=TempS->MaxBit;
			QAlignStart=TempS->QAlignStart;
			QAlignStop=TempS->QAlignStop;
			QLength=TempS->QLength;
			RelBit=TempS->RelBit;
			Start=TempS->Start;
			Stop=TempS->Stop;
		}
		return 0;
	}//close def.
	
	//Switch the representative for the query orf
	//Assumes that the first hit is at the top of the OtherHit list
	int SwitchRep2(){
		if(SecondaryHits.size()>1){
			Subject* TempS=SecondaryQ.top();
			SID=TempS->ID;
			Hit=TempS->Hit;
			HitID=TempS->HitID;
			HitOrg=TempS->HitOrg;
			HLength=TempS->HLength;
			Hypot=TempS->Hypot;
			ALength=TempS->ALength;
			Bit=TempS->Bit;
			EScore=TempS->EScore;
			EValue=TempS->EValue;
			HighBase=TempS->HighBase;
			LowBase=TempS->LowBase;
			MaxBit=TempS->MaxBit;
			QAlignStart=TempS->QAlignStart;
			QAlignStop=TempS->QAlignStop;
			QLength=TempS->QLength;
			RelBit=TempS->RelBit;
			Start=TempS->Start;
			Stop=TempS->Stop;
		}
		return 0;
	}//close def.

	//This function looks through the other hits/subjects for the orf to see
	//if there is a compatible orf and if so selects that one to be representative
	bool Incompatible(AARecord& Winner){
		priority_queue<Subject*,vector<Subject*>,MoreBit> CopyQ=PrimaryQ;//copy for refreshing PrimaryQ
		bool ToKO=true;
		//For each
		PrimaryQ.pop();//get rid of top hit
		while(!PrimaryQ.empty() && ToKO){//open while loop
			SwitchRep();//switch this Records Rep. hit
			int OverL=Overlap(Winner);
			if(OverL==0){//if the two do not overlap
				ToKO=false;
			}
			else{//else check if compatible
				ToKO=KnockOut(Winner,OverL);//check compatability
			}
			PrimaryQ.pop();//get rid of top hit
		}//close while loop

		PrimaryQ=CopyQ;//refresh PrimaryQ
		if(ToKO){//if ToKO is still true check secondary Q for compatible starts
			CopyQ=SecondaryQ;
			while(!SecondaryQ.empty() && ToKO){//open while loop
				SwitchRep2();//switch this Records Rep. hit
				int OverL=Overlap(Winner);
				if(OverL==0){//if the two do not overlap
					ToKO=false;
				}
				else{//else check if compatible
					ToKO=KnockOut(Winner,OverL);//check compatability
				}
				SecondaryQ.pop();//get rid of top hit
			}//close while loop
		}

		return ToKO;
	}//end def.
		



	//string ReportAcc ()const{//open definition
	//	return Accession;
	//}// close definition

//	string ReportSeq ()const{//open definition
//		return Sequence;
//	}// close definition

}; // close prototype

typedef multimap<long,AARecord*> RecordMap;
typedef map<string,AARecord*> IDMap;

#endif

