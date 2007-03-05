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
#include "CalcPack.h"
#define _USE_MATH_DEFINES
#include <math.h>
using std::vector;
using std::deque;
using std::priority_queue;
using std::set;

//container classs for holding maxbit and entropy calculations for a segment of sequence
class SeqCalc{//open prototype
public:
	double RawBit;
	double MaxBit;
	double Entropy;
	//default constructor
	SeqCalc(){
		RawBit=0;
		Entropy=0;
		MaxBit=0;
	}
	//paramaterized constructor
	SeqCalc(double RB, double MB, double E){
		RawBit=RB;
		Entropy=E;
		MaxBit=MB;
	}
	//copy constructor
	SeqCalc(const SeqCalc &Source){// open defintion
		RawBit=Source.RawBit;
		MaxBit=Source.MaxBit;
		Entropy=Source.Entropy;
	}
// 	//assignment operator
	SeqCalc& operator =(const SeqCalc &Source){// open defintion
		if (this!= &Source){// open non-self assignment consq.
			RawBit=Source.RawBit;
			MaxBit=Source.MaxBit;
			Entropy=Source.Entropy;
		}
	}
};//close prototype

typedef map<long,SeqCalc, std::greater<long> > SeqCalcMap;


class AARecord {//open prototype
	friend std::ostream& operator<<(std::ostream& ACOut, const AARecord& AC);
	friend std::ostream& operator<<(std::ostream& ACOut, AARecord* AC);

private:
//public:
	string Sequence;
	string ID; //unique for each record
	string Function;
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
	double MaxBit;//the maximum possible bit score for the query sequence
	double RelBit;//bit score * bit/maxbit
	double BitFrac;//bit Fraction Bit/MaxBit
	double HighScore;//the highest bit score out of all alignmnets for this query
	bool Hypot;//bool to tell whether description contains hypothetical
	string HitID;//id of the hit in db
	string HitOrg;//name of the organism in the db hit
	int SID;//identifier of subject that is current rep
	list<Subject> PrimaryHits;//other blast hits for this query orf
	list<Subject> SecondaryHits;//other blast hits for this query orf

	priority_queue<Subject*,vector<Subject*>,OrderSubject> SubjectQ;;
	map<string,Subject*> SubjectNames;//map of the subject names used to enforce unique subject addition in AddPrimary function
	//sequence specific calculations are stored here, according to increasing offset from the stop site
	 SeqCalcMap CalcMap;//for storing and retrieving RawBitValues and entropy values stored according to decreasing length
	double Entropy;
public:
	
	AARecord(){//default constructor
		Sequence="blah";
		ID="unassigned";
		Function="unassigned";
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
		Entropy=0;
	}

	//parameterized constructor
	AARecord( CalcPack& CP, string TID="unassigned", long St=0, long Sp=0, string HID="none", double B=0, string ES="none", long HL=0, long AL=0, long QASt=0, long QASp=0, string Func="none", string HOrg="none"){ // parameterized constructor1
		ID=TID;
		Function="unassigned";
		Start=St;
		Stop=Sp;
		HighBase=0;
		LowBase=0;
		Bit=0;
		EScore="Not Assigned";
		EValue=100000;
		HLength=0;
		Reverse =false;
		Blank =true;
		Defeated=false;;
		ALength=0;
		QAlignStart=0;
		QAlignStop=0;
		MaxBit=0;
		RelBit=0;
		HitID="none";
		HitOrg="none";
		SID=0;
		HighScore=0;//highest score so far 
		Blank=(B==0); //if the bit score is 0 then blank is true
		QLength=labs(St-Sp)+1;

		

		//This is a lazy addition. ToDo: Modify AArecord to use the values at the top of the BitQueue
		if(!Blank){
			AddPrimary(CP, St,Sp,HID,B,ES,HL,AL,QASt,QASp,Func,HOrg);//add Subject
		}
			
	}


		//Copy Constructor
	 AARecord(const AARecord &Source){// open defintion
		 Sequence=Source.Sequence;
		 ID=Source.ID; //unique for each record
		 Function=Source.Function;
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
		HighScore=Source.HighScore;
		Entropy=Source.Entropy;
		BitFrac=Source.BitFrac;
		CalcMap=Source.CalcMap;
		for(list<Subject>::iterator It=PrimaryHits.begin(); It!=PrimaryHits.end(); It++){
			SubjectQ.push(&(*It));
			SubjectNames.insert(map<string,Subject*>::value_type(It->GetID(),&(*It)));//insert pointer to Subject based on name
		}

	}// close definition

	 //Assignment Operator
	 AARecord& operator =(const AARecord &Source){// open defintion
		if (this!= &Source){// open non-self assignment consq.
			Sequence=Source.Sequence;
			ID=Source.ID; //unique for each record
			Function=Source.Function;
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
			Entropy=Source.Entropy;
			BitFrac=Source.BitFrac;
			HighScore=Source.HighScore;
			CalcMap=Source.CalcMap;
			for(list<Subject>::iterator It=PrimaryHits.begin(); It!=PrimaryHits.end(); It++){
				SubjectQ.push(&(*It));
				SubjectNames.insert(map<string,Subject*>::value_type(It->GetID(),&(*It)));//insert pointer to Subject based on name
			}
		}// close self assignment
		return *this;
	}// close definition
	
	 //destructor
	 ~AARecord(){
		//SubjectQ.clear();
		SubjectNames.clear();
	 }

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
				//return(Entropy>0.90||RHS.Entropy>0.90||QLength<300||RHS.QLength<300);
				return false;
			}
		}

		//double RPercentOLap=double(OverLen)/RHS.QLength;
		//double LPercentOLap=double(OverLen)/QLength;

		if(OverLen>OLapThreshold && OverLen>12){
			return true;
		}
		else return false;
	}//close definition
		
		

	//If there is no hit information for either ORF evaluate the entropy and length
	//LHS BetterThan RHS
	//Checks to see which is a more powerful discriminator length or entropy
	//and evaluates the LHS BetterThan RHS based on that discriminator
	bool operator ^(const AARecord& RHS)const{//open definition
		//return (Entropy<RHS.Entropy);
		//return (QLength>RHS.QLength);//LHS.QueryLength>RHS.QueryLength
		//return (Entropy<RHS.Entropy);
		double LengthFrac=0;
		double EntropyFrac=0;
		//create discriminator fractions by using (difference/Min(value1, value2))
		if(QLength>RHS.QLength){
			LengthFrac=double(QLength-RHS.QLength)/double(RHS.QLength);
		}
		else{
			LengthFrac=double(RHS.QLength-QLength)/double(QLength);
		}
		if(Entropy>RHS.Entropy){
			EntropyFrac=double(Entropy-RHS.Entropy)/double(RHS.Entropy);
		}
		else{
			EntropyFrac=double(RHS.Entropy-Entropy)/double(Entropy);
		}
		if (EntropyFrac>LengthFrac){
			return (Entropy<RHS.Entropy); //LHS is better than RHS if entropy is lower
		}
		else return (QLength>RHS.QLength); //LHS is better than RHS if length is greater
			
	}// close definition

	string ReportID(){
		return ID;
	}

	bool HasHit() const{//open defintion
		return !Blank;
	}//close definiton

	long ReportLowBase(){
		return LowBase;
	}

	long ReportHighBase(){
		return HighBase;
	}

	double ReportBitFrac(){
		return BitFrac;
	}

	double ReportEntropy(){
		return Entropy;
	}


	//Reporter function that tells whether a record has been knocked out
	bool Dead(){return Defeated;}

	//Tombstones this record
	int KnockOut(){
		Defeated=true;
		return 0;
	}

	//Adds Subjects to OtherHits and makes top RelBit value Subject to be Record Rep.
	int AddPrimary(CalcPack& CP, long St, long Sp, string& HID, double& B, string& ES, long& HL, long& AL, long& QASt, long& QASp, string& Func, string& HOrg){
		long LBase=0;
		long HBase=0;
		long OrigStart=St;
		if(HighScore<B){//record highest bit score for any alignment for this query
			HighScore=B;
		}
		if(St>Sp){
			Reverse=true;//set Reverse frame
			LBase=Sp;
			HBase=St;
		}
		else{
			Reverse=false;
			LBase=St;
			HBase=Sp;
		}

		while(CP.FindStarts(St,OrigStart,Sp,QAlignStart,Reverse)) {//find all start sites from aligned region back to original
			map<string,Subject*>::iterator FindIt;
			//Search for the subject ID
			FindIt=SubjectNames.find(HID);
			SeqCalc* CalcPointer= CalcSeqScore(CP,LBase,HBase,Reverse);
			if(FindIt!=SubjectNames.end()){//if the subject ID is found
				FindIt->second->AddAlign(St,Sp,B,ES,AL,QASt,QASp,CalcPointer->MaxBit);//add Alignment
			}
			else{//else add a new subject
				int TempID=PrimaryHits.size();
				PrimaryHits.push_back(Subject(TempID,St,Sp,HID,B,ES,HL,AL,QASt,QASp,CalcPointer->MaxBit,Func,HOrg));//add Subject
				SubjectNames.insert(map<string,Subject*>::value_type(HID,&PrimaryHits.back()));//insert pointer to Subject based on name
			}
		}
		//No need to add to primaryQ until each Subject has been scored based on HighScore
		//SubjectQ.push(&PrimaryHits.back());

		//if(BitQueue.top()->ID==OtherHits.back().ID && OtherHits.size()>0){//if the newest hit is the best
			//switch representative hit
		//}
		return 0;
	}//close def.



	//Calculate the entropy for the current sequence alignment representative
	int CalcEntropy(CalcPack& CP){
		Entropy=CP.GetEntropy(LowBase, HighBase,Reverse);
	}

	//Switch the representative for the query orf
	//Assumes SubjectQ is ordered in order of decreasing importantance
	int SwitchRep(){

		if(PrimaryHits.size()>1 && SID!=SubjectQ.top()->SubjectID){
			Subject* TopS=SubjectQ.top();
			TopS->GetInfo(SID, Function, HitID, HitOrg, HLength, Hypot, ALength, Bit, EScore, EValue, HighBase, LowBase, MaxBit, QAlignStart, QAlignStop, RelBit, Start, Stop);
		}
		return 0;
	}//close def.
	


	//This function looks through the other hits/subjects for the orf to see
	//if there is a compatible orf and if so selects that one to be representative
	bool Incompatible(AARecord& Winner){
		priority_queue<Subject*,vector<Subject*>,OrderSubject> CopyQ=SubjectQ;//copy for refreshing SubjectQ
		bool ToKO=true;
		//For each

		SubjectQ.pop();//get rid of top hit
		while(!SubjectQ.empty() && ToKO){//open while loop
			SwitchRep();//switch this Records Rep. hit
			int OverL=Overlap(Winner);
			if(OverL==0){//if the two do not overlap
				ToKO=false;
			}
			else{//else check if compatible
				ToKO=KnockOut(Winner,OverL);//check compatability
			}
			SubjectQ.pop();//get rid of top hit
		}//close while loop

		SubjectQ=CopyQ;//refresh SubjectQ

		return ToKO;
	}//end def.
		




	//This function calculates the maximum possible bit score for a given segment of the genome
	//Checks to see if the rawbit has been calculated for this LB HB Reverse Combo
	//Calculates the rawbit additively so that the same sequence is not iterated
	//over multiple times
	//Need to clean up this coordinate to string conversion +1 -1 stuff
	SeqCalc* CalcSeqScore(CalcPack& CP, const long& LowB, const long& HighB, const bool& Rev){
		double MxB=0;
		double RawBit=0;
		long LowerBound=0;//lower bound on calc raw bit
		long UpperBound=0;//upper bound on calc raw bit
		long Length=HighB-LowB+1;
		SeqCalcMap::iterator FindIt;
		FindIt=CalcMap.find(Length);//find according to the offset from stop which is always HighBase-LowBase

		//if the score has been found
		if(FindIt!=CalcMap.end()){
			return (&(FindIt->second));//return address of the Calculation container
		}
		//else calculate the score
		else{
			if(CalcMap.size()>0){//if the calc map has values
				//check to see if this Length is bigger than previous ones
				if(Length > CalcMap.begin()->first){ 
					//if it is, adjust Segment length for RawBit calc so that previous longest can be added to it
					if(Rev){
						LowerBound=LowB+CalcMap.begin()->first;
						UpperBound=HighB;
					}
					else{//forward orf
						UpperBound=HighB-CalcMap.begin()->first;
						LowerBound=LowB;
					}
					//Rawbit is the sum of both the subsequence and additional sequence
					RawBit=CP.CalcRawBit(LowerBound, UpperBound, Rev)+CalcMap.begin()->second.RawBit;
				}
				else{//not bigger so check to see if there is one small enough
					//(1)loop to check if there is a previous calc small enough
					FindIt=CalcMap.begin();
					while(FindIt!=CalcMap.end() && FindIt->first > Length){
						 FindIt++;
					}
					//if not small enough do self calc
					if(FindIt==CalcMap.end()){//if it walked off the end then there wasn't one small enough
						RawBit=CP.CalcRawBit(LowB,HighB,Rev);
					}
					//else small enough do additive calc
					else{
						if(Rev){
							LowerBound=LowB+FindIt->first;
							UpperBound=HighB;
						}
						else{//forward orf
							UpperBound=HighB-FindIt->first;
							LowerBound=LowB;
						}
						//Rawbit is the sum of both the subsequence and additional sequence
						RawBit=CP.CalcRawBit(LowerBound, UpperBound, Rev)+CalcMap.begin()->second.RawBit;
					}
				}
			}//close if the calc map has values
			else{
				RawBit=CP.CalcRawBit(LowB,HighB,Rev);
			}
		}
		MxB=((RawBit*CP.Lambda)-log(CP.K))/M_LN2;
		CalcMap.insert(map<long,SeqCalc>::value_type(Length,SeqCalc(RawBit,MxB,0)));//insert new SeqCalc based on this segment of sequence
		FindIt=CalcMap.find(Length);
		return (&(FindIt->second));//return pointer to SeqCalc structure that holds the scores
	}//close definition

	


}; // close prototype

typedef multimap<long,AARecord*> RecordMap;
typedef map<string,AARecord*> IDMap;

#endif

