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
	double EDR;//actually the entropy distance ratio
	int AACount[21];
	list<Alignment*> AlignV;//all the alignments for this Record for this segment of sequence
	//default constructor
	SeqCalc(){
		RawBit=0;
		EDR=0;
		MaxBit=0;
		for(int t=0; t<21; t++){
			AACount[t]=0;
		}
	}
	
	SeqCalc(double RB, double MB, double E){
		RawBit=RB;
		EDR=E;
		MaxBit=MB;
		for(int t=0; t<21; t++){
			AACount[t]=0;
		}
	}
	//copy constructor
	SeqCalc(const SeqCalc &Source){// open defintion
		RawBit=Source.RawBit;
		MaxBit=Source.MaxBit;
		EDR=Source.EDR;
		for(int t=0; t<21; t++){
			AACount[t]=Source.AACount[t];
		}
		AlignV=Source.AlignV;
	}
// 	//assignment operator
	SeqCalc& operator =(const SeqCalc &Source){// open defintion
		if (this!= &Source){// open non-self assignment consq.
			RawBit=Source.RawBit;
			MaxBit=Source.MaxBit;
			EDR=Source.EDR;
			AlignV=Source.AlignV;
			for(int t=0; t<21; t++){
				AACount[t]=Source.AACount[t];
			}
		}
		return *this;
	}
	//add the frequencies for the amino acids in Source to the frequencies in this SeqCalc object
	int AddCounts(const SeqCalc& Source){
		for(int t=0; t<21; t++){
			AACount[t]+=Source.AACount[t];
		}
	}
	//update the entropies for the alignments here
	int UpdateEntropy(){
		for(list<Alignment*>::iterator It=AlignV.begin(); It!=AlignV.end(); It++){
			(*It)->UpdateEDR(EDR);
		}
		return 0;
	}
};//close prototype

typedef map<long,SeqCalc, std::less<long> > SeqCalcMap;
typedef priority_queue<Subject*,vector<Subject*>,OrderSubject> PQSubject;
typedef set<Subject*> SubjectSet;
typedef map<GOFunction*,SubjectSet, OrderDepth> FuncToSubject;

//AARecord is a class representing the query ORF that is aligned in the BLAST search
//It contains a priority queue of Subjects to which the query ORF was aligned
//These subjects are ranked according to their alignments
//So each subject contains alignments. An alignment object is created
//for each possible start site suggested by the actual alignment.
//When a pop of the queue occurs the top (alignment,start) pair is removed
//from the queue until no more alignments to that subject exist
//then that subject is popped off the queue. This repeats untill all
//subjects are removed.

class AARecord {//open prototype
	friend std::ostream& operator<<(std::ostream& ACOut, const AARecord& AC);
	friend std::ostream& operator<<(std::ostream& ACOut, AARecord* AC);

private:
//public:
	string ID; //unique for each record
	string GenomeID;
	long Start;//these coordinates are stored here when the query orf has no hit
	long Stop;
	long LowBase;
	long HighBase;
	long Offset;
	bool Reverse; //Is it in a - frame
	bool Blank; //indicates whether the was a hit to this orf in the DB
	long QLength; //the length of the query orf
	long CurrentLength;
	bool Defeated; //marks whether this record has been knocked out(tombstone method)
	double HighScore;//the highest bit score out of all alignmnets for this query
	list<Subject> PrimaryHits;//other blast hits for this query orf
	PQSubject SubjectQ;;
	map<string,Subject*> SubjectNames;//map of the subject names used to enforce unique subject addition in AddPrimary function
	//sequence specific calculations are stored here, according to increasing offset from the stop site
	 SeqCalcMap CalcMap;//for storing and retrieving RawBitValues and entropy values stored according to decreasing length
	double EDR;//entropy distance ratio of (Coding/NonCoding)
	Subject* CurrentRep;//subject whose alignment is serving as the current representative of this orf
	FuncToSubject GOTerms;//maps the GO terms to the subjects from which they come (unique GO Terms and Subject pointers enforced)
	FuncToSubject ConsensusAnnot;//annotations indicated by multiple non-CurrentRep subjects currentRep is excluded becuase this would duplicate information

public:
	
	AARecord(){//default constructor

		ID="unassigned";
		GenomeID="NONE";
		Reverse =false;
		Blank =true;
		Defeated=false;
		QLength=0;
		CurrentLength=0;
		Offset=0;
		EDR=0;
		Start=Stop=LowBase=HighBase=0;
		CurrentRep=NULL;
	}

	//Initialize the Values for the record
	int InitRecord( CalcPack& CP, string TID="unassigned", long St=0, long Sp=0, string HID="none", double B=0, string ES="none", long HL=0, long AL=0, long QASt=0, long QASp=0, string Func="none", string HOrg="none"){ // parameterized constructor1
		ID=TID;
		ParseID();//parse the query ID
		//CP.SelectGenome(GenomeID);
		Start=St;
		Stop=Sp;
		Reverse =false;
		Blank =true;
		Defeated=false;;
		HighScore=0;//highest score so far 
		Blank=(B==0); //if the bit score is 0 then blank is true
		CurrentRep=NULL;
		CalcMap.clear();
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
		CurrentLength=QLength=HighBase-LowBase+1;
		

		//This is a lazy addition. ToDo: Modify AArecord to use the values at the top of the BitQueue
		if(!Blank){
			AddPrimary(CP, St,Sp,HID,B,ES,HL,AL,QASt,QASp,Func,HOrg);//add Subject
		}
		else{//it is blank set entropy
			CalcMap.insert(map<long,SeqCalc>::value_type(CurrentLength,SeqCalc()));//insert new SeqCalc based on this segment of sequence
			SeqCalcMap::iterator MarkIt=CalcMap.find(CurrentLength);
			CP.GetAACount(MarkIt->second.AACount,LowBase,HighBase,Reverse);
			EDR=CP.GetEntropy(MarkIt->second.AACount);
		}
		return 0;
	}


		//Copy Constructor
	 AARecord(const AARecord &Source){// open defintion
		ID=Source.ID; //unique for each record
		Start=Source.Start;
		Stop=Source.Stop;
		LowBase=Source.LowBase;
		HighBase=Source.HighBase;
		Defeated=Source.Defeated;
		Reverse=Source.Reverse; //Is it in a - frame
		Blank=Source.Blank;
		QLength=Source.QLength;//Orf length
		CurrentLength=Source.CurrentLength;
		PrimaryHits=Source.PrimaryHits;
		HighScore=Source.HighScore;
		EDR=Source.EDR;
		CalcMap=Source.CalcMap;
		GOTerms=Source.GOTerms;
		Offset=Source.Offset;
		ConsensusAnnot=Source.ConsensusAnnot;
		string TempName="none";
		for(list<Subject>::iterator It=PrimaryHits.begin(); It!=PrimaryHits.end(); It++){
			SubjectNames.insert(map<string,Subject*>::value_type(It->GetID(),&(*It)));//insert pointer to Subject based on name
		}
		PQSubject CopyQ=Source.SubjectQ;

		while(!CopyQ.empty()){
			TempName=CopyQ.top()->GetID();
			SubjectQ.push((SubjectNames.find(TempName)->second));//add each address for that was on the old queue
		}
		if(Source.CurrentRep!=NULL){
			TempName=Source.CurrentRep->GetID();
			CurrentRep=(SubjectNames.find(TempName)->second);
		}
		else {
			CurrentRep=NULL;
		}
	}// close definition

	 //Assignment Operator
	 AARecord& operator =(const AARecord &Source){// open defintion
		if (this!= &Source){// open non-self assignment consq.
			Start=Source.Start;
			Stop=Source.Stop;
			LowBase=Source.LowBase;
			HighBase=Source.HighBase;
			ID=Source.ID; //unique for each record
			Defeated=Source.Defeated;
			Offset=Source.Offset;
			Reverse=Source.Reverse; //Is it in a - frame
			Blank=Source.Blank;
			QLength=Source.QLength;//Orf length
			CurrentLength=Source.CurrentLength;
			PrimaryHits=Source.PrimaryHits;
			HighScore=Source.HighScore;
			EDR=Source.EDR;
			CalcMap=Source.CalcMap;
			string TempName="none";
			GOTerms=Source.GOTerms;
			ConsensusAnnot=Source.ConsensusAnnot;
			for(list<Subject>::iterator It=PrimaryHits.begin(); It!=PrimaryHits.end(); It++){
				SubjectNames.insert(map<string,Subject*>::value_type(It->GetID(),&(*It)));//insert pointer to Subject based on name
			}
			PQSubject CopyQ=Source.SubjectQ;
	
			while(!CopyQ.empty()){
				TempName=CopyQ.top()->GetID();
				SubjectQ.push((SubjectNames.find(TempName)->second));//add each address for that was on the old queue
			}
			if(Source.CurrentRep!=NULL){
				TempName=Source.CurrentRep->GetID();
				CurrentRep=(SubjectNames.find(TempName)->second);
			}
			else {
				CurrentRep=NULL;
			}
		}// close self assignment
		return *this;
	}// close definition
	
	//This function repalces all underscores with spaces
	int ParseID(){
		unsigned int ChPosition=ID.find("**");//look for '_' in ID indicating that there is a genome ID attached
		if(ChPosition!=string::npos){
			string TempID=ID;
			//TempID.replace(ChPosition,2," ");//replace '_' with a space
			//ChPosition=ID.find("**");
			//TempID.replace(ChPosition,2," ");//replace '_' with a space
			stringstream ParseSS;
			ParseSS<<TempID;
			getline(ParseSS,ID,'*');
			ParseSS.ignore();//ignore the next *
			//ParseSS>>ID;//pass orf id through
			getline(ParseSS,GenomeID,'*');
			//ParseSS>>GenomeID; //assign genome id
			ParseSS.ignore();//ignore next *
			ParseSS>>Offset;//assign offset value for contig coordinate conversioni
			ID+="_"+GenomeID;//reassign ORF ID to be orf_contig
		}
		else{
			GenomeID="NONE";
		}
		return 0;
	}

	 	//> OPERATOR overload
		//Both operators use the CurrentRep for the following two cases
		//Case1 One of the orfs have no alignment, in which case
		//there is no subject so currentRep==NULL
		//Case2 One of the orfs has alignments but they have all 
		//been removed in the dequing process which can only happen
		//in a comparison by comparison basis between AARecords
		//after which time the Queue is refreshed
	bool operator>(const AARecord& RHS)const{


		if(CurrentRep==NULL){
			return false;//LHS cannot be bigger if it has no score
		}
		else if(RHS.CurrentRep==NULL){
			return true;//if the RHS is has no score then LHS is bigger
		}
		else {
			return ((*CurrentRep)>(*RHS.CurrentRep));//else return who is bigger
		}
	}



		//< OPERATOR overload
	bool operator<(const AARecord& RHS)const{

	
		if(RHS.CurrentRep==NULL){
			return false;
		}
		else if(CurrentRep==NULL){
			return true;
		}
		else{
			return ((*CurrentRep)<(*RHS.CurrentRep));//else return who is bigger
		}

	}
	


	//Report Lowbase
	long ReportLowBase(){
		if(Blank){
			return LowBase;
		}
		else if(CurrentRep==NULL){
			cerr<<"LowBase is trying to be accessed when no representative\n";
			throw 20;
		}
		else{
			return CurrentRep->ReportLowBase();
		}
	}


	//Report HighBase
	long ReportHighBase(){
		if(Blank){
			return HighBase;
		}
		else if(CurrentRep==NULL){
			cerr<<"HighBase is trying to be accessed when no representative\n";
			throw 20;
		}
		else{
			return CurrentRep->ReportHighBase();
		}
	}


	//Update Coordinates
	//Because some orfs do not have alignments all coordinate information will
	//be evaluated at the AARecord level
	int UpdateCoord(){
		LowBase=ReportLowBase();
		HighBase=ReportHighBase();
		CurrentLength=HighBase-LowBase+1;
		if(Reverse){
			Start=HighBase;
			Stop=LowBase-3;
		}
		else{
			Start=LowBase;
			Stop=HighBase+3;
		}
		return 0;
	}



	//bool Overlap
	//returns the length of there being an overlap of ORFs
	//OR if two no_hit orfs are being compared it returns the distance between them if no overlap
	 int Overlap(AARecord& RHS){//open def
		int OverLen=0;
		 
		if (RHS.LowBase>=LowBase && RHS.LowBase <=HighBase){
			if(HighBase>=RHS.HighBase){OverLen=RHS.CurrentLength;}//if one frame encompasses the other
			else OverLen=HighBase-RHS.LowBase+1;
		}
		else if(RHS.LowBase <= LowBase && RHS.HighBase >= LowBase){
			 if(RHS.HighBase>=HighBase){OverLen=CurrentLength;}//if one frame encompasses the other
			 else OverLen=RHS.HighBase-LowBase+1;
		}
		//else if(RHS.Bit==0 && Bit==0){//neither have hits return distance between two orfs in possible intergenic region
		//	if(LowBase<RHS.LowBase){
		//		OverLen=HighBase-RHS.LowBase;//Negative overlap is distance between
		//	}
		//	else {
		//		OverLen=RHS.HighBase-LowBase;
		//	}
		//}
		return OverLen;
	 }// close defintion


	//bool CompatibleOverlap
	//meant to be run on Loser.CompatibleOverlap(Winner)
	//returns whether an overlap occurs on the three prime side of THIS ORF aka LHS
	//used to determine if there is any chance that adjusting the start site will result in compatible orfs
	//return true if threeprimeoverlap , return false if five prime overlap
	//if Loser is encompassed by winner return false
	//if Loser encompasses winner return true
	//otherwise return by case
	bool CompatibleOverlap(AARecord& Winner){
		if(Blank){
			return false;//if the loser has no alignment do not make an adjustment
		}

		if (Winner.LowBase>=LowBase && Winner.LowBase <=HighBase){
			if(HighBase>=Winner.HighBase){//if loser encompasses the winner
				return true;
			}
			else{//else overlap on high side of loser
				if (Reverse){//if its reversed the start site can be adjusted
					return true;
				}
				else{
					return false;
				}
			}
		}
		//if the overlap occurs on the low side of this loser orf
		else if(Winner.LowBase <= LowBase && Winner.HighBase >= LowBase){
			//if winner encompasses loser
			 if(Winner.HighBase>=HighBase){
				return false;//there is no hope
			}
			 else {//else overlap on low side of loser
				if(Reverse){
					return false;
				}
				else{
					return true;
				}
			}
		}
		
		//else why are you doing this check in the first place?
		else{
			cerr<<"Logic error in calling grc_overlap CompatibleOverlap\n";
		}
		return true;
	}






	//bool Knockout
	//returns whether two ORFs break overlapping threshold so that one should be knocked out
	//modified for different threshold overlaps based on blast score
	bool ToKnockOut(AARecord& RHS, int OverLen){//open def
		if(RHS.ID==ID){return true;}//if the same id then one has to go
		if(OverLen<=0){//else the overlap is negative which is distance between two no_hit orfs Possible intergenic region
				cerr<<"Logic error:ToKnockOut called on orfs that do not overlap.\n";
				throw 20;
		}

		bool NoScore =(RHS.Blank || Blank);//if either ORF does not have a score
		double MaxOverlap=45;
		float OEFactor=.05/3;//E-score to overlap threshold converstion factor 
		double OLapThreshold=15;//default overlap threshold is 12 nucl.

		if(!NoScore){
			OLapThreshold=CurrentRep->OverlapThreshold(*RHS.CurrentRep,MaxOverlap);
		}//close if both have score

		else if(RHS.Blank && Blank){
			if(RHS.EDR>EDR){
				OLapThreshold=MaxOverlap*(1-RHS.EDR);
			}
			else {
				OLapThreshold=MaxOverlap*(1-EDR);
			}
			
		}
		else if(RHS.Blank){
			OLapThreshold=MaxOverlap*(1-RHS.EDR);
		}
		else{//LHS has no hit
			OLapThreshold=MaxOverlap*(1-EDR);
		}


		if(OverLen>OLapThreshold || OverLen>45){
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
		//return (CurrentLength>RHS.CurrentLength);//LHS.QueryLength>RHS.QueryLength
		//return (Entropy<RHS.Entropy);

		double LengthFrac=0;
		double EntropyFrac=0;
		//create discriminator fractions by using (difference/Min(value1, value2))
		if(CurrentLength>RHS.CurrentLength){
			LengthFrac=double(CurrentLength-RHS.CurrentLength)/double(RHS.CurrentLength);
		}
		else{
			LengthFrac=double(RHS.CurrentLength-CurrentLength)/double(CurrentLength);
		}
		if(EDR>RHS.EDR){
			EntropyFrac=double(EDR-RHS.EDR)/double(RHS.EDR);
		}
		else{
			EntropyFrac=double(RHS.EDR-EDR)/double(EDR);
		}
		if (EntropyFrac*2>LengthFrac){//the difference in length has to 2 times the difference in entropy
			return (EDR<RHS.EDR); //LHS is better than RHS if entropy is lower
		}
		else return (CurrentLength>RHS.CurrentLength); //LHS is better than RHS if length is greater		
	}// close definition


	//Function to build inverted index of functions with pointers to the subjects from which
	//they come
	int BuildGOTerms(CalcPack& CP){
		vector<int> TempIDs;
		FuncToSubject::iterator FindIt;
		//for each Subject
		for(list<Subject>::iterator It= PrimaryHits.begin(); It!=PrimaryHits.end(); It++){
			//if the Subject has been annotated with GO terms
			if(It->ReportGO()){
				It->GOContent(TempIDs);//get go ids
				//for each go id
				for(vector<int>::iterator BuildIt=TempIDs.begin(); BuildIt!= TempIDs.end();BuildIt++){
					GOFunction* TempFunc=CP.GOAccess->Find(*BuildIt);//get pointer to go function
					if(TempFunc!=NULL){//if it exists in the GO ontology
						//check for previous existence
						FindIt=GOTerms.find(TempFunc);
						//if exists add subject pointer
						if(FindIt!=GOTerms.end()){
							FindIt->second.insert(&(*It));
						}
						//if not add GOID and pointer
						else{
							FuncToSubject::iterator MarkIt=(GOTerms.insert(FuncToSubject::value_type(TempFunc,SubjectSet()))).first;//insert Function
							MarkIt->second.insert(&(*It));
						}
					}//close if it exists
				}//close loop for these functions
				TempIDs.clear();//clear vector of functions for this subject
			}//close has GO
		}//close for loop all subjects
		return 0;
	}

	//Function for creating concensus annotations from GO terms in various subjects
	//Concensus annotations are those functions supported by multiple hits to this
	//protein but are not supported by the top hit (as returning these functions would
	//duplicate information already given in the annotation) Consensus annotations
	//must be of depth >=2 and have support from different subjects >=2
	//called as post processing for positives in main()
	int GOAnalysis(CalcPack& CP){
		ANCESTOR Family;//for retrieving ancestors of a GOID
		
		//if the query orf has hits
		if(!Blank){
			BuildGOTerms(CP);//create GO term subject map for this orf
			if(GOTerms.size()>0){//if there are GO terms
				//Build ancestor group to look for consensus annotations
				for(FuncToSubject::iterator AddIt=GOTerms.begin(); AddIt!=GOTerms.end(); AddIt++){
					//only create ancestor group for those do not come from CurrentRep
					//if(AddIt->second.find(CurrentRep)==AddIt->second.end()){
					CP.GOAccess->GetAllAncestors(AddIt->first->ReportID(), Family);//get the ancestors for this term
					//add any functions with support>2 that do not come from CurrentRep
					if(AddIt->second.size()>=2 && AddIt->second.find(CurrentRep)==AddIt->second.end()){//check for multiple support in the most specific annotations
						ConsensusAnnot.insert(*AddIt);//add this function as concensus annotation
					}
					//for each of the ancestors retrieved add it to the list
					for(ANCESTOR::iterator AncIt=Family.begin(); AncIt!=Family.end(); AncIt++){
						//if the depth>=2 assume top term depth 0
						if(AncIt->first->ReportDepth()>=2){
							//check if already exists
							FuncToSubject::iterator FindIt=ConsensusAnnot.find(AncIt->first);
							if(FindIt!=ConsensusAnnot.end()){ //if it does exist
								//add each subject to support set
								TransferSupport(FindIt->second, AddIt->second);
							}
							else{//if it does not exist
								//create new annotation and add support in the form of subjects with this function
								ConsensusAnnot.insert(FuncToSubject::value_type(AncIt->first,AddIt->second));
							}
						}//close if has depth >=2
					}//close for each ancestor
					//}//close if not from CurrentRep
					Family.clear();//clear the ancestors for this term
				}//close build ancestor group for loop
				CreateConsensus(CP);//filter concensus annotations
			}//close if there are GO terms
		}//close if has hit
		return 0;
	}//close definition


	//Transfer Support
	//this function is used by GOAnalysis to tranfer the support 
	//from a GO Term to one of its ancestors
	int TransferSupport(SubjectSet& Target, SubjectSet& Source){
		for(SubjectSet::iterator SIt=Source.begin(); SIt!=Source.end(); SIt++){
			Target.insert(*SIt);
		}
		return 0;
	}


	//Function for creating consensus annotations from a list of putative ones
	//takes all functions created as TempCon in GOAnalysis function
	int CreateConsensus(CalcPack& CP){
		ANCESTOR ConFamily;//ancestors of each consensus annotation
		set<GOFunction*> ToErase;
		//find those that do not have sufficient support
		//must have support of 2 alignments and one of those has to be at least 50% of the conservation potential for the query
		for(FuncToSubject::iterator It=ConsensusAnnot.begin(); It!=ConsensusAnnot.end(); It++){
			if(It->second.size()<2 || AnnotScore(It)<.50){
				ToErase.insert(It->first);
			}
		}
		//eliminate those that do not have sufficient support
		for(set<GOFunction*>::iterator EraseIt=ToErase.begin(); EraseIt!=ToErase.end(); EraseIt++){
			ConsensusAnnot.erase(*EraseIt);
		}
		ToErase.clear();//clear erase vector

		//Check for consensus conflicts by getting the ancestors of the consensus functions
		for(FuncToSubject::iterator AddIt=ConsensusAnnot.begin(); AddIt!=ConsensusAnnot.end(); AddIt++){
			int TempDepth=AddIt->first->ReportDepth();//this is just for debugging
			if(AddIt->second.find(CurrentRep)!=AddIt->second.end()){//if this is an annotation that will be made already
				ToErase.insert(AddIt->first);//then it should not be a consensus annotation
			}
			if(ToErase.find(AddIt->first)==ToErase.end()){//if it has not been marked for removal
				CP.GOAccess->GetAllAncestors(AddIt->first->ReportID(), ConFamily);//get the ancestors for this term
				for(ANCESTOR::iterator CheckIt=ConFamily.begin(); CheckIt!=ConFamily.end(); CheckIt++){
					if(CheckIt->first!=AddIt->first){//if the ancestor is not the current term
						FuncToSubject::iterator FindIt=ConsensusAnnot.find(CheckIt->first);//look for the ancestor in the concensus annotations
						//if one consensus annotation is the ancestor of another 
						//and it is not already marked for elimination then one must go
						if(FindIt!=ConsensusAnnot.end() && ToErase.find(FindIt->first)==ToErase.end()){
							//AddIt is the more specific function
							//LHS<RHS (more specific<ancestor)
							if(CompareAnnot(AddIt,FindIt)){//if AddIt loses stop checking its ancestors
								ToErase.insert(AddIt->first);
								ConFamily.clear();
								break;
							}
							else{//else the ancestor loses
								ToErase.insert(FindIt->first);
							}
						}
					}
				}//close for each ancestor
			}//close if not marked for removal
			ConFamily.clear();
		}//close for each consensus function
		//remove defeated annotations
		for(set<GOFunction*>::iterator EraseIt=ToErase.begin(); EraseIt!=ToErase.end(); EraseIt++){
			ConsensusAnnot.erase(*EraseIt);
		}
		return 0;
	}

	//Compares the annotations that represented by FuncToSubj iterators
	//compares the annotations based on depth, support, and conservation
	//LHS is expected to be the more specific item
	//This function is called by CreateConsensus
	//LHS<RHS for (Depth/LHS.Depth)*(Subjects/RHS.Subjects
	bool CompareAnnot(FuncToSubject::iterator& LHS, FuncToSubject::iterator&RHS){
		double LHScore=0;
		double RHScore=0;
		double RDepth=RHS->first->ReportDepth();
		double LDepth=LHS->first->ReportDepth();
		double LBitSum=0;
		double RBitSum=0;
		//sum the bit fractions for LHS
		for(SubjectSet::iterator LIt=LHS->second.begin(); LIt!=LHS->second.end(); LIt++){
			LBitSum+=(*LIt)->ReportTopBitFrac();
		}
		for(SubjectSet::iterator RIt=RHS->second.begin(); RIt!=RHS->second.end(); RIt++){
			RBitSum+=(*RIt)->ReportTopBitFrac();
		}
		LHScore=1*LBitSum/RBitSum;
		RHScore=(RDepth/LDepth)*((RBitSum-LBitSum)/RBitSum);
		if(RHScore==LHScore){//if its a tie go with the more specific function
			return false;//LHS is greater
		}
		else return (LHScore<RHScore);
	}

	//This function returns the TopScore (bitFrac) of all Subjects
	//that align to this query ORF that have the function Represented by the iterator
	double AnnotScore(FuncToSubject::iterator& Source){
		double MaxScore=0;
		double CurrentScore=0;
		for(SubjectSet::iterator It=Source->second.begin(); It!=Source->second.end(); It++){
			CurrentScore=(*It)->ReportTopBitFrac();
			if(CurrentScore>MaxScore){//if it is the highest score so far
				MaxScore=CurrentScore;
			}
		}
		return MaxScore;
	}


	//report ID of the record aka queryID
	string ReportID(){
		return ID;
	}

	bool HasHit() const{//open defintion
		return !Blank;
	}//close definiton



	//Reports the Bit/MaxBit value for the alignment currently
	//being used as rep
	double ReportBitFrac(){
		if(Blank){
			return 0;
		}
		else if(CurrentRep!=NULL){
			return CurrentRep->ReportBitFrac();
		}
		else {
			cerr<<"Logic error: Reporting BitFraction of orf with Hit but no Rep.\n";
			throw 20;
		}
		return 0;
	}



	double ReportEntropy(){
		return EDR;
	}


	//Reporter function that tells whether a record has been knocked out
	bool Dead(){return Defeated;}

	//Tombstones this record
	int KnockOut(){
		Defeated=true;
		if(!Blank){
			RefreshAll();//refresh all queues
		}
		return 0;
	}

	//Adds Subjects to OtherHits and makes top RelBit value Subject to be Record Rep.
	int AddPrimary(CalcPack& CP, long St, long Sp, string& HID, double& B, string& ES, long& HL, long& AL, long& QASt, long& QASp, string& Func, string& HOrg){
		//CP.SelectGenome(GenomeID);
		long LBase=0;
		long HBase=0;
		long OrigStart=St;
		double StartScore=0;
		Alignment* TempAlign=NULL;

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

		while(CP.FindStarts(St,OrigStart,Sp,QASt,Reverse, StartScore)) {//find all start sites from aligned region back to original
			//update low base and highbase
			if(Reverse){
				HBase=St;
			}
			else{
				LBase=St;
			}
			map<string,Subject*>::iterator FindIt;
			//Search for the subject ID
			FindIt=SubjectNames.find(HID);
			SeqCalcMap::iterator CalcIt=CalcSeqScore(CP,LBase,HBase,Reverse);
			if(FindIt!=SubjectNames.end()){//if the subject ID is found
				TempAlign=FindIt->second->AddAlign(St,Sp,B,ES,AL,QASt,QASp,CalcIt->second.MaxBit,StartScore, CalcIt->second.EDR);//add Alignment
				int SizeAlignV=CalcIt->second.AlignV.size();
				CalcIt->second.AlignV.push_back(TempAlign);
			}
			else{//else add a new subject
				int TempID=PrimaryHits.size();
				PrimaryHits.push_back(Subject(TempID,HID,HL,Func,HOrg));//add Subject
				TempAlign=PrimaryHits.back().AddAlign(St,Sp,B,ES,AL,QASt,QASp,CalcIt->second.MaxBit,StartScore, CalcIt->second.EDR);//add Alignment
				SubjectNames.insert(map<string,Subject*>::value_type(HID,&PrimaryHits.back()));//insert pointer to Subject based on name
				int SizeAlignV=CalcIt->second.AlignV.size();
				CalcIt->second.AlignV.push_back(TempAlign);
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
	int GetEntropy(){
		if(!Blank && CalcMap.find(CurrentLength)!=CalcMap.end()){
			EDR=CalcMap.find(CurrentLength)->second.EDR;
		}

		else if(!Blank){
			cerr<<"Logic Error:Could not find previously calculated entropy.";
			throw 20;
		}
		return 0;
	}

	//Switch the representative for the query orf
	//Assumes SubjectQ is ordered in order of decreasing importantance
	int SwitchRep(){
		if(SubjectQ.size()>0){
			CurrentRep=SubjectQ.top();
			UpdateCoord();
			GetEntropy();
		}
		else{
			CurrentRep=NULL;
		}
		return 0;
	}//close def.
	


	//This function looks through the other hits/subjects for the orf to see
	//if there is a compatible orf and if so selects that one to be representative
	//returns whether they are incompatible
	bool Incompatible(AARecord& Winner){

		if(!HasHit()){//if this loser does not have a hit
			return true;//no compatability can be found
		}
			
		if(!CompatibleOverlap(Winner)){
			return true;
		}

		bool ToKO=true;
		//For each
		while(!SubjectQ.empty() && ToKO){//open while loop
			PopTop();//get rid of top hit
			SwitchRep();//switch this Records Rep. hit
			if(CurrentRep!=NULL){
				int OverL=Overlap(Winner);
				if(OverL==0){//if the two do not overlap
					return false;
				}
				else{//else check if compatible
					ToKO=ToKnockOut(Winner,OverL);//check compatability
				}
			}
		}//close while loop

		if(ToKO){
			return true;
		}
		else return ToKO;
	}//end def.



	//This function removes the top alignment from the top subject
	//and then refreshes the subject priority Queue by emptying and reloading it
	int PopTop(){
		if(SubjectQ.size()>0){//if there is something to pop
			//pop the top alignment off the top
			SubjectQ.top()->PopTopAlign();
			RefreshSubjectQ();//Refresh the subject Q
		}
		return SubjectQ.size();
	}

	
	//This function refreshes the subjectQ adding only those subjects that still have alignments
	//in the alignment Q
	int RefreshSubjectQ(){
		while(!SubjectQ.empty()){
			SubjectQ.pop();
		}
		//for each subject in the list
		for(list<Subject>::iterator It=PrimaryHits.begin(); It!=PrimaryHits.end(); It++){
			//if the alignQ has alignments left in it
			if(It->ReportAlignQSize()>0){
				SubjectQ.push(&(*It));//add that subject to the Queue
			}
		}
		return SubjectQ.size();
	}



	//this function refreshes both subject and alignment queues
	int RefreshAll(){
		while(!SubjectQ.empty()){
			SubjectQ.pop();
		}
		for(list<Subject>::iterator It=PrimaryHits.begin(); It!=PrimaryHits.end(); It++){
			It->RefreshAlignQ();
			SubjectQ.push(&(*It));//put the subject back on the queue
		}
		SwitchRep();
		return 0;
	}



	//This function is designed to be run after all subjects/alignments have been added to the Record
	//It will pass in the highest BitScore out of all the alignments for this query so that
	//each (alignment,start site) pair can be scored according to its individual and relative characteristics
	int UpdateScores(){
	
		while(!SubjectQ.empty()){//clear the alignQ
			SubjectQ.pop();
		}
		//for each subject update the scores
		for(list<Subject>::iterator It=PrimaryHits.begin(); It!=PrimaryHits.end(); It++){
			It->UpdateScores(HighScore);
			SubjectQ.push(&(*It));//add to the subject Q
		}
		SwitchRep();
		return 0;
	}


	//Refresh entropies
	//this function refreshes the entropies of all alignments
	//and resets all the queues
	int RefreshEntropy(CalcPack& CP){
		if(Blank){
			EDR=CalcMap.begin()->second.EDR=CP.GetEntropy(CalcMap.begin()->second.AACount);
		}
		else{
			for(SeqCalcMap::iterator It=CalcMap.begin(); It!=CalcMap.end(); It++){
				It->second.EDR=CP.GetEntropy(It->second.AACount);
				It->second.UpdateEntropy();
			}
		}
		return 0;
	}

	//Reset Status 
	//reset as to whether this Record has been defeated
	int ResetStatus(){
		Defeated=false;
		return 0;
	}




	//This function calculates the maximum possible bit score for a given segment of the genome
	//Checks to see if the rawbit has been calculated for this LB HB Reverse Combo
	//Calculates the rawbit additively so that the same sequence is not iterated
	//over multiple times
	//Need to clean up this coordinate to string conversion +1 -1 stuff
	SeqCalcMap::iterator CalcSeqScore(CalcPack& CP, const long& LowB, const long& HighB, const bool& Rev){
		//CP.SelectGenome(GenomeID);
		long LowerBound=0;//lower bound on calc raw bit
		long UpperBound=0;//upper bound on calc raw bit
		long Length=HighB-LowB+1;
		SeqCalcMap::iterator FindIt;
		SeqCalcMap::iterator MarkIt;
		FindIt=CalcMap.find(Length);//find according to the offset from stop which is always HighBase-LowBase


		//if the score has been found
		if(FindIt!=CalcMap.end()){
			return ((FindIt));//return address of the Calculation container
		}
		//else calculate the score and create new SeqCalc object
		else{
			MarkIt=CalcMap.insert(map<long,SeqCalc>::value_type(Length,SeqCalc())).first;//insert new SeqCalc based on this segment of sequence
			int TempSize=MarkIt->second.AlignV.size();

			if(CalcMap.size()>1){//if the calc map has values other than the one just added
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
					MarkIt->second.RawBit=CP.CalcRawBit(LowerBound, UpperBound, Rev)+CalcMap.begin()->second.RawBit;
					//Calculate frequencies for the new segment
					CP.GetAACount(MarkIt->second.AACount,LowerBound,UpperBound,Rev);
					//Add AACounts of the rest of the orf to the current section's
					MarkIt->second.AddCounts(CalcMap.begin()->second);
				}
				else{//not bigger so check to see if there is one small enough
					//(1)loop to check if there is a previous calc small enough
					FindIt=CalcMap.begin();
					while(FindIt!=CalcMap.end() && FindIt->first >= Length){
						 FindIt++;
					}
					//if not small enough do self calc
					if(FindIt==CalcMap.end()){//if it walked off the end then there wasn't one small enough
						MarkIt->second.RawBit=CP.CalcRawBit(LowB,HighB,Rev);
						CP.GetAACount(MarkIt->second.AACount,LowB,HighB,Rev);
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
						MarkIt->second.RawBit=CP.CalcRawBit(LowerBound, UpperBound, Rev)+FindIt->second.RawBit;
						//Calculate frequencies for the new segment
						CP.GetAACount(MarkIt->second.AACount,LowerBound,UpperBound,Rev);
						//Add AACounts of the rest of the orf to the current section's
						MarkIt->second.AddCounts(FindIt->second);
					}
				}
			}//close if the calc map has values
			else{//no values
				MarkIt->second.RawBit=CP.CalcRawBit(LowB,HighB,Rev);
				CP.GetAACount(MarkIt->second.AACount,LowB,HighB,Rev);
			}
		}
		MarkIt->second.MaxBit=((MarkIt->second.RawBit*CP.Lambda)-log(CP.K))/M_LN2;
		MarkIt->second.EDR=CP.GetEntropy(MarkIt->second.AACount);
		return ((MarkIt));//return pointer to SeqCalc structure that holds the scores
	}//close definition


//Displays information in tab delimited format about AARecord
//and its current representatives
//this function violates information hiding
	int DisplayInfo(std::ostream& Out){
		double ConsensusConfidence=0;
		if(CurrentRep!=NULL){
			//print out any consensus GO annotations
			for(FuncToSubject::iterator It= ConsensusAnnot.begin(); It!= ConsensusAnnot.end(); It++){
				for(SubjectSet::iterator SIt=(It->second).begin(); SIt!=(It->second).end(); SIt++){
					if(ConsensusConfidence < (*SIt)->ReportTopBitFrac()){
						ConsensusConfidence= (*SIt)->ReportTopBitFrac();
					}
				}
				Out<<GO::IDToString(It->first->ReportID())<<" ( "<<ConsensusConfidence<<" ICA ) ";//these are all infered by consensus annotation
			}
			
			CurrentRep->DisplayInfo(Out);//display the information about the subject
			
		}
		else{
			cerr<<"Error in Displaying record information\n";
			throw 20;
		}
		return 0;
	}

//submits the nucleotide frequencies for creation of new EDP coding/ EDP non-coding
	int SubmitCount(CalcPack& CP){
		SeqCalcMap::iterator FindIt=CalcMap.find(CurrentLength);
		if (FindIt==CalcMap.end()){
			cerr<<"logic error: in finding AA frequencies in AARecord\n";
		}
		//submit the frequencies of bases, the length of the orf, and whether it was defeated
		CP.TotalEDP(FindIt->second.AACount, CurrentLength, Defeated);
		return 0;
	}

}; // close prototype

typedef multimap<long,AARecord*> RecordMap;
typedef map<string,AARecord*> IDMap;

#endif

