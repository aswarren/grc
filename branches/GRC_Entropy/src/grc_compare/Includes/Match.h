// Match.h


//Programmer: Andrew Warren
//email: anwarren@vt.edu
//Date of last modification: 04/xx/06      




#ifndef Match_H
#define Match_H

#include "Record.h"
#include "GO.h"
#include "GOMatch.h"
#include <set>
using std::set;

enum GOSTATUS{GRC, Ref, RefGRC, NoGO};//keeps of whether the match involves prediction/ref with GO terms assigned
//class Record;//forward declaration

	//NOTE The Match record is getting dual use
	//it could be improved to use inheritance from a base Match class for both uses
	//(1)Used for storing the pairs of putative orfs and the reference orf that they overlap with the most
	//(2)Used for storing pairs of orfs in the KnockAnalysis file which tracks the grc_overlap
	//results of which orf won and lost the contest between to overlapping orfs
class Match {//open prototype
		friend std::ostream& operator<<(std::ostream& Out, const Match& M);
public:
	double RefOLapPercent; //Subject Sequence percent overlap
	double GRCOLapPercent; //Query sequence percent overlap
	double OverLen; //length of overlap
	double RefWordMatchPercent; //setubal word percent
	double GRCWordMatchPercent; //grc word percent
	double CombinedOLap; //keep track of combined overlap
	result StatScore; //tracks the result for this match
	list<string> MTerms; //Terms that match
	Record* GRCRecord;//GRC prediction ORF
	Record* RefRecord;//Ref ORF
	int SmallRef; //tracks if the Reference ORF is less than 300 bp
	set<GOMatch> Confirmed;//GO annotations that are confirmed (ancestor of a MostSpecificAnnotation)
	FunctionMap NotCompatible;//Annotations that are not compatible (do not contain a msa on path to root)
	set<GOMatch> Compatible; //annotations that have a msa as an ancestor
	GOSTATUS GOStat;//reports of go status of the match
	


	Match(){//default constructor
		RefOLapPercent=0;
		GRCOLapPercent=0;
		CombinedOLap=0;
		OverLen=0;
		GRCWordMatchPercent=0;
		RefWordMatchPercent=0;
		GRCRecord=NULL;
		RefRecord=NULL;
		StatScore=TN;
		SmallRef=0;
		GOStat=NoGO;
	}// close default

	~Match(){//default destructor
		GRCRecord=NULL;
		RefRecord=NULL;
	}//close destructor



			//Copy Constructor
	 Match(const Match &Source){// open defintion
		 RefOLapPercent=Source.RefOLapPercent;
		 GRCOLapPercent=Source.GRCOLapPercent;
		 GRCWordMatchPercent=Source.GRCWordMatchPercent;
		 RefWordMatchPercent=Source.RefWordMatchPercent;
		 StatScore=Source.StatScore;
		 CombinedOLap=Source.CombinedOLap;
		 OverLen=Source.OverLen;
		 MTerms=Source.MTerms;
		 GRCRecord=Source.GRCRecord;
		 RefRecord=Source.RefRecord;
		 SmallRef=Source.SmallRef;
		 Compatible=Source.Compatible;
		 NotCompatible=Source.NotCompatible;
		 Confirmed=Source.Confirmed;
		 GOStat=Source.GOStat;
	 }// close definition


	 //Assignment Operator
	 Match& operator =(const Match &Source){// open defintion
		 if (this!= &Source){// open non-self assignment consq.
			 RefOLapPercent=Source.RefOLapPercent;
			GRCOLapPercent=Source.GRCOLapPercent;
			GRCWordMatchPercent=Source.GRCWordMatchPercent;
			RefWordMatchPercent=Source.RefWordMatchPercent;
			StatScore=Source.StatScore;
			CombinedOLap=Source.CombinedOLap;
			OverLen=Source.OverLen;
			MTerms=Source.MTerms;
			GRCRecord=Source.GRCRecord;
			RefRecord=Source.RefRecord;
			SmallRef=Source.SmallRef;//assign the small reference boolean state
			Compatible=Source.Compatible;
			NotCompatible=Source.NotCompatible;
			Confirmed=Source.Confirmed;
			GOStat=Source.GOStat;
			}
		 return *this;
	 }// close definition

	//The Match record is getting dual use
	//it could be improved to use inheritance from a base Match class for both uses
	//(1)Used for storing the pairs of putative orfs and the reference orf that they overlap with the most
	//(2)Used for storing pairs of orfs in the KnockAnalysis file which tracks the grc_overlap
	//results of which orf won and lost the contest between to overlapping orfs

	 Match(Record* BP, Record* CP, const double& OLap, bool PosOrNeg, int NumSmall=0, GO* GOAccess=NULL){//parameterized constructor
		GRCRecord=BP;
		RefRecord=CP;
		OverLen=OLap;
		GOStat=NoGO;
		
		//Set RefEval Pointer
		if(RefRecord->Ref){//if the refernece record  really is a reference record and not a GRC record for FNanalysis
			RefOLapPercent=OverLen/RefRecord->OLength; //reference annotation overlap
			GRCOLapPercent=OverLen/GRCRecord->OLength; //grc annotation overlap
			StatScore=GRCRecord->Evaluation;//get the result of the comparison
			if(StatScore ==TP || StatScore==FN){//if the reference record is worth tracking
				GRCRecord->RefEval=RefRecord;//set RefEval to be used in FNAnalysis
				RefRecord->RefMatched=true;//set RefMatched to keep FP and TN from using the reference
			}
		}
		 //if GRCRecord is not null which is the case if this match object is being used to record a knockout for entropy
		if(GRCRecord!=NULL){
			if(StatScore!=TN){
				FindMatchTerms();//Find which terms match between the two records.
			}//close if not TN
	
			SmallRef=NumSmall;//keep track of the number of Ref. ORFs <300bp that the GRC ORF overlaps with
			double WMatches=MTerms.size();
			double GW=BP->Hit.size();
			double SW=CP->Hit.size();
			if(GW==0)GRCWordMatchPercent=0;
			else GRCWordMatchPercent=WMatches/GW;
			if(SW==0)RefWordMatchPercent=0;	
			else RefWordMatchPercent=WMatches/SW;
	
			CombinedOLap=RefOLapPercent*GRCOLapPercent; //keep track of combined overlap
			if(RefRecord->MaxOLap<CombinedOLap){//if the maximum overlap for the reference orf is less than the current update it
				RefRecord->MaxOLap=CombinedOLap;
				RefRecord->MaxTerms=RefWordMatchPercent;
			}
			if(GRCRecord->MaxOLap<CombinedOLap){//if the overlap with respect to the grc set orf is the greatest so far update it
				GRCRecord->MaxOLap=CombinedOLap;//update the maximum overlap for this GRC predicted orf
				GRCRecord->MaxTerms=GRCWordMatchPercent;//update the percentage of terms that cooresponds to the maxolap
				
			}
			if(GOAccess!=NULL && StatScore==TP){//if the ontology is available
				FindGOStat();//check the status of the match
				if(GOStat==RefGRC){//if the annotations can be verified
					CheckGO(GOAccess);//run checkGO for each of the TP orfs with GO terms
				}//close if verifiable
			}//close ontology available
		}
	 }//close definition


	 //This function checks for the precesence of GO terms in ref. or GRC initializes
	 //the GOStatus variable GOStat for the match and if appropriate
	 int FindGOStat(){
		if (GRCRecord->HasGO){//if the GRC prediction has GO term
			if(RefRecord->HasGO){
				GOStat=RefGRC;
			}
			else {
				GOStat=GRC;
			}
		}//close has GRC prediction
		else if(RefRecord->HasGO){//reference has GO
			GOStat=Ref;
		}
		else{
			GOStat=NoGO;
		}
		return 0;
	 }//close defintion


	 //This function performs checks to see if annotation is Confirmed, Compatible, or NotCompatible
	 int CheckGO(GO* GOAccess){//open def.
		 ANCESTOR RefAncestors;
		 set<int> Verified;//all the grc annotations that are confirmed or compatible
		 //get all reference annotation ancestors
		 for(FunctionMap::iterator It=RefRecord->GOTerms.begin(); It!=RefRecord->GOTerms.end(); It++){//open for
			 GOAccess->GetAllAncestors(It->first,RefAncestors);
		 }
		 //check for confirmation
		 //get the ancestors of the annotations of reference and check if any of the predicted functions
		 //exist in the ancestor set
		 for(ANCESTOR::iterator AncIt=RefAncestors.begin(); AncIt!=RefAncestors.end();AncIt++){
			 if(AncIt->first!=NULL){//NULL check
				 FunctionMap::iterator GRCIt=GRCRecord->GOTerms.find(AncIt->first->ID);//if any of GRC terms are ancestors of any Ref term
				 if(GRCIt!=GRCRecord->GOTerms.end()){//if its an ancestor
					 GOFunction* OrigRef=GOAccess->Find(AncIt->first->DistID);//find the original Ref. that created this ancestor
					 FunctionMap::iterator RefIt=RefRecord->GOTerms.find(AncIt->first->DistID);//find the original in the refernce
					 if (OrigRef !=NULL && RefIt!=RefRecord->GOTerms.end()){
						 Verified.insert(GRCIt->first);//insert id as being verified
						 Confirmed.insert(GOMatch(OrigRef->ID, RefIt->second, OrigRef->Depth, GRCIt->first, GRCIt->second, (OrigRef->Depth)-(AncIt->first->Distance), AncIt->first->Distance));
					 }//close if OrigRef Found
					 else {//else original refernce not found
						 cerr<<"WARNING: Error in GOCheck Analysis at "<<AncIt->first->DistID<<'\n';
					 }
				 }//close if ancestor
			 }//close null check
		 }//close confirmation for loop
		 RefAncestors.clear();

		 //check for compatibility
		 //get the ancestors of the predictions and see if any of the reference annotations 
		 //are contained therein

		 for(FunctionMap::iterator GRCIt=GRCRecord->GOTerms.begin(); GRCIt!=GRCRecord->GOTerms.end(); GRCIt++){
			 ANCESTOR GRCAncestors;
			 GRCAncestors=GOAccess->GetAncestors(GRCIt->first);//get the ancestors of the GRC id
			 for(ANCESTOR::iterator AncIt=GRCAncestors.begin(); AncIt!=GRCAncestors.end();AncIt++){
				 if(AncIt->first!=NULL){//NULL check
					 FunctionMap::iterator RefIt=RefRecord->GOTerms.find(AncIt->first->ID);//if any of Ref terms are ancestors of any GRC term
					 if(RefIt!=RefRecord->GOTerms.end()){//if its an ancestor
						 GOFunction* OrigGRC=GOAccess->Find(AncIt->first->DistID);//find the original GRC. that created this ancestor
						 if (OrigGRC !=NULL){//if the original exists in hierarchy
							 Verified.insert(GRCIt->first);//insert id as being verified
							 Compatible.insert(GOMatch(RefIt->first, RefIt->second, (OrigGRC->Depth)-(AncIt->first->Distance), GRCIt->first, GRCIt->second, OrigGRC->Depth, AncIt->first->Distance));
						 }
						 else{
							 cerr<<"WARNING: Error in GOCheck Analysis at "<<AncIt->first->DistID<<'\n';
						 }
					 }//close if ancestor
				 }//close null check
			 }//close ancestor loop
		 }//close compatibility loop

         
		 //else they are NotCompatible
		 for(FunctionMap::iterator FindIt=GRCRecord->GOTerms.begin(); FindIt!=GRCRecord->GOTerms.end(); FindIt++){
			 if(Verified.find(FindIt->first)==Verified.end()){//if the term is not found in the verified set
				 NotCompatible.insert(*FindIt);
			 }
		 }//close for every GRC GOTerm
		 return 0;
	 }//close def.


	//This function is for listing the functional description terms of both the GRC and reference record
	int ListTerms(ostream& Out){
		Out<<"------------------------------------------------------------------\n";
		Out<<"GRC pgene:\t";
		for(list<string>::iterator It1=GRCRecord->Description.begin(); It1!=GRCRecord->Description.end(); It1++){
			Out<<(*It1)<<"\t";
		}
		Out <<"\n";//end of GRC word list
		
		Out<<"Ref gene:\t";
		for(list<string>::iterator It2=RefRecord->Description.begin(); It2!=RefRecord->Description.end(); It2++){
			Out<<(*It2)<<"\t";
		}
		Out<<"\n";//end of Ref word list
	return 1;
	}//close definition

	//Function to print Match and associated stats
	int MatchOut(std::ostream& Out, string Delim){
		Out<<"Stat\t";
		Out<<PrintEval(StatScore)<<"\t";


		for (list<string>::iterator It1 =MTerms.begin(); It1!=MTerms.end(); It1++ ){
			Out<<*It1<<" ";
		}

		Out<<OverLen<<"\n";
		GRCRecord->RecordOut(Out);
		RefRecord->RecordOut(Out);
		for(set<GOMatch>::iterator Con=Confirmed.begin(); Con!=Confirmed.end(); Con++){
			Out<<"Confirmed\t"<<Con->GOMatchString()<<"\n";
		}
		for(set<GOMatch>::iterator Com=Compatible.begin(); Com!=Compatible.end(); Com++){
			Out<<"Compatible\t"<<Com->GOMatchString()<<"\n";
		}
		for(FunctionMap::iterator NCom=NotCompatible.begin(); NCom!=NotCompatible.end(); NCom++){
			Out<<"Incompatible\t"<<GO::IDToString(NCom->first);
			for(set<string>::iterator GIt=NCom->second.begin(); GIt!=NCom->second.end(); GIt++){//open for loop
				if(GIt!=NCom->second.begin()){
					Out<<" "<<*GIt;
				}
				else{
					Out<<"\t"<<*GIt;
				}
			}//close inner for loop
		}//close outer for loop
		Out<<Delim<<"\n";

		return 0;
	}


//Function to print FN-Match and associated stats
	int FNOut(std::ostream& Out, string Delim){
		Out<<"Stat\t";


		for (list<string>::iterator It1 =MTerms.begin(); It1!=MTerms.end(); It1++ ){
			Out<<*It1<<" ";
		}
		if(MTerms.size()==0){
			Out<<'\t';
		}

		Out<<OverLen<<'\t'<<RefOLapPercent<<'\n';
		if(GRCRecord==NULL){
			Out<<"ENTROPY\n";
		}
		else{
			Out<<PrintEval(GRCRecord->Evaluation)<<"\t";
			GRCRecord->RecordOut(Out);
		}
		Out<<PrintEval(RefRecord->Evaluation)<<"\t";
		RefRecord->RecordOut(Out);
		Out<<"NA\t";
		if(RefRecord->RefEval != NULL){
			RefRecord->RefEval->RecordOut(Out);
		}
		else{
			Out<<"Error: Reference record not found\n";
		}
		Out<<Delim<<"\n";

		return 0;
	}



	//This function looks for matching terms between the two records
	int FindMatchTerms(){
		for (list<string>::const_iterator It1 =GRCRecord->Description.begin(); It1!=GRCRecord->Description.end(); It1++ ){
				for (list<string>::const_iterator It2 =RefRecord->Description.begin(); It2!=RefRecord->Description.end(); It2++ ){
					if(*It1==*It2) {MTerms.push_back(*It1);}
				}
		//(Description.npos!=Description.find("hypothetical")
		}
		return 0;
	}//close def



	//This function is for printing out the result of a comparison TP,FP etc.
	string PrintEval(const result& StatScore){
		switch(StatScore){//depending on Result add to statistic map
			case TN:
				return "TN";
				break;
			case FN:
				return "FN";
				break;
			case TP:
				return "TP";
				break;
			case FP:
				return "FP";
				break;
			case NRP:
				return "NRP";
				break;
			case NRN:
				return "NRN";
				break;
			default: break;
		}
		return "NA";
	}//close definition

	int GOConfirmed(){
		return (Confirmed.size());
	}
	int GOCompatible(){
		return (Compatible.size());
	}
	int GONotCompat(){
		return (NotCompatible.size());
	}

	//Sum up the Depth and Distance Stats for each Confirmed Annotation
	int SumConfirmedStats(double& Depth, double& Distance){
		for(set<GOMatch>::iterator It=Confirmed.begin(); It!=Confirmed.end(); It++){
			Depth=Depth+It->GRCDepth;
			Distance=Distance+It->Distance;
		}
		return 0;
	}

	//Sum up the Depth and Distance Stats for each Compatible Annotation
	int SumCompatibleStats(double& Depth, double& Distance){
		for(set<GOMatch>::iterator It=Compatible.begin(); It!=Compatible.end(); It++){
			Depth=Depth+It->GRCDepth;
			Distance=Distance+It->Distance;
		}
		return 0;
	}

};//close prototype





#endif
