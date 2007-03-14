/*GRC.Overlap
This component of the GRC program finds the ORFs that overlap from the BLAST results file
and eliminates all but the one with the best BLAST hit.*/

/*********************************NOTE
REMOVING LOWCOMPLEXITY READ SO THAT CAN USE OLDER BLAST FILES
REMEMBER TO ADD THIS BACK IN BEFORE RUNNING TESTS
**********************************/

#include <list>
#define _USE_MATH_DEFINES
#include <math.h>
#include "AARecord.h"
#include "Compete.h"
#include "CalcPack.h"
#include <map>
#include <utility>



using std::list;
using std::cout;
using std::ios;
using std::fstream;
using std::map;
using std::pair;




int DumpList(list<AARecord>& InitList);
int DumpList(list<AARecord*>& InitList, string PosName);
int Nulify(list<AARecord*>& InitList, list<AARecord*>& InitList2);
int Compare(RecordMap& PositionMap, list<AARecord*>& WinnerList, list<AARecord*>& LoserList, CompeteMap& KOMap);
int EntropyFilter(list<AARecord*>& WinnerList, list<AARecord*>& LoserList, CompeteMap& KOMap, double EntCutoff);//removes orfs that are "winners" that have high entropy
int PrintLosers(list<AARecord*>& InitList, string NegName);
void DisplayKO(ostream& Out, CompeteMap& KOMap);




//run as GRC_overlap -i BlastResults.txt
int main (int argc, char* argv[]) {   //  Main is open
	char* BlastFile = argv[1]; //get the name of the blast test results file
	char* GenomeName= argv[2];//the name of the target genome
	string GenomeFile= argv[3];//the name of the fna file
	string Matrix=argv[4];//the matrix used for blast
	string TransFile=argv[5];//file used for translating sequences in entropy calculations



	list<AARecord> RecordList;//storage for all of the records
	RecordMap PositionMap; //the initial map	 for the aa records ordered by HighBase of the orf
	//RecordMap::_Pairib MapIt;//for keeping track of insertions into the map

	map<string,AARecord*> HitList; //the list of AARecord pointers to ORFs aka pGenes that have hits (hash by ID)
	RecordMap::iterator MapIt;//pair for keeping track of insertions into the record map

	list<AARecord*> LoserList; //the list of the eliminated records
	list<AARecord*> WinnerList; //the final list of highest non overlaping blast hits
	CompeteMap KOMap;//create map for tracking win/lose relationships


	CalcPack InfoPack(Matrix, GenomeFile, TransFile);//create information package
	

	string GName;//for storing input and creating output name from it
	GName=GenomeName;//convert input to string

	//Create the names of the Positive and negative output files
	string Negatives=GName;
	string Positives=GName;
	string Neg=".Neg";
	string Pos=".Pos";
	Negatives+=Neg;//create false/true negatives filename
	Positives+=Pos;//create fasle/treu positives filename






	ifstream BlastIn; //input for the blast results
	BlastIn.open(BlastFile); //open the blast input

	string OldID="none";
	long OldQAS=-1;
	double OldBit=0;
	string OldOrg="none";
	string OldHID="none";
	int OldCount=0;

	map<string,double>::iterator FindMB;
	while (BlastIn){//read in the input file
		string ID;//
		long Stop;
		long Start;
		long QAlignStart;
		long QAlignStop;
		long SAlignStart;
		long SAlignStop;
		int GapOpens;
		int Mismatches;
		long SubjectLen;
		string Function;
		long Length;
		double Bit;
		double PercentIdent;
		string ES;//escore
		long ALength; //the length of the alignment
		string HitOrg;//organism from which db hit comes
		string HitID;//ID from which db hit comes
		

		BlastIn >>ID; //read in the ID
		
		if (!BlastIn){//if the input hit the EOF
			break;//break from the loop
		}

		BlastIn >>Start; //read in the Start position
		BlastIn >>Stop; //read in the Stop position
		getline(BlastIn,HitID,'\t'); //skip next tab
		getline(BlastIn,HitID,'\t'); //get line for hit ID/No_hits
		long LowBase;
		long HighBase;
		double LowComplexity=0;//fraction of low complexity AA's as determined by SEG in fsa-blast
		
		bool Rev=false;//indicates if ORF is in reverse reading frame
		if (Start<Stop){
			LowBase=Start;
			HighBase=Stop;
		}
		else {
			LowBase=Stop;
			HighBase=Start;
			Rev=true;
		}

		if (HitID =="No_hits"){//open consq.
			//In>>ES; //read in the delimiter
			//Insert Record into Initial RecordMap
			BlastIn>>LowComplexity;
			RecordList.push_back(AARecord(InfoPack, ID, Start, Stop, HitID));//add record to the list
			//EditList.push_back(&((MapIt.first)->second));//add a pointer to the location of the record
		}//close consq.
		else {//there is a hit
			//string TempHack=HitID;//SWITCHED TO ACCOMODATE INVERSED ID AND DESCRIPTION IN OLD OUTPUT********
			getline(BlastIn,Function,'\t'); //get line for hit description 
			//HitID=Hit;//SWITCHED TO ACCOMODATE INVERSED ID AND DESCRIPTION IN OLD OUTPUT********
			//Hit=TempHack;//SWITCHED TO ACCOMODATE INVERSED ID AND DESCRIPTION IN OLD OUTPUT********

			getline(BlastIn,HitOrg,'\t'); //organism from which the hit comes
			BlastIn>>PercentIdent;
			BlastIn>>ALength; //alignment length
			BlastIn>>SubjectLen; //subject length
			BlastIn>>Mismatches;
			BlastIn>>GapOpens;
			BlastIn>>QAlignStart;
			BlastIn>>QAlignStop;
			BlastIn>>SAlignStart;
			BlastIn>>SAlignStop;
			BlastIn>>ES;
			BlastIn>>Bit;
			BlastIn>>LowComplexity;//ONLY needs to be entered once per query ID since all is query sequence dependent
			long OrigStart=Start;//for start searching purposes
			//boolean value to determine whether any new information is being contributed
			//TO DO: Add org and HitID check to bool Old and update structure for storing
			//same alignment regions with multiple organisms,functions,and hit_ID's
			
			//if it doesn't have align start further in or as high a bit score
			bool Old;
			Old=(ID==OldID&&HitID==OldHID&&Bit==OldBit&&QAlignStart==OldQAS);

			if(!Old){
				OldID=ID;//update old values
				OldOrg=HitOrg;
				OldHID=HitID;
				OldBit=Bit;
				OldQAS=QAlignStart;
			}
			

			//See if start can be adjusted then set lowbase
			//To Do:read in all hits for same query ID then run adjust start on each
			//After this is done a CalcMaxBit function that calculates each MaxBit
			//for each subject in 1 move across the genome
			if(!Old){
				IDMap::iterator FindID;
				FindID=HitList.find(ID);
				if(FindID!=HitList.end()){//if the orf already has a hit 
				//Add Subject orf to Record
					FindID->second->AddPrimary(InfoPack, Start, Stop, HitID, Bit, ES, SubjectLen, ALength, QAlignStart, QAlignStop, Function, HitOrg);
				}
				else{//else there has not been a hit for this orf so far. so create a new AARecord
					RecordList.push_back(AARecord(InfoPack, ID, Start, Stop, HitID, Bit, ES, SubjectLen, ALength, QAlignStart, QAlignStop, Function, HitOrg));
					HitList.insert(IDMap::value_type(ID,&RecordList.back()));//add pointer to the record to the hit list
				}
			}

			else{
				OldCount++;
			}
			//In>>ES; //read in delimiter
		}//close hit
		//OldID=ID;//remember this id in next iteration
	}// close read input
	BlastIn.close();



	//create position map
	//and switch the start site to the one with the highest conservation
	//and calculate entropy

	//this is not accounted for yet
	for(list<AARecord>::iterator PosIt=RecordList.begin(); PosIt!=RecordList.end(); PosIt++){

		if(PosIt->HasHit()){//if the query has a hit
			PosIt->UpdateScores();
		}
		PositionMap.insert(RecordMap::value_type(PosIt->ReportLowBase(), &(*PosIt))); //Add to position map
	}//close for loop


	InfoPack.Genome.clear();//clear the genome (no longer needed)
	
	ofstream ChkOut;

	//DumpList(InitList);//print out the orfs from initlist
	


	//compare the ORFs and remove the ones that conflict due to overlap
	Compare(PositionMap,WinnerList,LoserList, KOMap);
	
	//tally the number of winners with hits and
	//calculate avg. entropy of winners with hits
	double ConservCutoff=.90;//Bit fraction cutoff for creating the entropy cutoff
	double AvgEntropy=0;//The average entropy
	int NumWinHits=0;
	int NumForAvg=0;
	for (list<AARecord*>::iterator AvgIt =WinnerList.begin(); AvgIt!=WinnerList.end(); AvgIt++ ){
		if((*AvgIt)->HasHit()){//if this orf has a hit
			NumWinHits++;
			if((*AvgIt)->ReportBitFrac()>ConservCutoff){
				NumForAvg++;
				AvgEntropy+=(*AvgIt)->ReportEntropy();
			}
		}
	}
	AvgEntropy=AvgEntropy/double(NumForAvg);//finish avg entropy calc. by dividing by number of orfs with hits

	//Calculate the std deviation
	double Variance=0;
	for (list<AARecord*>::iterator EntIt =WinnerList.begin(); EntIt!=WinnerList.end(); EntIt++ ){
		if((*EntIt)->HasHit() && (*EntIt)->ReportBitFrac()>ConservCutoff){//if this orf has a hit
			double Diff=(*EntIt)->ReportEntropy()-AvgEntropy;//get difference between mean and current
			Variance+=(Diff*Diff);//square the difference of mean and current and add to total variance
		}
	}
	double EntropyDev=sqrt((Variance/(NumWinHits-1)));//calculate std deviation

	int NumFiltered=EntropyFilter(WinnerList, LoserList, KOMap, AvgEntropy+EntropyDev);//filter the orfs based on entropy
	//Tally results
	cout<<"Average entropy is\t"<<AvgEntropy<<"\n";
	cout<<"Std Dev of entropy is\t"<<EntropyDev<<"\n";
	cout<<"Number of orfs filtered from entropy\t"<<NumFiltered<<"\n";

	cout<<"******************GRC v0.01******************"<<"\n";
	cout<<"Total # of initial ORFS created:   "<<PositionMap.size()<<"\n";
	cout<<"Total # non-overlapping orfs:   "<<WinnerList.size()<<"\n";
	cout<<"Total # annotated:   "<<NumWinHits<<"\n\n";
	

	DumpList(WinnerList, Positives);
	//DumpList(InitList);
	PrintLosers(LoserList, Negatives);

	Nulify(WinnerList, LoserList);

	ofstream KOut;
	KOut.open("KnockList.txt");
	DisplayKO(KOut,KOMap);
	KOut.close();





	return 0;
}


//function to print out master list of ORFS
int DumpList(RecordMap& PositionMap){//open definition
	ofstream Out;
	Out.open("x2blastparse.txt");
		for (RecordMap::iterator It1 =PositionMap.begin(); It1!=PositionMap.end(); It1++ ){
			Out<<(It1->second);//print out the Records
		}
		Out.close();
		return 0;
}

// NULL OUT lists of AARecord pointers
int Nulify(list<AARecord*>& InitList, list<AARecord*>& InitList2){
		for (list<AARecord*>::iterator It1 =InitList.begin(); It1!=InitList.end(); It1++ ){
			*It1=NULL;//set to null
		}
		for (list<AARecord*>::iterator It2 =InitList2.begin(); It2!=InitList2.end(); It2++ ){
			*It2=NULL;//set to null
		}
		return 0;
}


//function to compare orfs and populate final list
//basically goes through the list of orfs and nullifys the ones 
//that overlap another orf but have a inferior BLAST score
//ORFs that have a tieing score the longer will be kept
//MODIFIED TO keep track of who knocks out who in the WList
//Compare(PositionMap,WinnerList,LoserList, KOMap)
int Compare(RecordMap& PositionMap, list<AARecord*>& WinnerList, list<AARecord*>& LoserList, CompeteMap& KOMap){//open defintion
	bool FirstVictim=true;
	for (RecordMap::iterator It1 =PositionMap.begin(); It1!=PositionMap.end(); It1++ ){
		//FirstVictim=true;
		int OverLen=-5555; //overlapping length of the two orfs
		string Spot="nada";//Key for finding the Competition Record
		if(!(It1->second)->Dead()){// if It1 isn't on a dead record
			for (RecordMap::iterator It2 =It1; It2!=PositionMap.end()&& OverLen!=0; It2++ ){//inner loop only needs to go through list from position It1
				OverLen=-5555;

				if(It1!=It2 && !(It2->second->Dead()) && !(It1->second->Dead())){//if not the same
					
					OverLen=It1->second->Overlap(*(It2->second));//get the amount the two overlap if at all
				//	if(((It1->second.ID=="T2432")||(It2->second.ID=="T2432"))){
				//		cout<<"Whats going on?";
				//	}

					bool KO=false;//default is to not knock one out
					if(OverLen!=0){//if they overlap
						KO=It1->second->ToKnockOut(*(It2->second), OverLen);//determine if one should be knocked out
					}
					if(KO){//if one of the ORFs should be knocked out

						if (*(It1->second)>*(It2->second)){//If record 1 has a higher score
							KO=(It2->second->Incompatible(*(It1->second)));//KO if its incompatible
							if(KO){//if KO is still true
								Spot=It1->second->ReportID();
								CompeteMap::iterator TempC=KOMap.find(Spot);
								if (TempC!=KOMap.end()){
									(TempC->second).AddLoser((It2->second));
								}
								else {
									KOMap.insert(CompeteMap::value_type(Spot, Compete((It1->second),(It2->second))));
								}
									
								LoserList.push_back((It2->second));//loser
								It2->second->KnockOut();//then record 2 is out
							}
						}
						else if(*(It1->second)<*(It2->second)){//if record 1 has a lower score
							KO=It1->second->Incompatible(*(It2->second));//KO if its incompatible
							if(KO){
								Spot=((It2)->second->ReportID());
								CompeteMap::iterator TempC=KOMap.find(Spot);
								if (TempC!=KOMap.end()){
									TempC->second.AddLoser((It1->second));
								}
								else {
									KOMap.insert(CompeteMap::value_type(Spot, Compete((It2->second), (It1->second))));
								}
									
								LoserList.push_back((It1->second));//loser
								It1->second->KnockOut();//then record 1 is out
								break;//break from the inner loop since record 1 is dead
							}
						}
						else if(*(It1->second)^*(It2->second)){//else if It1 better than It2
							KO= It2->second->Incompatible(*(It1->second));
							if(KO){
								Spot=((It1->second)->ReportID());
								CompeteMap::iterator TempC=KOMap.find(Spot);
								if (TempC!=KOMap.end()){
									(TempC->second).AddLoser((It2->second));
								}
								else {
									KOMap.insert(CompeteMap::value_type(Spot, Compete((It1->second),(It2->second))));
								}
								LoserList.push_back((It2->second));//loser
								It2->second->KnockOut(); //record 2 is out
							}
						}
						else if(*(It2->second)^*(It1->second)){//else It2 is better
							KO=It1->second->Incompatible(*(It2->second));
							if(KO){
								Spot=((It2->second)->ReportID());
								CompeteMap::iterator TempC=KOMap.find(Spot);
								if (TempC!=KOMap.end()){
									(TempC->second).AddLoser((It1->second));
								}
								else {
									KOMap.insert(CompeteMap::value_type(Spot, Compete((It2->second), (It1->second))));
								}
								LoserList.push_back((It1->second));
								It1->second->KnockOut();
								break; //break from inner loop
							}
						}
						else{//else its a tie
						}
					}//close KO orfs
				}//close if not the same	
			}//close inner loop
			if (!It1->second->Dead()){//if the record still hasn't been eliminated
				WinnerList.push_back((It1->second));//add it to the final list
			}//close if the record still hasn't been elim
		}//close initial It1 Null check
	}//close outer loop
	return 0;
}// close function


//function to print out master list of ORFS
int DumpList(list<AARecord*>& InitList, string PosName){//open definition
		ofstream ChkOut;
		ChkOut.open(PosName.c_str());
		ChkOut<<"ID\tStart\tStop\tLength(nt)\tStrand\tDBFunc.\tDBID\t DBOrg\tBit\tEScore\tHitLength\t%QueryAligned\t%HSPAligned\n";
	
		for (list<AARecord*>::iterator It1 =InitList.begin(); It1!=InitList.end(); It1++ ){
			if(*It1!=NULL){
				//cout<<**It1;//print out the Records
				ChkOut<<*It1;//print to minout
			}
			else ChkOut<<"NULL record: grc_overlap error\n";
		}
		ChkOut.close();
		return 0;
}


//function to print the comparisons that have been made
int PrintLosers(list<AARecord*>& InitList, string NegName){//open definition
	ofstream Out;
	Out.open(NegName.c_str());
	int count=0; //for printing delimeters for pair
	//list<compete>::iterator It2 =TR.begin();
	Out<<"ID\tStart\tStop\tLength(nt)\tStrand\tDBFunc.\tDBID\t DBOrg\tBit\tEScore\tHitLength\t%QueryAligned\t%HSPAligned\n";
	for (list<AARecord*>::iterator It1 =InitList.begin(); It1!=InitList.end(); It1++ ){
		Out<<*It1;
		/*switch(*It2){//depending on win loss print appropriate message
			case WIN:
				Out<<"Winner:\t";
				break;
			case LOSE:
				Out<<"Loser:\t";
				break;
			case TIE:
				Out<<"Tie:\t";
				break;
			default: break;
		}
		*/
		//Out<<(**It1);//print out the Records
			

		//It2++; //increment win lose list iterator
		count++; //increment print +++ delimeter counter
	}
	Out.close();
	return 0;
}

/*std::ostream& operator<<(std::ostream& ACOut, const AARecord& AC){
	ACOut<<"PPCG**"<<"\n";
	ACOut<<"ID:\t"<<AC.ID<<"\n";
	ACOut<<"ORF Start:\t"<<AC.Start<<"\n";
	ACOut<<"ORF Stop:\t"<<AC.Stop<<"\n";
	ACOut<<"ORF Length(nt):\t"<<AC.QLength<<"\n";
	if (AC.Reverse){
		ACOut<<"ORF Strand:\t"<<'-'<<"\n";
	}
	else 	ACOut<<"ORF Strand:\t"<<'+'<<"\n";

	ACOut<<"DB Hit:\t"<<AC.Function<<"\n"; //print hit description or no hit line
	if(!AC.Blank){//if there is a hit
		double ALength=AC.ALength;
		double OLength=AC.QLength;
		double HLength=AC.HLength;
		ACOut<<"Hit ID:\t"<<AC.HitID<<"\n";
		ACOut<<"Hit Org:\t"<<AC.HitOrg<<"\n";
		ACOut<<"Bit Score:\t"<<AC.Bit<<"\n";
		ACOut<<"E Score:\t"<<AC.EScore<<"\n";
		ACOut<<"HSP Length(aa):\t"<<AC.HLength<<"\n";
		ACOut<<"% Query Aligned:\t"<<(ALength/(OLength/3))*100<<"\n";
		ACOut<<"% HSP Aligned:\t"<<(ALength/HLength)*100<<"\n";
	}
	//ACOut<<"Blank:     "<<AC.Blank<<"\n";
	//ACOut<<"OLength:   "<<AC.QLength<<"\n";
	ACOut<<"**PPCG"<<"\n";

	return ACOut;
}*/


//ofstream overload for Compete class
std::ostream& operator<<(std::ostream& Out, const Compete& C){
	if(C.Winner==NULL){
		Out<<"Entropy";
	}
	else Out<<(C.Winner)->ReportID();
	for (list<AARecord*>::const_iterator It1 =C.Losers.begin(); It1!=C.Losers.end(); It1++ ){
		Out<<"\t"<<(*It1)->ReportID();
	}//close for loop
	Out<<"\n";
	return Out;
}//close definition





//output overload for AARecord
std::ostream& operator<<(std::ostream& ChkOut, AARecord* AC){
	//ChkOut<<"PPCG**"<<"\t";
	ChkOut<<AC->ID<<"\t";
	ChkOut<<AC->Start<<"\t";
	ChkOut<<AC->Stop<<"\t";
	ChkOut<<AC->CurrentLength<<"\t";
	
	if (AC->Reverse){
		ChkOut<<'-'<<"\t";
	}
	else 	ChkOut<<'+'<<"\t";

	ChkOut<<AC->EDR<<"\t";


	if(!AC->Blank){//if there is a hit
		AC->DisplayInfo(ChkOut);
	}

	else{
		ChkOut<<"No_hits"<<"\t";
		ChkOut<<"-\t";//HitID
		ChkOut<<"-\t";//HitOrg
		ChkOut<<"-\t";//Bit
		ChkOut<<"-\t";//EScore
		ChkOut<<"-\t";//HLength
		ChkOut<<"-\t";//% Query aligned
		ChkOut<<"-\n";//%Subject aligned
	}
	//ChkOut<<"Blank:     "<<AC.Blank<<"\n";
	//ChkOut<<"OLength:   "<<AC.QLength<<"\n";
	//ChkOut<<"**PPCG"<<"\n";
	return ChkOut;
}

//Display function for the Direct Hash<Compete>
void DisplayKO(ostream& Out, CompeteMap& KOMap){
	Out<<"Winner\tLosers...\n";
	for (CompeteMap::iterator KOIt= KOMap.begin(); KOIt != KOMap.end(); KOIt++){//open for loop
		Out<<(KOIt->second);
	}//close for loop
}//close definition






//this function removes ORFs that do not have hits from the WinnersList and places them in the losers list
//it also creates a KO record that the orf was removed due to entropy
 int EntropyFilter(list<AARecord*>& WinnerList, list<AARecord*>& LoserList, CompeteMap& KOMap, double EntCutoff){
	list<AARecord*> NuWinners;
	int OrigSize=WinnerList.size();
	//check each winner with no hit for high entropy
	for (list<AARecord*>::iterator CheckIt=WinnerList.begin(); CheckIt!=WinnerList.end(); CheckIt++){
		if((*CheckIt)->HasHit()||(*CheckIt)->ReportEntropy()<=EntCutoff){//if the record has a hit OR its entropy is under cutoff
			NuWinners.push_back((*CheckIt));//add it to the final list
		}
		else{
			CompeteMap::iterator TempC=KOMap.find("Entropy");//look for the entropy entry in the KOMap
			if(TempC!=KOMap.end()){//if its found add another loser
				(TempC->second).AddLoser(*CheckIt);
			}
			else{//if not found add a new entry
				KOMap.insert(CompeteMap::value_type("Entropy", Compete(NULL,(*CheckIt))));
			}
			LoserList.push_back((*CheckIt));//add to the loser list
			
		}
		*CheckIt=NULL;
	}
	WinnerList=NuWinners;
	NuWinners.clear();
	return OrigSize-WinnerList.size();//return number of orfs removed
}

			
	









