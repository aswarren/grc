/*GRC.Overlap
This component of the GRC program finds the ORFs that overlap from the BLAST results file
and eliminates all but the one with the best BLAST hit.*/



#include <list>
#define _USE_MATH_DEFINES
#include <math.h>
#include "AARecord.h"
#include "Compete.h"
#include <map>
#include <utility>
using std::list;
using std::cout;
using std::ios;
using std::fstream;
using std::map;
using std::pair;

string ltos(long i)	// convert long to string
	{
		stringstream s;
		s << i;
		return s.str();
	}
string btos(bool i)	// convert long to string
	{
		stringstream s;
		s << i;
		return s.str();
	}



int DumpList(list<AARecord>& InitList);
int DumpList(list<AARecord*>& InitList, string PosName);
int Nulify(list<AARecord*>& InitList, list<AARecord*>& InitList2);
int Compare(RecordMap& EList, list<AARecord*>& FList, list<AARecord*>& CList, CompeteMap& WList);
int PrintLosers(list<AARecord*>& InitList, string NegName);
void DisplayKO(ostream& Out, CompeteMap& KOMap);
int LowerTheCase(string& Seq);
int AdjustStart1(long& St, const long& Sp, const long& QAS, const bool& Reverse, const string& Genome);
bool AdjustStart2(long& St, const long& OSt, const long& Sp, const long& QAS, const bool& Reverse, const string& Genome);
bool ReverseStart(const string& Forward);
bool ForwardStart(const string& Codon);
double CalcSS(const string& Codon, const double& Travel);
string ReverseComp(const string& Codon);//returns the reverse complement of a codon
char Complement(const char& Base);//returns the complement of a nucl.base
int InitCodes(double& Lambda, double& K, const string& Matrix, map <string,int>& ConsValue);
double CalcMaxBit(const long& LB, const long& HB, const string& Genome, map <string,int>& ConsValue, const bool& Reverse, const double& Lambda, const double& K);
double GetEntropy(const long& LB, const long& HB, const string& Genome, const bool& Reverse);

//run as GRC_overlap -i BlastResults.txt
int main (int argc, char* argv[]) {   //  Main is open
	char* InFile = argv[1]; //get the name of the blast test results file
	char* GenomeName= argv[2];//the name of the target genome
	char* GenomeFile= argv[3];//the name of the fna file
	string Matrix=argv[4];//the matrix used for blast
	string Genome; //for storing the genome
	int Status=0;//status variable

	ifstream In; //input for the blast results
	In.open(InFile); //open the blast input

	list<AARecord> RecordList;//storage for all of the records
	RecordMap PositionMap; //the initial map	 for the aa records ordered by HighBase of the orf
	//RecordMap::_Pairib MapIt;//for keeping track of insertions into the map

	map<string,AARecord*> HitList; //the list of AARecord pointers to ORFs aka pGenes that have hits (hash by ID)
	RecordMap::iterator MapIt;//pair for keeping track of insertions into the record map

	list<AARecord*> LoserList; //the list of the eliminated records
	list<AARecord*> WinnerList; //the final list of highest non overlaping blast hits
	CompeteMap KOMap;//create map for tracking win/lose relationships

	double Lambda;//constant in MaxBit calculation
	double K;//constant in MaxBit calculation
	map <string,int> ConsValue;//map for storing the max value of a conserved amino acid
	Status= InitCodes(Lambda, K, Matrix, ConsValue);//read in the values.
	if(Status!=0){
		cout<<"\nexiting\n";
		return 1;
	}
	

	string GName;//for storing input and creating output name from it
	GName=GenomeName;//convert input to string

	//Create the names of the Positive and negative output files
	string Negatives=GName;
	string Positives=GName;
	string Neg=".Neg";
	string Pos=".Pos";
	Negatives+=Neg;//create false/true negatives filename
	Positives+=Pos;//create fasle/treu positives filename



	ifstream In2;//ofstream operator for reading in the genomic sequence
	In2.open(GenomeFile);//open up the translated file
	
	string GenomeID;//for reading in the id
	getline(In2,GenomeID); //get line for description
	string Seq;//for reading in the sequence
	Genome="";//Initialize to empty string

	while(In2){//read in the genome file

		In2>>Seq; //read in the genomic sequence
		LowerTheCase(Seq);
		Genome+=Seq;//concat. each line to Genome

		/*FindIt=HitList.find(ID.substr(1,(ID.length())-1));

		if(FindIt!=HitList.end()){//if its found then its a hit
			FindIt->second->Sequence=Seq;//assign the sequence
		}*/
	}
	In2.close();//close the input stream

	string OldID="none";
	long OldQAS=-1;
	double OldBit=0;
	string OldOrg="none";
	string OldHID="none";
	int OldCount=0;
	map<string,double> MaxBitMap;//for storing and retrieving MaxBitValues
	map<string,double>::iterator FindMB;
	while (In){//read in the input file
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
		string Hit;
		long Length;
		double Bit;
		double PercentIdent;
		string ES;//escore
		long ALength; //the length of the alignment
		string HitOrg;//organism from which db hit comes
		string HitID;//ID from which db hit comes
		

		In >>ID; //read in the ID
		
		if (!In){//if the input hit the EOF
			break;//break from the loop
		}

		In >>Start; //read in the Start position
		In >>Stop; //read in the Stop position
		getline(In,HitID,'\t'); //skip next tab
		getline(In,HitID,'\t'); //get line for hit ID/No_hits
		long LowBase;
		long HighBase;
		double MaxBit;//the maximum bit score possible for a given sequence
		double LowComplexity=0;//fraction of low complexity AA's as determined by SEG in fsa-blast

		if (HitID =="No_hits"){//open consq.
			//In>>ES; //read in the delimiter
			//Insert Record into Initial RecordMap
			In>>LowComplexity;
			if (Start<Stop){
				LowBase=Start;
			}
			else {
				LowBase=Stop;
			}
			RecordList.push_back(AARecord(ID, Start, Stop, HitID, LowComplexity));//add record to the list
			//PositionMap.insert(RecordMap::value_type(LowBase, &RecordList.back()));
			//EditList.push_back(&((MapIt.first)->second));//add a pointer to the location of the record
		}//close consq.
		else {//there is a hit
			getline(In,Hit,'\t'); //get line for hit description
			getline(In,HitOrg,'\t'); //organism from which the hit comes
			In>>PercentIdent;
			In>>ALength; //alignment length
			In>>SubjectLen; //subject length
			In>>Mismatches;
			In>>GapOpens;
			In>>QAlignStart;
			In>>QAlignStop;
			In>>SAlignStart;
			In>>SAlignStop;
			In>>ES;
			In>>Bit;
			In>>LowComplexity;//ONLY needs to be entered once per query ID since all is query sequence dependent
			bool Rev=Start>Stop;
			long OrigStart=Start;//for start searching purposes
			//boolean value to determine whether any new information is being contributed
			//TO DO: Add org and HitID check to bool Old and update structure for storing
			//same alignment regions with multiple organisms,functions,and hit_ID's
			
			//if it doesn't have align start further in or as high a bit score
			bool Old;
			if(ID==OldID){
				Old=(Bit==OldBit&&QAlignStart==OldQAS);
			}
			else{
				Old=false;
				MaxBitMap.clear();//no more max bit calculations for that stop
			}
			OldID=ID;//update old values
			OldOrg=HitOrg;
			OldHID=HitID;
			OldBit=Bit;
			OldQAS=QAlignStart;
	
			

			//See if start can be adjusted then set lowbase
			//To Do:read in all hits for same query ID then run adjust start on each
			//After this is done a CalcMaxBit function that calculates each MaxBit
			//for each subject in 1 move across the genome
			if(!Old){
				AdjustStart1(Start,Stop,QAlignStart,Rev,Genome);//attempt adjust start site to alignment
				if (!Rev){
					LowBase=Start;
					HighBase=Stop;
				}
				else {
					LowBase=Stop;
					HighBase=Start;
				}
				//Add the original start site
				string LowHigh=ltos(LowBase)+","+ltos(HighBase)+","+btos(Rev);
				FindMB=MaxBitMap.find(LowHigh);//see if the maxbit has been previously calc.
				if(FindMB==MaxBitMap.end()){
					MaxBit=CalcMaxBit(LowBase,HighBase,Genome,ConsValue,Rev, Lambda, K);
					MaxBitMap.insert(map<string,double>::value_type(LowHigh,MaxBit));
				}
				else{//else its found
					MaxBit=FindMB->second;
				}
				
				IDMap::iterator FindID;
				FindID=HitList.find(ID);
				if(FindID!=HitList.end()){//if the orf already has a hit
					//Add Subject orf to Record
					FindID->second->AddPrimary(Start, Stop, Hit, Bit, ES, SubjectLen, ALength, QAlignStart, QAlignStop, MaxBit, HitID, HitOrg);
				}
				else{//else there has not been a hit for this orf so far
					RecordList.push_back(AARecord(ID, Start, Stop, Hit, LowComplexity, Bit, ES, SubjectLen, ALength, QAlignStart, QAlignStop, MaxBit, HitID, HitOrg));
					HitList.insert(IDMap::value_type(ID,&RecordList.back()));//add pointer to the record to the hit list
				}
				while(AdjustStart2(Start,OrigStart,Stop,QAlignStart,Rev,Genome)){//find all start sites from aligned region back to original
					if (!Rev){
						LowBase=Start;
						HighBase=Stop;
					}
					else {
						LowBase=Stop;
						HighBase=Start;
					}
					LowHigh=ltos(LowBase)+","+ltos(HighBase)+","+btos(Rev);
					FindMB=MaxBitMap.find(LowHigh);//see if the maxbit has been previously calc.
					if(FindMB==MaxBitMap.end()){
						MaxBit=CalcMaxBit(LowBase,HighBase,Genome,ConsValue,Rev, Lambda, K);
						MaxBitMap.insert(map<string,double>::value_type(LowHigh,MaxBit));
					}
					else{//else its found
						MaxBit=FindMB->second;
					}
					
					FindID=HitList.find(ID);
					if(FindID!=HitList.end()){//if the orf already has a hit
						//Add Subject orf to Record
						FindID->second->AddSecondary(Start, Stop, Hit, Bit, ES, SubjectLen, ALength, QAlignStart, QAlignStop, MaxBit, HitID, HitOrg);
					}
					else{//else there has not been a hit for this orf so far
						RecordList.push_back(AARecord(ID, Start, Stop, Hit, LowComplexity, Bit, ES, SubjectLen, ALength, QAlignStart, QAlignStop, MaxBit, HitID, HitOrg));
						HitList.insert(IDMap::value_type(ID,&RecordList.back()));//add pointer to the record to the hit list
					}
				}//close find start sites
			}//close if not Old
			else{
				OldCount++;
			}
			//In>>ES; //read in delimiter
		}//close hit
		//OldID=ID;//remember this id in next iteration
	}// close read input
	In.close();

	Genome.clear();//clear the genome (no longer needed)

	//create position map
	for(list<AARecord>::iterator PosIt=RecordList.begin(); PosIt!=RecordList.end(); PosIt++){
		PosIt->SwitchRep();
		PositionMap.insert(RecordMap::value_type(PosIt->LowBase, &(*PosIt))); //Add to position map
	}//close for loop

	
	
	ofstream ChkOut;
	//ChkOut.open(Positives.c_str());
	//ChkOut.close();
	//DumpList(InitList);//print out the orfs from initlist

	//Adjust starts according to aligned regions


	//compare the ORFs
	Compare(PositionMap,WinnerList,LoserList, KOMap);
	
	//Tally results
	int NumOrfs=0;
	for (list<AARecord*>::iterator It1 =WinnerList.begin(); It1!=WinnerList.end(); It1++ ){
		if((**It1).ReportHit()){ NumOrfs++;}
	}
	cout<<"******************GRC v0.01******************"<<"\n";
	cout<<"Total # of initial ORFS created:   "<<PositionMap.size()<<"\n";
	cout<<"Total # non-overlapping orfs:   "<<WinnerList.size()<<"\n";
	cout<<"Total # annotated:   "<<NumOrfs<<"\n\n";
	

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
int Compare(RecordMap& EList, list<AARecord*>& FList, list<AARecord*>& CList, CompeteMap& KOMap){//open defintion
	bool FirstVictim=true;
	for (RecordMap::iterator It1 =EList.begin(); It1!=EList.end(); It1++ ){
		//FirstVictim=true;
		int OverLen=-5555; //overlapping length of the two orfs
		string Spot="nada";//Key for finding the Competition Record
		if(!(It1->second)->Dead()){// if It1 isn't on a dead record
			for (RecordMap::iterator It2 =It1; It2!=EList.end()&& OverLen!=0; It2++ ){//inner loop only needs to go through list from position It1
				OverLen=-5555;

				if(It1!=It2 && !(It2->second->Dead()) && !(It1->second->Dead())){//if not the same
					
					OverLen=It1->second->Overlap(*(It2->second));//get the amount the two overlap if at all
				//	if(((It1->second.ID=="T2432")||(It2->second.ID=="T2432"))){
				//		cout<<"Whats going on?";
				//	}

					bool KO=false;//default is to not knock one out
					if(OverLen!=0){//if they overlap
						KO=It1->second->KnockOut(*(It2->second), OverLen);//determine if one should be knocked out
					}
					if(KO){//if one of the ORFs should be knocked out

						if (*(It1->second)>*(It2->second) && It2->second->Incompatible(*(It1->second))){//If record 1 has a higher score
							Spot=It1->second->ID;
							CompeteMap::iterator TempC=KOMap.find(Spot);
							if (TempC!=KOMap.end()){
								(TempC->second).AddLoser((It2->second));
							}
							else {
								KOMap.insert(CompeteMap::value_type(Spot, Compete((It1->second),(It2->second))));
							}
								
							CList.push_back((It2->second));//loser
							It2->second->KnockOut();//then record 2 is out
						}
						else if(*(It1->second)<*(It2->second) && It1->second->Incompatible(*(It2->second))){//if record 1 has a lower score
							Spot=((It2)->second->ID);
							CompeteMap::iterator TempC=KOMap.find(Spot);
							if (TempC!=KOMap.end()){
								TempC->second.AddLoser((It1->second));
							}
							else {
								KOMap.insert(CompeteMap::value_type(Spot, Compete((It2->second), (It1->second))));
							}
								
							CList.push_back((It1->second));//loser
							It1->second->KnockOut();//then record 1 is out
							break;//break from the inner loop since record 1 is dead
						}
						else if(*(It1->second)^*(It2->second) && It2->second->Incompatible(*(It1->second))){//else if It1 is longer than It2
							Spot=((It1->second)->ID);
							CompeteMap::iterator TempC=KOMap.find(Spot);
							if (TempC!=KOMap.end()){
								(TempC->second).AddLoser((It2->second));
							}
							else {
								KOMap.insert(CompeteMap::value_type(Spot, Compete((It1->second),(It2->second))));
							}
							CList.push_back((It2->second));//loser
							It2->second->KnockOut(); //record 2 is out
						}
						else if(*(It2->second)^*(It1->second)  && It1->second->Incompatible(*(It2->second))){//else It2 is longer
							Spot=((It2->second)->ID);
							CompeteMap::iterator TempC=KOMap.find(Spot);
							if (TempC!=KOMap.end()){
								(TempC->second).AddLoser((It1->second));
							}
							else {
								KOMap.insert(CompeteMap::value_type(Spot, Compete((It2->second), (It1->second))));
							}
							CList.push_back((It1->second));
							It1->second->KnockOut();
							break; //break from inner loop
						}
						else{//else its a tie
						}
					}//close KO orfs
				}//close if not the same	
			}//close inner loop
			if (!It1->second->Dead()){//if the record still hasn't been eliminated
				FList.push_back((It1->second));//add it to the final list
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

std::ostream& operator<<(std::ostream& ACOut, const AARecord& AC){
	ACOut<<"PPCG**"<<"\n";
	ACOut<<"ID:\t"<<AC.ID<<"\n";
	ACOut<<"ORF Start:\t"<<AC.Start<<"\n";
	ACOut<<"ORF Stop:\t"<<AC.Stop<<"\n";
	ACOut<<"ORF Length(nt):\t"<<AC.QLength<<"\n";
	if (AC.Reverse){
		ACOut<<"ORF Strand:\t"<<'-'<<"\n";
	}
	else 	ACOut<<"ORF Strand:\t"<<'+'<<"\n";

	ACOut<<"DB Hit:\t"<<AC.Hit<<"\n"; //print hit description or no hit line
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
}


//ofstream overload for Compete class
std::ostream& operator<<(std::ostream& Out, const Compete& C){
	Out<<(C.Winner)->ID;
	for (list<AARecord*>::const_iterator It1 =C.Losers.begin(); It1!=C.Losers.end(); It1++ ){
		Out<<"\t"<<(*It1)->ID;
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
	ChkOut<<AC->QLength<<"\t";
	if (AC->Reverse){
		ChkOut<<'-'<<"\t";
	}
	else 	ChkOut<<'+'<<"\t";

	ChkOut<<AC->Hit<<"\t"; //print hit description or no hit line
	if(!AC->Blank){//if there is a hit
		double ALength=AC->ALength;
		double OLength=AC->QLength;
		double HLength=AC->HLength;
		ChkOut<<AC->HitID<<"\t";
		ChkOut<<AC->HitOrg<<"\t";
		ChkOut<<AC->Bit<<"\t";
		ChkOut<<AC->EScore<<"\t";
		ChkOut<<AC->HLength<<"\t";
		ChkOut<<(ALength/(OLength/3))*100<<"\t";
		ChkOut<<(ALength/HLength)*100<<"\n";
	}

	else{
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

//function to lower the case of all characters in a string
int LowerTheCase(string & Seq){
	for(int i=0; i<Seq.length(); i++){
		Seq[i]=tolower(Seq[i]);
	}
	return 1;
}//close definition
	
//Function to Adjust the start site based on where the alignment begins
//returns bool stop searching or not
bool AdjustStart2(long& St, const long& OSt, const long& Sp, const long& QAS, const bool& Reverse, const string& Genome){//open definition
	long Start=St;
	long Stop=Sp;
	long OrigStart=OSt;
	long QAlignStart=QAS;
	int StartSearch;
	double StartScore=0;
	double MaxStartScore=0;
	double Travel=0;
	double NuclDist=(3*QAlignStart);
	int Halt=0;//defines when the search has gone back to original

	if(QAlignStart<2){//if there is no room to search for a start between the aligned region and current start
		return false;//stop the search
	}
	else if(Reverse){//if the pGene is reversed
		if(Start==OrigStart){
			StartSearch=(Start-1)-(3*QAlignStart);// start search position
			Halt=(QAlignStart*3);
		}
		else{//start from the position of the last start found
			StartSearch=Start+2;//next codon 
			Halt=OrigStart-Start-3;//room left between orig and current search position
		}

		for (int s=0; s<=Halt; s=s+3){//search codons in the upstream direction
			if(s%3==0){//if its the next codon
				if(ReverseStart(Genome.substr(StartSearch+s-2, 3))){//if its a start
					St=StartSearch+s+1;
					return St!=OrigStart;//if back to original start, stop searching
					//start score calc. from codon and distance traveled to orig. start
					//Travel=(1-(s/NuclDist));
					//StartScore=CalcSS(Genome.substr(StartSearch+s-2, 3), Travel);
					//if (StartScore>MaxStartScore){
                        			//St=StartSearch+s+1;
						//MaxStartScore=StartScore;
					//}//close max start score
				}//close is start
			}//close next codon
		}//close search codons
	}//close reverse pGene

	else{//else its not reversed
		if(Start==OrigStart){
			StartSearch=(Start-1)+(3*QAlignStart);
			Halt=(QAlignStart*3);
		}
		else{
			StartSearch=Start-4;//start at the next codon
			Halt=Start-OrigStart-3;//room left between orig and current search position
		}

		for (int s=0; s<=Halt; s=s+3){//search 3 codons in the upstream direction
			if(s%3==0){//if its the next codon
				if(ForwardStart(Genome.substr(StartSearch-s, 3))){//if its a start
					St=StartSearch-s+1;
					return St!=OrigStart;
					//Travel=(1-(s/NuclDist));
					//StartScore=CalcSS(Genome.substr(StartSearch-s, 3), Travel);
					//if (StartScore>MaxStartScore){
                        			//St=StartSearch-s+1;
						//MaxStartScore=StartScore;
					//}//close if max score
				}//close if start
			}//close if next codon
		}//close search next codons
	}//close not reversed
	return false;
}//close defintion


//This function is designed to check if submitted string is reverse codon
bool ReverseStart(const string& Codon){//open definition
	if(Codon=="cat"||Codon=="cac"||Codon=="caa"){
		return true;
	}
	else return false;
}

//This function is designed to check if submitted string is reverse codon
bool ForwardStart(const string& Codon){//open definition
	if(Codon=="atg"||Codon=="gtg"||Codon=="ttg"){
		return true;
	}
	else return false;
}

//This function is intended to give a likelihood of correctness relative to other starts
//this is done by considering the different start codons to have different inherent probab.
//this is then impacted by mult. times the 1- percent distance away from alignment
double CalcSS(const string& Codon, const double& Travel){//open definition
	if (Codon=="atg"||Codon=="cat"){
		return (.77*Travel);
	}
	else if (Codon=="gtg"||Codon=="cac"){
		return (.14*Travel);
	}
	else if (Codon=="ttg"||Codon=="caa"){
		return (.08*Travel);
	}
	return 0;
}

//This function returns the complement of a nucleotide passed as a parameter
char Complement(const char& Base){//open definition
	switch(Base){
		case 'a':
			return 't';
			
		case 't':
			return 'a';
			
		case 'c':
			return 'g';
			
		case 'g':
			return 'c';
		default: return '*';
	}
}//close definition

//This function returns the reverse complement of a given sequence
string ReverseComp(const string& Forward){

	string Comp="";
	for(int s= int(Forward.length())-1; s>=0; s--){
		Comp+=Complement(Forward[s]);
	}
	return Comp;
}//close definition



	int InitCodes(double& Lambda, double& K, const string& Matrix, map <string,int>& ConsValue){//open definition
	//Get the Genetic Codes for translation
		ifstream InCode;
		InCode.open(Matrix.c_str());//open matrix specified
		if(!InCode.is_open()){
			cout<<"Error opening MaxMatrix file";
			return 1;
		}
		string Delim;
		string OpCase="";
		int AA; //Max value for conservation of the amino acid
		//InCode>>Delim;

		while(InCode){
			InCode>>Delim;
			if(isupper(Delim[0])){//create lower case in map
				OpCase=tolower(Delim[0]);
				for(int t=1; t<Delim.size(); t++){
					OpCase+=tolower(Delim[t]);
				}
			}
			else {
				OpCase=Delim;
			}

			InCode>>AA; //Read in the AA conservation value
			//CodeArray[TableNum].insert(map<string,char>::value_type(Delim,AA));
			ConsValue.insert(map<string,int>::value_type(OpCase,AA));
			
		}//close while loop
		InCode.close();

		//Set constants
		if(string::npos!=Matrix.find("MaxB62",0)){//if the matrix is BLOSUM62
			Lambda=0.267;
			K=0.041;
		}

		else if(string::npos!=Matrix.find("MaxB80",0)){//if the matrix is BLOSUM80
			Lambda=0.299;
			K=0.071;
		}

		else if(string::npos!=Matrix.find("MaxB45",0)){//if the matrix is BLOSUM45
			Lambda=0.195;
			K=0.032;
		}

		else if(string::npos!=Matrix.find("MaxP30",0)){//if the matrix is PAM30
			Lambda=0.309;
			K=0.15;
		}

		else if(string::npos!=Matrix.find("MaxP70",0)){//if the matrix is PAM30
			Lambda=0.270;
			K=0.060;
		}

		else{//default B62
			Lambda=0.267;
			K=0.041;
		}


		return 0;
	}//close definition

	//Function to Adjust the start site based on where the alignment begins
int AdjustStart1(long& St, const long& Sp, const long& QAS, const bool& Reverse, const string& Genome){//open definition
	long Start=St;
	long Stop=Sp;
	long QAlignStart=QAS;
	int StartSearch;
	double StartScore=0;
	double MaxStartScore=0;
	double Travel=0;
	double NuclDist=(3*QAlignStart);

	if(QAlignStart<2){
		return 1;
	}
	else if(Reverse){//if the pGene is reversed
		StartSearch=(Start-1)-(3*QAlignStart);// start search position

		for (int s=0; s<9; s++){//search codons in the upstream direction
			if(s%3==0){//if its the next codon
				if(ReverseStart(Genome.substr(StartSearch+s-2, 3))){//if its a start
					St=StartSearch+s+1;
					return 1;
					//}//close max start score
				}//close is start
			}//close next codon
		}//close search codons
	}//close reverse pGene

	else{//else its not reversed
		StartSearch=(Start-1)+(3*QAlignStart);

		for (int s=0; s<(9); s++){//search 3 codons in the upstream direction
			if(s%3==0){//if its the next codon
				if(ForwardStart(Genome.substr(StartSearch-s, 3))){//if its a start
					St=StartSearch-s+1;
					return 1;
				}//close if start
			}//close if next codon
		}//close search next codons
	}//close not reversed
	return 1;
}//close defintion


//This function calculates the entropy for the section of the genome specified
double GetEntropy(const long& LB, const long& HB, const string& Genome, const bool& Reverse){//open definition
	long Length=HB-LB+1;
	long Begin=LB-1;
	string Run="./entropy/entropy-score ";
	
	if(Reverse){//if the pGene is reversed
		Run=Run+ReverseComp(Genome.substr(Begin, Length));
		return system(Run.c_str());//get entropy reverse complement
	}
	else{//not reverse complement
		Run=Run+ReverseComp(Genome.substr(Begin, Length));
		return system(Run.c_str());
	}
}
	


	//This function calculates the maximum possible bit score for a given segment of the genome
	//Takes parameters LowBase, HighBase, Genome, Map of conservation values, lambda, and K blast values
	//Need to clean up this coordinate to string conversion +1 -1 stuff
double CalcMaxBit(const long& LB, const long& HB, const string& Genome, map <string,int>& ConsValue, const bool& Reverse, const double & Lambda, const double & K){
		double MaxBit=0;
		double RawBit=0;
		long StartSearch=LB-1;//Subtract one to convert to string coordinates
		long Length=HB-LB;
		string Codon;
		map<string,int>::iterator FindIt;
		//even though its in reverse go forward through the sequence and get reverse complement
		if(Reverse){
			for (long s=0; s<Length+1; s=s+3){//calc max possible score
				if(s%3==0){//if its the next codon
					Codon=ReverseComp(Genome.substr(StartSearch+s, 3));//get reverse complement
					FindIt=ConsValue.find(Codon);//find the score of this codon
					if(FindIt!=ConsValue.end()){
						RawBit=RawBit+FindIt->second;//add up score
					}
				}
			}//close max loop
		}//close if Reverse

		else {
			for (long s=0; s<Length+1; s=s+3){//calc max possible score
				if(s%3==0){//if its the next codon
					Codon=Genome.substr(StartSearch+s, 3);//codon
					FindIt=ConsValue.find(Codon);//find the score of this codon
					if(FindIt!=ConsValue.end()){
						RawBit=RawBit+(FindIt->second);//add up score
					}
				}
			}//close max loop
		}//close not Reverse

		MaxBit=((RawBit*Lambda)-log(K))/M_LN2;
		return MaxBit;
	}//close definition





