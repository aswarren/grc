/*GRC.Overlap
 * This component of the GRC program finds the ORFs that overlap from the BLAST results file
 * and eliminates all but the one with the best BLAST hit.*/


#include <list>
#define _USE_MATH_DEFINES
#include <math.h>
#include <map>
#include <utility>
#include "AARecord.h"
#include "Compete.h"
#include "CalcPack.h"
#include "GO.h"
#include "FastaRead.h"




using std::list;
using std::cout;
using std::ios;
using std::fstream;
using std::map;
using std::pair;


typedef map<string, string> SSMap;

int DumpList(list<AARecord>& InitList);
int DumpList(list<AARecord*>& InitList, string PosName, const int& GFMin);
int Nulify(list<AARecord*>& InitList, list<AARecord*>& InitList2);
int Compare(RecordMap& PositionMap, list<AARecord*>& WinnerList, list<AARecord*>& LoserList, CompeteMap& KOMap);
int EntropyFilter(list<AARecord>& RecordList, list<AARecord*>& LoserList, CompeteMap& KOMap, const double& EntCutoff, const double& BitCutoff);//removes orfs that are "winners" that have high entropy
int SmallFilter(list<AARecord>& RecordList, list<AARecord*>& LoserList, CompeteMap& KOMap, const double& EntCutoff, const double& BitCutoff);
void DisplayKO(ostream& Out, CompeteMap& KOMap, const int& GFMin);
int RefreshRecords(list<AARecord>& RecordList, CalcPack& CP);
int TrainEDP(list<AARecord*>& WinnerList, list<AARecord*>& LoserList, CalcPack& CP);
int ProcessID(string& ID, long& Start, long& Stop, long& Offset);
int GetBlastResults(const string BlastFile, list<AARecord>& RecordList, map<string, AARecord*> & HitList, CalcPack& InfoPack);
SSMap ParseCommandLine(const int& ac, char* const av[]);
int SetEFilter(StringSet& ECodeFilter, const string& ECodeTxt);
int FastaPrint(list<AARecord*>& RecList, const string& FilePrefix, const int& GFMin, CalcPack& InfoPack, const string& OptFasta);
string ConstructHeader(const string& prefix, const string& delim, vector<string>& Terms);


//run as GRC_overlap
int main(int argc, char* argv[]) {   //  Main is open
    
    bool ICAon=false;
    SSMap Options;
    Options=ParseCommandLine(argc, argv);
    if(Options.size()==0){
        cerr<<"Usage: grc_annotate -b [blast results file] -o [output name] -g [genome file] -m [blast matrix file] -t [translation tables file] -n [trans. table num.] -s [start codon file] -l [min. gene length] OPTIONAL -y [Gene Ontology file] -a [Use Ontology MF, BP, CC (e.g. 'mbc')] -f [Filter evidence codes (e.g. 'IEA ND')] -p [N or A (nucleotide or amino acid fasta)] -x [% subject aligned cutoff] -c (create consensus annotations) -d [minimum depth of GO term]\n";
        return -1;
    }
    string BlastFile = Options.find("-b")->second; //get the name of the blast test results file
    string GenomeName= Options.find("-o")->second;//the name of the target genome
    string GenomeFile= Options.find("-g")->second;//the name of the fna file
    string Matrix=Options.find("-m")->second;//the matrix used for blast
    string TransFile=Options.find("-t")->second;//file used for translating sequences in entropy calculations
    int TableNum=stoi(Options.find("-n")->second);
    string StartFile=Options.find("-s")->second;//file used to specify which stops are used
    string GOFile="none";
    stringstream Convert;
    int GFMinLength=300;
    Convert<<string(Options.find("-l")->second);
    Convert>>GFMinLength;
    GO Ontology;//ontology object for storing Gene Ontology
    CalcPack InfoPack(Matrix, GenomeFile, TransFile, TableNum);//create information package
    InfoPack.SetStarts(StartFile);
    StringSet ECodeFilter;
    int MinDepth=-1;//minimum depth for GO term filtering
    double sbj_cutoff=0;
    
    SSMap::iterator Oit=Options.end();
    if(Options.find("-y")!=Options.end()){//if GO.obo specified
        GOFile=Options.find("-y")->second;
        bool MFunc, BProc, CComp;
        MFunc=BProc=CComp=true;
        Oit=Options.find("-a");
        if(Oit!=Options.end()){
            Ontology.FindWarningOff();
            MFunc=BProc=CComp=false;
            string aOptions=Oit->second;
            if(aOptions.find("m",0)!=string::npos){
                MFunc=true;
            }
            if(aOptions.find("b",0)!=string::npos){
                BProc=true;
            }
            if(aOptions.find("c",0)!=string::npos){
                CComp=true;
            }
        }
        Oit=Options.find("-f");
        if(Oit!=Options.end()){//if there is a filter option
            SetEFilter(ECodeFilter, Oit->second);
        }
        Oit=Options.find("-d");
        if(Oit!=Options.end()){
            MinDepth=stoi(Oit->second);
        }
        ifstream GOIn;
        GOIn.open(GOFile.c_str());
        Ontology.ReadOBO(&GOIn, MFunc, BProc, CComp);
        GOIn.close();
        InfoPack.SetGOAccess(&Ontology);//set go access pointer
        if(Options.find("-c")!=Options.end()){//if GO.obo specified
            ICAon=true;
        }
    }
    //check what type of fasta output
    string OptFasta="x";
    Oit=Options.find("-p");
    if(Oit!=Options.end()){
        OptFasta=Oit->second;
    }
    
    Oit=Options.find("-x");
    if(Oit!=Options.end()){
        sbj_cutoff=stod(Oit->second);
    }
    
    
    list<AARecord> RecordList;//storage for all of the records
    RecordMap PositionMap; //the initial map for the aa records ordered by HighBase of the orf
    //RecordMap::_Pairib MapIt;//for keeping track of insertions into the map
    
    map<string, AARecord*> HitList; //the list of AARecord pointers to ORFs aka pGenes that have hits (hash by ID)
    RecordMap::iterator MapIt;//pair for keeping track of insertions into the record map
    
    list<AARecord*> TrainingLosers;//list used to store losers from the initial comparison and train new EDP
    list<AARecord*> TrainingWinners;
    CompeteMap TrainingKO;//create map for tracking win/lose relationships
    
    list<AARecord*> LoserList; //the list of the eliminated records
    list<AARecord*> WinnerList; //the final list of highest non overlaping blast hits
    CompeteMap KOMap;//create map for tracking win/lose relationships
    
       
    
    
    //Create the names of the Positive and negative output files
    string Negatives=GenomeName;
    string Positives=GenomeName;
    string Neg=".Neg";
    string Pos=".Pos";
    Negatives+=Neg;//create false/true negatives filename
    Positives+=Pos;//create fasle/treu positives filename
    
    
    
    GetBlastResults(BlastFile, RecordList, HitList, InfoPack );
    
    //create position map
    //and switch the start site to the one with the highest conservation
    //and calculate entropy
    for(list<AARecord>::iterator PosIt=RecordList.begin(); PosIt!=RecordList.end(); PosIt++){
        
        if(PosIt->HasHit()){//if the query has a hit
            PosIt->UpdateScores();//now that all blast scores have been read in for each record create priority queues for start sites
        }
        PositionMap.insert(RecordMap::value_type(PosIt->ReportLowBase(), &(*PosIt))); //Add to position map
    }//close for loop
    
    //DumpList(InitList);//print out the orfs from initlist
    
    //compare the ORFs and develop a set of training orfs for creating
    //new entropy density profiles for coding and non coding genes
    Compare(PositionMap, TrainingWinners, TrainingLosers, TrainingKO);
    //activate the use of small orf EDP if sufficient number
    
    //train new EDP's
    TrainEDP(TrainingWinners, TrainingLosers, InfoPack);
    //refresh all records using new EDP
    RefreshRecords(RecordList, InfoPack);
    
    
    
    //calculate avg. entropy of winners with hits from training data using new entropies
    
    double AvgEntropy=0;//The average entropy
    int NumForAvg=0;
    for (list<AARecord*>::iterator AvgIt =TrainingWinners.begin(); AvgIt!=TrainingWinners.end(); AvgIt++ ){
        if((*AvgIt)->ReportLength()>300){
            NumForAvg++;
            AvgEntropy+=(*AvgIt)->ReportEntropy();
        }
    }
    AvgEntropy=AvgEntropy/double(NumForAvg);//finish avg entropy calc. by dividing by number of orfs with hits
    
    //Calculate the std deviation
    double Variance=0;
    for (list<AARecord*>::iterator EntIt =TrainingWinners.begin(); EntIt!=TrainingWinners.end(); EntIt++ ){
        //if((*EntIt)->HasHit() && (*EntIt)->ReportBitFrac()>ConservCutoff){//if this orf has a hit
        if((*EntIt)->ReportLength()>300){
            double Diff=(*EntIt)->ReportEntropy()-AvgEntropy;//get difference between mean and current
            Variance+=(Diff*Diff);//square the difference of mean and current and add to total variance
        }
        //}
    }
    double EntropyDev=sqrt((Variance/double(NumForAvg)));//calculate std deviation
    double EDRCutoff=1.0;
    double EDRStrict=1.0;
    //if for some reason the EDRCutoff is too lenient
    //that 1.0 is > Avg+2*Std Dev.
    if(EDRCutoff>AvgEntropy+(2*EntropyDev)){
        EDRCutoff=AvgEntropy+(2*EntropyDev);
    }
    if(EDRStrict>AvgEntropy+(EntropyDev)){
        EDRStrict=AvgEntropy+(EntropyDev);
    }
    double BitCutoff=.50;//Bit fraction cutoff for being evaluated by entropy Bit/MaxBit
    double BitStrict=.80;//bit fraction necessary for short alignments
    int NumFiltered=EntropyFilter(RecordList, LoserList, KOMap, EDRCutoff, BitCutoff);//filter the orfs based on entropy
    NumFiltered+=SmallFilter(RecordList, LoserList, KOMap, EDRStrict, BitStrict);
    
    //compare the ORFs and remove the ones that conflict due to overlap
    Compare(PositionMap, WinnerList, LoserList, KOMap);
    
    /*if(GOFile!="none"){//if there is a GO file specified
        //cout<<"Reading in the ontology file "<<GOFile<<"\n";
        ifstream GOIn;
        GOIn.open(GOFile.c_str());
        Ontology.ReadOBO(&GOIn, true, true, true);
        GOIn.close();
        InfoPack.SetGOAccess(&Ontology);//set go access pointer
    }//close if there is a obo file*/
    
    //Tally results
    //cout<<"Average entropy is\t"<<AvgEntropy<<"\n";
    //cout<<"Std Dev of entropy is\t"<<EntropyDev<<"\n";
    int NumWinHits=0;
    if(InfoPack.GOAccess!=NULL){
        for (list<AARecord*>::iterator CountIt =WinnerList.begin(); CountIt!=WinnerList.end(); CountIt++ ){
            (*CountIt)->BuildGOTerms(InfoPack, ECodeFilter, MinDepth);//post process all GO information
            if(ICAon){
                (*CountIt)->GOAnalysis(InfoPack);//finalize GO annotations
            }
        }
    }
    
    for (list<AARecord*>::iterator CountIt =WinnerList.begin(); CountIt!=WinnerList.end(); CountIt++ ){
        if(sbj_cutoff!=0){
            (*CountIt)->ScreenFunction(sbj_cutoff);
        }
        if((*CountIt)->HasHit()){//if this orf has a hit
            NumWinHits++;
        }
    }
    
    
    cout<<"******************GRC Annotate******************"<<"\n";
    cout<<"Total # of initial ORFS created:   "<<PositionMap.size()<<"\n";
    cout<<"Total # non-overlapping orfs:   "<<WinnerList.size()<<"\n";
    cout<<"Total # annotated:   "<<NumWinHits<<"\n\n";
    cout<<"Number of orfs filtered from entropy\t"<<NumFiltered<<"\n";
    
    FastaPrint(WinnerList, GenomeName, GFMinLength, InfoPack, OptFasta);
    DumpList(WinnerList, Positives, GFMinLength);
    //DumpList(InitList);
    DumpList(LoserList, Negatives, GFMinLength);
    
    Nulify(WinnerList, LoserList);
    
    ofstream KOut;
    KOut.open("KnockList.txt");
    DisplayKO(KOut, KOMap, GFMinLength);
    KOut.close();
    InfoPack.ClearGenome();//clear the genome (no longer needed)
    
    return 0;
}


/*Parses the command line and returns a dictionary of strings initializing each
 * option-tag to its corresponding value
 * cerr<<"Usage: grc_overlap -b [blast results file] -o [output name] -g [genome file] -m [blast matrix file] -t [translation tables file] -l [min. gene length] OPTIONAL -b [Gene Ontology file] -a [Use Ontology MF, BP, CC (e.g. 'mbc')] -f [Filter evidence codes (e.g. 'IEA,ND') \n";
 */
SSMap ParseCommandLine(const int& ac, char* const av[]){
    SSMap Result;
    SSMap::iterator MarkIt;
    MarkIt=Result.end();
    for(int t=1; t<ac; t++){
        string Part=av[t];
        //if it is an option-tag
        if(Part[0]=='-'){
            MarkIt=Result.insert(SSMap::value_type(Part, "")).first;
        }
        else if(MarkIt==Result.end()){
            Result.clear();
            return(Result);
        }
        else if(MarkIt->second.size()==0){
            MarkIt->second=Part;
        }
        else{
            MarkIt->second+=" "+Part;
        }
    }
    if(Result.size()==0){
        return(Result);
    }
    if(Result.find("-b")==Result.end()){
        cerr<<"Missing parameter -b\n";
        Result.clear();
    }
    if(Result.find("-o")==Result.end()){
        cerr<<"Missing parameter -o\n";
        Result.clear();
    }
    if(Result.find("-g")==Result.end()){
        cerr<<"Missing parameter -g\n";
        Result.clear();
    }
    if(Result.find("-m")==Result.end()){
        cerr<<"Missing parameter -m\n";
        Result.clear();
    }
    if(Result.find("-t")==Result.end()){
        cerr<<"Missing parameter -t\n";
        Result.clear();
    }
    if(Result.find("-l")==Result.end()){
        cerr<<"Missing parameter -l\n";
        Result.clear();
    }
    if(Result.find("-s")==Result.end()){
        cerr<<"Missing parameter -s\n";
        Result.clear();
    }
    if(Result.find("-n")==Result.end()){
        cerr<<"Missing parameter -n\n";
        Result.clear();
    }
    return(Result);
}

//Reads in the BLAST results file and initializes:
//RecordList stores Records of ORFs and their blast results
//HitList hashes Pointers to those Records that have hits according to ORF ID.
//This function uses InfoPack which contains the genomic sequence and other useful functions
int GetBlastResults(string BlastFile, list<AARecord>& RecordList, map<string, AARecord*> & HitList, CalcPack& InfoPack){
    ifstream BlastIn; //input for the blast results
    BlastIn.open(BlastFile.c_str()); //open the blast input
    
    string OldID="none";
    long OldQAS=-1;
    double OldBit=0;
    string OldOrg="none";
    string OldHID="none";
    int OldCount=0;
    
    map<string, double>::iterator FindMB;
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
        long Offset=0;//the offset is the total number of bases for replicons that occur before the current replicon
        string ES;//escore
        long ALength; //the length of the alignment
        string HitOrg;//organism from which db hit comes
        string HitID;//ID from which db hit comes
        string HitGeneName;//If ptt is used this is the gene name from the ptt
        string HitSynonym;//synonym name from the ptt
        
        
        BlastIn >>ID; //read in the ID
        
        if (!BlastIn){//if the input hit the EOF
            break;//break from the loop
        }
        
        BlastIn >>Start; //read in the Start position
        BlastIn >>Stop; //read in the Stop position
        
        ProcessID(ID, Start, Stop, Offset); //process the ID to see if it indicates multiple replicons
        
        
        getline(BlastIn, HitID, '\t'); //skip next tab
        //getline(BlastIn,HitID,'\t'); //get line for hit ID/No_hits
        BlastIn>>HitID;
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
            //BlastIn>>LowComplexity;
            RecordList.push_back(AARecord());//add record to the list
            RecordList.back().InitRecord(InfoPack, ID, Start, Stop, HitID, Offset);
            //EditList.push_back(&((MapIt.first)->second));//add a pointer to the location of the record
        }//close consq.
        else {//there is a hit
            string HitID2="";
            getline(BlastIn, HitID2, '\t');
            HitID=HitID+HitID2;
            getline(BlastIn, HitGeneName, '\t');
            getline(BlastIn, HitSynonym, '\t');
            getline(BlastIn, Function, '\t'); //get line for hit description
            
            getline(BlastIn, HitOrg, '\t'); //organism from which the hit comes
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
            //BlastIn>>LowComplexity;//ONLY needs to be entered once per query ID since all is query sequence dependent
            long OrigStart=Start;//for start searching purposes
            //boolean value to determine whether any new information is being contributed
            //TO DO: Add org and HitID check to bool Old and update structure for storing
            //same alignment regions with multiple organisms,functions,and hit_ID's
            
            //if it doesn't have align start further into the interior and does not have the high bit score for the same (query, subject) pair
            //then it does not contain new information and is classified as Old
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
                    FindID->second->AddPrimary(InfoPack, Start, Stop, HitID, HitGeneName, HitSynonym, Bit, ES, SubjectLen, ALength, QAlignStart, QAlignStop, Function, HitOrg);
                }
                else{//else there has not been a hit for this orf so far. so create a new AARecord
                    RecordList.push_back(AARecord());
                    RecordList.back().InitRecord(InfoPack, ID, Start, Stop, HitID, Offset, HitGeneName, HitSynonym, Bit, ES, SubjectLen, ALength, QAlignStart, QAlignStop, Function, HitOrg);
                    HitList.insert(IDMap::value_type(ID, &RecordList.back()));//add pointer to the record to the hit list
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
                                    KOMap.insert(CompeteMap::value_type(Spot, Compete((It1->second), (It2->second))));
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
                                    KOMap.insert(CompeteMap::value_type(Spot, Compete((It1->second), (It2->second))));
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
int DumpList(list<AARecord*>& InitList, string PosName, const int& GFMin){//open definition
    ofstream ChkOut;
    ChkOut.open(PosName.c_str());
    ChkOut<<"ID\tStart\tStop\tLength(nt)\tStrand\tEDR\tDBID1\tDBID2\tDBID3\tDBOrg\tGOTerms(Conf)\tDescription(Conf)\tBit\tEScore\tHitLength\t%QueryAligned\t%HSPAligned\n";
    
    for (list<AARecord*>::iterator It1 =InitList.begin(); It1!=InitList.end(); It1++ ){
        if(*It1!=NULL){
            //cout<<**It1;//print out the Records
            if((*It1)->ReportLength()>=GFMin){//only print those records over the minimum gene length
                ChkOut<<*It1;//print to minout
            }
        }
        else ChkOut<<"NULL record: grc_annotate error\n";
    }
    ChkOut.close();
    return 0;
}


//function to print out master list of ORFS
int FastaPrint(list<AARecord*>& RecList, const string& FilePrefix, const int& GFMin, CalcPack& InfoPack, const string &OptFasta){//open definition
    string Delimiter=" ";
    string Prefix="lcl|";
    if(OptFasta.find("N",0)!=string::npos){
        string FileName=FilePrefix+".ffn";
        ofstream ChkOut;
        ChkOut.open(FileName.c_str());

        for (list<AARecord*>::iterator It =RecList.begin(); It!=RecList.end(); It++ ){
            if(*It!=NULL){
                //cout<<**It1;//print out the Records
                if((*It)->ReportLength()>=GFMin){//only print those records over the minimum gene length
                    vector<string> components=(*It)->GetBasicInfo();
                    ChkOut<<ConstructHeader(Prefix, Delimiter, components)<<"\n";
                    FastaRead::OutputSeq((*It)->GetNuclSeq(InfoPack), ChkOut);
                }
            }
            else ChkOut<<"NULL record: grc_annotate error\n";
        }
        ChkOut.close();
    }
    if(OptFasta.find("A",0)!=string::npos){
        string FileName=FilePrefix+".faa";
        ofstream ChkOut;
        ChkOut.open(FileName.c_str());

        for (list<AARecord*>::iterator It =RecList.begin(); It!=RecList.end(); It++ ){
            if(*It!=NULL){
                //cout<<**It1;//print out the Records
                if((*It)->ReportLength()>=GFMin){//only print those records over the minimum gene length
                    vector<string> components=(*It)->GetBasicInfo();
                    ChkOut<<ConstructHeader(Prefix, Delimiter, components)<<"\n";
                    FastaRead::OutputSeq((*It)->GetTrans(InfoPack), ChkOut);
                }
            }
            else ChkOut<<"NULL record: grc_annotate error\n";
        }
        ChkOut.close();
    }
    if(OptFasta.find("T",0)!=string::npos){
        string FileName=FilePrefix+".pep";
        string TblName=FilePrefix+".tbl";
        string FsaName=FilePrefix+".fsa";
        ofstream FsaOut;
        FsaOut.open(FsaName.c_str());
        InfoPack.WriteGenome(FsaOut);
        FsaOut.close();
        ofstream TblOut;
        TblOut.open(TblName.c_str());
        ofstream ChkOut;
        ChkOut.open(FileName.c_str());
        TblOut<<">feature "<<InfoPack.CurrentGenomeID<<"\n";

        for (list<AARecord*>::iterator It =RecList.begin(); It!=RecList.end(); It++ ){
            if(*It!=NULL){
                //cout<<**It1;//print out the Records
                if((*It)->ReportLength()>=GFMin){//only print those records over the minimum gene length
                    vector<string> components=(*It)->GetBasicInfo();
                    ChkOut<<ConstructHeader(Prefix, Delimiter, components)<<"\n";
                    FastaRead::OutputSeq((*It)->GetTrans(InfoPack), ChkOut);
                    (*It)->WriteTBL(TblOut, InfoPack, Prefix);
                }
            }
            else ChkOut<<"NULL record: grc_annotate error\n";
        }
        ChkOut.close();
        TblOut.close();
    }
    return 0;
}


//Take header components, prefix, and delimiter to construct header line
string ConstructHeader(const string& prefix, const string& delim, vector<string>& Terms){
    string result=">"+prefix;
    for(vector<string>::iterator It=Terms.begin(); It!=Terms.end(); It++){
        if(It==Terms.begin()){
            result+=(*It);
        }
        else{
            result+=delim+(*It);
        }
    }
    return result;
}




/*std::ostream& operator<<(std::ostream& ACOut, const AARecord& AC){
 * ACOut<<"PPCG**"<<"\n";
 * ACOut<<"ID:\t"<<AC.ID<<"\n";
 * ACOut<<"ORF Start:\t"<<AC.Start<<"\n";
 * ACOut<<"ORF Stop:\t"<<AC.Stop<<"\n";
 * ACOut<<"ORF Length(nt):\t"<<AC.QLength<<"\n";
 * if (AC.Reverse){
 * ACOut<<"ORF Strand:\t"<<'-'<<"\n";
 * }
 * else 	ACOut<<"ORF Strand:\t"<<'+'<<"\n";
 *
 * ACOut<<"DB Hit:\t"<<AC.Function<<"\n"; //print hit description or no hit line
 * if(!AC.Blank){//if there is a hit
 * double ALength=AC.ALength;
 * double OLength=AC.QLength;
 * double HLength=AC.HLength;
 * ACOut<<"Hit ID:\t"<<AC.HitID<<"\n";
 * ACOut<<"Hit Org:\t"<<AC.HitOrg<<"\n";
 * ACOut<<"Bit Score:\t"<<AC.Bit<<"\n";
 * ACOut<<"E Score:\t"<<AC.EScore<<"\n";
 * ACOut<<"HSP Length(aa):\t"<<AC.HLength<<"\n";
 * ACOut<<"% Query Aligned:\t"<<(ALength/(OLength/3))*100<<"\n";
 * ACOut<<"% HSP Aligned:\t"<<(ALength/HLength)*100<<"\n";
 * }
 * //ACOut<<"Blank:     "<<AC.Blank<<"\n";
 * //ACOut<<"OLength:   "<<AC.QLength<<"\n";
 * ACOut<<"**PPCG"<<"\n";
 *
 * return ACOut;
 * }*/


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
//"ID\tStart\tStop\tLength(nt)\tStrand\tEDR\tDBID1\tDBID2\tDBID3\tDBOrg\tGOTerms(Conf)\tDescription(Conf)\tBit\tEScore\tHitLength\t%QueryAligned\t%HSPAligned\n"
std::ostream& operator<<(std::ostream& ChkOut, AARecord* AC){
    //ChkOut<<"PPCG**"<<"\t";
    ChkOut<<AC->ID<<"\t";
    ChkOut<<(AC->ReportStart())<<"\t";
    ChkOut<<(AC->ReportStop())<<"\t";
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
        ChkOut<<"No_hits"<<"\t";//DBID1
        ChkOut<<"-\t";//DBID2
        ChkOut<<"-\t";//DBID3
        ChkOut<<"-\t";//DBOrg
        ChkOut<<"-\t";//GOTerms
        ChkOut<<"-\t";//Description
        ChkOut<<"-\t";//Bit
        ChkOut<<"-\t";//EScore
        ChkOut<<"-\t";//HLength
        ChkOut<<"-\t";//% Query aligned
        ChkOut<<"-";//%Subject aligned
    }
    ChkOut<<"\n";
    //ChkOut<<"Blank:     "<<AC.Blank<<"\n";
    //ChkOut<<"OLength:   "<<AC.QLength<<"\n";
    //ChkOut<<"**PPCG"<<"\n";
    return ChkOut;
}

//Display function for the Direct Hash<Compete>
void DisplayKO(ostream& Out, CompeteMap& KOMap, const int& GFMin){
    Out<<"Winner\tLosers...\n";
    for (CompeteMap::iterator KOIt= KOMap.begin(); KOIt != KOMap.end(); KOIt++){//open for loop
        KOIt->second.DisplayKO(Out, GFMin);
    }//close for loop
}//close definition






//this function removes ORFs that do not have hits from the WinnersList and places them in the losers list
//it also creates a KO record that the orf was removed due to entropy
int EntropyFilter(list<AARecord>& RecordList, list<AARecord*>& LoserList, CompeteMap& KOMap, const double& EntCutoff, const double& BitCutoff){
    
    int NumFiltered=0;
    //check each orf for high EDR
    for (list<AARecord>::iterator CheckIt=RecordList.begin(); CheckIt!=RecordList.end(); CheckIt++){
        if(CheckIt->ReportEntropy()>=EntCutoff && (!CheckIt->HasHit() || CheckIt->ReportBitFrac()<BitCutoff)){//if the record has a no/poor hit and its entropy is over cutoff
            NumFiltered++;
            CheckIt->KnockOut();//this orf is out of contention
            LoserList.push_back(&(*CheckIt));//add it to the loser list
            CompeteMap::iterator TempC=KOMap.find("Entropy");//look for the entropy entry in the KOMap
            if(TempC!=KOMap.end()){//if its found add another loser
                (TempC->second).AddLoser(&(*CheckIt));
            }
            else{//if not found add a new entry
                KOMap.insert(CompeteMap::value_type("Entropy", Compete(NULL, (&(*CheckIt)))));
            }
            
        }
    }
    
    return NumFiltered;//return number of orfs removed
}

//this function subjects ORFs that are <300bp to a second round of filtering
int SmallFilter(list<AARecord>& RecordList, list<AARecord*>& LoserList, CompeteMap& KOMap, const double& EntCutoff, const double& BitCutoff){
    
    int NumFiltered=0;
    //check each orf for high EDR
    for (list<AARecord>::iterator CheckIt=RecordList.begin(); CheckIt!=RecordList.end(); CheckIt++){
        if(!CheckIt->Dead() && CheckIt->ReportLength()<300 && CheckIt->ReportEntropy()>=EntCutoff && (!CheckIt->HasHit() || CheckIt->ReportBitFrac()<BitCutoff)){//if the record has a no/poor hit and its entropy is over cutoff
            NumFiltered++;
            CheckIt->KnockOut();//this orf is out of contention
            LoserList.push_back(&(*CheckIt));//add it to the loser list
            CompeteMap::iterator TempC=KOMap.find("Entropy");//look for the entropy entry in the KOMap
            if(TempC!=KOMap.end()){//if its found add another loser
                (TempC->second).AddLoser(&(*CheckIt));
            }
            else{//if not found add a new entry
                KOMap.insert(CompeteMap::value_type("Entropy", Compete(NULL, (&(*CheckIt)))));
            }
            
        }
    }
    
    return NumFiltered;//return number of orfs removed
}



//Train new EDPs for coding and noncoding genes based on the
//training winners and losers
int TrainEDP(list<AARecord*>& WinnerList, list<AARecord*>& LoserList, CalcPack& CP){
    //submit counts to Calc Pack for averaging
    list<AARecord*>::iterator It;
    for(It=WinnerList.begin(); It!=WinnerList.end(); It++){
        (*It)->SubmitCount(CP);
    }
    for(It=LoserList.begin(); It!=LoserList.end(); It++){
        (*It)->SubmitCount(CP);
    }
    CP.CreateOrgEDPs();//averages entropy profiles for coding and non-coding based on losers and winners
    return 0;
}


//RefreshRecords
//refresh the EDR after creating new organism specific EDPs
//and refresh all the queues in all of the records
int RefreshRecords(list<AARecord>& RecordList, CalcPack& CP){
    //for every record
    for(list<AARecord>::iterator It=RecordList.begin(); It!=RecordList.end(); It++){
        It->RefreshEntropy(CP);
        It->RefreshAll();
        It->ResetStatus();
    }
    return 0;
}

//This function parses the incoming ID incase and adjusts coord. based on offset of multiple genomes
int ProcessID(string& ID, long& Start, long& Stop, long& Offset){
    
    unsigned int ChPosition=ID.find("|REPLICON|");//look for '_' in ID indicating that there is a genome ID attached
    unsigned int OPosition=ID.find("|OFFSET|");
    string GenomeID;
    string Junk;
    if(ChPosition!=string::npos && OPosition!=string::npos){
        string TempID=ID;
        //TempID.replace(ChPosition,2," ");//replace '_' with a space
        //ChPosition=ID.find("**");
        //TempID.replace(ChPosition,2," ");//replace '_' with a space
        stringstream ParseSS;
        ParseSS<<TempID;
        getline(ParseSS, ID, '|');
        //ParseSS.ignore();//ignore the next |
        getline(ParseSS, Junk, '|');
        //ParseSS>>ID;//pass orf id through
        getline(ParseSS, GenomeID, '|');
        getline(ParseSS, Junk, '|');
        //ParseSS>>GenomeID; //assign genome id
        ParseSS>>Offset;//assign offset value for contig coordinate conversioni
        Start=Start+Offset;//adjust start/stop positions for multiple contigs. so that concatenated genome sequence can be used
        Stop=Stop+Offset;
        ID+="_"+GenomeID;//reassign ORF ID to be orf_contig
    }
    else{
        GenomeID="NONE";
    }
    return 0;
}

//Parses ECodes from the command line into a set
//and makes sure GRC recognizes it.
int SetEFilter(StringSet& ECodeFilter, const string& ECodeTxt){
    stringstream ReadSS;
    ReadSS<<ECodeTxt;
    string TempS="";
    while(ReadSS>>TempS){
        if(GO::IsECode(TempS)){
            ECodeFilter.insert(TempS);
        }
        else{
            cerr<<"GRC does not recognize Ecode "<<TempS<<". if this is valid"\
                    <<" you may want to add it to GO.cpp::IsECode and recompile.";
        }
    }
    return 0;
}







