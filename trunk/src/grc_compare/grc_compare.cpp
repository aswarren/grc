/*GRC. COMPARE
 *
 * USAGE: grc_compare input1 input2
 *
 * This component of the GRC program compares GRC performance to a given reference file
 * NOTE: There are lots of inefficiencies in this code.  For instance when checking for overlaps
 * the orfs should be orded by low base similar to grc_overlap.  This is not done here.
 * Also there are some polymorphic issues with Record functioning as the class for
 * a Reference record and a GRC prediction.  Also Match pairs up (GRC, Ref) as well as (GRC, GRC) for
 * FN analysis.  Some inheritance would be nice here.
 */




#include "Record.h"
#include "Match.h"
#include "MatchStats.h"
#include "DirectHash.h"
using std::multimap;
using std::vector;

typedef map<string, string> SSMap;

struct OLapID{
    string ID1;
    string ID2;
};
//the variables in the *.Pos and *.Neg input
struct RecVar{
    string ID;
    long Start;
    long Stop;
    double Entropy;
    string DBID;
    string DBOrg;
    string GOString;
    string Hit;
    double Bit;
    string ES;
    long HLength;
    double QPercent;
    double HPercent;
};

/*long stol(const string& target){
    stringstream TempSS;
    TempSS<<target;
    long result;
    TempSS>>result;
    return result;
}*/



int DumpList(list<Record>& InitList);
int DumpList(list<Record*>& InitList);
int Nulify(list<Record*>& InitList, list<Record*>& InitList2);
int Compare(list<Record*>& GList, list<Record*>& RList, list<Record*>& NList, list<Match>& MList, int& NumNO, const int& GFMin, GO* GOAccess);
int PrintCompare(list<Match>& ML);//function for printing out comparison chart of terms and sequence overlap
int LenVSOLap(list<Match>& ML); //function for printing out data for sequence overlap with respect to length of ORF
int PrintSeqMatch(list<Match>& ML, char* RInput);//function for printing distribution chart of combined sequence statistic over reference annotations
int PrintPositive(list<Match>& ML, int& NumNO, double NumRefG, const int& NumNeg, const int& NumSmallDE, const int& NumSmallRef, const int& GFMin, GO* GOAccess, const bool& DumpStats, const string& DumpName);//fuction for printing summary statistics based on FP, TP, FN, TN
int GetKnock(DirectHash<Record*>& Putatives, list<Match>& MList, string KnockFile);//function to get data for analysis
int KnockAnalysis(list<Match>& KML);//function to perform Knockout Analysis
int CountGOStat(Match& TempM, int& GR, int& R, int& G, int& N, int& V);//update the total number of cases for all the grc,ref pairs
int PrintCountDist(int Array[], int size, ostream& Out);//prints out array of values
int PrintFracDist(double Array[], int size, ostream& Out);
SSMap ParseCommandLine(const int& ac, char* const av[]);
int SetColumnKeys(vector<string>& PositionKey, map<string, string>& VariableKey, const char& Separator, const string& HeaderLine);//Figures out the order that columns are in for *.Pos and *.Neg
int ReadVariables(RecVar& TempVar, vector<string>& PositionKey, map<string, string> VariableKey, const char& Separator, const string& CurrentLine);
//run as GRC_compare
int main(int argc, char* argv[]) {   //  Main is open
    
    string RefFile; //get the reference annotation
    string PosFile; //get the name of the parsed test results files
    string NegFile;// the negatives from grc
    string GOFile="none"; //the name of the obo file if there is one
    string KnockFile="none";//specifies the name of the knocklist file that records who knocked out what
    int GFMinLength=300; //the minimum gene finding length
    stringstream Convert;//for converting command line parameters
    bool DumpStats=false;
    vector<string> PositionKey;//for allocating array that stores position of values in input
    map<string, string> VariableKey;
    const char Separator='\t';
    
            
    SSMap Options;
    Options=ParseCommandLine(argc, argv);
    if(Options.size()==0){
        cerr<<"Usage: grc_compare -r [reference annotation] -p [grc results *.Pos] -n [grc results *.Neg] -k [knocklist] -l [min. gene length] OPTIONAL -y [Gene Ontology file] -d (dumps stats to *.Pos.stats)\n";
        return -1;
    }
    
    RefFile = Options.find("-r")->second; //get the reference annotation
    PosFile = Options.find("-p")->second; //get the name of the parsed test results files
    NegFile= Options.find("-n")->second;// the negatives from grc
    KnockFile= Options.find("-k")->second;
    Convert<< Options.find("-l")->second; //the minimum gene finding length
    Convert>>GFMinLength;
    
    if(Options.find("-y")!=Options.end()){//if GO.obo specified
        GOFile=Options.find("-y")->second;
    }
    
    if(Options.find("-d")!=Options.end()){
        DumpStats=true;
    }
    
    
    ifstream In; //input for the setubal results
    ifstream In2; //input for the GRC results
    In.open(RefFile.c_str()); //open the blast input
    
    list<Record> orfRefList; //the Reference List of AA records
    list<Record*> RefList; //the Reference annotations pointer list
    list<Record> PositiveList; //the GRC list for the aa records
    list<Record*> PosList; //the GRC list pointers that will be edited
    list<Match> MList; //list for matches pairs up a putative orf and orf from reference file
    list<Match> FNMList;//the FN knockout list pairs up two orfs involved in knockout
    list<Record> NegativeList; //the negative list
    list<Record*> NegList; //the negative list pointers that will be edited
    long ORFLength=0;
    int NumLessML=0; //the number of Reference orfs less than the minimum finding length
    DirectHash<Record*> IDHash(20000, NULL);//Tracks the predicted orfs based on ID
    GO Ontology;
    GO* GOAccess=NULL;//pointer passed to functions for GO access
    
    
    if(GOFile!="none"){//if there is a GO file specified
        cerr<<"Reading in the ontology file "<<GOFile<<"\n";
        ifstream GOIn;
        GOIn.open(GOFile.c_str());
        Ontology.ReadOBO(&GOIn, true, true, true);
        GOIn.close();
        GOAccess=&Ontology;//set go access pointer
    }//close if there is a obo file
    
    string ID;
    long Stop;
    long Start;
    long Length=0;;
    double Bit=0;
    string ES="nothing";
    long ALength=0;//the length of the query
    string Tag; //for reading in Tags
    string Skip ="nothing"; //string for skipping lines
    double QPercent=0; //the Query percent aligned
    double HPercent=0; //the HSP percent aligned
    string DBID;
    string DBOrg;
    double Entropy=0;
    
    while (In){//read in the Reference input file
        string Hit;
        ID=Hit=ES="nothing";
        Stop=Start=Length=ALength=0;
        Bit=0;
        
        
        
        In >>ID; //read in the ID
        if (!In){//if the input hit the EOF
            break;//break from the loop
        }
        In >>Start; //read in the Start position
        In >>Stop; //read in the Stop position
        
        
        ORFLength=labs(Start-Stop)+1;
        if(ORFLength<GFMinLength){
            NumLessML++;
        }//close if less than min length
        
        getline(In, Hit, '\t');//skip next tab
        getline(In, Hit, '\n');//get the function
        
        orfRefList.push_back(Record(ID, Start, Stop, Hit, 9999, true, GOAccess)); //add a hit record
        RefList.push_back(&(orfRefList.back()));//add pointer to the record
    }// close read input
    In.close();
    
    In2.open(PosFile.c_str()); //open the second input
    
    
    //skip to the first header
    //while (Skip !="PPCG**"){ //read in terms
    //	In2>>Skip;
    //}
    
    
    
    //cout<<"Have Read in Refernce\n";
    getline(In2,Skip);
    SetColumnKeys(PositionKey, VariableKey, Separator, Skip);
    RecVar TempVar;//struct for reading in variables
    string CurrentLine="";
    while (getline(In2, CurrentLine)){//read in the GRC input file
        ReadVariables(TempVar,PositionKey,VariableKey,Separator,CurrentLine);
        
        if(TempVar.DBID=="No_hits"){//no hit
            PositiveList.push_back(Record(TempVar.ID, TempVar.Start, TempVar.Stop, TempVar.DBID, TempVar.Entropy, false)); //add a no-hit record
            PosList.push_back(&(PositiveList.back()));//add pointer to the record
        }//end of no hit
        
        else {
            PositiveList.push_back(Record(TempVar.ID, TempVar.Start, TempVar.Stop, TempVar.GOString+" "+TempVar.Hit, TempVar.Entropy, false, GOAccess, TempVar.Bit, TempVar.ES, TempVar.HLength, TempVar.QPercent, TempVar.HPercent, TempVar.DBID, TempVar.DBOrg)); //add a hit record
            PosList.push_back(&(PositiveList.back()));//add pointer to the record
        }
        //add a pointer to hash table that uses ID to generate key
        unsigned int Position = IDHash.HashingKey(TempVar.ID);
        IDHash.InsertKey(Position, &(PositiveList.back()));
        CurrentLine="";
    }//close while loop for second input
    
    In2.close();
    
    ifstream In3;
    In3.open(NegFile.c_str());
    //cout<<"Have Read In Positives\n";
    getline(In3, Skip);//skip the header line
    
    while (getline(In3, CurrentLine)){//read in the Negatives input file
        ReadVariables(TempVar,PositionKey,VariableKey,Separator,CurrentLine);

        if(TempVar.DBID=="No_hits"){//no hit
            NegativeList.push_back(Record(TempVar.ID, TempVar.Start, TempVar.Stop, TempVar.DBID, TempVar.Entropy, false)); //add a hit record
            NegList.push_back(&(NegativeList.back()));//add pointer to the record
            unsigned int Position = IDHash.HashingKey(TempVar.ID);
            IDHash.InsertKey(Position, &(NegativeList.back()));
        }//end of no hit
        
        else {
            NegativeList.push_back(Record(TempVar.ID, TempVar.Start, TempVar.Stop, TempVar.Hit, TempVar.Entropy, false, NULL, TempVar.Bit, TempVar.ES, TempVar.HLength, TempVar.QPercent, TempVar.HPercent, TempVar.DBID, TempVar.DBOrg)); //add a hit record
            NegList.push_back(&(NegativeList.back()));//add pointer to the record
            unsigned int Position = IDHash.HashingKey(TempVar.ID);
            IDHash.InsertKey(Position, &(NegativeList.back()));
        }
        CurrentLine="";
    }//close while loop for negatives input
    
    //cout<<"have read in negatives\n";
    
    
    //DumpList(orfRefList);//print out the orfs from initlist
//cout<<"have made it past refernce statistics\n";
    //compare the ORFs
    int NumNO=0;
    int NumMatch=Compare(PosList, RefList, NegList, MList, NumNO, GFMinLength, GOAccess); //compare the orfs
    //LenVSOLap(MList); //print out the length versus overlap data
    
    
    
    //Determine the number of overlaps inherent in the Refernce ORF list
    int TotalRefOLaps=0;
    int SameDirRefOLap=0;
    double SumOLap=0;//keeps track of overlap
    double MaxOLap=0;
    double NextMaxOLap=0;
    double MinOLap=1000;
    double CurOLap=0;
    double SumOfSquares=0;
    int NumSmall=0;
    long RefORFLength=0;
    double SumPercentOLap=0;
    double MaxPercent=0;
    string MaxRefOLapID1;
    string MaxRefOLapID2;
    multimap <double, OLapID> Top10OLap;
    int NumDExist=0;//keeps track of the number of reference orfs whose stop and frame dont exist in the pool of orfs generated by long_orfs
    int NumDExistSmall=0;//keeps track of ref. that small and don't exist
    
    for (list<Record>::iterator It11 =orfRefList.begin(); It11!=orfRefList.end(); It11++ ){//outer loop
        RefORFLength=It11->OLength;
        if(RefORFLength<GFMinLength){NumSmall++;}//count the number of reference ORFs less than 300
        if(!It11->RefMatched){
            NumDExist++;
            if (RefORFLength<GFMinLength){NumDExistSmall++;}
        }
        for (list<Record>::iterator It22 =It11; It22!=orfRefList.end(); It22++){ //inner loop
            CurOLap=0;//reset the current Overlap value
            bool SameDirec=false;
            if((It11)!=(It22) && It11->JustOverlap(*It22, SameDirec, CurOLap)){
                SumPercentOLap+=((CurOLap/RefORFLength)+(CurOLap/It22->OLength));
                TotalRefOLaps++;//increase the number of overlapping Reference ORFs
                SumOLap+=CurOLap;//add current overlap to the running total
                SumOfSquares+=(CurOLap*CurOLap);//tally the sum of squares
                OLapID Dummy;
                Dummy.ID1=It11->ID;
                Dummy.ID2=It22->ID;
                Top10OLap.insert(multimap<double, OLapID>::value_type(CurOLap, Dummy));
                if (Top10OLap.size()>10){
                    Top10OLap.erase(Top10OLap.begin());
                }
                if(CurOLap<MinOLap){MinOLap=CurOLap;}//update the minimum
                if(CurOLap>MaxOLap){
                    MaxRefOLapID1=It11->ID;
                    MaxRefOLapID2=It22->ID;
                    NextMaxOLap=MaxOLap;
                    MaxOLap=CurOLap;
                    MaxPercent=CurOLap/RefORFLength;//record percentage as well
                }//update the maximum
                if(SameDirec){
                    SameDirRefOLap++;
                }
            }//close if overlap
        }//close inner loop
    }//close outer loop
    
    
    
    double AvgOLap=SumOLap/TotalRefOLaps;//calculate average
    double AvgPercentOLap=SumPercentOLap/TotalRefOLaps;
    if(MaxOLap>1.1*AvgOLap){
        AvgOLap=(SumOLap-MaxOLap)/(TotalRefOLaps-1);//recalculate average
        AvgPercentOLap=(SumPercentOLap-MaxPercent)/(TotalRefOLaps-1);//recalculate average percent
    }
    
    double StdDev=sqrt((double(TotalRefOLaps) * SumOfSquares - SumOLap * SumOLap)/(double(TotalRefOLaps)*double(TotalRefOLaps-1))); // calculate std dev of overlaps
    
    //Tally results
    //int NumOrfs=0;
    //for (list<Record*>::iterator It1 =FinalList.begin(); It1!=FinalList.end(); It1++ ){
    //	if((**It1).ReportHit()){ NumOrfs++;}
    //}
    cout<<"******************GRC_COMPARISON v0.01******************"<<"\n";
    cout<<"Total # Ref ORFs:\t"<<orfRefList.size()<<"\n";
    cout<<"Total # Ref ORFs<MinLength:\t"<<NumLessML<<"\n";
    cout<<"Total # GRC ORFs:\t"<<PositiveList.size()<<"\n";
    cout<<"Total # Ref/GRC overlaps:\t"<<NumMatch<<"\n\n";
    cout<<"Reference Statistics:\n";
    cout<<"Total # of overlaps within Ref ORFs\t"<<TotalRefOLaps<<"\n";
    cout<<"Total # same direction overlaps in Ref\t"<<SameDirRefOLap<<"\n";
    cout<<"Average length of Reference overlaps\t"<<AvgOLap<<"\n";
    cout<<"Average percent of overlapping ORF length\t"<<AvgPercentOLap<<"\n";
    cout<<"Minimum Reference overlap\t"<<MinOLap<<"\n";
    cout<<"Top 10 Reference overlaps:\n";
    cout<<"Length(nt)\tID1\tID2\n";
    for(multimap<double, OLapID>::iterator PrintIt=Top10OLap.begin(); PrintIt!=Top10OLap.end(); PrintIt++){
        cout<<PrintIt->first<<"\t"<<PrintIt->second.ID1<<"\t"<<PrintIt->second.ID2<<"\n";
    }
    cout<<"Total Reference ORFs <"<<GFMinLength<<"bp\t"<<NumSmall<<"\n";
    cout<<"Total Reference ORFs that DE in orf pool\t"<<NumDExist<<"\n";
    cout<<"Total Reference ORFs that DE <"<<GFMinLength<<"bp\t"<<NumDExistSmall<<"\n";
    //PrintCompare(MList);
    //PrintSeqMatch(MList, RetroIn);//print out distribution chart based on combined statistic
    
    GetKnock(IDHash, FNMList, KnockFile); //read in the knock out list of who did what to whom
    KnockAnalysis(FNMList);
    
    PrintPositive(MList, NumNO, orfRefList.size(), NegList.size(), NumDExistSmall, NumSmall, GFMinLength, GOAccess, DumpStats, PosFile);//Print Summary Statistics
    //print out matches
    cout<<"RESULTS:\n";
    cout<<"Stat\tResult\tMatching_terms\tOverlap\n";
    cout<<"GRC\tID\tStart\tStop\tLength\tStrand\tAnnotation\tBit\te-value\tHSP_length\tDB_ID\tDB_Org\n";
    cout<<"Ref\tID\tStart\tStop\tLength\tStrand\tAnnotation\n";
    if(GOAccess!=NULL){
        cout<<"Confirmed\tRef_GOID\tEvidenceCode\tDepth\tGRC_GOID\tEvidenceCode\tConfidence\tDepth\tDistance\n";
        cout<<"Compatible\tRef_GOID\tEvidenceCode\tDepth\tGRC_GOID\tEvidenceCode\tConfidence\tDepth\tDistance\n";
        cout<<"Incompatibe\tGRC_GOID\tEvidenceCode\tConfidence\n";
    }
    int DCount=0;
    string Delim="-------------------------------------------------------------";
    cout<<Delim<<'\n';
    for (list<Match>::iterator It1 =MList.begin(); It1!=MList.end(); It1++ ){
        result SS = It1->StatScore;
        if(SS==FP || SS==TP || SS==FN){ //only print out the overlap if its a GRC prediction
            /*if(DCount%2==0){
             * Delim="**";
             * }
             * else {
             * Delim="++";
             * }*/
            It1->MatchOut(cout, Delim);
            DCount++;
        }
    }//close for loop
    
    
    
    
    
    //DumpList(orfRefList);
    //PrintCompare(CompareList, WinList);
    Nulify(RefList, PosList);
    
    return 0;
}


//function to print out master list of ORFS
int DumpList(list<Record>& InitList){//open definition
    ofstream Out;
    Out.open("x2blastparse.txt");
    for (list<Record>::iterator It1 =InitList.begin(); It1!=InitList.end(); It1++ ){
        Out<<*It1;//print out the Records
    }
    Out.close();
    return 0;
}

//function to print out data for doing analysis of ORF's that overlap with respect to length
int LenVSOLap(list<Match>& ML){//open definition
    string Percents[10]; Percents[9]="x-90%"; Percents[8]="89-80%"; Percents[7]="79-70%";
    Percents[6]="69-60%"; Percents[5]="59-50%"; Percents[4]="49-40%"; Percents[3]="39-30%";
    Percents[2]="29-20%"; Percents[1]="19-10%"; Percents[0]="09-00%";
    ofstream Out;
    Out.open("OLapvsLengthRef.txt");
    Out<<"TruePositive_ORFS\n"<<"Overlap vs. GRC_ORF_Length(nt)\n";
    Out<<"Overlap= (OLap/GRCLength) * (OLap/RefLength)\n";
    int TPDistORFS[4][10];//2dimensional array to show distribution length versus overlap
    int FPDistORFS[4][10];//2dimensional array to show distribution length versus overlap
    int TNDistORFS[4][10];//2dimensional array to show distribution length versus overlap
    int FNDistORFS[4][10];//2dimensional array to show distribution length versus overlap
    int TPDistORFS2[4][10];//2dimensional array to show distribution length versus overlap
    int FPDistORFS2[4][10];//2dimensional array to show distribution length versus overlap
    int TNDistORFS2[4][10];//2dimensional array to show distribution length versus overlap
    int FNDistORFS2[4][10];//2dimensional array to show distribution length versus overlap
    int FPDistStop[4][10];//for looking at difference between stop sites
    result PorN;
    int x=0;
    int t=0;
    //initalize counter array
    for (x=0; x<4; x++){//outer loop
        for(t=0; t<10; t++){//inner loop
            TPDistORFS[x][t]=0;
            TPDistORFS2[x][t]=0;
            FPDistORFS[x][t]=0;
            FPDistORFS2[x][t]=0;
            FNDistORFS[x][t]=0;
            FNDistORFS2[x][t]=0;
            TNDistORFS[x][t]=0;
            TNDistORFS2[x][t]=0;
            FPDistStop[x][t]=0;
        }//close inner looop
    }//close outer loop
    
    
    int Coord1=0;//for storing the calculated
    int Coord2=0;
    int Coord12=0;
    int CoordStop=0;
    for (list<Match>::iterator It1 =ML.begin(); It1!=ML.end(); It1++ ){
        PorN=It1->StatScore;// get the result of the overlap
        Coord2=It1->CombinedOLap /.1;
        Coord1=It1->GRCRecord->OLength /300;//print out the Records
        Coord12=It1->RefRecord->OLength /300;
        
        if(Coord2>9){Coord2=9;}
        if(Coord1>3){Coord1=3;}
        if(Coord2<0) {cout<<"Error less than zero in OLapvsLength\n";}
        if(Coord1<0) {cout<<"Error less than zero in OLapvsLength\n";}
        if(Coord12>3){Coord12=3;}
        if(Coord12<0) {cout<<"Error less than zero in OLapvsLength\n";}
        
        
        
        switch(PorN){//depending on Result add to statistic map
            case TN:
                TNDistORFS[Coord1][Coord2]++;//increment the count for that bin
                TNDistORFS2[Coord12][Coord2]++;//increment the count for that bin
                break;
            case FN:
                FNDistORFS[Coord1][Coord2]++;//increment the count for that bin
                FNDistORFS2[Coord12][Coord2]++;//increment the count for that bin
                
                break;
            case TP:
                TPDistORFS[Coord1][Coord2]++;//increment the count for that bin
                TPDistORFS2[Coord12][Coord2]++;//increment the count for that bin
                break;
            case FP:
                FPDistORFS[Coord1][Coord2]++;//increment the count for that bin
                FPDistORFS2[Coord12][Coord2]++;//increment the count for that bin
                CoordStop=(((It1->GRCRecord->Stop)-(It1->RefRecord->Stop))/10);
                if(CoordStop<0) {CoordStop=CoordStop*-1;}
                if(CoordStop>3){CoordStop=3;}
                FPDistStop[CoordStop][Coord2]++;//increment the count for that bin
                break;
            default: break;
        }
        
    }
    //print out distribution
    Out<<"Length\tx-900\t900-600\t600-300\t300-0\n";
    Out<<"OvrLap\n";
    for (x=9; x>-1; x--){//outer loop
        Out<<Percents[x]<<"\t";
        for(t=3; t>-1; t--){//inner loop
            Out<<TPDistORFS[t][x]<<"\t";
        }//close inner looop
        Out<<"\n";
    }//close outer loop
    Out<<"\n\n\n";
    
    
    Out<<"FalsePositive_ORFS\n"<<"Overlap vs. ORF_Length(nt)\n";
    Out<<"Overlap= (OLap/GRCLength) * (OLap/RefLength)\n";
    
    //print out distribution
    Out<<"Length\tx-900\t900-600\t600-300\t300-0\n";
    Out<<"OvrLap\n";
    for (x=9; x>-1; x--){//outer loop
        Out<<Percents[x]<<"\t";
        for(t=3; t>-1; t--){//inner loop
            Out<<FPDistORFS[t][x]<<"\t";
        }//close inner looop
        Out<<"\n";
    }//close outer loop
    Out<<"\n\n\n";
    
    Out<<"TrueNegative_ORFS\n"<<"Overlap vs. ORF_Length(nt)\n";
    Out<<"Overlap= (OLap/GRCLength) * (OLap/RefLength)\n";
    //print out distribution
    Out<<"Length\tx-900\t900-600\t600-300\t300-0\n";
    Out<<"OvrLap\n";
    for (x=9; x>-1; x--){//outer loop
        Out<<Percents[x]<<"\t";
        for(t=3; t>-1; t--){//inner loop
            Out<<TNDistORFS[t][x]<<"\t";
        }//close inner looop
        Out<<"\n";
    }//close outer loop
    Out<<"\n\n\n";
    
    Out<<"FalseNegative_ORFS\n"<<"Overlap vs. ORF_Length(nt)\n";
    Out<<"Overlap= (OLap/GRCLength) * (OLap/RefLength)\n";
    //print out distribution
    Out<<"Length\tx-900\t900-600\t600-300\t300-0\n";
    Out<<"OvrLap\n";
    for (x=9; x>-1; x--){//outer loop
        Out<<Percents[x]<<"\t";
        for(t=3; t>-1; t--){//inner loop
            Out<<FNDistORFS[t][x]<<"\t";
        }//close inner looop
        Out<<"\n";
    }//close outer loop
    Out<<"\n\n\n";
    
    
    Out<<"********THIS SECTION BASED ON REFERENCE LENGTH\n";
    Out<<"TruePositive_ORFS\n"<<"Overlap vs. Ref_ORF_Length(nt)\n";
    Out<<"Overlap= (OLap/GRCLength) * (OLap/RefLength)\n";
    //print out distribution
    Out<<"Length\tx-900\t900-600\t600-300\t300-0\n";
    Out<<"OvrLap\n";
    for (x=9; x>-1; x--){//outer loop
        Out<<Percents[x]<<"\t";
        for(t=3; t>-1; t--){//inner loop
            Out<<TPDistORFS2[t][x]<<"\t";
        }//close inner looop
        Out<<"\n";
    }//close outer loop
    Out<<"\n\n\n";
    
    
    Out<<"FalsePositive_ORFS\n"<<"Overlap vs. ORF_Length(nt)\n";
    Out<<"Overlap= (OLap/GRCLength) * (OLap/RefLength)\n";
    
    //print out distribution
    Out<<"Length\tx-900\t900-600\t600-300\t300-0\n";
    Out<<"OvrLap\n";
    for (x=9; x>-1; x--){//outer loop
        Out<<Percents[x]<<"\t";
        for(t=3; t>-1; t--){//inner loop
            Out<<FPDistORFS2[t][x]<<"\t";
        }//close inner looop
        Out<<"\n";
    }//close outer loop
    Out<<"\n\n\n";
    
    Out<<"TrueNegative_ORFS\n"<<"Overlap vs. ORF_Length(nt)\n";
    Out<<"Overlap= (OLap/GRCLength) * (OLap/RefLength)\n";
    //print out distribution
    Out<<"Length\tx-900\t900-600\t600-300\t300-0\n";
    Out<<"OvrLap\n";
    for (x=9; x>-1; x--){//outer loop
        Out<<Percents[x]<<"\t";
        for(t=3; t>-1; t--){//inner loop
            Out<<TNDistORFS2[t][x]<<"\t";
        }//close inner looop
        Out<<"\n";
    }//close outer loop
    Out<<"\n\n\n";
    
    Out<<"FalseNegative_ORFS\n"<<"Overlap vs. ORF_Length(nt)\n";
    Out<<"Overlap= (OLap/GRCLength) * (OLap/RefLength)\n";
    //print out distribution
    Out<<"Length\tx-900\t900-600\t600-300\t300-0\n";
    Out<<"OvrLap\n";
    for (x=9; x>-1; x--){//outer loop
        Out<<Percents[x]<<"\t";
        for(t=3; t>-1; t--){//inner loop
            Out<<FNDistORFS2[t][x]<<"\t";
        }//close inner looop
        Out<<"\n";
    }//close outer loop
    
    Out<<"\n\n\n";
    
    Out<<"FalsePostive_ORFS\n"<<"Overlap vs. StopDistance(nt)\n";
    Out<<"Overlap= (OLap/GRCLength) * (OLap/RefLength)\n";
    //print out distribution
    Out<<"Length\tx-40\t39-20\t19-10\t9-0\n";
    Out<<"OvrLap\n";
    for (x=9; x>-1; x--){//outer loop
        Out<<Percents[x]<<"\t";
        for(t=3; t>-1; t--){//inner loop
            Out<<FPDistStop[t][x]<<"\t";
        }//close inner looop
        Out<<"\n";
    }//close outer loop
    Out.close();
    return 0; //return zero function complete
}

/*Parses the command line and returns a dictionary of strings initializing each
 * option-tag to its corresponding value
 * "Usage: grc_compare -r [reference annotation] -p [grc results *.Pos] -n [grc results *.Neg] -k [knocklist] -l [min. gene length] OPTIONAL -y [Gene Ontology file]
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
        else{
            MarkIt->second=Part;
        }
    }
    if(Result.size()==0){
        return(Result);
    }
    if(Result.find("-r")==Result.end()){
        cerr<<"Missing parameter -r\n";
        Result.clear();
    }
    if(Result.find("-p")==Result.end()){
        cerr<<"Missing parameter -p\n";
        Result.clear();
    }
    if(Result.find("-n")==Result.end()){
        cerr<<"Missing parameter -n\n";
        Result.clear();
    }
    if(Result.find("-k")==Result.end()){
        cerr<<"Missing parameter -k\n";
        Result.clear();
    }
    if(Result.find("-l")==Result.end()){
        cerr<<"Missing parameter -l\n";
        Result.clear();
    }
    return(Result);
}

// NULL OUT lists of Record pointers
int Nulify(list<Record*>& InitList, list<Record*>& InitList2){
    for (list<Record*>::iterator It1 =InitList.begin(); It1!=InitList.end(); It1++ ){
        *It1=NULL;//set to null
    }
    for (list<Record*>::iterator It2 =InitList2.begin(); It2!=InitList2.end(); It2++ ){
        *It2=NULL;//set to null
    }
    return 0;
}


//Compares the GRC orfs versus the Reference ORFs to see if they are overlapping
//The program can be simplified by having this function put Matches for TP, FP, FN, TN on separate lists**********************************
// THIS FUNCTION NO Longer tracks all overlaps with matches It only keeps track of those overlaps that meet the criteria for TP, FP, TN, FN
int Compare(list<Record*>& GList, list<Record*>& RList, list<Record*>& NList, list<Match>& MList, int& NumNO, const int& GFMin, GO* GOAccess){//open defintion
    int count=0;
    bool TruePos=false;//keep track of whether there has been a TruePositive for a GRC ORF
    double MaxPercent=0;
    double HighSS;
    double HighGS;
    double NRPEntropy=0;//Temp variable for outputing average no reference entropy
    double HighOLap;
    int NumSmall;//for keeping track of the number of Ref orfs <300bp that a GRC orf overlaps with
    Record* HighFP=NULL;//for keeping track of the highest overlapping false positive
    for (list<Record*>::iterator It1 =GList.begin(); It1!=GList.end(); It1++ ){
        if(*It1 !=NULL){// if It1 isn't on a null spot
            TruePos=false;
            MaxPercent=0;
            HighFP=NULL;
            HighSS=0;
            HighGS=0;
            HighOLap=0;
            NumSmall=0;
            double OLap=0;
            for (list<Record*>::iterator It2 =RList.begin(); It2!=RList.end(); It2++ ){//inner loop only needs to go through list from position It1
                if(*It2 !=NULL && *It1 !=NULL){//if its ok to compare
                    double GRCOLapPercent;
                    double RefOLapPercent;
                    int WordNum;
                    OLap=0;
                    if((*It1)->Overlap(*(*It2), OLap)){//if the two ORFS overlap
                        count++;
                        
                        if((*It1)->SameFrame(*(*It2))){//if its a true positive
                            if((*It2)->OLength<GFMin) NumSmall++;//increase the number of Ref<Minbp interactions for this prediction
                            (*It1)->Evaluation=TP;
                            MList.push_back(Match(*It1, *It2, OLap, true, NumSmall, GOAccess));//add to match list
                            TruePos=true;
                        }//close true positive
                        else {//IF the prediction is a TP then one is generated if its not then a FP is generated at the end of the loop for the Reference record with which the prediction overlaps the most
                            GRCOLapPercent=OLap/(*It1)->OLength;
                            RefOLapPercent=OLap/(*It2)->OLength;
                            if(MaxPercent<(GRCOLapPercent*RefOLapPercent)){
                                MaxPercent=GRCOLapPercent*RefOLapPercent;//set the new maxpercent
                                HighFP=*It2;//set to the reference ORF
                                HighGS=GRCOLapPercent;
                                HighSS=RefOLapPercent;
                                HighOLap=OLap;
                            }//close maxpercent
                        }//close else false positive
                    }//close overlaping orfs
                }//close ok to compare
            }//close inner loop
            if(!TruePos){//if no true positive has been found
                if(HighFP!=NULL){//if the grc orf overlaps with something
                    (*It1)->Evaluation=FP;
                    if(HighFP->OLength<GFMin) NumSmall++;//increase the number of Ref<Minbp interactions for this prediction
                    MList.push_back(Match(*It1, HighFP, HighOLap, true, NumSmall, NULL));//add to match list
                    HighFP=NULL;
                }
                else {//else the positive doesn't overlap anything
                    (*It1)->Evaluation=NRP;
                    NRPEntropy+=(*It1)->Entropy;
                    NumNO++;
                }
            }
        }//close initial It1 Null check
    }//close outer loop
    
    cout<<"NRP AVERAGE ENTROPY:\t"<<NRPEntropy/double(NumNO)<<"\n";
    bool FalseNeg=false;
    Record* HighTN=NULL;
    //for finding the True and false Negatives
    for (list<Record*>::iterator It3 =NList.begin(); It3!=NList.end(); It3++ ){
        if(*It3 !=NULL){// if It3 isn't on a null spot
            MaxPercent=0;
            FalseNeg=false;
            HighTN=NULL;
            HighGS=0;
            HighSS=0;
            HighOLap=0;
            NumSmall=0;
            double OLap=0;
            for (list<Record*>::iterator It4 =RList.begin(); It4!=RList.end(); It4++ ){//inner loop
                if(*It4 !=NULL && *It3 !=NULL){//if its ok to compare
                    double GRCOLapPercent=0;
                    double RefOLapPercent=0;
                    OLap=0;
                    int WordNum;
                    if((*It3)->Overlap(*(*It4), OLap)){//if the two ORFS overlap
                        
                        if((*It3)->SameFrame(*(*It4)) && !(*It4)->RefMatched){//if its a False Negative
                            if((*It4)->OLength<GFMin) NumSmall++;//increase the number of Ref<300bp interactions for this prediction
                            //count++;
                            (*It3)->Evaluation=FN;//set false negative
                            MList.push_back(Match(*It3, *It4, OLap, false, NumSmall, NULL));//add to match list
                            FalseNeg=true;
                        }//close FalseNegative
                        else{
                            GRCOLapPercent=OLap/(*It3)->OLength;
                            RefOLapPercent=OLap/(*It4)->OLength;
                            if(MaxPercent<(GRCOLapPercent*RefOLapPercent)){
                                MaxPercent=GRCOLapPercent*RefOLapPercent;//set the new maxpercent
                                HighTN=*It4;//set to the reference orf
                                HighSS=RefOLapPercent;
                                HighGS=GRCOLapPercent;
                                HighOLap=OLap;
                            }//close new maxpercent
                        }
                    }//close overlaping orfs
                }//close ok to compare
            }//close inner loop
            if(!FalseNeg){
                if(HighTN!=NULL){
                    (*It3)->Evaluation=TN;
                    MList.push_back(Match(*It3, HighTN, HighOLap, false, NumSmall, NULL));//add to match list
                    HighTN=NULL;
                }
                else{
                    (*It3)->Evaluation=NRN;
                }
            }//close no False negative
        }//close initial It3 Null check
    }//close outer loop
    return count;
}// close function


//function to print out master list of ORFS
int DumpList(list<Record*>& InitList){//open definition
    for (list<Record*>::iterator It1 =InitList.begin(); It1!=InitList.end(); It1++ ){
        if(*It1!=NULL){
            cout<<**It1;//print out the Records
        }
        else cout<<"NULL\n";
    }
    return 0;
}


//function to print the comparisons that have been made
int PrintCompare(list<Match>& ML){//open definition
    
    MatchStats GRCStatMap[10][10];
    MatchStats OtherStatMap[10][10];
    for (list<Match>::iterator It1 =ML.begin(); It1!=ML.end(); It1++ ){
        int GRCI1=floor(It1->GRCOLapPercent/.1);
        int GRCI2=floor(It1->GRCWordMatchPercent/.1);
        int OI1=floor(It1->RefOLapPercent/.1);
        int OI2=floor(It1->RefWordMatchPercent/.1);
        
        
        if(GRCI1>9)GRCI1=9;
        if(GRCI2>9)GRCI2=9;
        if(OI1>9)OI1=9;
        if(OI2>9)OI2=9;
        
        if(GRCI1<0||GRCI2<0||OI1<0||OI2<0){
            
            cout<<"Mapping error check PrintCompare function\n";
            cout<<"GRCI2 "<<It1->GRCWordMatchPercent<<"\n";
            cout<<GRCI1<<"\t"<<GRCI2<<"\t"<<OI1<<"\t"<<OI2<<"\n";
            cout<<*It1;
            return 0;
        }
        
        GRCStatMap[GRCI1][GRCI2].MultiMatch.push_back(&(*It1));
        OtherStatMap[OI1][OI2].MultiMatch.push_back(&(*It1));
    }
    
    
    
    string Range[10]; Range[9]="x-90%"; Range[8]="89-80%"; Range[7]="79-70%";
    Range[6]="69-60%"; Range[5]="59-50%"; Range[4]="49-40%"; Range[3]="39-30%";
    Range[2]="29-20%"; Range[1]="19-10%"; Range[0]="09-00%";
    
    cout<<"Distribution of annotation comparisons (GRC):\n";
    int x=9;
    cout<<setw(7)<<"Terms";
    for(x=9; x>-1; x--){
        cout<<setw(7)<<Range[x];
    }
    cout<<"\n";
    cout<<setw(7)<<"ORFOLap"<<"\n";
    int i=9;
    int t=9;
    for (i=9; i>-1; i--){
        cout<<setw(7)<<Range[i];
        for (t=9; t>-1; t--){
            cout<<setw(7)<<GRCStatMap[i][t].MultiMatch.size();
        }
        cout<<"\n\n\n";
    }
    
    cout<<"Distribution of annotation comparisons (Ref):\n";
    cout<<setw(7)<<"Terms";
    for(x=9; x>-1; x--){
        cout<<setw(7)<<Range[x];
    }
    cout<<"\n";
    cout<<setw(7)<<"ORFOLap"<<"\n";
    for (i=9; i>-1; i--){
        cout<<setw(7)<<Range[i];
        for (t=9; t>-1; t--){
            cout<<setw(7)<<OtherStatMap[i][t].MultiMatch.size();
        }
        cout<<"\n\n\n";
    }
    
    
    return 0;
}


//Print SequenceMatch Score

//function to print the comparisons that have been made
//uses a combined percentage of overlap from the GRC and reference perspective to produce a
//single distribution chart over the number of annotation terms from the reference
int PrintSeqMatch(list<Match>& ML, char* RInput){//open definition
    
    MatchStats GRCStatMap[10][10];
    MatchStats OtherStatMap[10][10];
    for (list<Match>::iterator It1 =ML.begin(); It1!=ML.end(); It1++ ){
        double CombSeqPercent= (It1->GRCOLapPercent) * (It1->RefOLapPercent);
        //int GRCI1=floor(It1->GRCOLapPercent/.1);
        //int GRCI2=floor(It1->GRCWordMatchPercent/.1);
        int OI1=floor(CombSeqPercent/.1);
        int OI2=floor(It1->RefWordMatchPercent/.1);
        
        if(OI1>9)OI1=9;
        if(OI2>9)OI2=9;
        
        if(OI1<0||OI2<0){
            
            cout<<"Mapping error check PrintSeqMatch function\n";
            
            cout<<OI1<<"\t"<<OI2<<"\n";
            cout<<*It1;
            return 0;
        }
        if(OI1==0){cout<<"bleep\n";}
        OtherStatMap[OI1][OI2].MultiMatch.push_back(&(*It1));
    }
    
    
    
    string Range[10]; Range[9]="x-90%"; Range[8]="89-80%"; Range[7]="79-70%";
    Range[6]="69-60%"; Range[5]="59-50%"; Range[4]="49-40%"; Range[3]="39-30%";
    Range[2]="29-20%"; Range[1]="19-10%"; Range[0]="09-00%";
    
    cout<<"Distribution of annotation comparisons CombinedSequenceStatistic (Ref):\n";
    cout<<"ORFOLap = (Overlap/GRCLength)*(Overlap/RefLength)\n";
    int x=9;
    int i=9;
    int t=9;
    
    cout<<setw(7)<<"Terms";
    for(x=9; x>-1; x--){
        cout<<setw(7)<<Range[x];
    }
    cout<<"\n";
    cout<<setw(7)<<"ORFOLap"<<"\n";
    for (i=9; i>-1; i--){
        cout<<setw(7)<<Range[i];
        for (t=9; t>-1; t--){
            cout<<setw(7)<<OtherStatMap[i][t].MultiMatch.size();
        }
        cout<<"\n\n\n";
    }
    
    return 0;
}



//outstream overload for record

std::ostream& operator<<(std::ostream& ACOut, const Record& AC){
    if (AC.Ref){ACOut<<"Reference_PPCG**"<<"\n";}
    else {ACOut<<"GRC_PPCG**"<<"\n";}
    ACOut<<"Result:\t"<<AC.Evaluation<<"\n";
    ACOut<<"ID:\t"<<AC.ID<<"\n";
    ACOut<<"ORF Start:\t"<<AC.Start<<"\n";
    ACOut<<"ORF Stop:\t"<<AC.Stop<<"\n";
    ACOut<<"ORF Length(nt):\t"<<AC.OLength<<"\n";
    if (AC.Reverse){
        ACOut<<"ORF Strand:\t"<<'-'<<"\n";
    }
    else 	ACOut<<"ORF Strand:\t"<<'+'<<"\n";
    
    ACOut<<"Annotation:\t";
    for (list<string>::const_iterator It1 =AC.Description.begin(); It1!=AC.Description.end(); It1++ ){
        ACOut<<*It1<<" "; //output annotation
    }
    ACOut<<"\n";
    //ACOut<<"DB Hit:\t"<<AC.Hit<<"\n"; //print hit description or no hit line
    if(!AC.Blank){//if there is a hit
        //double ALength=AC.ALength;
        double OLength=AC.OLength;
        double HLength=AC.HLength;
        ACOut<<"Bit Score:\t"<<AC.Bit<<"\n";
        ACOut<<"E Score:\t"<<AC.EScore<<"\n";
        ACOut<<"E Value:\t"<<AC.EValue<<"\n";
        ACOut<<"HSP Length(aa):\t"<<AC.HLength<<"\n";
        //ACOut<<"% Query Aligned:\t"<<(ALength/(OLength/3))*100<<"\n";
        //ACOut<<"% HSP Aligned:\t"<<(ALength/HLength)*100<<"\n";
    }
    //ACOut<<"Blank:     "<<AC.Blank<<"\n";
    //ACOut<<"OLength:   "<<AC.OLength<<"\n";
    if (AC.Ref){ACOut<<"**Reference_PPCG"<<"\n";}
    else {ACOut<<"**GRC_PPCG"<<"\n";}
    return ACOut;
}


//outstream overload for the match class
std::ostream& operator<<(std::ostream& Out, const Match& M){
    Out<<"Match**\n";
    switch(M.StatScore){//depending on Result add to statistic map
        case TN:
            Out<<"True Negative\n";
            break;
        case FN:
            Out<<"False Negative\n";
            break;
        case TP:
            Out<<"True Positive\n";
            break;
        case FP:
            Out<<"False Positive\n";
            break;
        default: break;
    }
    Out<<"# Matching annotation terms:\t"<<M.MTerms.size()<<"\n";
    Out<<"MATCHING TERMS: ";
    for (list<string>::const_iterator It1 =M.MTerms.begin(); It1!=M.MTerms.end(); It1++ ){
        Out<<*It1<<"\t";
    }
    Out<<"\n";
    Out<<"% Ref ORF overlap:\t"<<M.RefOLapPercent<<"\n";
    Out<<"% GRC ORF overlap:\t"<<M.GRCOLapPercent<<"\n";
    Out<<"% Ref Annotation match:\t"<<M.RefWordMatchPercent<<"\n";
    Out<<"% GRC Annotation match:\t"<<M.GRCWordMatchPercent<<"\n";
    Out<<*(M.GRCRecord);
    Out<<*(M.RefRecord);
    Out<<"**Match\n";
    return Out;
}

//PrintPositive Function for printing out statistics in relation to the false positive,
//false negative, true positive, and true negative numbers
//GRCRecord is the GRC orf RefRecord is the reference orf
int PrintPositive(list<Match>& ML, int& NumNO, double NumRefG, const int& NumNeg, const int& NumSmallDE, const int& NumSmallRef, const int& GFMin, GO* GOAccess, const bool& DumpStats, const string& DumpName){//open definition
//cout<<"made it into PrintPositive\n";
    int NumFP=0;//counters to keep track of positives and negatives
    int NumTP=0;
    int NumFN=0;
    int NumTN=0;
//counters for tracking interactions with ref orfs<300bp
    int NumSmallTP=0;
    int NumSmallFP=0;
    int NumSmallFN=0;
    double TPPosSumStartDif=0;//keeping track of avg. difference in start position
    double TPNegSumStartDif=0;
    int StartDif=0;
    int TPNumPosStart=0;
    int TPNumNegStart=0;
    double FNPosSumStartDif=0;//keeping track of avg. difference in start position
    double FNNegSumStartDif=0;
    int FNNumPosStart=0;
    int FNNumNegStart=0;
    int FNExactStart=0;
    int TPExactStart=0;
    int FNNumwHit=0;
    int NumTPMatchTerm=0;
    int NumTPMatch50=0;
    int TPGORefGRC=0;//the number of TP with GO for Ref and GRC
    int TPGOGRC=0;//just GRC
    int TPGORef=0;//just Ref
    int TPGONone=0;//neither
    int TPGOVerify=0;//counter that keeps track of the number of TP where at least one GO cat used is the same as GO cat used for Ref
    int TPConfirmed=0;//number of TP with confirmed annotations
    int TPCompatible=0;//number of TP with compatible annotations
    int TPNotCompat=0;//number with not compatible
    int TPTotalConfirm=0;
    int TPTotalCompat=0;
    int TPTotalNotCompat=0;
    int FNGORefGRC=0;//the number of TP with GO for Ref and GRC
    int FNGOGRC=0;//just GRC
    int FNGORef=0;//just Ref
    int FNGONone=0;//neither
    double TPAvgConDist=0;
    double TPAvgConDepth=0;//for GRC annotations what is the average confirmed depth
    double TPAvgComDist=0;
    double TPAvgComDepth=0;//for GRC annotations what is the average confirmed depth
    double FNEntropy=0;
    double FPEntropy=0;
    double TPEntropy=0;
    double TNEntropy=0;
    double HitPEntropy=0;
    double NumHitPos=0;
    int CountConfirm[10];//keeps track of counts in bins of .10 for confidence of confirmed GO terms
    int CountCompat[10];
    int CountIncompat[10];
    double PercentCompat[10]; //The fraction of (compatibleRANGE/allcompatible) all the compatible terms in that bin out of all terms in that bin
    double PercentConfirm[10];
    double PercentIncompat[10];
    double PercentTotal[10];
    ofstream WordOut;//test ofstream
  
    
    for(int t=0; t<10; t++){
        //initialize counter arrays
        CountConfirm[t]=CountCompat[t]=CountIncompat[t]=0;
    }
    
    
    result RefCompare;
    for (list<Match>::iterator It1 =ML.begin(); It1!=ML.end(); It1++ ){
        double CombSeqPercent= (It1->CombinedOLap);
        //int GRCI1=floor(It1->GRCOLapPercent/.1);
        //int GRCI2=floor(It1->GRCWordMatchPercent/.1);
        
        
        
        RefCompare=It1->StatScore;//get the result of the overlap
        
        
        //if(OI1==0){cout<<"Print Positive function index is ZERO!\n";}
        switch(RefCompare){//depending on Result add to statistic map
            case TN:
                NumTN++;
                TNEntropy+=It1->GRCRecord->Entropy;
                break;
            case FN:
                NumFN++;
                StartDif=(It1->GRCRecord->Start-It1->RefRecord->Start);
                FNEntropy+=It1->GRCRecord->Entropy;;
                if(!It1->GRCRecord->Reverse){StartDif=StartDif*-1;}
                if(StartDif>0){
                    FNNumPosStart++;
                    FNPosSumStartDif+=StartDif;
                }
                else if(StartDif<0){
                    FNNumNegStart++;
                    FNNegSumStartDif+=abs(StartDif);
                    //cout<<*(It1->GRCRecord);
                }
                else {FNExactStart++;}
                if(It1->SmallRef>0){//if the reference overlaps an orf less than 300
                    //NumSmallFN=NumSmallFN+It1->SmallRef;
                    NumSmallFN++;
                }
                if(!It1->GRCRecord->Blank){
                    FNNumwHit++;
                }
                
                break;
            case TP:
                if(GOAccess!=NULL){//if GO analysis is turned on
                    CountGOStat(*It1, TPGORefGRC, TPGORef, TPGOGRC, TPGONone, TPGOVerify);//update go statistics
                    if(It1->GOCompatible()>0){//if there are compatible terms
                        TPCompatible++;
                        TPTotalCompat=TPTotalCompat+It1->GOCompatible();
                        It1->SumCompatibleStats(TPAvgComDepth, TPAvgComDist);
                        It1->GOCompatDist(CountCompat);
                    }
                    if(It1->GOConfirmed()>0){//if there are confirmed terms
                        TPConfirmed++;
                        TPTotalConfirm=TPTotalConfirm+It1->GOConfirmed();
                        It1->SumConfirmedStats(TPAvgConDepth, TPAvgConDist);
                        It1->GOConfirmDist(CountConfirm);
                    }
                    if(It1->GONotCompat()>0){
                        TPNotCompat++;
                        TPTotalNotCompat=TPTotalNotCompat+It1->GONotCompat();
                        It1->GOIncompatDist(CountIncompat);
                    }
                }//close if GO Analysis
                
                NumTP++;
                TPEntropy+=It1->GRCRecord->Entropy;
                StartDif=(It1->GRCRecord->Start-It1->RefRecord->Start);
                if(!It1->GRCRecord->Reverse){StartDif=StartDif*-1;}
                if(StartDif>0){
                    TPNumPosStart++;
                    TPPosSumStartDif+=StartDif;
                }
                else if(StartDif<0){
                    TPNumNegStart++;
                    TPNegSumStartDif+=abs(StartDif);
                }
                else{TPExactStart++;}
                
                if(It1->SmallRef>0){//if the reference overlaps an orf less than 300
                    
                    NumSmallTP++;
                }
                if(It1->MTerms.size()>0){
                    NumTPMatchTerm++;
                }
                else {//if TP but no terms match list respective annotations in TPnoterm.txt
                    //It1->ListTerms(WordOut);
                }
                if(It1->GRCWordMatchPercent>.50){
                    NumTPMatch50++;
                }
                if(!It1->GRCRecord->Blank){
                    HitPEntropy+=It1->GRCRecord->Entropy;
                    NumHitPos++;
                }
                break;
            case FP:
                NumFP++;
                FPEntropy+=It1->GRCRecord->Entropy;
                if(It1->SmallRef>0){//if the reference overlaps an orf less than 300
                    //NumSmallFP=NumSmallFP+It1->SmallRef;
                    NumSmallFP++;
                }
                if(!It1->GRCRecord->Blank){
                    HitPEntropy+=It1->GRCRecord->Entropy;
                    NumHitPos++;
                }
                break;
            default: break;
        }
    }//close Match for loop
    
    int AdjustedTP=NumTP-NumSmallTP;//don't count any TP that have a ref. under the minimum gene length
    int AdjustedFN=NumFN-NumSmallFN;
    double AdjustedRef=NumRefG-NumSmallRef;
    //WordOut.close();
    //output summary statistics
    //cout<<"Made it to summary statisitcs\n";
    cout<<"\nSummary Statistics:\n"<<"TP:\t"<<AdjustedTP<<"\n"<<"FP:\t"<<NumFP<<"\t\tNRP:\t"<<NumNO<<"\n"<<"TN:\t"<<NumNeg-NumFN<<"\n";
    cout<<"FN:\t"<<AdjustedFN<<"\tDE:\t"<<int(AdjustedRef-(AdjustedFN+AdjustedTP))<<"\n\n"<<"\nPrecision:\t"<<double(AdjustedTP)/double(AdjustedTP+NumFP)<<"\nSensitivity:\t"<<double(AdjustedTP)/AdjustedRef<<"\n\n";
    cout<<"FN w/ hits:\t"<<FNNumwHit<<"\n\n";
    cout<<"Statistics as to whether a prediction (TP, FP, TN)) overlaps with a Reference ORF that is <"<<GFMin<<"bp\n";
    cout<<"TP w/ ref. <"<<GFMin<<":\t"<<NumSmallTP<<"\n";
    cout<<"FP w/ ref. <"<<GFMin<<":\t"<<NumSmallFP<<"\n";
    cout<<"FN w/ ref. <"<<GFMin<<":\t"<<NumSmallFN<<"\n\n";
    
    cout<<"Start Site statistics:\n";
    cout<<"True Positives:\n";
    cout<<"NumLong\t"<<TPNumPosStart<<"\n";
    if(TPNumPosStart>0){
        cout<<"AvgLong(nt)\t"<<TPPosSumStartDif/TPNumPosStart<<"\n";
    }
    else {
        cout<<"AvgLong(nt)\t"<<"NA\n";
    }
    cout<<"NumShort\t"<<TPNumNegStart<<"\n";
    if(TPNumNegStart>0){
        cout<<"AvgShort(nt)\t"<<TPNegSumStartDif/TPNumNegStart<<"\n";
    }
    else {
        cout<<"AvgShort(nt)\t"<<"NA\n";
    }
    
    cout<<"Exactly Right\t"<<TPExactStart<<"\n\n";
    cout<<"False Negatives:"<<"\n";
    cout<<"NumLong\t"<<FNNumPosStart<<"\n";
    if(FNNumPosStart>0){
        cout<<"AvgLong(nt)\t"<<FNPosSumStartDif/FNNumPosStart<<"\n";
    }
    else {
        cout<<"AvgLong(nt)\t"<<"NA\n";
    }
    cout<<"NumShort\t"<<FNNumNegStart<<"\n";
    if(FNNumNegStart>0){
        cout<<"AvgShort(nt)\t"<<FNNegSumStartDif/FNNumNegStart<<"\n";
    }
    else {
        cout<<"AvgShort(nt)\t"<<"NA\n";
    }
    
    cout<<"ExactlyRight\t"<<FNExactStart<<"\n\n\n";
    
    //output information with no labels
    if(DumpStats){//print Name, MinGeneLength, TP, FP, NRP, TN, FN, DE, Precision, Sensitvity, Exact Start Correct, Exact Start Fraction
        string DumpFile=DumpName+".stats";
        WordOut.open(DumpFile.c_str());//the ofstream operator for records with terms not equivalent
        if(GOAccess==NULL){
            WordOut<<DumpName<<"\t"<<GFMin<<"\t"<<AdjustedTP<<"\t"<<NumFP<<"\t"<<NumNO<<"\t"<<NumNeg-NumFN<<"\t"\
            <<AdjustedFN<<"\t"<<int(AdjustedRef-(AdjustedFN+AdjustedTP))<<"\t"<<double(AdjustedTP)/double(AdjustedTP+NumFP)<<"\t"<<double(AdjustedTP)/AdjustedRef<<"\t"<<TPExactStart<<"\t"<<double(TPExactStart)/double(NumTP)<<"\n";
        }
    }
    double AvgHitEntropy=HitPEntropy/NumHitPos;
    double Variance=0;
    for (list<Match>::iterator EntIt =ML.begin(); EntIt!=ML.end(); EntIt++ ){
        if((EntIt->StatScore==TP || EntIt->StatScore==FP) && !EntIt->GRCRecord->Blank){
            double Diff=EntIt->GRCRecord->Entropy-AvgHitEntropy;
            Variance+=(Diff*Diff);//square the difference of average and actual and add to total variance
        }
    }
    double EntropyDev=sqrt((Variance/(NumHitPos-1)));//calculate std deviation
    
    cout<<"Entropy Statistics:\n";
    cout<<"Avg. Entropy\n";
    //cout<<"AllPositives\t"<<(TPEntropy+FPEntropy+NRPEntropy)/double((NumTP+NumFP+NumNO))<<"\n";
    cout<<"RefPositives\t"<<(TPEntropy+FPEntropy)/double((NumTP+NumFP))<<"\n";
    cout<<"HitPositives\t"<<AvgHitEntropy<<"\n";
    cout<<"StdDevHitPos:\t"<<EntropyDev<<"\n";
    cout<<"TP\t"<<TPEntropy/double(NumTP)<<"\n";
    cout<<"FP\t"<<FPEntropy/double(NumFP)<<"\n";
    cout<<"Negatives\t"<<(TNEntropy+FNEntropy)/double((NumFN+NumTN))<<"\n";
    cout<<"TN\t"<<TNEntropy/double(NumTN)<<"\n";
    cout<<"FN\t"<<FNEntropy/double(NumFN)<<"\n\n\n";
    
    
    cout<<"Annotation Statistics:\n";
    cout<<"TP with a matching term\t"<<NumTPMatchTerm<<"\n";
    //cout<<"TP >.50 matching terms\t"<<NumTPMatch50<<"\n\n\n";
    if(GOAccess!=NULL){//if there is GO analysis
        cout<<"TP Gene Ontology Stats:\n";
        cout<<"Number of TP w/ both GRC and Ref. GO Terms\t"<<TPGORefGRC<<'\n';
        cout<<"Number of TP w/ same category GO Terms(verifiable)\t"<<TPGOVerify<<'\n';
        cout<<"Number of TP w/ just GRC GO Terms\t"<<TPGOGRC<<'\n';
        cout<<"Number of TP w/ just Ref. GO Terms\t"<<TPGORef<<'\n';
        cout<<"Number of TP w/ no GO Terms\t"<<TPGONone<<'\n'<<'\n';
        cout<<"Verifiable TP GO Stats:\n";
        cout<<"Number of TP w/ confirmed annotations\t"<<TPConfirmed<<'\n';
        cout<<"Number of TP w/ compatible annotations\t"<<TPCompatible<<'\n';
        cout<<"Number of TP w/ incompatible annotations\t"<<TPNotCompat<<'\n';
        cout<<"Total number confirmed annotations\t"<<TPTotalConfirm<<'\n';
        cout<<"Total number compatible annotations\t"<<TPTotalCompat<<'\n';
        cout<<"Total number incompatible annotations\t"<<TPTotalNotCompat<<'\n';
        TPAvgConDist=TPAvgConDist/TPTotalConfirm;
        TPAvgConDepth=TPAvgConDepth/TPTotalConfirm;
        TPAvgComDist=TPAvgComDist/TPTotalCompat;
        TPAvgComDepth=TPAvgComDepth/TPTotalCompat;
        cout<<"Average distance for confirmed annotations\t"<<TPAvgConDist<<'\n';
        cout<<"Average depth for confirmed annotations\t"<<TPAvgConDepth<<'\n';
        cout<<"Average distance for compatible annotations\t"<<TPAvgComDist<<'\n';
        cout<<"Average depth for compatible annotations\t"<<TPAvgComDepth<<'\n';
        
        cout<<"Distribution of confidence values:\n";
        string Percents[10]; Percents[9]="1.0-.90"; Percents[8]=".89-.80"; Percents[7]=".79-.70";
        Percents[6]=".69-.60"; Percents[5]=".59-.50"; Percents[4]=".49-.40"; Percents[3]=".39-.30";
        Percents[2]=".29-.20"; Percents[1]=".19-.10"; Percents[0]="0.9-0";
        cout<<"\t";//spacer for row label
        //print labels and calc. fraction distribution
        for(int x=0; x<10; x++){
            cout<<Percents[x]<<"\t";//print labels
            if(TPTotalCompat==0 || CountCompat[x]==0){
                PercentCompat[x]=0;
            }
            else{
                PercentCompat[x]=(double(CountCompat[x])/double(TPTotalCompat));
            }
            if(TPTotalConfirm==0 || CountConfirm[x]==0){
                PercentConfirm[x]=0;
            }
            else{
                PercentConfirm[x]=(double(CountConfirm[x])/double(TPTotalConfirm));
            }
            
            if(TPTotalNotCompat==0 || CountIncompat[x]==0){
                PercentIncompat[x]=0;
            }
            else{
                PercentIncompat[x]=(double(CountIncompat[x])/double(TPTotalNotCompat));
            }
            PercentTotal[x]=double(CountCompat[x]+CountIncompat[x]+CountConfirm[x])/double(TPTotalConfirm+TPTotalCompat+TPTotalNotCompat);
        }
        //print count distributions and fractions
        cout<<"\nCountConfirmed\t";
        PrintCountDist(CountConfirm, 10, cout);
        cout<<"\nFracConfirmed\t";
        PrintFracDist(PercentConfirm, 10, cout);
        cout<<"\nCountCompatible\t";
        PrintCountDist(CountCompat, 10, cout);
        cout<<"\nFracCompatible\t";
        PrintFracDist(PercentCompat, 10, cout);
        cout<<"\nCountIncompatible\t";
        PrintCountDist(CountIncompat, 10, cout);
        cout<<"\nFracIncompatible\t";
        PrintFracDist(PercentIncompat, 10, cout);
        cout<<"\nFracTotal\t";
        PrintFracDist(PercentTotal, 10, cout);
        if(DumpStats){
            WordOut<<DumpName<<"\t"<<GFMin<<"\t"<<NumTP<<"\t"<<TPGOVerify<<"\t"<<TPGOGRC<<"\t"<<TPGORef<<"\t"<<TPGONone<<"\t"<<TPConfirmed<<"\t"<<\
                    TPCompatible<<"\t"<<TPNotCompat<<"\t"<<TPTotalConfirm<<"\t"<<TPTotalCompat<<"\t"<<TPTotalNotCompat<<"\t"<<\
                    TPAvgConDist<<"\t"<<TPAvgConDepth<<"\t"<<TPAvgComDist<<"\t"<<TPAvgComDepth<<"\n";
            /*PrintCountDist(CountConfirm, 10, WordOut);
            WordOut<<"\n";
            PrintFracDist(PercentConfirm, 10, WordOut);
            WordOut<<"\n";
            PrintCountDist(CountCompat, 10, WordOut);
            WordOut<<"\n";
            PrintFracDist(PercentCompat, 10, WordOut);
            WordOut<<"\n";
            PrintCountDist(CountIncompat, 10, WordOut);
            WordOut<<"\n";
            PrintFracDist(PercentIncompat, 10, WordOut);
            WordOut<<"\n";
            PrintFracDist(PercentTotal, 10, WordOut);*/
        }
    }//close if GO
    cout<<"\n\n";
    
    
    return 0;
}//close definition

int PrintCountDist(int Array[], int size, ostream& Out){
    for(int t=0; t<size; t++){
        Out<<Array[t]<<"\t";
    }
    return 0;
}

int PrintFracDist(double Array[], int size, ostream &Out){
    for(int t=0; t<size; t++){
        Out<<Array[t]<<"\t";
    }
    return 0;
}

//Knock Analysis function reads in the KnockList from GRC_Overlap that
//specifies who knocked out who and does analysis find
//the characteristics of what caused the elimination
int GetKnock(DirectHash<Record*>& Putatives, list<Match>& MList, string KnockFile){//open def.
    ifstream GetK;
    string Tyson="nothing";//the ID of the one doing the knocking out
    string Target="nothing";//the ID of the thing knocked out
    GetK.open(KnockFile.c_str());
    GetK>>Tyson;//read in the delimiter
    unsigned int Ty=0;
    unsigned int Tg=0;
    Record** TarPP=NULL;
    Record* TarP=NULL;
    Record* TyP=NULL;
    Record** TyPP=NULL;
    double OLapLen=0;//the length of an overlap
    bool SameDirec=false;//bool value if the two orf's that overlap do so in the same direction
    string Header;
    getline(GetK, Header, '\n');
    while(GetK){//get the knocklist
        string Line;
        getline(GetK, Line);//get the line
        stringstream ss(Line);
        ss>>Tyson;//set KOer
        Ty=Putatives.HashingKey(Tyson);
        if(Tyson!="Entropy"){//if not entropy knockout
            
            while(ss>>Target){//read in the rest of the line until delimeter
                Tg=Putatives.HashingKey(Target);//hash target
                TarPP=(Putatives.FindKey(Tg));//find the putative in the hash
                if(TarPP!=NULL){//if it exists
                    TarP=*TarPP;
                    if((TarP->Evaluation)==FN){//a putative got knocked out that was not supposed to
                        OLapLen=0;
                        TyPP=(Putatives.FindKey(Ty));//get pointer to knock out'er
                        if(TyPP!=NULL){//if the two orfs in the knockout relationship are found
                            TyP=*TyPP;
                            TarP->JustOverlap(*TyP, SameDirec, OLapLen);
                            MList.push_back(Match(TyP, TarP, OLapLen, false));//add to match list
                        }
                        else{
                            //Eliminated by entropy
                            cerr<<Ty<<" does not exist in KnockOut hashTable\n";
                        }
                    }//close putative knocked out
                }//close if TarP exists
                else{cerr<<Tg<<" does not exist in KnockOut hashTable\n";}
            }//close read in rest of line
        }//close if not entropy
        else{//knocked out by entropy
            while(ss>>Target){//read in the rest of the line until delimeter
                Tg=Putatives.HashingKey(Target);//hash target
                TarPP=(Putatives.FindKey(Tg));//find the putative in the hash
                if(TarPP!=NULL){//if it exists
                    TarP=*TarPP;
                    if((TarP->Evaluation)==FN){//a putative got knocked out that was not supposed to
                        OLapLen=0;
                        MList.push_back(Match(NULL, TarP, OLapLen, false));//add to match list
                    }//close putative knocked out
                }//close if TarP exists
                else{cerr<<Tg<<" does not exist in hashTable\n";}
            }//close read in rest of line
        }//close knocked out by entropy
    }//close get the knocklist
    
    return 0;
}//close defintion


//Perform knock analysis
int KnockAnalysis(list<Match>& KML){//
    int TotalOLaps=0;
    int SameDirOLap=0;
    double SumOLap=0;//keeps track of overlap
    double MaxOLap=0;
    double NextMaxOLap=0;
    double MinOLap=1000;
    double CurOLap=0;
    double SumOfSquares=0;
    int NumSmall=0;
    long ORFLength=0;
    double SumPercentOLap=0;
    double MaxPercent=0;
    double SumEScore=0;
    double SumBitScore=0;
    int CountTP=0;
    int CountFP=0;
    int CountTN=0;
    int CountFN=0;
    int CountEN=0;
    int FNNoHit=0;
    int CountNRN=0;
    int CountNRP=0;
    double MLSize=KML.size();
    double NumOL=KML.size();
    ofstream Out;
    Out.open("FNAnalysis.txt");
    for (list<Match>::iterator It11 =KML.begin(); It11!=KML.end(); It11++ ){//outer loop
        if(It11->OverLen==0){
            NumOL--;
        }
        SumOLap+=(It11->OverLen);
        SumPercentOLap+=(It11->RefOLapPercent);//with respect to the FN
        double TempBit=(It11->RefRecord->Bit);
        SumBitScore+=TempBit;
        if (TempBit!=0){
            SumEScore+=(It11->RefRecord->EValue);
        }
        else{	FNNoHit++;}
        if(It11->GRCRecord==NULL){//then it was eliminated by entropy
            CountEN++;//increase entropy count
        }
        else{
            switch((It11->GRCRecord->Evaluation)){
                case TP:
                    CountTP++;
                    break;
                case FP:
                    CountFP++;
                    break;
                case TN:
                    CountTN++;
                    break;
                case FN:
                    CountFN++;
                    break;
                case NRP:
                    CountNRP++;
                    break;
                case NRN:
                    CountNRN++;
                    break;
                default: break;
            }
        }
    }//close for loop
    double AvgEScore=(SumEScore/(MLSize-FNNoHit));
    double AvgBitScore=(SumBitScore/(MLSize-FNNoHit));
    double AvgOLap=(SumOLap/NumOL);
    double AvgPerOLap=(SumPercentOLap/NumOL);
    Out<<"FN analysis\n";
    Out<<"Number of FN:"<<MLSize<<"\n";
    Out<<"FN with Hits:"<<(MLSize-FNNoHit)<<"\n";
    Out<<"Avg. EScore for FN:"<<AvgEScore<<"\n";
    Out<<"Avg. BitScore for FN:"<<AvgBitScore<<"\n";
    Out<<"Avg. Overlap (nt):"<<AvgOLap<<"\n";
    Out<<"Avg. Percent Overlap(FN-nt):"<<AvgPerOLap<<"\n\n";
    Out<<"Total FN Eliminated by TP:"<<CountTP<<"\n";
    Out<<"Total FN Eliminated by FP:"<<CountFP<<"\n";
    Out<<"Total FN Eliminated by NRP:"<<CountNRP<<"\n";
    Out<<"Total FN Eliminated by TN:"<<CountTN<<"\n";
    Out<<"Total FN Eliminated by FN:"<<CountFN<<"\n";
    Out<<"Total FN Eliminated by NRN:"<<CountNRN<<"\n";
    Out<<"Total FN Eliminated by Entropy:"<<CountEN<<"\n\n";
    //Out.flush();
    string Delim="-------------------------------------------------------------";
    //print out knock analysis
    Out<<"Stat\tMatchTerms\tOverlap(nt)\tFN%OLap\n";
    Out<<"Result\tGRC\tID\tStart\tStop\tLength\tStrand\tAnnotation\tBit\te-value\tHSP_length\tDB_ID\tDB_Org\n";
    Out<<"Result\tGRC\tID\tStart\tStop\tLength\tStrand\tAnnotation\tBit\te-value\tHSP_length\tDB_ID\tDB_Org\n";
    Out<<"Result\tRef\tID\tStart\tStop\tLength\tStrand\tAnnotation\n";
    Out<<Delim<<'\n';
    for (list<Match>::iterator It1 =KML.begin(); It1!=KML.end(); It1++ ){
        It1->FNOut(Out, Delim);
    }
    Out.close();
    return 0;
}//close definition

//This function updates a running count of how many Matches have
//GO terms for both the GRC and Reference Records
int CountGOStat(Match& TempM, int& GR, int& R, int& G, int& N, int& V){//open definition
    switch(TempM.GOStat){
        case NoGO:
            N++;
            break;
        case Ref:
            R++;
            break;
        case GRC:
            G++;
            break;
        case RefGRC:
            GR++;
            if(TempM.SameGO){
                V++;
            }
            break;
        default: break;
    }//close switch
    return 0;
}//close def.

//This function provides a mapping of the column names to variable names so that
//if the order of the input columns changes it will handle it automatically
//This is the function to change if the column header names change.
int SetColumnKeys(vector<string> &PositionKey, map<string, string>& VariableKey, const char& Separator,const string& HeaderLine){
    stringstream HeaderSS;
    HeaderSS<<HeaderLine;
    string Temp;
    while(getline(HeaderSS, Temp, Separator)){
        PositionKey.push_back(Temp);
    }
    if(PositionKey.size()<=0){
        cerr<<"Error parsing command line for grc_compare";
        exit(-1);
    }

    //ID\tStart\tStop\tLength(nt)\tStrand\tEDR\tDBID1\tDBID2\tDBID3\tDBOrg\tGOTerms(Conf)\tDescription(Conf)\tBit\tEScore\tHitLength\t%QueryAligned\t%HSPAligned\n";
    //set mapping from column names to variable names format map<ColName, VarName>
    VariableKey.insert(map<string,string>::value_type("ID", "ID"));
    VariableKey.insert(map<string,string>::value_type("Start", "Start"));
    VariableKey.insert(map<string,string>::value_type("Stop", "Stop"));
    VariableKey.insert(map<string,string>::value_type("Length(nt)", "Nothing"));
    VariableKey.insert(map<string,string>::value_type("Strand", "Nothing"));
    VariableKey.insert(map<string,string>::value_type("EDR", "Entropy"));
    VariableKey.insert(map<string,string>::value_type("DBID1", "DBID"));
    VariableKey.insert(map<string,string>::value_type("DBID2", "Nothing"));
    VariableKey.insert(map<string,string>::value_type("DBID3", "Nothing"));
    VariableKey.insert(map<string,string>::value_type("DBOrg", "DBOrg"));
    VariableKey.insert(map<string,string>::value_type("GOTerms(Conf)", "GOString"));
    VariableKey.insert(map<string,string>::value_type("Description(Conf)", "Hit"));
    VariableKey.insert(map<string,string>::value_type("Bit", "Bit"));
    VariableKey.insert(map<string,string>::value_type("EScore", "ES"));
    VariableKey.insert(map<string,string>::value_type("HitLength", "HLength"));
    VariableKey.insert(map<string,string>::value_type("%QueryAligned", "QPercent"));
    VariableKey.insert(map<string,string>::value_type("%HSPAligned", "HPercent"));
    return 0;
}

//This function takes a line of *.Pos or *.Neg input and parses it assigning it to the appropriate variables
int ReadVariables(RecVar& TempVar, vector<string>& PositionKey, map<string, string> VariableKey, const char& Separator, const string& CurrentLine){
    //set default values
    string blank="-";
    TempVar.ID=blank;
    TempVar.Start=0;
    TempVar.Stop=0;
    TempVar.Entropy=0;
    TempVar.DBID=blank;
    TempVar.DBOrg=blank;
    TempVar.GOString=blank;
    TempVar.Hit=blank;
    TempVar.Bit=0;
    TempVar.ES=blank;
    TempVar.HLength=0;
    TempVar.QPercent=0;
    TempVar.HPercent=0;
    stringstream TempSS;
    string ReadStr="";
    TempSS<<CurrentLine;
    int position=0;
    map<string,string>::iterator FindIt=VariableKey.end();
    
    
    while(getline(TempSS, ReadStr, Separator)){
        if(position >=PositionKey.size()){
            cerr<<"Error parsing *.Pos or *.Neg: number of headers does not match number of columns\n";
            exit(-1);
        }
        FindIt=VariableKey.find(PositionKey.at(position));
        if(FindIt!=VariableKey.end() && ReadStr!=blank){
            if(FindIt->second=="ID"){
                TempVar.ID=ReadStr;
            }
            else if (FindIt->second=="Start"){
                TempVar.Start=stol(ReadStr);
            }
            else if (FindIt->second=="Stop"){
                TempVar.Stop=stol(ReadStr);
            }
            else if (FindIt->second=="Entropy"){
                TempVar.Entropy= stod(ReadStr);
            }
            else if (FindIt->second=="DBID"){
                TempVar.DBID=ReadStr;                
            }
            else if (FindIt->second=="DBOrg"){
                TempVar.DBOrg=ReadStr;
            }
            else if (FindIt->second=="GOString"){
                TempVar.GOString=ReadStr;
            }
            else if (FindIt->second=="Hit"){
                TempVar.Hit=ReadStr;
            }
            else if (FindIt->second=="Bit"){
                TempVar.Bit=stod(ReadStr);
            }
            else if (FindIt->second=="ES"){
                TempVar.ES=ReadStr;
            }
            else if (FindIt->second=="HLength"){
                TempVar.HLength=stol(ReadStr);
            }
            else if (FindIt->second=="QPercent"){
                TempVar.QPercent=stod(ReadStr);
            }
            else if (FindIt->second=="HPercent"){
                TempVar.HPercent=stod(ReadStr);
            }
        }
        else if(FindIt==VariableKey.end()){
            cerr<<"Error parsing *.Pos or *.Neg "<<PositionKey[position]<<" not found\n";
            exit(-1);
        }
        position++;
    }
    
    if(position<PositionKey.size()){
        cerr<<"Not enough columns to match header for ID "<<TempVar.ID<<" Headers:"<<PositionKey.size()<<" Columns:"<<position<<"\n";
        exit(-1);
    }
    return 0;
}







