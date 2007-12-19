/*GRC.ORFS
This component of the GRC program finds the ORFs in a genome or contig file.*/




#include "FastaRead.h"
#include "OrfFinder.h"

#include <map>
using std::multimap;
using std::map;


string itos(int i){	// convert int to string
	stringstream s;
	s << i;
	return s.str();
}




int main (int argc, char* argv[]) {   //  Main is open
	if(argc<4){//if enough parameters
		cerr<<"Usage: grc_orfs [sequence file]  [minimum gene length] [output file]\n";
		return -1;
		
	}
	char* SeqFile = argv[1]; //get the name of the sequences file
	char* ML = argv[2];
	char* OutFile = argv[3];

	int MinLength = atoi(ML);
	FastaRead Reader;//Object designed for reading in fasta files
	
	ifstream In; //input for the 
	In.open(SeqFile); //open the input file
	Reader.SetInput(&In);//set the input file in the reader
	
	OrfFinder Finder(MinLength);//create OrfFinder with specified minimum length
	CalcPack Calculator;//for calculating sequence related things
	int OrfCounter=0;//for constructing ID's of orfs found
	multimap<int,string> StartStopMap;//stores orf information according to LowBase in increasing order
	map<string,string> IDToSeq;//map that provides access to nucleotide sequence based on ID ASSUMES unique FASTA style ID
	string Header="";
	string Sequence="";
	int Offset=0;//the number of bases processed so far (from previous contigs/genomes)
	int TotalSeq=0;
	//for each nucleotide sequence in the file
	while(Reader.ReadFasta(Header,Sequence)){
		//find +1 orfs
		
		int Start=-1;
		int Stop=-1;
		string SeqTag=Reader.HeaderToID(Header);
		IDToSeq.insert(map<string,string>::value_type(SeqTag,Sequence));
		while(Finder.OrfsForwardFrame(Sequence,Start,Stop,1)){
			//construct string with orf information
			string OrfInfo="";
			OrfInfo+=SeqTag+"\t"+itos(Start)+"\t"+itos(Stop)+"\t"+itos(Offset);//string will be SeqID\tStart\tStop\tOffset
			StartStopMap.insert(multimap<int,string>::value_type(Start+Offset,OrfInfo));
		}
		while(Finder.OrfsForwardFrame(Sequence,Start,Stop,2)){
			//construct string with orf information
			string OrfInfo="";
			OrfInfo+=SeqTag+"\t"+itos(Start)+"\t"+itos(Stop)+"\t"+itos(Offset);//string will be SeqID\tStart\tStop\tOffset
			StartStopMap.insert(multimap<int,string>::value_type(Start+Offset,OrfInfo));
		}
		while(Finder.OrfsForwardFrame(Sequence,Start,Stop,3)){
			//construct string with orf information
			string OrfInfo="";
			OrfInfo+=SeqTag+"\t"+itos(Start)+"\t"+itos(Stop)+"\t"+itos(Offset);//string will be SeqID\tStart\tStop\tOffset
			StartStopMap.insert(multimap<int,string>::value_type(Start+Offset,OrfInfo));
		}
		while(Finder.OrfsReverseFrame(Sequence,Start,Stop,1)){
			//construct string with orf information
			string OrfInfo="";
			OrfInfo+=SeqTag+"\t"+itos(Start)+"\t"+itos(Stop)+"\t"+itos(Offset);//string will be SeqID\tStart\tStop\tOffset
			StartStopMap.insert(multimap<int,string>::value_type(Stop+Offset,OrfInfo));
		}
		while(Finder.OrfsReverseFrame(Sequence,Start,Stop,2)){
			//construct string with orf information
			string OrfInfo="";
			OrfInfo+=SeqTag+"\t"+itos(Start)+"\t"+itos(Stop)+"\t"+itos(Offset);//string will be SeqID\tStart\tStop\tOffset
			StartStopMap.insert(multimap<int,string>::value_type(Stop+Offset,OrfInfo));
		}
		while(Finder.OrfsReverseFrame(Sequence,Start,Stop,3)){
			//construct string with orf information
			string OrfInfo="";
			OrfInfo+=SeqTag+"\t"+itos(Start)+"\t"+itos(Stop)+"\t"+itos(Offset);//string will be SeqID\tStart\tStop\tOffset
			StartStopMap.insert(multimap<int,string>::value_type(Stop+Offset,OrfInfo));
		}
		Offset+=Sequence.size();//set offset so the start stop coordinates can be adjusted for a concatenated genome
	}

	ofstream Out;
	Out.open(OutFile);
	map<string, string>::iterator SeqFinder;
	//for each orf found look up its "genomic" sequence that it is from and output the sequence of the gene
	for(multimap<int,string>::iterator It=StartStopMap.begin(); It!= StartStopMap.end(); It++){
		stringstream InfoSS;
		InfoSS<<It->second;//read in the orf info
		bool Reverse=false;
		int LB=-1;
		int HB=-1;
		int StartOut=-1;
		int StopOut=-1;
		int OffsetOut=-1;
		string LookupID="";
		InfoSS>>LookupID;//get the lookup ID
		InfoSS>>StartOut;
		InfoSS>>StopOut;//get the stop coordinate
		InfoSS>>OffsetOut;
		
		if(StartOut>StopOut){//if the gene is reverse frame
			LB=StopOut;
			HB=StartOut;
			Reverse=true;
		}
		else{
			LB=StartOut;
			HB=StopOut;
			Reverse=false;
		}
		if(IDToSeq.size()>1){//if there is more than one genome
			SeqFinder=IDToSeq.find(LookupID);
		}
		else{
			SeqFinder=IDToSeq.begin();
		}
		if(SeqFinder!=IDToSeq.end()){//if the sequence is found
			TotalSeq+=(HB-LB+1);
			OrfCounter++;
			string GeneSeq=Calculator.GeneSequence(SeqFinder->second,LB,HB,Reverse);
			string OutID="";
			if(IDToSeq.size()>1){//if there is more than one genome adjust IDs to indicate
				Out<<">T"<<OrfCounter<<"|REPLICON|"<<LookupID<<"|OFFSET|"<<OffsetOut<<"\t"<<StartOut<<"\t"<<StopOut<<"\n";
			}
			else {//else no need
				Out<<">T"<<OrfCounter<<"\t"<<StartOut<<"\t"<<StopOut<<"\n";
			}
			Reader.OutputSeq(GeneSeq, Out);//output the sequence
			//cout<<GeneSeq<<"\n";
		}
		else{//Logic error could not find the name of the replicon
			cerr<<"grc_orfs: internal logic error, cannot find ID in sequence lookup.\n";
			return -1;
		}
	}//close for each orf
	
	cout<<"grc_orfs:\n";
	cout<<"orfs found\t"<<StartStopMap.size()<<"\n";
	cout<<"total seq. length\t"<<TotalSeq<<"\n";
	return 0;
}//close main

