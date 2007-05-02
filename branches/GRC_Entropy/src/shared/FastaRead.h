//FastaRead.h
//The following is a class for reading a generic fasta file
//Expects header lines to begin with > and the all lines in between header lines are
//sequence lines.


#ifndef FastaRead_H
#define FastaRead_H

		
#include <iostream>
#include <string>
#include <fstream>
#include <map>
#include <sstream>

		
using std::cout;
using std::cerr;
using std::string;
using std::ostream;
using std::ofstream;
using std::ifstream;
using std::istream;
using std::map;
using std::stringstream;

class FastaRead{//open prototype
	istream* Input;//pointer to current input source
	string Line;
	bool MoreStatus;
public:
		
	
	FastaRead(){
		Input=NULL;
		Line="NONE";
		MoreStatus=true;
	}
	
	//Copy Constructor
	FastaRead(const FastaRead &Source){// open defintion
		Input=Source.Input;
		Line=Source.Line;
		MoreStatus=Source.MoreStatus;
	}
	
	//Assignment operator
	FastaRead& operator =(const FastaRead &Source){// open defintion
		if(this!=&Source){
			Input=Source.Input;
			Line=Source.Line;
			MoreStatus=Source.MoreStatus;
		}
		return *this;
	}
	
	//this function must be used to set the input file that
	//will be read from
	int SetInput(istream* In){
		Input = In;
		Line="NONE";
		MoreStatus=true;
		return 0;
	}
	
	//This function takes two parameters
	//1. The string that ID will be assigned to
	//2. The string that the sequence will be assigned to
	//The function expects to encounter >ID\n  Sequence
	//Returns true if still values returned are good.
	//Returns false if an error is encountered or end of the file is reached
	bool ReadFasta(string& ID, string& Seq){
		ID ="NONE";
		Seq ="";
		while(MoreStatus && (Line == "NONE" || Line=="")){//if its the first line or a blank line
			MoreStatus=getline(*Input,Line); //get line for ID
		}
		if(MoreStatus && Line[0]=='>'){
			ID = Line; // assign ID value
			while(MoreStatus=getline(*Input,Line)){//while not yet at the end of file, get next line
				if(Line[0]=='>'){
					break; //break from loop if at next ID
				}
				else{
					LowerTheCase(Line);//add line to sequence
					Seq+=Line;
				}
			}
			return true;//these values are good
		}
		else {
			MoreStatus=false;//else there is a parsing error
		}
		return MoreStatus; //return whether there are more records to be read
	}
	
	//This function is for creating an ID from the Header line in a fasta file
	//It returns the ID immediately following the > character
	string HeaderToID(string Header){
		stringstream ParseStream;
		ParseStream<<Header;
		string SeqTag="";//Tag created from genome ID
		ParseStream>>SeqTag;
		if(SeqTag[0]=='>'){
			SeqTag.erase(0,1);//remove the '>' character
		}
		return SeqTag;
	}
	
	//function for outputing sequence in fasta format
	int OutputSeq(const string& Seq, ostream& Out){
		int LineLength=70;
		for (int s=0; s<Seq.size(); s+=LineLength){
			if((Seq.size()-s)>=LineLength){//if its the end of the next codon
				Out<<Seq.substr(s,LineLength)<<"\n";
			}
			else{
				Out<<Seq.substr(s,Seq.size()-s)<<"\n";
			}
		}
		return 0;
	}
	
		//function to lower the case of all characters in a string
	int LowerTheCase(string & Seq){
		for(int i=0; i<Seq.length(); i++){
			Seq[i]=tolower(Seq[i]);
		}
		return 1;
	}//close definition
	
}; //close prototype

#endif
