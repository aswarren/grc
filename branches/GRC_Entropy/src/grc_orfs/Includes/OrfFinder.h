//OrfFinder.h
//The following is a class for taking a nucleotide sequence and
//finding orfs of a specified minimum length

#ifndef OrfFinder_H
#define OrfFinder_H

		
#include <iostream>
#include <string>
#include <fstream>
#include <map>
#include <sstream>
#include <set>
#include "CalcPack.h"
		
		
using std::cout;
using std::cerr;
using std::string;
using std::ostream;
using std::ofstream;
using std::ifstream;
using std::istream;
using std::map;
using std::stringstream;
using std::set;

class OrfFinder: public CalcPack{//open prototype
	int Position;
	bool MoreStatus;
	int MinLength;
	string* CurrentSeq;//keeps track of the current sequence
	int CurrentFrame;
	
	public:
		
	//Default constructor
	OrfFinder(){
		Position=0;
		MinLength=300;
		MoreStatus=true;
		CurrentSeq=NULL;
		CurrentFrame=0;
	}
	//Parameterized constructor
	OrfFinder(const int& MinL){
		MinLength=MinL;
		Position=0;
		MoreStatus=true;
		CurrentSeq=NULL;
		CurrentFrame=0;
	}
	//Copy Constructor
	OrfFinder(const OrfFinder &Source){// open defintion
		Position=Source.Position;
		MinLength=Source.MinLength;
		MoreStatus=Source.MoreStatus;
		CurrentSeq=Source.CurrentSeq;
		CurrentFrame=Source.CurrentFrame;
	}
	
	//Assignment operator
	OrfFinder& operator =(const OrfFinder &Source){// open defintion
		if(this!=&Source){
			Position=Source.Position;
			MinLength=Source.MinLength;
			MoreStatus=Source.MoreStatus;
			CurrentSeq=Source.CurrentSeq;
			CurrentFrame=Source.CurrentFrame;
		}
		return *this;
	}

	
	//This function takes three parameters
	//1. The sequence in which to find ORFS
	//2. The int that the Start Position will be assigned to
	//3. The int that the Stop Position will be assigned to
	//The function expects to encounter >ID\n  Sequence
	//Returns true if still values returned are good.
	//Returns false if an error is encountered or end of the sequence is reached
	//expects frame 1, 2, or 3
	bool OrfsForwardFrame(string& Seq, int& Start, int& Stop, const int& Frame){
		Start=Stop=-1;
		MoreStatus=true;
		if(&Seq != CurrentSeq){
			Position=CurrentFrame;
			CurrentSeq=&Seq;//remember what sequence currently operating on
		}
		if(Frame-1!=CurrentFrame){
			Position=CurrentFrame=Frame-1;
		}
		//
		for(int s=Position; s<=(Seq.size()-3); s+=3){
			if(Start==-1 && ForwardStart(Seq.substr(s,3))){//if its a start codon and not yet assigned
				Start=s+1;//set start coordinate
			}
			else if(ForwardStop(Seq.substr(s,3))){//if its a stop codon
				if(Start!=-1){// start site already found
					Stop=s+3;//stop includes the stop codon
					if((Stop-Start)+1-3>=MinLength){
						Position=s+3;//advance to the next codon
						return true;
					}
					else{//else reset
						Stop=-1;
						Start=-1;
					}
				}
			}
		}
		
		//if we have walked off the end of the sequence and there is an orf not yet returned
		if(Stop!=-1 && Start!=-1 && (Stop-Start)+1>=MinLength){
			Position=Seq.size();
			return true;
		}
		else{
			CurrentFrame=-1;//made it all the way through this frame reset CurrentFrame value
			return false; //return whether there are more records to be read
		}
	}
	
	//This function takes three parameters
	//1. The sequence in which to find ORFS
	//2. The int that the Start Position will be assigned to
	//3. The int that the Stop Position will be assigned to
	//The function expects to encounter >ID\n  Sequence
	//Returns true if still values returned are good.
	//Returns false if an error is encountered or end of the sequence is reached
	//expects frame 1, 2, or 3
	bool OrfsReverseFrame(string& Seq, int& Start, int& Stop, const int& Frame){
		Start=Stop=-1;//-1 is the unassigned value
		MoreStatus=true;
		if(&Seq != CurrentSeq){
			Position=CurrentFrame;
			CurrentSeq=&Seq;//remember what sequence currently operating on
		}
		if(Frame-1!=CurrentFrame){
			Position=CurrentFrame=Frame-1;
		}
		//
		for(int s=Position; s<=(Seq.size()-3); s+=3){


			if(ReverseStart(Seq.substr(s,3))){//if its a start codon
				Start=s+1+2;//set start coordinate +1 to get to genome coordinate and +2 because first base of the reverse start
			}
			else if(ReverseStop(Seq.substr(s,3))){//if its a stop codon
				if(Start!=-1 && Stop!=-1){// start site and stop already found (at next in frame stop)
					if((Start-Stop)+1-3>=MinLength){//minus 3 since the stop coordinate includes the stop codon
						Position=s;//begin here to pick up this stop
						return true;
					}
					else{//else reset start record stop and keep going
						Stop=s+1;
						Start=-1;
					}
				}
				else if(Start==-1){//if start is unassigned then no start found before this next in frame stop OR this is first stop found
					Stop=s+1;//stop is currently includes the stop codon +1 to genome coord.
				}
				else{//Start!=-1 && Stop==-1 invalid start with this stop
					Stop=s+1;
					Start=-1;
				}
			}
		}
	
		//if we have walked off the end of the sequence and there is an orf not yet returned
		if(Stop!=-1 && Start!=-1 && (Start-Stop)+1>=MinLength){
			Position=Seq.size();//set position to end of sequence
			return true;
		}
		else{
			CurrentFrame=-1;//made it all the way through this frame reset CurrentFrame value
			return false; //return whether there are more records to be read
		}
	}
	
	
};//close prototype

#endif
