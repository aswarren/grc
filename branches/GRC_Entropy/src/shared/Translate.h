// Translate.h


//Programmer: Andrew Warren
//email: anwarren@vt.edu
//Date of last modification: 02/xx/07




#ifndef Translate_H
#define Translate_H

#include <map>
#include <set>
#include <iostream>
#include <string>
#include <fstream>
#include <stdio.h>

using std::map;
using std::set;
using std::cout;
using std::cerr;
using std::string;
using std::stringstream;
using std::ostream;
using std::ofstream;
using std::ifstream;

class Translate {//open prototype

private:
	int TransCode;
	int TransCodePosition;
	map <string,char> TTables[23]; //translation tables
	string InputFile;
public:
	//default constructor
	Translate(){
		TransCode=11;
		TransCodePosition=TransCode-1;
		InputFile="none";
	}

	//parameterized constructor
	Translate(int TC, string IF){
		InputFile=IF;
		TransCode=TC;
		TransCodePosition=TransCode-1;
		InitCodes();//read in translation tables
	}

	//copy constructor
	Translate(const Translate &Source){// open defintion
		TransCode=Source.TransCode;
		InputFile=Source.InputFile;
		TransCodePosition=Source.TransCodePosition;
		for(int t=0; t<23; t++){
			TTables[t]=Source.TTables[t];
		}
	}

	//assignment operator
	 Translate& operator =(const Translate &Source){// open defintion
		if (this!= &Source){// open non-self assignment consq.
			TransCode=Source.TransCode;
			InputFile=Source.InputFile;
			TransCodePosition=Source.TransCodePosition;
			for(int t=0; t<23; t++){
				TTables[t]=Source.TTables[t];
			}
		}
		return *this;
	}// close definition


	//Set Input File
	int SetTransFile(const int& TN, const string& IF){
		TransCode=TN;
		TransCodePosition=TransCode-1;
                InputFile=IF;
		return 0;
	}

/* Name:InitCodes
**	@param1: Array of 23 maps that are to be initialized to contain translation information
**	Function: This function initializes all of the genetic codes to be used in the translation of the input file
*/
	int InitCodes(){//open definition
	//Get the Genetic Codes for translation
		ifstream InCode;
		InCode.open(InputFile.c_str());
		int TableNum;// for reading in the table number
		string Delim;
		string OpCase;
		char AA; //character for storing the amino acid
		//InCode>>Delim;

		while(InCode){
			InCode>>Delim;
			if(Delim=="++"){
				InCode>>TableNum;
				TableNum--;
			}
			else{
				if(isupper(Delim[0])){//create opposite case upper/lower in map
					OpCase=tolower(Delim[0]);
					for(int t=1; t<Delim.size(); t++){
						OpCase+=tolower(Delim[t]);
					}
				}
				else {
					OpCase=toupper(Delim[0]);
					for(int x=1; x<Delim.size(); x++){
						OpCase+=toupper(Delim[x]);
					}
				}

				InCode>>AA; //Read in the AA translation
				TTables[TableNum].insert(map<string,char>::value_type(Delim,AA));
				TTables[TableNum].insert(map<string,char>::value_type(OpCase,AA));
			}
		}//close while loop
		
		InCode.close();
		return 0;
	}//close definition
	
        //function to lower the case of all characters in a string
	int LowerTheCase(string & Seq){
		for(int i=0; i<Seq.length(); i++){
			Seq[i]=tolower(Seq[i]);
		}
		return 1;
	}//close definition
        
        //Return the codons that are marked as stops in the translation table
       set <string> GetStops(){
            set<string> TempSet;
            if(TTables[TransCodePosition].size()==0){
                cerr<<"\nNo translation table for determining stop codons\n";
                throw 20;
            }
            else{
                for(map<string,char>::iterator It=TTables[TransCodePosition].begin(); It!=TTables[TransCodePosition].end(); It++){
                    if(It->second=='*'){
                        string Codon=It->first;
                        LowerTheCase(Codon);
                        TempSet.insert(Codon);
                    }
                }
            }
            return TempSet;
        }


//TranslateSeq
//function for translating provided sequence
	string TranslateSeq(const string& Sequence){
		int len=3;
		char AA ='X';
		string Translation="";
		string SubSeq;
		map<string,char>::iterator FindIt;
		for (int s=0; s<Sequence.size(); s=s+3){
			if((s)%3==0){//if its the end of the next codon
				SubSeq=Sequence.substr(s,len);				
				FindIt=TTables[TransCodePosition].find(SubSeq);//look for the codon
				if(FindIt!=TTables[TransCodePosition].end()){//if found
					AA=FindIt->second;
				}
				Translation+=AA;//append
				AA='X';
			}
		}
		//enforce that the first codon in an orf gives Met
		if(Translation.size()>0 && Translation[0]!='M'){
			Translation[0]='M';
		}
		return Translation;
	}

};//close prototype

#endif
