/*GRC.Translate
This component of the GRC program translates the orfs generated from the long orfs component into amino acid sequences based on hard-coded translation maps.
The default is to use the bacteria map but a command line parameter is available to change to Universal/Standard. The information for these tables was obtained from the Taxonomy division of NCBI http://ncbi.nih.gov/Taxonomy/Utils/wprintgc.cgi?mode=t#SG1
*/



#include <map>
#include <iostream>
#include <string>
#include <fstream>
#include <stdio.h>

using std::map;
using std::cout;
using std::cerr;
using std::string;
using std::stringstream;
using std::ostream;
using std::ofstream;
using std::ifstream;


int InitCodes(map <string,char>* CodeArray);

//run as grc_translate extract.out <tablenum> 
int main (int argc, char* argv[]) {   //  Main is open
	char* InFile = argv[1]; //get the name of the input file to translate
	char* TCode = "11"; //the number of the translation to use Default 11
	int TransCode;
	if(argc>2){
		TCode = argv[2];//assign the translation code to use
	}
		

	ifstream In; //input for the 
	In.open(InFile); //open the input file

	TransCode =atoi(TCode); //convert the translation parameter to an int
	TransCode--;//decrease to match indicies

	map <string,char> TTables[23]; //translation tables
	InitCodes(TTables);


	//Output written to cout
	//string coutName=InFile;
	//OutName+=".transl";
	//ofstream Out;
	//cout.open(OutName.c_str());// open the output file

	//Begin Read/Write/Translation Loop
	string Line;
	string Sequence;
	string SubSeq;
	const int LineLength=70;//number of AA per line
	const int len=3; //length of codon
	char AA='*';
	bool EndLine=true;
	map<string,char>::iterator FindIt;//iterator for finding amino acid
	while(In){
		Line="!!!";
		In >>Line;
		if(!In){
			break;
		}
		if(Line[0]=='>'){ //If its a header line
			cout<<Line;
			getline(In,Line,'\n');
			cout << Line << "\n";
		}
		else { //get the Sequence
			Sequence=Line;
			
			for (int s=0; s<Sequence.size(); s=s+3){
				if((s)%3==0){//if its the end of the next codon
					SubSeq=Sequence.substr(s,len);
					if(s==0){//if its the first codon set to M
						AA='M';
					}
					else{
						FindIt=TTables[TransCode].find(SubSeq);//look for the codon
						if(FindIt!=TTables[TransCode].end()){//if found
							AA=FindIt->second;
						}
					}
					if(((s+3)/3)%LineLength ==0){//if there have been the specified number of aa displayed
						cout<<AA<<'\n';
						EndLine=false;
					}
					else{
						cout<<AA;
						EndLine=true;
					}
					AA='*';
				}
			}
			if(EndLine){
				cout<<"\n";//next line
			}
		}//close get the sequence
	}//close while

	//Out.close();
	In.close();

	return 0;

} //close main



/* Name:InitCodes
**	@param1: Array of 23 maps that are to be initialized to contain translation information
**	Function: This function initializes all of the genetic codes to be used in the translation of the input file
*/
	int InitCodes(map <string,char>* CodeArray){//open definition
	//Get the Genetic Codes for translation
		ifstream InCode;
		InCode.open("GCode.txt");
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
				CodeArray[TableNum].insert(map<string,char>::value_type(Delim,AA));
				CodeArray[TableNum].insert(map<string,char>::value_type(OpCase,AA));
			}
		}//close while loop
		
		InCode.close();

		return 0;
	}//close definition
			

	


