/*GRC.Translate
This component of the GRC program translates the orfs generated from the long orfs component into amino acid sequences based on hard-coded translation maps.
The default is to use the bacteria map but a command line parameter is available to change to Universal/Standard. The information for these tables was obtained from the Taxonomy division of NCBI http://ncbi.nih.gov/Taxonomy/Utils/wprintgc.cgi?mode=t#SG1
*/



#include <map>
#include <iostream>
#include <string>
#include <fstream>
#include <stdio.h>
#include "FastaRead.h"
#include "Translate.h"

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
	FastaRead Reader;//object for reading in fasta files
	
	if(argc<5){
		cerr<<"Usage: grc_translate [translation tables file] [table number] [nucleotide file (fasta format)] [output file]\n";
		cerr<<"This program uses NCBI translation tables specified in the format seen in GCode.txt\n";
	}
	string InCode= argv[1];//translation tables
	char* TCode = argv[2]; //the number of the translation to use Default 11
	char* InFile = argv[3]; //get the name of the input file to translate
	char* OutFile =argv[4]; //output file

	

	int TransCode =atoi(TCode); //convert the translation parameter to an int
	Translate Translator(TransCode,InCode);//object for translating nucleotide sequences to AA
	ifstream In; //input for the 
	In.open(InFile); //open the input file
	ofstream Out;
	Out.open (OutFile);

	Reader.SetInput(&In);//set the input file in the reader

	


	//Output written to cout
	//string coutName=InFile;
	//OutName+=".transl";
	//ofstream Out;
	//cout.open(OutName.c_str());// open the output file

	//Begin Read/Write/Translation Loop
	string ID;
	string Sequence;
	string Translation;
	string SubSeq;
	const int LineLength=70;//number of AA per line
	const int len=3; //length of codon
	char AA='X';
	bool EndLine=true;
	map<string,char>::iterator FindIt;//iterator for finding amino acid
	while(Reader.ReadFasta(ID,Sequence)){
		/*Line="!!!";
		In >>Line;
		if(!In){
			break;
		}
		if(Line[0]=='>'){ //If its a header line
			cout<<Line;
			getline(In,Line,'\n');
			cout << Line << "\n";
		}
		else { //get the Sequence*/
			Out << ID <<"\n";
			Translation=Translator.TranslateSeq(Sequence);
			Reader.OutputSeq(Translation, Out);
	}//close while

	//Out.close();
	In.close();

	return 0;

} //close main




			

	


