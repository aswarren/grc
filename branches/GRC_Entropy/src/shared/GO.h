//GO.h


//Programmer: Andrew Warren
//email: anwarren@vt.edu
//Date of last modification: 08/xx/06      

//Gene Ontology lib for processing GO.obo files and doing some analyses


#ifndef GO_H
#define GO_H

#include "PointHash.h"
#include "GOFunction.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <list>
using std::ios;
using std::cout;
using std::cerr;
using std::string;
using std::ostream;
using std::ofstream;
using std::ifstream;
using std::istream;
using std::stringstream;
using std::list;
using std::map;

typedef map<GOFunction*,string,LessAncestor> ANCESTOR;//for storing ancestors of a given node

typedef map<int,int> ANCESTORID;//for storing the ID's of ancestors


class GO{
	friend class Match;
private:
	list<GOFunction*> Storage; //the place where functions actually exist
	PointHash<GOFunction*> GOPoint; //access to functions by GOID Internal is GOFunction**
	GOFunction* PHead; //pointer to the head of the ontology
	GOFunction* CHead; //pointer to the head of the ontology
	GOFunction* FHead; //pointer to the head of the ontology
	vector<GOFunction*> PHeads; //potential head of bprocess hiearchy
	vector<GOFunction*> FHeads; //potential head of mfunction hiearchy
	vector<GOFunction*> CHeads; //potential head of ccomponent hiearchy
	int BuildOntology();//build hash table for ID's and initialize relationship pointers
	bool MFBool;
	bool BPBool;
	bool CCBool;
	int GatherA(ANCESTOR* AP, GOFunction* Start, int Dist, const int& OrigID, const int& UpdateID);//recursive function for gathering ancestors
	int AllAnCounter;//this counter is used in the all ancestors function to keep the origID variable the  same for a group of terms.
	set<int> NotFound;//keep track of whether there has previously been a message about not finding a GO ID
        set<string> UnknownTag;//keep track of whether the does not recognize a tag
        bool WarnFound;
public:
	
	GO();//default constructor
	GO(const GO& Source);//copy constructor
	GO& operator =(const GO& Source); //assignment operator
	~GO();//default destructor
	int ReadOBO(istream *In, const bool& MFunc, const bool& BProc, const bool& CComp);
	ANCESTOR GetAncestors(const int& GOID);//return the ancestors for a given GO ID
	GOFunction* Find(const int& FindMe);//function for finding goid
	GOFunction* Find(const string& ToFind);//functin for finding goid string
	ANCESTORID GetAncestorIDs(const int& GOID);//return the GOID's for an ancestor
	static string IDToString(const int& TempID);//converts xxxx to GO:000xxxx
	static int StringToID(const string& TempS);//converts GO:000xxxx to xxxx
	static bool StringIsGO(const string& TempS);//returns whether string is a GO Term
	static bool IsECode(const string& TempS);//returns bool string is evidence code
	int GetAllAncestors(const int& GOID, ANCESTOR& Family);
        int ExpandAltID(GOFunction** CF);
        int FindWarningOff();
	

};//close prototype for GO

#endif



