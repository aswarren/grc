// GOMatch.h
//GOMatch stores the confirmation
//of a GRCAnnotation by a Reference annotation
//It pairs the two together and provides relevant distance/depth information


//Programmer: Andrew Warren
//email: anwarren@vt.edu
//Date of last modification: 10/xx/06    

#ifndef GOMatch_H
#define GOMatch_H

#include "GO.h"
#include <set>
#include <string>
using std::string;
using std::set;

string itos(int i)	// convert int to string
	{
		stringstream s;
		s << i;
		return s.str();
	}


class GOMatch {//open prototype
public:
	int RefGOID;
	int GRCGOID;
	set<string> RefECodes;
	set<string> GRCECodes;
	int Distance;//distance from one another
	int GRCDepth;
	int RefDepth;
	GOFunction* RefTerm;
	GOFunction* GRCTerm;


	GOMatch(){//default constructor
		RefGOID=-1;
		GRCGOID=-1;
		Distance=-1;
		GRCDepth=-1;
		RefDepth=-1;
		RefTerm=NULL;
		GRCTerm=NULL;
	}//close def

	GOMatch(const GOMatch &Source){// copy constructor
		RefGOID=Source.RefGOID;
		GRCGOID=Source.GRCGOID;
		RefECodes=Source.RefECodes;
		GRCECodes=Source.GRCECodes;
		Distance=Source.Distance;
		GRCDepth=Source.GRCDepth;
		RefDepth=Source.RefDepth;
		RefTerm=Source.RefTerm;
		GRCTerm=Source.GRCTerm;
	}//close def

	 //Assignment Operator
	 GOMatch& operator =(const GOMatch &Source){// open defintion
		 if (this!= &Source){// open non-self assignment consq
				RefGOID=Source.RefGOID;
				GRCGOID=Source.GRCGOID;
				RefECodes=Source.RefECodes;
				GRCECodes=Source.GRCECodes;
				Distance=Source.Distance;
				GRCDepth=Source.GRCDepth;
				RefDepth=Source.RefDepth;
				RefTerm=Source.RefTerm;
				GRCTerm=Source.GRCTerm;
			}
		 return *this;
	 }//close def.

	 //Parameterized constructor
	 GOMatch(int RID, set<string> RECodes, int RDepth, int GID, set<string> GECodes, int GDepth, int Dist){
		RefGOID=RID;
		GRCGOID=GID;
		RefECodes=RECodes;
		GRCECodes=GECodes;
		Distance=Dist;
		GRCDepth=GDepth;
		RefDepth=RDepth;
		RefTerm=NULL;
		GRCTerm=NULL;
	 }//close constructor
 
	 //> OPERATOR overload
	bool operator>(const GOMatch& RHS)const{
		return(GRCGOID>RHS.GRCGOID);
	}

		//< OPERATOR overload
	bool operator<(const GOMatch& RHS)const{
		return(GRCGOID<RHS.GRCGOID);
	}

	 	//== OPERATOR overload
	bool operator==(const GOMatch& RHS)const{
		return(GRCGOID==RHS.GRCGOID);
	}

		//!= OPERATOR overload
	bool operator!=(const GOMatch& RHS)const{
		return(GRCGOID!=RHS.GRCGOID);
	}
	
	//Assemble string for the GOMatch
	//Ref_GOID\tEvidenceCode\tDepth\tGRC_GOID\tEvidenceCode\tDepth\tDistance\n
	string GOMatchString()const{
		string Info="";
		string Convert;
		Info=Info+GO::IDToString(RefGOID);
		for(set<string>::const_iterator RIt=RefECodes.begin(); RIt!=RefECodes.end(); RIt++){//open for loop
			if(RIt!=RefECodes.begin()){
				Info=Info+" "+*RIt;
			}
			else{
				Info=Info+"\t"+*RIt;
			}
		}
		Convert=itos(RefDepth);
		Info=Info+"\t"+Convert+"\t"+GO::IDToString(GRCGOID);
		for(set<string>::const_iterator GIt=GRCECodes.begin(); GIt!=GRCECodes.end(); GIt++){//open for loop
			if(GIt!=RefECodes.begin()){
				Info=Info+" "+*GIt;
			}
			else{
				Info=Info+"\t"+*GIt;
			}
		}
		Convert=itos(GRCDepth);
		Info=Info+"\t"+Convert+"\t";
		Convert=itos(Distance);
		Info=Info+Convert;

		return Info;
	}//close def.

};//close prototype


#endif


