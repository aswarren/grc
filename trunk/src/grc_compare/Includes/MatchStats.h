// MatchStats.h


//Programmer: Andrew Warren
//email: anwarren@vt.edu
//Date of last modification: 01/xx/05      




#ifndef MatchStats_H
#define MatchStats_H

#include "Match.h"


class MatchStats {//open prototype

public:
	list<Match*> MultiMatch;

	MatchStats(){//default constructor

	}// close default

	~MatchStats(){//default destructor
		for (list<Match*>::iterator It1 =MultiMatch.begin(); It1!=MultiMatch.end(); It1++ ){
			*It1=NULL;
		}
	}//close destructor



			//Copy Constructor
	 MatchStats(const MatchStats &Source){// open defintion
		 MultiMatch=Source.MultiMatch;
	 }// close definition


	 //Assignment Operator
	 MatchStats& operator =(const MatchStats &Source){// open defintion
		 if (this!= &Source){// open non-self assignment consq.
			MultiMatch=Source.MultiMatch;
			}
		 return *this;
	 }// close definition


};//close protoptype

#endif