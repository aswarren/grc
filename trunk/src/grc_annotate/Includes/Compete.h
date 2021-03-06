 // Compete.h


//Programmer: Andrew Warren
//email: anwarren@vt.edu
//Date of last modification: 01/xx/05      




#ifndef Compete_H
#define Compete_H

#include "AARecord.h"
using std::map;



class Compete {//open prototype
		friend std::ostream& operator<<(std::ostream& Out, const Compete& C);

private:
	list<AARecord*> Losers;
	AARecord* Winner;
public:
	Compete(){//default constructor
		Winner=NULL;
	}

	//parameterized constructor
	Compete(AARecord* Win, AARecord* Lose){ // parameterized constructor1
		Winner=Win;//asssign the winner
		Losers.push_back(Lose);//add the loser
	}


		//Copy Constructor
	 Compete(const Compete &Source){// open defintion
		Winner=Source.Winner;
		Losers=Source.Losers;
	}// close definition


	
	 //Assignment Operator
	Compete& operator =(const Compete &Source){// open defintion
		if (this!= &Source){// open non-self assignment consq.
			Winner=Source.Winner;
			Losers=Source.Losers;
		}// close self assignment
		return *this;
	}// close definition

	 //Equality Operator
	bool operator ==(const Compete &Source){// open defintion
		return Winner==Source.Winner;
	}// close definition

	 //Inequality Operator
	bool operator !=(const Compete &Source){// open defintion
		return Winner!=Source.Winner;
	}// close definition

	int AddLoser(AARecord* Lose){//open definition
		Losers.push_back(Lose);
		return 0;
	}//close definition

	int DisplayKO(ostream& Out, const int& GFMin){
		if(Winner==NULL){
			Out<<"Entropy";
		}
		else if(Winner->ReportLength()>=GFMin){
			Out<<(Winner)->ReportID();
		}
		else {
			return 0;//don't report any KOs for this
		}
		for (list<AARecord*>::const_iterator It1 =Losers.begin(); It1!=Losers.end(); It1++ ){
			if((*It1)->ReportLength()>=GFMin){
				Out<<"\t"<<(*It1)->ReportID();
			}
		}//close for loop
		Out<<"\n";
		return 0;
	}

};//close class Compete

typedef map<string,Compete> CompeteMap;


#endif
