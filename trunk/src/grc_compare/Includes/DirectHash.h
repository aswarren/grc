// DirectHash.h


//Programmer: Andrew Warren
//email: anwarren@vt.edu
//Date of last modification: 04/xx/06      

//DirectHash assumes simple unique hashing of a key string that contains a unique number


#ifndef DirectHash_H
#define DirectHash_H


#include <vector>
using std::vector;
#include <string>
using std::string;
#include <iostream>
using std::ostream;
using std::cerr;

#include <cctype> // for converting characters
  using std::isalnum;
  using std::isdigit;
  using std::isupper;
  using std::tolower;


template <typename T> class DirectHash{


private:
	vector<T>* HTable;//table where everything will hash to
	typename vector<T>::iterator HTIt;//iterator for the vector
	T DefaultValue;//Default value for empty slots in the table
	

public:


	//constructor
	DirectHash(unsigned int StartSize, T DV=T()){//constructor must specify a default value and size
		DefaultValue=DV;//set default value
		HTable =new vector<T>(StartSize,DefaultValue);//allocate size vector of default values
	}//close constructor

	//default destructor
	~DirectHash(){
		HTable->clear();//clear vector
		delete HTable;//deallocate vector
	}

	//Copy Constructor
	 DirectHash(const DirectHash &Source){// open defintion
		 DefaultValue=Source.DefaultValue;
		 HTable=new vector<T>(*(Source.HTable));
	}// close definition


	
	 //Assignment Operator
	 DirectHash& operator =(const DirectHash &Source){// open defintion
		if (this!= &Source){// open non-self assignment consq.
			DefaultValue=Source.DefaultValue;
			HTable=new vector<T>(*(Source.HTable));
		}// close self assignment
		return *this;
	}// close definition

	int InsertKey(unsigned int Key, T Element){//open definition
			if (Key>(HTable->size()-1)){//if the key is greater than the available slots
				if(Key<(HTable->max_size()-1)/2){//if the key is less than half max size
					HTable->resize(Key+1000,T());//resize the hash table
				}
				else if(Key<(HTable->max_size()-1)){//else resize to max size
					HTable->resize(HTable->max_size(),T());
				}
				else{//else report Hash Table overflow
					cerr<<"An integer key is too large for DirectHash\n";
				}
			}//close need resize

			if(HTable->at(Key)==DefaultValue){//if the slot is empty
				HTable->at(Key)=Element;//COPIES the value of the element
			}
			else{//else report a non-unique key
				//cerr<<"An insert of Key "<<Key<<" was attempted twice.\n";
				return -1;
			}

            return 0;
	}//close definition

	//FindKey returns pointer of type T to Key NULL if not there
	T* FindKey(unsigned int Key){
		unsigned int CurSize=HTable->size();
	
		if((Key> CurSize) || (HTable->at(Key))==DefaultValue){
			return NULL;
		}

		else return &(HTable->at(Key));
	}//close FindKey




	unsigned int HashingKey(string ReadIn){//remove all non-numeric characters from string
		unsigned int Result=0;
		//remove leading non-numeric
		while((ReadIn.length()!=0) && (!isdigit(ReadIn[0])) && (ReadIn[0]!='-')){
			ReadIn.erase(0,1);
		}//close while


		for(int t=0; t<ReadIn.size()&& isdigit(ReadIn[t]); t++){//convert to final result
			Result=(Result*10)+(ReadIn[t]-'0');
		}//close convert
		return Result;
	}//close definiton

	
	//friend std::ostream& operator<< <>( std::ostream& Out, const DirectHash& D );
	void Display(ostream& Out){
		for (HTIt= HTable->begin(); HTIt!= HTable->end(); HTIt++){//open for loop
			if((*(HTIt))!=DefaultValue){
				Out<<(*(HTIt));
			}
		}//close for loop
	}//close definition


};//end DirectHash

#endif



