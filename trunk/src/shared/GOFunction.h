//GOFunction.h


//Programmer: Andrew Warren
//email: anwarren@vt.edu
//Date of last modification: 08/xx/06      

//Gene Ontology function stores category, id, and function name

#ifndef GOFunction_H
#define GOFunction_H

#include <map>
#include <set>
#include <string>
using std::map;
using std::string;
using std::set;



class GOFunction{
	friend class GO;
	friend struct LessAncestor;
public:
	int ID; //the GOID
	vector<string> AltIDs; //alternate ID's for this GO Term
	string Name; //the function name
	string Category; //the gene ontology Biological process 'b', molecular function 'f', cellular component 'c'
	map<int,string> Parents;
	vector<int> TempPart;//vector storing part_of
	vector<int> TempIs;//vector storing is_a
	map<int,string> Children;
	int Depth;
	string Definition;
	string Comment;
	vector<string> Subset;
	bool Obsolete;
	vector<string> Syns;
	vector<string> XRef;
	int Distance;
	int DistID;
	vector<int> UseInstead;//vector for use_term id's
        vector<int> ConsiderInstead;
	int VisitNum;
	int UpdateID;
public:
	//default constructor
	GOFunction(){
		ID=-1;
		Name="none";
		Category="none";
		Depth =-1;
		Definition="none";
		Comment="none";
		Obsolete=false;
		Distance=-1;
		DistID=-1;
		VisitNum=0;
		UpdateID=-1;
	}//close definition

	//destructor
	~GOFunction(){
		Parents.clear();
		Children.clear();
	}

	//parameterized constructor
	GOFunction(const int& I, const string& N, const string& C){
		ID=I;
		Name=N;
		Category=C;
		Depth =-1;
		Definition="none";
		Comment="none";
		Obsolete=false;
		Distance=-1;
		DistID=-1;
		VisitNum=0;
		UpdateID=-1;
	}//close def
	
	//copy constructor
	GOFunction(const GOFunction& Source){
		ID=Source.ID;
		AltIDs=Source.AltIDs;
		Name=Source.Name;
		Category=Source.Category;
		Parents=Source.Parents;
		Children=Source.Children;
		Depth =Source.Depth;
		Definition=Source.Definition;
		Comment=Source.Comment;
		Subset=Source.Subset;
		Obsolete=Source.Obsolete;
		Syns=Source.Syns;
		XRef=Source.XRef;
		TempPart=Source.TempPart;
		TempIs=Source.TempIs;
		Distance=Source.Distance;
		DistID=Source.DistID;
		UseInstead=Source.UseInstead;
		VisitNum=Source.VisitNum;
		UpdateID=Source.UpdateID;
                ConsiderInstead=Source.ConsiderInstead;
	}//close def.

	//assignment operator
	GOFunction& operator =(const GOFunction& Source){
		if (this != &Source){
			ID=Source.ID;
			AltIDs=Source.AltIDs;
			Name=Source.Name;
			Category=Source.Category;
			Parents=Source.Parents;
			Children=Source.Children;
			Depth=Source.Depth;
			Definition=Source.Definition;
			Comment=Source.Comment;
			Subset=Source.Subset;
			Obsolete=Source.Obsolete;
			Syns=Source.Syns;
			XRef=Source.XRef;
			TempPart=Source.TempPart;
			TempIs=Source.TempIs;
			Distance=Source.Distance;
			DistID=Source.DistID;
			UseInstead=Source.UseInstead;
			VisitNum=Source.VisitNum;
			UpdateID=Source.UpdateID;
                        ConsiderInstead=Source.ConsiderInstead;
		}
		return *this;
	}//close definition
        
        //This function returns a copy of the Alternate IDs for this Function
        vector<string> GetAltID(){
            return AltIDs;
        }
	
	//returns a string indicating the ontology type
	string OntoType(){
		if(Category=="p" || Category=="b"){
			return "BP";
		}
		else if(Category=="c"){
			return "CC";
		}
		else if(Category=="f" || Category=="m"){
			return "MF";
		}
		else return Category;
	}

	//Function to update distance relative to descendant
	//Also provides the functionality of a cycle check by updating visitnum
	//returns true if cycle is suspected
	bool UpdateDist(const int& D, const int& OrigID, const int& UID){
		if (UID==UpdateID){
			VisitNum++;
			if(VisitNum>10){
				return true;
			}
			if(D<Distance){
				Distance=D;
				DistID=OrigID;//set who the distance is from
			}
		}
		else {
			VisitNum=0;
			Distance=D;
			UpdateID=UID;
			DistID=OrigID;
		}
		if(Parents.size()==0){//if its a root node then the visitnum must be reset here
			VisitNum=0;
		}
		return false;
	}//close definition


		//Function to update depth relative to descendant
	//Takes parameter MaxDistance
	int UpdateDepth(const int& MD, const int& UID){
		VisitNum=0;//reset visit number
		if (UID==UpdateID){
			if((MD-Distance)<Depth ||Depth==-1){
				Depth=MD-Distance;
			}
		}
		else {//else the current up/down traversal is not behaving appropriately
			Depth=(MD-Distance);//based on order of calls this shouldn't happen
			cerr<<"WARNING: error in depth logic\n";
		}
		return 0;
	}//close definition

	//Report Depth
	int ReportDepth(){
		return Depth;
	}

	//Report ID
	int ReportID(){
		return ID;
	}


}; //end GOFunction



struct LessAncestor {
	bool operator()(GOFunction* a1, GOFunction* a2){
		if(a1==NULL){
			return false;
		}
		else if(a2==NULL){
			return false;
		}
		else{
			if(a1->ID==a2->ID){//don't allow equivalent Functons
				return false;
			}
			else if(a1->Distance==a1->Distance){//else if the distances are equivalent but the ID's are not
				return (a1->ID<a2->ID);//ordering depends on ID
			}
			else {//compare distances
				return (a1->Distance<a2->Distance);
			}
		}
	}//close def
};

struct OrderDepth {
	bool operator()(GOFunction* a1, GOFunction* a2){
		if(a1==NULL){
			return false;
		}
		else if(a2==NULL){
			return false;
		}
		else{
			if(a1->ID==a2->ID){//don't allow equivalent Functons
				return false;
			}
			else if(a1->Depth==a1->Depth){//else if the distances are equivalent but the ID's are not
				return (a1->ID<a2->ID);//ordering depends on ID
			}
			else {//compare distances
				return (a1->Depth<a2->Depth);
			}
		}
	}//close def
};

struct FuncReport {
	int ID; //the GOID
	vector<int> AltIDs; //alternate ID's for this GO Term
	string Name; //the function name
	string Category; //the gene ontology Biological process 'b', molecular function 'f', cellular component 'c'
	map<int,string> Parents;
	int Depth;
	string Definition;
	string Comment;
	vector<string> Subset;
	bool Obsolete;
	vector<string> Syns;
	vector<string> XRef;
	int Distance;
};


typedef std::map<int, string> RELATE;//for storing parents id's and edge type





#endif

