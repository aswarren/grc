#include "GO.h"


GO::GO(){//default constructor
		PHead=NULL;
		FHead=NULL;
		CHead=NULL;
		BPBool=false;
		MFBool=false;
		CCBool=false;
		AllAnCounter=-1;

	}//close definition

GO::GO(const GO& Source){//copy constructor
		Storage=Source.Storage;
		GOPoint=Source.GOPoint;
		PHead=Source.PHead;
		PHeads=Source.PHeads;
		FHead=Source.FHead;
		FHeads=Source.FHeads;
		CHead=Source.CHead;
		CHeads=Source.CHeads;
		BPBool=Source.BPBool;
		MFBool=Source.MFBool;
		CCBool=Source.CCBool;
		AllAnCounter=Source.AllAnCounter;
		NotFound=Source.NotFound;
	}

GO& GO::operator =(const GO& Source){//assignment operator
		if (this != &Source){
			Storage=Source.Storage;
			GOPoint=Source.GOPoint;
			PHead=Source.PHead;
			PHeads=Source.PHeads;
			FHead=Source.FHead;
			FHeads=Source.FHeads;
			CHead=Source.CHead;
			CHeads=Source.CHeads;
			BPBool=Source.BPBool;
			MFBool=Source.MFBool;
			CCBool=Source.CCBool;
			AllAnCounter=Source.AllAnCounter;
			NotFound=Source.NotFound;
		}
		return *this;
}//close definition

GO::~GO(){//default destructor
	PHead =NULL;
	PHeads.clear();
	FHead =NULL;
	FHeads.clear();
	CHead =NULL;
	CHeads.clear();
	for(list<GOFunction*>::iterator i=Storage.begin();i!=Storage.end();i++){
		delete (*i);
	}
	Storage.clear();
}

//ReadOBO is for parsing obo file and initializing the various ontologies
//this function is meant to parse the GO obo file and initialize GO
//structures provided by reference
//Parsing structure obtained from GO library in GAIN by Murali et. al.
int GO::ReadOBO(istream *In, const bool& MFunc, const bool& BProc, const bool& CComp){
	MFBool=MFunc;//boolean variables to determine whether to process ontologies
	BPBool=BProc;
	CCBool=CComp;
	int MaxID=0;

		//parse the gene ontology file
	//precondition: gene ontology is correctlt formatted
	string line;
	string Namespace;
	bool NewTerm=false;
	GOFunction  *CurrentFunc=NULL;
	int linenum=0;
	while(getline(*In,line))
	{
		stringstream ss;
		ss<<line;
		
		linenum++;		//count lines
		if(line[0]=='!')	//handle comment
			continue;

		//get initial tag
		string tag;
		getline(ss,tag,':');
		ss.ignore(80,' ');
		
		if(!NewTerm)
		{
			//process [Term] header
			if(line=="[Term]")
			{
				NewTerm=true;
				CurrentFunc=new GOFunction();
//				CurrentFunc->name_space=default_namespace;
				CurrentFunc->Category=Namespace;
			}
			//process default-namespace: header stanza tag
			else if(tag=="default-namespace")
				getline(ss,Namespace);		
				
			continue;
		}
		//handle blank line
		if(line=="")
		{
			if(NewTerm)
			{
				NewTerm=false;
				if(CurrentFunc->Category=="p" && BPBool){
					Storage.push_back(CurrentFunc);
					GOPoint.InsertKey(CurrentFunc->ID, &(Storage.back()));//insert pointer into hash table
					if(!CurrentFunc->Obsolete && !(CurrentFunc->Parents.size()>0)){
						PHead=CurrentFunc;
						PHeads.push_back(CurrentFunc);
						PHead->Depth=0;
					}
				}
				else if(CurrentFunc->Category=="f" && MFBool){
					Storage.push_back(CurrentFunc);
					GOPoint.InsertKey(CurrentFunc->ID, &(Storage.back()));//insert pointer into hash table
					if(!CurrentFunc->Obsolete && !(CurrentFunc->Parents.size()>0)){
						FHead=CurrentFunc;
						FHeads.push_back(CurrentFunc);
						FHead->Depth=0;
					}
				}
				else if(CurrentFunc->Category=="c" && CCBool){
					Storage.push_back(CurrentFunc);
					GOPoint.InsertKey(CurrentFunc->ID, &(Storage.back()));//insert pointer into hash table
					if(!CurrentFunc->Obsolete && !(CurrentFunc->Parents.size()>0)){
						CHead=CurrentFunc;
						CHeads.push_back(CurrentFunc);
						CHead->Depth=0;
					}
				}
				else if(CurrentFunc !=NULL) {
					delete CurrentFunc;
				}
				CurrentFunc=NULL;
			}
			continue;
		}

		//required tags
		if(tag=="id")
		{
			//Currently, assumes all ID's are "GO:" prepended
			ss.ignore(80,':');
			ss>>CurrentFunc->ID;
		}
		else if(tag=="name")
			getline(ss,CurrentFunc->Name,'\n');
		
		//supported, optional tags
		else if(tag=="alt_id")
		{
			//Currently, assumes all ID's are "GO:" prepended
			//ss.ignore(80,':');
			string alt;
			ss>>alt;
			CurrentFunc->AltIDs.push_back(alt);
		}
		else if(tag=="namespace")
				{
//			getline(ss,CurrentFunc->name_space,'\n');
					getline(ss,CurrentFunc->Category,'\n');
					// HACKHACHKHACK. convert namespace to abbreviation.
					if ("biological_process" == CurrentFunc->Category){
						CurrentFunc->Category = "p";
					}
					else if ("molecular_function" == CurrentFunc->Category){
						CurrentFunc->Category = "f";
					}
					else{
						CurrentFunc->Category = "c";
					}
				}
	            
		else if(tag=="def")
		{
			ss.ignore(80,'"');
			getline(ss,CurrentFunc->Definition,'"');
			//TODO: Handle dbxref's here
			//"followed by a dbxref list containing dbxrefs that describe the origin of this definition"
		}
		else if(tag=="comment")
			getline(ss,CurrentFunc->Comment);
		else if(tag=="subset")
		{
			
			string val;
			ss>>val;
			CurrentFunc->Subset.push_back(val);
		}
		else if(tag=="is_obsolete")
		{
			string value;
			getline(ss,value,'\n');
			if(value=="true")
				CurrentFunc->Obsolete=true;
			if(value=="false")
				CurrentFunc->Obsolete=false;
		}
		else if(tag=="use_term"){//Currently, assumes all ID's are "GO:" prepended
			ss.ignore(80,':');
			int Instead;
			ss>>Instead;
			CurrentFunc->UseInstead.push_back(Instead);
		}

		else if(tag=="synonym" || tag=="related_synonym" || tag=="exact_synonym" || tag=="broad_synonym" || tag=="narrow_synonym")
		{
			ss.ignore(80,'"');
			string Synonym;
			getline(ss,Synonym,'"');
			CurrentFunc->Syns.push_back(Synonym);
			//TODO: Handle dbxref's here
			//"followed by a dbxref list containing dbxrefs that describe the origin of this definition"
			//TODO: Handle distinct synonym types
			//http://www.geneontology.org/GO.synonyms.shtml
		}
		else if(tag=="xref_analog")
		{
			string Ref;
			getline(ss,Ref);
			CurrentFunc->XRef.push_back(Ref);
		}
		else if(tag=="xref")
		{
			string Ref;
			getline(ss,Ref);
			CurrentFunc->XRef.push_back(Ref);
		}
		else if(tag=="is_a")
		{
			ss.ignore(80,':');
			int Val;
			ss>>Val;
			CurrentFunc->Parents.insert(RELATE::value_type(Val,tag));
		}
		else if(tag=="relationship")
		{
			string rel;
			getline(ss,rel,' ');
			ss.ignore(80,':');
			int Val;
			ss>>Val;

			if(rel=="is_a"){
				CurrentFunc->Parents.insert(RELATE::value_type(Val,tag));
			}
			else if(rel=="part_of")
			{
				CurrentFunc->Parents.insert(RELATE::value_type(Val,tag));
			}
			else
			{
				cerr<<"WARNING: Unknown Relationship \""<<rel<<"\" Found. Line:"<<linenum<<"\n";
				cerr<<" >> "<<line<<"\n";
			}
		}		
		//unsupported alternate tags
		//these tags are not currently correctly handled by this library
		//they are currently ignored but fair warning should be given to the user
		else if(tag=="is_cyclic")	{
			cerr<<"WARNING: Cyclic Tag Found. Line:"<<linenum<<"\n";
			cerr<<" >> "<<line<<"\n";
		}
		else if(tag=="is_transitive")	{
			cerr<<"WARNING: Transitive Tag Found. Line:"<<linenum<<"\n";
			cerr<<" >> "<<line<<"\n";
		}

		else if(tag=="is_symmetric")	{
			cerr<<"WARNING: Symmetric Tag Found. Line:"<<linenum<<"\n";
			cerr<<" >> "<<line<<"\n";
		}
		else if(tag=="range")	{
			cerr<<"WARNING: Range Tag Found. Line:"<<linenum<<"\n";
			cerr<<" >> "<<line<<"\n";
		}
		else if(tag=="domain")	{
			cerr<<"WARNING: Domain Tag Found. Line:"<<linenum<<"\n";
			cerr<<" >> "<<line<<"\n";
		}
		else if(tag=="use_tag")	{
			cerr<<"WARNING: Use_tag Tag Found. Line:"<<linenum<<"\n";
			cerr<<" >> "<<line<<"\n";
		}
		else if(tag=="xref_unkown")	{
			cerr<<"WARNING: Xref_unknown Tag Found. Line:"<<linenum<<"\n";
			cerr<<" >> "<<line<<"\n";
		}
		else	{
			cerr<<"WARNING: Unknown Tag \""<<tag<<"\" Found. Line:"<<linenum<<"\n";
			cerr<<" >> "<<line<<"\n";
		}
	}
//	cerr<<"done parsing"<<"\n";
	BuildOntology();

	return 0;
	}//close def.


	//  Precondition, read in all GOFunctions to Storage
	//Right now this just initializes the pointers to the top node of each ontology
	//TO DO: Is it useful to initialize the Parent/Child pointers since have direct hashing?
	int GO::BuildOntology(){
		if(BPBool){
			if (PHeads.size()>1){
				cerr<<"WARNING: more than one head functon\n";
				 PHead=PHeads.front();
			}
		}
		if(MFBool){
			if (FHeads.size()>1){
				cerr<<"WARNING: more than one head functon\n";
				FHead=FHeads.front();
			}
		}
		if(CCBool){
			if (CHeads.size()>1){
				cerr<<"WARNING: more than one head functon\n";
				CHead=CHeads.front();
			}
		}

		return 0;
	}//close definition

	//This function returns all the ancestors of a given GOID
	ANCESTOR GO::GetAncestors(const int& GOID){//open definition
		GOFunction* FuncP=Find(GOID);
		ANCESTOR Relatives;


		if (FuncP!=NULL){
			//GOFunction* FuncP=(*FindP);//get pointer to function
			FuncP->UpdateDist(0,GOID,GOID);
			GatherA(&Relatives,FuncP,0,GOID,GOID);//recursive function to gather all pointers into Relatives	
			FuncP=NULL;
		}
		return Relatives;
	}//close definition


	//This function simply returns the ID's of all the ancestors in a set
	ANCESTORID GO::GetAncestorIDs(const int& GOID){//open definition
		GOFunction* FuncP=Find(GOID);
		ANCESTOR Relatives;
		ANCESTORID RelativeID;

		if (FuncP!=NULL){
			//GOFunction* FuncP=(*FindP);//get pointer to function
			FuncP->UpdateDist(0,GOID,GOID);
			GatherA(&Relatives,FuncP,0,GOID,GOID);//recursive function to gather all pointers into Relatives	
			FuncP=NULL;
		}
		int t=0;
		
		for(ANCESTOR::iterator It=Relatives.begin(); It!=Relatives.end(); It++){//open for loop
			if(It->first!=NULL){
				RelativeID.insert(ANCESTORID::value_type(It->first->ID, It->first->Distance));
			}
		}
		return RelativeID;
	}//close defintion

	//GetAllAncestors function gets all the ancestor functons for
	//a group of functons with each distance initialized to
	//the shortest distance from its respective starting point
	int GO::GetAllAncestors(const int& GOID, ANCESTOR& Family){
		if(Family.size()==0){//this counter serves to group id's together when new Family is submitted
			AllAnCounter--;
		}
		GOFunction* FuncP=Find(GOID);
		//ANCESTOR Relatives;

		if (FuncP!=NULL){
			//GOFunction* FuncP=(*FindP);//get pointer to function
			FuncP->UpdateDist(0,GOID,AllAnCounter);
			Family.insert(ANCESTOR::value_type(FuncP,"self"));//insert self UNIQUE to AllAncestors function
			GatherA(&Family,FuncP,0,GOID,AllAnCounter);//recursive function to gather all pointers into Relatives	
			FuncP=NULL;
		}
		
		return 0;
	}


	//recursive function that gathers ancestors into the ANCESTOR map
	//Uses Start Function Pointer to gather parent pointers into ANCESTOR map
	//and then calls GatherA on each Parent
	int GO::GatherA(ANCESTOR* AP, GOFunction* Start,int Dist, const int& OrigID, const int& UpdateID){
		int MaxD=0;
		if(Start->Parents.size()==0){
			MaxD=Dist;//set max distance
		}
		else{//else its got parents

			Dist++;//update distance from start function
			for(RELATE::iterator It=Start->Parents.begin(); It!=Start->Parents.end(); It++){
				int GID=It->first;//get the ID of the parent
				string ParentType=It->second;//get the parent type
				GOFunction* ParentPP=Find(GID);
				if(ParentPP!=NULL){
					bool Cycle=false;
					Cycle=ParentPP->UpdateDist(Dist,OrigID,UpdateID);//update distance from orignal node (Going up)
					if(Cycle){//cycle check
						cerr<<"WARNING: potential cycle detected in GO at "<<GID<<"\n";
						return MaxD;
					}
					AP->insert(ANCESTOR::value_type(ParentPP,ParentType));
					MaxD=GatherA(AP,ParentPP, Dist, OrigID, UpdateID);//recursive call to gather all ancestors
					Start->UpdateDepth(MaxD,UpdateID);//update depth (Going down)
				}
				ParentPP=NULL;
			}//close for loop
		}//close else
		return MaxD;
	}//close def.

	//Function for returning a GOFunction* when given an int ID
	//If the term is obsolete it will return a pointer to a use_term functoin
	//If no use_term ID exists then return pointer to obsolete function
	GOFunction* GO::Find(const int& FindMe){
		GOFunction** FindP=GOPoint.FindKey(FindMe);
		if(FindP!=NULL){
			if(!(*FindP)->Obsolete){//if the term is not obsolete
				return *FindP;
			}
			else if((*FindP)->UseInstead.size()>0){//else if there is a replacement term
				int Alternate=(*FindP)->UseInstead.at(0);
				cerr<<"WARNING: "<<FindMe<<" is obsolete using "<<Alternate<<" instead\n";
				return Find(Alternate);//run Find on first term
			}
			else return (*FindP);
		}
		else {
			if(NotFound.find(FindMe)==NotFound.end()){//if no error message for this ID
				cerr<<"WARNING: "<<FindMe<<" could not be found in the ontology.\n";
				NotFound.insert(FindMe);
			}
			return NULL;//else could not find the function
		}
	}//close definition

	//Function for finding string version of GOID
	GOFunction* GO::Find(const string& ToFind){
			stringstream Convert;
			int TempID;
			Convert<<ToFind;
			Convert.ignore(80,':');
			Convert>>TempID;
			return Find(TempID);
	}

	//This function returns the integer verion of the GO ID string
	int GO::StringToID(const string& TempS){//open definition
		int TempID;
		stringstream SS(TempS);
		SS.ignore(80,':');
		SS>>TempID;
		return TempID;
	}//close defintion

	//Returns bool whether string is a GOID
	bool GO::StringIsGO(const string& TempS){//open defintion
		return TempS.substr(0,3)=="GO:";
	}//close definition

	//Converts an integer ID to a GO formatted string
	string GO::IDToString(const int& TempID){//open definition
		stringstream SS;
		SS<<TempID;//read in the ID
		string TempS;
		SS>>TempS;
		int NumZero=7-TempS.size();//get the number of zero's
		for(int t=0; t<NumZero; t++){
			TempS="0"+TempS;
		}
		TempS="GO:"+TempS;
		return TempS;
	}//close definition

	//returns whether a string is an evidence code
	bool GO::IsECode(const string & TempS){//open definition
		return (TempS=="IEA" || TempS=="IDA" || TempS=="IMP" || TempS=="IC"
			|| TempS=="IEP" || TempS=="IGI" || TempS=="IPI" || TempS=="ISS"
			|| TempS=="NAS" || TempS=="ND" || TempS=="RCA" || TempS=="TAS"
			|| TempS=="NR");
	}




		
