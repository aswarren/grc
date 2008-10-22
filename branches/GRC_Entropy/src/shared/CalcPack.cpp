//CalckPack.cpp
//This file is for calculating various sequence and genome statistics
#include "CalcPack.h"

	CalcPack::CalcPack(){//default constructor
		Lambda=0;
		K=0;
		NumSmallC=NumLargeC=NumLargeNC=NumSmallNC=0;
		TransFile="";
		GenomeFile="";
		OStartCount=0;
		UseCProfile=DefaultCProfile;
		UseNCProfile=DefaultNCProfile;
		UseSmallProf=false;
		GOAccess=NULL;
		CurrentGenomeID="NONE";
		GenomeSize=0;
		for(int t=0; t<20; t++){
			DefaultCProfile[t]=DefaultNCProfile[t]=CProfile[t]=NCProfile[t]\
			=SmallNCProfile[t]=SmallCProfile[t]=0;
		}
	}

	CalcPack::CalcPack(string Matx, string GF, string TF, int TN){//parameterized constructor
		TransFile=TF;
		Matrix=Matx;
		OStartCount=0;
		NumSmallC=NumLargeC=NumLargeNC=NumSmallNC=0;
		UseSmallProf=false;
		UseCProfile=DefaultCProfile;
		UseNCProfile=DefaultNCProfile;
		GOAccess=NULL;
		CurrentGenomeID="NONE";
		GenomeSize=0;
		int Status=InitCodes();//read in the values.
		for(int t=0; t<20; t++){
			DefaultCProfile[t]=DefaultNCProfile[t]=CProfile[t]=NCProfile[t]\
			=SmallNCProfile[t]=SmallCProfile[t]=0;
		}
		if(Status!=0){
			cerr<<"\ngrc_overlap unable to initialize BLAST parameters:exiting...\n";
			throw 20;
		}
		SetDefaultEDP();
                SetupTrans(TN, TransFile);
		Translator.InitCodes();
		GenomeFile=GF;
		GetGenome();
	}

	//Copy Constructor
	//Should not be USED!! pack contains GENOME
	 CalcPack::CalcPack(const CalcPack &Source){// open defintion
		cerr<<"Error Trying to copy CalcPack object";
		throw 15;//throw exception
	}

	//Assignment operator
	CalcPack& CalcPack::operator =(const CalcPack &Source){// open defintion
		if (this!= &Source){// open non-self assignment consq.
			cerr<<"Error trying to copy calcpack object";
			throw 15;//throw exception
		}
		return *this;
	}

	//Function for setting ontological access
	int CalcPack::SetGOAccess(GO* Access){
		GOAccess=Access;
	}


	//Function to initialize parameters based on blast matrix used
	int CalcPack::InitCodes(){//open definition
		ifstream InCode;
		InCode.open(Matrix.c_str());//open matrix specified
		if(!InCode.is_open()){
			cout<<"Error opening MaxMatrix file";
			return 1;
		}
		string Delim;
		string OpCase="";
		int AA; //Max value for conservation of the amino acid
		//InCode>>Delim;

		while(InCode){
			InCode>>Delim;
			if(isupper(Delim[0])){//create lower case in map
				OpCase=tolower(Delim[0]);
				for(int t=1; t<Delim.size(); t++){
					OpCase+=tolower(Delim[t]);
				}
			}
			else {
				OpCase=Delim;
			}

			InCode>>AA; //Read in the AA conservation value
			//CodeArray[TableNum].insert(map<string,char>::value_type(Delim,AA));
			ConsValue.insert(map<string,int>::value_type(OpCase,AA));
			
		}//close while loop
		InCode.close();

		//Set constants Values obtained from FSABLAST
		if(string::npos!=Matrix.find("MaxB62",0)){//if the matrix is BLOSUM62
			Lambda=0.267;
			K=0.041;
		}

		else if(string::npos!=Matrix.find("MaxB80",0)){//if the matrix is BLOSUM80
			Lambda=0.299;
			K=0.071;
		}

		else if(string::npos!=Matrix.find("MaxB45",0)){//if the matrix is BLOSUM45
			Lambda=0.195;
			K=0.032;
		}

		else if(string::npos!=Matrix.find("MaxP30",0)){//if the matrix is PAM30
			Lambda=0.309;
			K=0.15;
		}

		else if(string::npos!=Matrix.find("MaxP70",0)){//if the matrix is PAM30
			Lambda=0.270;
			K=0.060;
		}

		else{//default B62
			Lambda=0.267;
			K=0.041;
		}

		return 0;
	}//close definition
        
        //this function is for reading in the start codon's and their weight
        int CalcPack::SetStarts(const string & filename){
            ifstream In;
            string line="";
            In.open(filename.c_str());
            while(getline(In, line)){
                if(line.find("#", 0 )==string::npos && line!=""){
                    stringstream ss;
                    string Codon;
                    string rCodon;
                    double Weight;
                    ss <<line;
                    ss >>Codon;
                    ss >>Weight;
                    LowerTheCase(Codon);
                    rCodon=ReverseComp(Codon);
                    FStartCodons.insert(map<string,double>::value_type(Codon,Weight));
                    RStartCodons.insert(map<string,double>::value_type(rCodon,Weight));
                }
                line="";
            }
            In.close();
            return 0;
        }
        
        //This function initializes that stop codon sets for forward and reverse
        //frames. 
        //TODO: set the functions for checking codon status for f and r stops
        int CalcPack::SetStops(){
            //get stops from translator
            FStopCodons=Translator.GetStops();
            //set the stops for the resverse frame
            for(set<string>::iterator It=FStopCodons.begin(); It!=FStopCodons.end(); It++){
                string Codon;
                Codon=*It;
                RStopCodons.insert(ReverseComp(Codon));
            }
            return 0;
        }
        
        //Setup the translator and get the stop codons from it
        //Still need to set the Translate to take a Table Number
        int CalcPack::SetupTrans(const int & TN, const string & TF){
            Translator.SetTransFile(TN, TF);
            Translator.InitCodes();
            SetStops();
        }
        
        bool CalcPack::CheckStarts(){
            int FSize=FStartCodons.size();
            int RSize=RStartCodons.size();
            if(FSize==0 || RSize==0){
                cerr<<"error: starts not initialized\n";
		throw 20;
                return false;
            }
            else return true;
        }
        bool CalcPack::CheckStops(){
            if(FStopCodons.size()==0 || RStopCodons.size()==0){
                cerr<<"error: stops not initialized\n";
		throw 20;
                return false;
            }
            else return true;
        }


	int CalcPack::SetDefaultEDP(){
	//Annoying
	//no list initialization of non-static data member....
		DefaultCProfile[0]=0.08468;
		DefaultCProfile[1]=0.01606;
		DefaultCProfile[2]=0.05739;
		DefaultCProfile[3]=0.05752;
		DefaultCProfile[4]=0.04328;
		DefaultCProfile[5]=0.07042;
		DefaultCProfile[6]=0.02942;
		DefaultCProfile[7]=0.05624;
		DefaultCProfile[8]=0.04442;
		DefaultCProfile[9]=0.05620;
		DefaultCProfile[10]=0.03029;
		DefaultCProfile[11]=0.03975;
		DefaultCProfile[12]=0.05116;
		DefaultCProfile[13]=0.04098;
		DefaultCProfile[14]=0.05989;
		DefaultCProfile[15]=0.08224;
		DefaultCProfile[16]=0.05660;
		DefaultCProfile[17]=0.06991;
		DefaultCProfile[18]=0.02044;
		DefaultCProfile[19]=0.03310;

		//setup Default noncoding profile
		DefaultNCProfile[0]=0.07434;
		DefaultNCProfile[1]=0.03035;
		DefaultNCProfile[2]=0.05936;
		DefaultNCProfile[3]=0.04729;
		DefaultNCProfile[4]=0.05662;
		DefaultNCProfile[5]=0.07704;
		DefaultNCProfile[6]=0.05777;
		DefaultNCProfile[7]=0.05328;
		DefaultNCProfile[8]=0.03360;
		DefaultNCProfile[9]=0.05581;
		DefaultNCProfile[10]=0.01457;
		DefaultNCProfile[11]=0.03718;
		DefaultNCProfile[12]=0.04594;
		DefaultNCProfile[13]=0.05977;
		DefaultNCProfile[14]=0.08489;
		DefaultNCProfile[15]=0.05990;
		DefaultNCProfile[16]=0.04978;
		DefaultNCProfile[17]=0.07227;
		DefaultNCProfile[18]=0.01050;
		DefaultNCProfile[19]=0.01974;
		return 0;
	}


	//this function returns the coordinate of the amino acid in the array
	//expects that any non-AA coding seq. enountered will be translated to '*'
	int CalcPack::MapAA(char AA){
		int Value=20;
		if(AA=='*'){
			return Value;
		}
		else{
			Value=AA-'A';
		}
		//26 letters in the alphabet but only 20 AA convert to correct indecies
		if(Value>23){
			Value-=5;
		}
		else if(Value>20){
			Value-=4;
		}
		else if(Value>14){
			Value-=3;
		}
		else if(Value>9){
			Value-=2;
		}
		else if(Value>1){
			Value-=1;
		}
		
		return Value;
	}



	//This function calculates the entropy distance for the section of the genome specified
	double CalcPack::GetEntropy(int AACount[]){//open definition

		double Entropy=-1;
		double TempFreq[20];
		CountToEDP(AACount,TempFreq);
		Entropy=GetEDR(TempFreq);
	
		/*if ( ( TempF = popen ( Command.c_str(), "r") ) != NULL ){
			fgets (TempC, sizeof(double), TempF);
			pclose(TempF); 
			stringstream StoD;
			StoD<<TempC;
			StoD>>Entropy;
			return Entropy;
		}
		else{
			cerr<<"Could not run command with popen in GetEntropy.";
			return Entropy;
		}*/
		return Entropy;
	}//close function



	//create EDP from AACounts
	int CalcPack::CountToEDP(int AACount[], double ResultEDP[]){
		int TotalCount=0;
		int t=0;
		double TotalEntropy=0;
		for(t=0; t<20; t++){
			TotalCount+=AACount[t];
		}
		//convert counts to frequencies
		for(t=0; t<20; t++){
			ResultEDP[t]=double(AACount[t])/(double(TotalCount));
		}
		//convert frequencies to entropies
		for(t=0; t<20; t++){
			if(ResultEDP[t]!=0){
				ResultEDP[t]=-1*ResultEDP[t]*log10(ResultEDP[t]);
				TotalEntropy+=ResultEDP[t];
			}
		}
		for(t=0;t<20;t++){
			ResultEDP[t]=ResultEDP[t]/TotalEntropy;
		}
		return 0;
	}


	//Create new EDP for coding and non coding genes specific to this organism
	int CalcPack::CreateOrgEDPs(){
		//check to see whether there is sufficient data to retrain coding profiles

		//check to see whether to use small profile
		if(NumSmallC>20 && NumSmallNC>20){
			UseSmallProf=true;//default is false
		}
		else{//do not use small profile
			if (NumSmallC>0){
				for(int t=0; t<20; t++){//if not using small profile incorporate that information into coding profile
					CProfile[t]+=SmallCProfile[t];
				}
				NumLargeC+=NumSmallC;
			}
			if (NumSmallNC>0){
				for(int t=0; t<20; t++){
					NCProfile[t]+=SmallNCProfile[t];
				}
				NumLargeNC+=NumSmallNC;
			}
		}//close else do not use small profile
				
		if(UseSmallProf){//if small profiles are needed
			
			for(int t=0; t<20; t++){
				SmallCProfile[t]/=double(NumSmallC);
				SmallNCProfile[t]/=double(NumSmallNC);
				CProfile[t]/=double(NumLargeC);
				NCProfile[t]/=double(NumLargeNC);
			}
		}
		else{//calculate EDP's for just Large CProfile and NCProfile
			for(int t=0; t<20; t++){
				CProfile[t]/=double(NumLargeC);
				NCProfile[t]/=double(NumLargeNC);
			}
		}
		//set the pointers to use the new profiles
		if(NumLargeC>100){//make sure there is sufficient data for retraining
			UseCProfile=CProfile;
		}
		if(NumLargeNC>100){
			UseNCProfile=NCProfile;
		}
		return 0;
	}
        
        //retrieves a specific subsequence from the current genome
        string CalcPack::GenomeSubseq(const bool& Reverse, const long& LB, const long& HB){
            long Length=HB-LB+1;
            long Begin=LB-1;
            string result="";
            if(Reverse){//if in the reverse frame
                    result=ReverseComp((CurrentGenome->second).substr(Begin,Length));	
            }
            else{//not reverse complement
                    result=(CurrentGenome->second).substr(Begin, Length);
            }
            return result;
        }
        
        string CalcPack::GetTrans(const bool& Reverse, const long& LB, const long& HB){
            string Translation="";
            Translation=Translator.TranslateSeq(GeneSequence(CurrentGenome->second, LB, HB, Reverse));
            return Translation;
        }

	//function for finding setting the frequencies of the amino acids in a sequence
	int CalcPack::GetAACount(int AACount[],const long& LB, const long& HB, const bool& Reverse){
		string Translation="";
                Translation=Translator.TranslateSeq(GeneSequence(CurrentGenome->second, LB, HB, Reverse));
		//for each amino acid in the sequence
		for (int t=0; t<Translation.size(); t++){
			int Coord=MapAA(Translation[t]);
			if(Coord>=0&&Coord<=20){
				AACount[Coord]++;
			}
			else{
				cerr<<"logic error in grc_overlap translation\n";
			}
		}
		return 0;
	}
	
	//function for retrieving the Gene's sequence from the genome given the LowBase HighBase and Reverse status
	string CalcPack::GeneSequence(const string& GenomeSeq, const long& LB, const long& HB, const bool& Reverse){
		long Length=HB-LB+1;
		long Begin=LB-1;
		if(Reverse){
			return(ReverseComp(GenomeSeq.substr(Begin,Length)));
		}
		else{
			return(GenomeSeq.substr(Begin, Length));
		}
	}


	//Get the Entropy Distance Ratio
	double CalcPack::GetEDR(double EntropyProfile[]){
		int x=0;
		double CodingDist=0;
		double NonCodingDist=0;
		for (x=0; x<20; x++){
			double CodingDif=EntropyProfile[x]-UseCProfile[x];
			CodingDist+=(CodingDif*CodingDif);
			double NonCodingDif=EntropyProfile[x]-UseNCProfile[x];
			NonCodingDist+=(NonCodingDif*NonCodingDif);
		}
		CodingDist=sqrt(CodingDist);
		NonCodingDist=sqrt(NonCodingDist);
		return CodingDist/NonCodingDist;
	}

	//total up the EDPs from all the orfs
	//to come up with a new EDP for coding and noncoding
	int CalcPack::TotalEDP(int Counts[], const long& Length, const bool& Defeated){

		double TempEDP[20];
		CountToEDP(Counts, TempEDP);

		if (Length<300){
			if (Defeated){//if its a loser add it to the non-coding profile
				NumSmallNC++;
				AddEDP(TempEDP, SmallNCProfile);
			}
			else{//if its a winner add it to the coding profile
				NumSmallC++;
				AddEDP(TempEDP, SmallCProfile);
			}
		}
		else {
			if (Defeated){//if its a loser add it to the non-coding profile
				NumLargeNC++;
				AddEDP(TempEDP, NCProfile);
			}
			else{//if its a winner add it to the coding profile
				NumLargeC++;
				AddEDP(TempEDP, CProfile);
			}
		}
		return 0;
	}


//Add the EDP to the profile submitted so that it can be averaged
	int CalcPack::AddEDP(double TempEDP[], double Profile[]){
			
		for(int t=0; t<20; t++){
			Profile[t]+=TempEDP[t];
		}
		return 0;
	}


	//function to lower the case of all characters in a string
	int CalcPack::LowerTheCase(string & Seq){
		for(int i=0; i<Seq.length(); i++){
			Seq[i]=tolower(Seq[i]);
		}
		return 1;
	}//close definition
		

	
	//This function is designed to check if submitted string is reverse codon
	bool CalcPack::ReverseStart(const string& Codon){//open definition
            CheckStarts();
            return (RStartCodons.find(Codon)!=RStartCodons.end());
	}

	//This function is designed to check if submitted string is forward start
	bool CalcPack::ForwardStart(const string& Codon){//open definition
            CheckStarts();
            return (FStartCodons.find(Codon)!=FStartCodons.end());
	}	

	//This function is designed to check if submitted string a forward stop
	bool CalcPack::ForwardStop(const string& Codon){//open definition
            CheckStops();
            return (FStopCodons.find(Codon)!=FStopCodons.end());
	}
	//This function is designed to check if submitted string a forward stop
	bool CalcPack::ReverseStop(const string& Codon){//open definition
            CheckStops();
            return (RStopCodons.find(Codon)!=RStopCodons.end());
	}
	

	//This function is intended to give a likelihood of correctness relative to other starts
	double CalcPack::CalcSS(const string& Codon){//open definition
            CheckStarts();
            map<string,double>::iterator FindIt;
            FindIt=FStartCodons.end();
            FindIt=FStartCodons.find(Codon);
            if (FindIt!=FStartCodons.end()){
                return (FindIt->second);
            }
            else{
                FindIt=RStartCodons.find(Codon);
                if (FindIt!=RStartCodons.end()){
                    return (FindIt->second);
                }
                else{
                    return .08;
                }
            }
	}//close definition



	//This function returns the complement of a nucleotide passed as a parameter
	char CalcPack::Complement(const char& Base){//open definition
		switch(Base){
			case 'a':
				return 't';
				
			case 't':
				return 'a';
				
			case 'c':
				return 'g';
				
			case 'g':
				return 'c';
			case 'w':
				return 'w';
			case 's':
				return 's';
			case 'k':
				return 'm';
			case 'm':
				return 'k';
			case 'y':
				return 'r';
			case 'r':
				return 'y';
			default: return 'n';
		}
	}//close definition

	//This function gets the genome based on the genome file name
	int CalcPack::GetGenome(){
		ifstream In2;//ofstream operator for reading in the genomic sequence
		In2.open(GenomeFile.c_str());//open up the translated file
		Reader.SetInput(&In2);
		string GenomeID;//for reading in the id
		string Seq;//for reading in the sequence
		string GenomeSeq="";//Initialize to empty string
		string TempID="";
	
		while(Reader.ReadFasta(GenomeID, GenomeSeq)){//read in the genome file
	
			GenomeID=Reader.HeaderToID(GenomeID);//Parse the unessary information from the GenomeID
			
			//used to do individual lookup of each contig/genome based on GenomeID
			//but that is not possible since the overlap function needs relative coordinates
			//from each ORF to correctly judge overlap
			if(Genomes.size()==0){
				Genomes.insert(map<string,string>::value_type(GenomeID, GenomeSeq));//keep track of all the genomes
			}
			else {
				(Genomes.begin()->second)+=GenomeSeq;//cat to create genome
			}
	
			/*FindIt=HitList.find(ID.substr(1,(ID.length())-1));
	
			if(FindIt!=HitList.end()){//if its found then its a hit
				FindIt->second->Sequence=Seq;//assign the sequence
			}*/
		}
		In2.close();//close the input stream
		GenomeSize=(Genomes.begin()->second).size();
		SelectGenome(CurrentGenomeID);
		return 0;
	}

	//This function returns the reverse complement of a given sequence
	string CalcPack::ReverseComp(const string& Forward){
	
		string Comp="";
		for(int s= int(Forward.length())-1; s>=0; s--){
			Comp+=Complement(Forward[s]);
		}
		return Comp;
	}//close definition


	//Calculate the RawBit score from sequence coordinates
	double CalcPack::CalcRawBit(const long& LowB, const long& HighB, const bool& Rev){
		long LB=LowB;
		long HB=HighB;
		long StartSearch=LB-1;//Subtract one to convert to string coordinates
		long Length=HB-LB+1;
		string Codon="";
		double RawBit=0;
		map <string,int>::iterator FindIt;		
	//even though its in reverse go forward through the sequence and get reverse complement
		if(Rev){
			for (long s=0; s<Length; s=s+3){//calc max possible score
				if(s%3==0){//if its the next codon
					Codon=ReverseComp((CurrentGenome->second).substr(StartSearch+s, 3));//get reverse complement
					FindIt=ConsValue.find(Codon);//find the score of this codon
					if(FindIt!=ConsValue.end()){
						RawBit=RawBit+FindIt->second;//add up score
					}
				}
			}//close max loop
		}//close if Reverse

		else {
			for (long s=0; s<Length; s=s+3){//calc max possible score
				if(s%3==0){//if its the next codon
					Codon=(CurrentGenome->second).substr(StartSearch+s, 3);//codon
					FindIt=ConsValue.find(Codon);//find the score of this codon
					if(FindIt!=ConsValue.end()){
						RawBit=RawBit+(FindIt->second);//add up score
					}
				}
			}//close max loop
		}//close not Reverse
		return RawBit;
	}//close definition

	//Function for setting the current Genome to be used based on that genomes fasta ID
	int CalcPack::SelectGenome(string& SeqID){
		if(SeqID==CurrentGenomeID && SeqID!="NONE"){
			return 0;
		}
		if(Genomes.size()>0){
			if(SeqID=="NONE"){
				CurrentGenome=Genomes.begin();
				CurrentGenomeID=Genomes.begin()->first;
			}
			else{
				CurrentGenome=Genomes.find(SeqID);
				if(CurrentGenome==Genomes.end()){
					cerr<<"error in setting genome to be used in calcpack\n";
					throw 20;
				}
				else{
					CurrentGenomeID=CurrentGenome->first;
				}
			}
		}
		else{
			cerr<<"Error: trying to select a genome when none exists\n";
			throw 20;
		}
		return 0;
	}
	
	//Free memory associate with genome
	//Clear Genome free memory associated with genomes
	int CalcPack::ClearGenome(){
		Genomes.clear();
		CurrentGenome=Genomes.end();
		return 0;
	}

	//This function continually adjusts the start site until its back at the original
	//assumes there is a query alignment offset to start at and an original start site to come back to
	bool CalcPack::FindStarts(long& St, const long& OSt, const long& Sp, const long& QAS, const bool& Reverse, double& StartScore){//open definition
                CheckStarts();
                long Start=St;
		long Stop=Sp;
		long OrigStart=OSt;
		long QAlignStart=QAS;
		int StartSearch=0;
		double MaxStartScore=0;
		double Travel=0;
		double NuclDist=(3*QAlignStart);
		int Halt=0;//defines when the search has gone back to original

		//if there is no room to search for a start between the aligned region and current start
		//increase orig start count and return
		//if have returned to Original start (OStartCount>0) then return false
		if(Start==OrigStart){
			if(OStartCount>0){	
				OStartCount=0;//reset for next use
				return false;
			}
			else if(QAlignStart<2){//if there isn't any room for finding another start return true this time and false next time
				OStartCount++;
				if(Reverse){
					StartScore=CalcSS((CurrentGenome->second).substr(Start-3, 3));//get the probability of the codon being start site
				}
				else{
					StartScore=CalcSS((CurrentGenome->second).substr(Start-1, 3));//get the probability of the codon being start site
				}
				return true;
			} 
		}
		if(Reverse){//if the pGene is reversed LOOKING FOR reverse starts in THE FORWARD DIRECTION!!!
			if(Start==OrigStart){
				OStartCount++;//increase count that originalstart has been processed
				StartSearch=(Start-3-1)-(3*QAlignStart);// start search position //minus 3 to account for subtraction difference of(+1*3) and -1 to convert to string coordinates
				Halt=(QAlignStart*3)+3;
			}
			else{//start from the position of the last start found
				StartSearch=Start+2;//next codon 
				Halt=OrigStart-Start;//room left between orig and current search position
			}
	
			for (int s=0; s<=Halt; s=s+3){//search codons in the upstream direction
				if(s%3==0){//if its the next codon
					if(ReverseStart((CurrentGenome->second).substr(StartSearch+s-2, 3))){//if its a start -2 because looking in forward direction
						StartScore=CalcSS((CurrentGenome->second).substr(StartSearch+s-2, 3));//get the probability of the codon being start site
						St=StartSearch+s+1;
						return true;//if back to original start, stop searching
						//}//close max start score
					}//close is start
				}//close next codon
			}//close search codons
		}//close reverse pGene
	
		else{//else its not reversed
			if(Start==OrigStart){
				OStartCount++;//increase the count that original start has been processed
				StartSearch=(Start-1)+(3*QAlignStart);
				Halt=(QAlignStart*3);
			}
			else{
				StartSearch=Start-4;//start at the next codon
				Halt=Start-OrigStart-3;//room left between orig and current search position
			}
	
			for (int s=0; s<=Halt; s=s+3){//search 3 codons in the upstream direction
				if(s%3==0){//if its the next codon
					if(ForwardStart((CurrentGenome->second).substr(StartSearch-s, 3))){//if its a start
						StartScore=CalcSS((CurrentGenome->second).substr(StartSearch-s, 3));//get the probability of the codon being start site
						St=StartSearch-s+1;
						return true;
						//}//close if max score
					}//close if start
				}//close if next codon
			}//close search next codons
		}//close not reversed
		//Default just return original start
		St=OrigStart;
		OStartCount++;
		return true;
	}//close defintion
        
        
        //Writes the genomic sequence to a file
        int CalcPack::WriteGenome(std::ostream& Out){
            for(map<string,string>::iterator It=Genomes.begin(); It!=Genomes.end(); It++){
                Out<<">"<<It->first<<" [mol_type=genomic DNA] [gcode="<<Translator.GetTransCode()<<"]\n";
                FastaRead::OutputSeq(It->second, Out);
            }
            return 0;
        }


/*	//Function to Adjust the start site based on where the alignment begins
int AdjustStart1(long& St, const long& Sp, const long& QAS, const bool& Reverse, const string& Genome){//open definition
	long Start=St;
	long Stop=Sp;
	long QAlignStart=QAS;
	int StartSearch;
	double StartScore=0;
	double MaxStartScore=0;
	double Travel=0;
	int NuclDist=(3*QAlignStart);


	if(QAlignStart<2){
		return 1;
	}
	else if(Reverse){//if the pGene is reversed
		StartSearch=(Start-1)-NuclDist;// start search position

		for (int s=0; s<NuclDist; s+=3){//search codons in the upstream direction
			if(s%3==0){//if its the next codon
				if(ReverseStart(Genome.substr(StartSearch+s-2, 3))){//if its a start
					St=StartSearch+s+1;
					return 1;
					//}//close max start score
				}//close is start
			}//close next codon
		}//close search codons
	}//close reverse pGene

	else{//else its not reversed
		StartSearch=(Start-1)+NuclDist;

		for (int s=0; s<NuclDist; s+=3){//search 3 codons in the upstream direction
			if(s%3==0){//if its the next codon
				if(ForwardStart(Genome.substr(StartSearch-s, 3))){//if its a start
					St=StartSearch-s+1;
					return 1;
				}//close if start
			}//close if next codon
		}//close search next codons
	}//close not reversed
	return 1;
}//close defintion*/ 
