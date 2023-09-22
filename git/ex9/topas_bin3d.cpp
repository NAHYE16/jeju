/*
#2015.12.20: by Jungwook Shin
#purpose: to reduce number of scripts in mcauto-tools
	//2. Operation: sum, avg, weightby (for let), diff 

	#3. operation
	sum: input1 a.bin b.bin c.bin input2 1 2 3 -> 1*a.bin +  if no input2: all 1
	avg: input1 a.bin b.bin c.bin input2 1 2 3 -> sum/sum()
	weightby : input1 & input2 are all distribution 
	sqweightby : sqrt(input1) & input2 are all distribution 
    diff: input1 & input2 are all distribution 
	//2. Operation: sum, avg, weightby (for let), diff 
	#3. operation
	sum: input1 a.bin b.bin c.bin input2 1 2 3 -> 1*a.bin +  if no input2: all 1
	avg: input1 a.bin b.bin c.bin input2 1 2 3 -> sum/sum()
	weightby : input1 & input2 are all distribution 
	sqweightby : sqrt(input1) & input2 are all distribution 
    diff: input1 & input2 are all distribution 

*/

#include <array>
#include <iostream>
#include <valarray>
#include <string>
#include <vector>
#include <map>
#include <sstream>
#include <fstream>
#include <limits>
//for load_input1
#include <sys/stat.h>
//#include <errno.h>
#include <fcntl.h>
#include <unistd.h>
#include <cstdlib>
#include <sys/mman.h>

#ifdef __APPLE__
//Available in OS X v10.0 and later
#include <machine/endian.h>
#include <libkern/OSByteOrder.h>

#define bswap_16 OSSwapConstInt16
#define bswap_32 OSSwapConstInt32 
#define bswap_64 OSSwapConstInt64

#else
#include <byteswap.h>
#include <endian.h>
//maybe 
#define bswap_16 __bswap_16
#define bswap_32 __bswap_32

#endif

using namespace std;

typedef struct {
    string    output_type ; //0:topas, 1: geometry, 2:set
    string    geometry    ; //only for .geometry header output
	vector<string> contents; //additional info to be written in header
    string  X; string Y; string Z; //name
    int    nX; int   nY; int   nZ; //number of voxels
	float  dX; float dY; float dZ; //voxel dimension
    string uX; string uY; string uZ; //dim units
} DoseGeometry;


void print_help(char* s);

std::valarray<double>* load_input1(string filename);

void median_filter(std::valarray<double>& target, int n_pixels, int nX, int nY);
void flip(std::valarray<double>& target, string type, int nX, int nY);

void calculate_std(int n,
                   std::valarray<double>& data, 
                   std::valarray<double>& mean,
                   std::valarray<double>& sq_mean
                   );

void finalize_std(int n,
                std::valarray<double>& mu,
                std::valarray<double>& s);

template<class T>
void write_out(std::valarray<T>& output, const string filename);

template<class T>
void write_out(std::valarray<T>& output, const string filename, DoseGeometry header);


template<class T>
void double2int(std::valarray<T>& out, std::valarray<double>& in, DoseGeometry& geo );

void double2float(std::valarray<float>& out, std::valarray<double>& in, DoseGeometry& geo );

template<class T>
void check_overflow(std::valarray<T>& out, double min, double max);

DoseGeometry read_header(const string filename, vector<string>& header_parameters);
void         write_header(const string filename, DoseGeometry geo);

int main(int argc, char** argv){
	//1. Build a option table
	map<string, vector<string> > parameters = { 
		//option name,  {parameters}
		{"--header"    ,{} },  //header type: topas (default) or geometry
		{"--flip"      ,{} },  //
		{"--filter"    ,{} },  //
		{"--cut"       ,{} },  //(below this value, make all given value, requires two), not yet implemented
		{"--scale"     ,{} },  
		{"--bits"      ,{} }, 
		{"--operation" ,{} }, 
		{"--input1"    ,{} }, 
		{"--input2"    ,{} }, 
		{"--output1"   ,{} }, 
		{"--output2"   ,{} },
		{"--std"       ,{} }, //caculate mean std (%) and print out
	};  
	std::cout<<"# of arguments: "<< argc<<endl;
	if (argc ==1) {print_help(argv[0]); exit(1);};

	for(int i=1; i < argc ; ++i){
		//Find a parameter
		auto it = parameters.find(argv[i]);
		if( it != parameters.end()){
			//Accumulate argv until it meets next option		
			int j = i+1;
			do{
				if(string(argv[j]).compare(0,2,"--") == 0) break;
				it->second.push_back(argv[j])	;
				j++;
			}while( j < argc );
			//Print out options 
			std::cout<<it->first<<" : " ;
			for(auto parm:it->second) std::cout<< parm <<" ";
			std::cout<<std::endl;
		}
	}
    //Defaults for header
	if (parameters["--header"].size() == 0){
        parameters["--header"].push_back("topas"); //default
	}

	bool   mean_std = false;
	//string std;     
	if( parameters["--std"].size() == 1 ){
		//parsing type & endian
		if(parameters["--std"][0] == "1"){
			mean_std = true;
		}else if(parameters["--std"][0] == "0"){
			mean_std = false;
		}else{
			std::cerr<<"--std: should be '1' for yes or '0' no. "<< std::endl; 
            exit(1);
		}
	}else if(parameters["--std"].size() > 1 ){
        std::cerr<<"--std should have 1 or 0. "<< std::endl;
		exit(1);
	}else{
        //do nothing
        mean_std = false;
    }

    bool   median_filter_flag = false;
	if( parameters["--filter"].size() == 1 ){
		//parsing type & endian
		if(parameters["--filter"][0] == "median"){
			median_filter_flag = true;
		}else{
			std::cerr<<"--filter: should be 'median'  or nothing no. "<< std::endl;
            exit(1);
		}
	}

	
	//----- Result distribution -----//
	DoseGeometry geometry  = read_header( parameters["--input1"][0]+"header", parameters["--header"]);
	std::cout<<"\t --- X,Y,Z:"<< geometry.X <<", "<< geometry.Y <<", "<< geometry.Z<<std::endl;
	std::cout<<"\t ---     n:"<< geometry.nX <<", "<< geometry.nY <<", "<< geometry.nZ<<std::endl;
	std::cout<<"\t ---     d:"<< geometry.dX <<", "<< geometry.dY <<", "<< geometry.dZ<<std::endl;
	std::cout<<"\t ---     s:"<< geometry.uX <<", "<< geometry.uY <<", "<< geometry.uZ<<std::endl;
	
    std::valarray<double> result     (0.0, (geometry.nX*geometry.nY*geometry.nZ))  ; //by default, zero filled

    //To calculate std
    std::valarray<double> result_mean (0.0, (geometry.nX*geometry.nY*geometry.nZ))  ; //by default, zero filled
    std::valarray<double> result_std  (0.0, (geometry.nX*geometry.nY*geometry.nZ))  ; //by default, zero filled

	const string op_type = parameters["--operation"][0];
	if (op_type == "sum" || op_type == "avg"){
		std::vector< double > input2 ;

		//2.2 input2: setup numbers 
		if( parameters["--input2"].size() >0){
			if( parameters["--input2"].size() != parameters["--input1"].size() ){
				std::cerr<<"Number of parameters do not match in input1 & input2. "<< std::endl;
				exit(1);
			}
			for(auto f: parameters["--input2"] ) input2.push_back( stod(f) );
		}else{
			input2.resize( parameters["--input1"].size() );
			std::fill(input2.begin(), input2.end(), 1.0);	
		}
		
		std::valarray< double > input2val( input2.data(), input2.size() ) ;

		//2.3 distribution * number

		for(std::size_t i = 0; i < parameters["--input1"].size(); ++i){
		    std::valarray< double >* in = load_input1( parameters["--input1"][i]);
            result += (  (*in) * input2val[i] );

            //std dev
            std::cout<<"weights_sum:"<<input2val[i]<< ", sum:"<< result.sum() <<std::endl;
            if (mean_std) calculate_std(i, (*in), result_mean, result_std );
            delete in;
	    }

		//2.4 Divide by weight numbers	
		if( op_type == "avg"){
			std::cout<<"avg_sum:"<<input2val.sum()<<std::endl;
			result /= input2val.sum();
		}

		std::cout<<"Final sum:"<< result.sum() <<std::endl;
        
						
	}else if( op_type == "weightby"){
		//2.3 distribution1 * distribution2
		//For Dose weighted LET
        std::valarray<double> denum(0.0, (geometry.nX*geometry.nY*geometry.nZ)); 

        for(std::size_t i = 0; i < parameters["--input1"].size(); ++i){
		    std::valarray< double >* in1 = load_input1( parameters["--input1"][i]);
		    std::valarray< double >* in2 = load_input1( parameters["--input2"][i]);
            std::valarray<double> d =  (*in1) * (*in2);
		    result += (d);
			denum  += (*in2);

			std::cout<<"weightby numerator   sum:"<<result.sum()<< std::endl;
			std::cout<<"weightby denumerator sum:"<<denum.sum() << std::endl;

            //std dev
            if(mean_std) calculate_std(i, d, result_mean, result_std );
            
            delete in1;
            delete in2;
	    }

		//2.4 sum(distribution1 * distribution2)/sum(distribution2)
        //result[denum < denum.max()*0.001] = 0.0; //to prevent devide by 0 
		result /= denum; 

	}else if( op_type == "sqweightby"){
        // (sqrt(distribution1) * distribution2)^2
        //std::tuple<>
		// return 
		
        // sum(sqrt(distribution1) * distribution2)/sum(distribution2)
        //result[denum < denum.max()*0.001] = 0.0; //to prevent devide by 0 
		
        std::valarray<double> denum(0.0, (geometry.nX*geometry.nY*geometry.nZ)); 

        for(std::size_t i = 0; i < parameters["--input1"].size(); ++i){
		    std::valarray< double >* in1 = load_input1( parameters["--input1"][i]);
		    std::valarray< double >* in2 = load_input1( parameters["--input2"][i]);
            std::valarray<double> d =  sqrt(*in1) * (*in2);

		    result += d;
			denum  += (*in2);

			std::cout<<"sqweightby numerator   sum:"<<result.sum()<< std::endl;
			std::cout<<"sqweightby denumerator sum:"<<denum.sum() << std::endl;

            //std dev
            if(mean_std) calculate_std(i, d, result_mean, result_std );

            delete in1;
            delete in2;
	    }

        result /= denum; 
        result *= result;                //do we need? ask Maria
        
       
    }else if( op_type == "multi"){
		//2.3 sum(distribution1 * distribution2)
        for(std::size_t i = 0; i < parameters["--input1"].size(); ++i){
		    std::valarray< double >* in1 = load_input1( parameters["--input1"][i]);
		    std::valarray< double >* in2 = load_input1( parameters["--input2"][i]);
            std::valarray<double> d     = (*in1)*(*in2);
		    result += d;
			std::cout<<"multiplication  sum:"<<result.sum()<< std::endl;

            //std dev
            if (mean_std) calculate_std(i, d, result_mean, result_std );

            delete in1;
            delete in2;
	    }

	}else if( op_type == "diff"){
		//2.3 sum(input1 - input2)
		//For Dose weighted LET
		
        for(std::size_t i = 0; i < parameters["--input1"].size(); ++i){
		    std::valarray< double >* in1 = load_input1( parameters["--input1"][i]);
		    std::valarray< double >* in2 = load_input1( parameters["--input2"][i]);
            std::valarray<double> d = (*in1)-(*in2);
		    result += d;
			std::cout<<"multiplication  sum:"<<result.sum()<< std::endl;

            //std dev
            if (mean_std) calculate_std(i,  d, result_mean, result_std );

            delete in1;
            delete in2;
	    }

	} else if( op_type == "mcnamara"){
		//2.3 
		//input1 Dose_phy_fx LETd
        //input2 a/b   Fx
	    const double              ab = stod( parameters["--input2"][0]);
		const double              fx = stod( parameters["--input2"][1]);
        //std::cout<<"alpha/beta: "<<	ab <<", fx size:"<< fx << std::endl;
		std::valarray< double > Dp   = (*load_input1( parameters["--input1"][0]))/fx;
		std::valarray< double > LETd = (*load_input1( parameters["--input1"][1]));
	    result   = 1.0/(2*Dp) * (
                  sqrt( (ab*ab) + 4.0*ab*Dp*(0.999064+(0.35605/ab)*LETd)  + 
                        4.0*Dp*Dp*pow((1.1012-0.0038703*sqrt(ab)*LETd),2)) -
                  ab);
       

	} else{
		std::cerr<<"Operation types should be one of 'sum', 'avg', 'weightby', 'diff' "<< std::endl;
		exit(1);
	}
        
    result[result != result ] = 0.0;
    
    if ( mean_std && parameters["--input1"].size() > 1 ) finalize_std( parameters["--input1"].size(), result_mean, result_std);
	
	//3. apply scale if it is
	if (parameters["--scale"].size() > 1){
		std::cerr<<"--scale should have one or no. "<< std::endl;
		exit(1);
	}
	const double scale = (parameters["--scale"].size()==1) ? stod(parameters["--scale"][0]) : 1.0 ;
	result *= scale;
	std::cout<<"Final sum, max, min :"<< result.sum() <<", " << result.max()
             <<", "<<result.min() << " after scaled by "<<scale<< std::endl;

    //3.a median filter(?)
    if (median_filter_flag){

        median_filter(result, 1, geometry.nX, geometry.nY);

        std::cout<< parameters["--filter"][0] <<" filter is applied... "<<std::endl;

    }

	//4. flip: combination of flipped axis
    //e.g) XYZ -> flipped X & Y & Z
    //     X   -> flipped X only
	switch(parameters["--flip"].size()){
	case 1:
		{
			//we need to read header at leat
			//like awk '$5 == "bins" || $5 == "bin" {print $2, $4, $7, $8}'  XiOProtonLET.binheader
			const int nX = geometry.nX;
			const int nY = geometry.nY;
    		flip(result, parameters["--flip"][0], nX, nY);
    		flip(result_std, parameters["--flip"][0], nX, nY);
			std::cout<<"flipping done... "<< parameters["--flip"][0]<<std::endl;
		}
		break;
		//call 
	case 0:
		//do nothing
		break;
	default:
		std::cerr<<"--flip should have one or no. "<< std::endl;
		exit(1);
	}

    /*
	# 5. bits (f64 to int32, uint32, int16, uint16) | little or big
	# f64l (default),  f64b 
	*/
	string bits;
	bool   little_endian = true;
	if( parameters["--bits"].size() == 1 ){
		//parsing type & endian
		bits  		  = parameters["--bits"][0].substr(0,3);
		little_endian = true;
		if(parameters["--bits"][0].substr(3,4) == "l"){
			little_endian = true;
		}else if(parameters["--bits"][0].substr(3,4) == "b"){
			little_endian = false;
		}else{
			std::cerr<<"--bits: endian should be 'l' or 'b'. "<< std::endl; //'b' for swap , 'l' no action
			exit(1);
		}//if (l or b)
	}else if(parameters["--bits"].size() == 0 ){
		bits  		  = "f64";
		little_endian = true;
	}else{
		std::cerr<<"--bits should have one or no. "<< std::endl;
		exit(1);
	}

	//6. Header filling
	for(auto p: parameters){
		geometry.contents.push_back(string("# "+p.first));
		for(auto t: p.second){
			geometry.contents.push_back(string("#\t"+t));
		}
	}

	//endian conversion & output
	if(bits == "u16"){
		std::valarray<uint16_t> out(result.size()); 
		double2int(out, result, geometry); 
		//not sure. to apply same swap function unsigned and signed.
		if(!little_endian){ for(std::size_t i = 0 ; i < out.size(); ++i  )out[i] = bswap_16(out[i]); }
		if(parameters["--output1"].size() > 0 )write_out(out, parameters["--output1"][0], geometry);

	}else if(bits == "s16"){
		std::valarray<int16_t> out(result.size()); 
		double2int(out, result, geometry);
		if(!little_endian){ for(std::size_t i = 0 ; i < out.size(); ++i  )out[i] = bswap_16(out[i]); }
		if(parameters["--output1"].size() > 0 )write_out(out, parameters["--output1"][0], geometry);

	}else if(bits == "u32"){
		std::valarray<uint32_t> out(result.size()); 
		double2int(out, result,geometry);
		if(!little_endian){ for(std::size_t i = 0 ; i < out.size(); ++i  )out[i] = bswap_32(out[i]); }
		// for(auto f: out) f = bswap_32(f); <-!!! don't use auto, it massed up your data
		if(parameters["--output1"].size() > 0 )write_out(out, parameters["--output1"][0], geometry);
	}else if(bits == "s32"){
		std::valarray<int32_t> out(result.size()); 
		double2int(out, result,geometry);
		if(!little_endian){ for(std::size_t i = 0 ; i < out.size(); ++i  )out[i] = bswap_32(out[i]); }
		if(parameters["--output1"].size() > 0 )write_out(out, parameters["--output1"][0], geometry);
	}else if(bits == "f32"){
		std::valarray<float> out(result.size()); 
		double2float(out,result,geometry);
		if(parameters["--output1"].size() > 0 )write_out( out, parameters["--output1"][0], geometry);
	}else if(bits == "f64"){
		//if(!little_endian){ for(int i = 0 ; i < out.size(); ++i  )out[i] = bswap_32(out[i]);? }
		if(parameters["--output1"].size() > 0 ) write_out( result, parameters["--output1"][0], geometry);
	}else{
		//do nothing, because default is f64
	}

    if(mean_std){
        write_out( result_mean, parameters["--output1"][0]+"_mean", geometry);
        write_out( result_std, parameters["--output1"][0]+"_std", geometry);
    }
    
    return 0;
}

void print_help(char* s){
  cout<<"Usage:   "<<s<<" [-option] [argument]"<<endl;
  cout<<"option:  "<<"--bits [u,s]32[b,l], [u,s]16[b,l], or f32[b,l] "<<endl;
  cout<<"         "<<"--flip XYZ"<<endl;
  cout<<"         "<<"--scale value"<<endl;
  cout<<"         "<<"--operation [avg | sum | diff | weightby]"<<endl;
  cout<<"         "<<"--input1 [file1 file2 file3 ...] "<<endl;
  cout<<"         "<<"--input2 [file1 file2 file3 ...] or [num1 num2 num3 ...]"<<endl;
  cout<<"         "<<"--output1 outputfile"<<endl;
  cout<<"         "<<"--output2 "<<endl;
}

//----- load_input1
std::valarray<double>* load_input1(string filename){
	std::cout<<"Reading ... "<<filename;
	int fd = open(filename.c_str(), O_RDONLY);
	struct stat fds;

	if(fstat(fd,&fds) < 0 ){
		std::cerr<<"Error: "<< filename << std::endl;
		exit(1);
	}
	long total_voxels = fds.st_size/sizeof(double);	
	double* dose  =  (double*) mmap( (void*)0, fds.st_size, PROT_READ,MAP_SHARED, fd,0 );
	std::valarray<double>* cube = new std::valarray<double>(dose, total_voxels );
    munmap(dose, total_voxels*sizeof(double));
    std::cout<<" cube size: "<<cube->size()<<std::endl;

     close(fd);
	 return cube; 

}

//----- flip array
DoseGeometry read_header(const string filename, vector<string>& output_params){

        DoseGeometry headerInfo; 
        
        headerInfo.output_type = output_params[0];
        for(size_t i = 1 ; i<output_params.size() ;++i){
            headerInfo.geometry += output_params[i];
            if(i<output_params.size()-1){
                headerInfo.geometry +=", "     ;
            }  
        }

        // Load header & read
        ifstream fileHeader;
		std::cout<<"Checking file header: "<<filename<<std::endl;	
        fileHeader.open(filename.c_str(), ios::in);
        if (!fileHeader) {
                cout << "ConvertDose ERROR: CANNOT OPEN DOSE HEADER " <<filename << endl;
                exit(EXIT_FAILURE);
        }

        string hline;
        string dummy;
		int    dim_count = 0;
		while(getline(fileHeader, hline)){
        	istringstream inputX(hline);
			vector<string> parsed ;
			string word ;
			headerInfo.contents.push_back(hline);
			while(inputX >> word ){ parsed.push_back(word); }
			//for(auto t: parsed) cout<<"parsed: "<<t<<endl;
			switch( dim_count ){
				case 0: 
				{
					if( parsed[1] == "X" || parsed[1] == "R"  ){
						headerInfo.X = parsed[1];
						headerInfo.nX =stoi(parsed[3]);
						headerInfo.dX =stof(parsed[6]);
						headerInfo.uX =parsed[7];
						++dim_count;
					}
					break;
				}
				case 1:
				{
					if( parsed[1] == "Y" || parsed[1] == "Phi"  ){
						headerInfo.Y = parsed[1];
						headerInfo.nY =stoi(parsed[3]);
						headerInfo.dY =stof(parsed[6]);
						headerInfo.uY =parsed[7];
						++dim_count;
					}
					break;
				}
				case 2:
				{
					if( parsed[1] == "Z" || parsed[1] == "Theta"  ){
						headerInfo.Z = parsed[1];
						headerInfo.nZ =stoi(parsed[3]);
						headerInfo.dZ =stof(parsed[6]);
						headerInfo.uZ =parsed[7];
						++dim_count;
					}
					break;
				}
				default: //do nothing
					break;
			}//switch

		};//while(line)
        
        fileHeader.close();

        return headerInfo;

}

//----- write TOPAS header? mdose header? or dicom? 
void write_header(const string filename, DoseGeometry header ){

    if ( header.output_type == "topas" ){ //file+header
	    ofstream file1( filename + "header", ios::out);
	    for(auto f: header.contents){
		    file1 << f<<std::endl;
	    }
        file1.close();
    }else if(header.output_type == "geometry"){//file+".geometry"
        //this should be updated later
	    ofstream file1( filename + ".geometry", ios::out);
        //for(auto f; header.){
        //      << header.nX <<","<<header.nY<<","<<header.nZ<<std::endl;
        file1 << header.geometry<<std::endl;
        file1.close();
    }else if(header.output_type == "set"){
	    ofstream file1( filename + ".set", ios::out);
        file1 <<"calc-vol \""  << header.geometry <<"\""<<std::endl;
        file1 <<"ct-dir \"\""  <<std::endl;
        string filename_only = "";
        auto p = filename.find_last_of("/");
        if  (p == std::string::npos ){
            p = -1;
        }
        filename_only = filename.substr(p+1);

        file1 <<"dose     \""  <<  filename_only <<"\""<<std::endl; //file name w/o path information
        file1 <<"Ref-dose \""  <<  filename_only <<"\""<<std::endl;
        file1.close();
    }else{
        std::cerr<<"The header type:"<< header.output_type << " is not supported."<<std::endl;
        exit(1);
    }

}

void median_filter(std::valarray<double>& target, int n_pixels, int nX, int nY){
    std::valarray<double> tmp0(target.size());
	tmp0 = target; //copy for temporally

    long int id_from = 0;
	long int id_to   = 0;
	const int nZ     = target.size()/(nX*nY) ;

	//TODO: edge pixels
	//

	//for(int k = n_pixels ; k < nZ-n_pixels ; ++k){
	for(int k = 0 ; k < nZ ; ++k){ //2D
		for(int j=n_pixels ; j < nY-n_pixels ; ++j){
			for(int i=n_pixels ; i < nX-n_pixels ; ++i){
                id_to    = k*nX*nY + j*nX + i;

                //iteration over filter kernel
			    //std::vector<double> median_kernel( (2*n_pixels+1)*(2*n_pixels+1)*(2*n_pixels+1),0 );
			    std::vector<double> median_kernel ; //( (2*n_pixels+1)*(2*n_pixels+1),0 );
			    //for(int kk = -n_pixels; kk <= n_pixels; ++kk){
			        for(int jj = -n_pixels; jj <= n_pixels; ++jj){
			            for(int ii = -n_pixels; ii <= n_pixels; ++ii){
                            //id_from = (k+kk)*nX*nY + (j+jj)*nX + (i+ii);
                            id_from = (k)*nX*nY + (j+jj)*nX + (i+ii);
                            median_kernel.push_back( tmp0[ id_from ]);
			            }
			        }
			    //}
			    std::sort(median_kernel.begin(), median_kernel.end());
			    target[id_to] = median_kernel[ median_kernel.size()/2 ]; //13th element for kernel of 1 n_pixels, i.e., = 27/2
			}//x
		}//y
	}//z


}


//----- flip array
void flip(std::valarray<double>& target, string type, int nX, int nY){
    std::valarray<double> tmp0(target.size()); 
	tmp0 = target; //copy for temporally
	long int id_from = 0;
	long int id_to   = 0;
	const int nZ     = target.size()/(nX*nY) ;
    bool x_flip = false;
    bool y_flip = false;
    bool z_flip = false;

    //For XiO, YZ flipped
    if( type.find("X") != std::string::npos  ) x_flip = true ;
    if( type.find("Y") != std::string::npos  ) y_flip = true ;
    if( type.find("Z") != std::string::npos  ) z_flip = true ;
    
    long int idx = 0;
    long int idy = 0;
    long int idz = 0;
	for(int k =0 ; k < nZ ; ++k){
		for(int j=0 ; j < nY ; ++j){
			for(int i=0 ; i < nX ; ++i){
                idx = (x_flip) ? (nX-1-i)       : i;
                idy = (y_flip) ? (nY-1-j)*nX    : j*nX;
                idz = (z_flip) ? (nZ-1-k)*nX*nY : k*nX*nY;
				id_to    = k*nX*nY + j*nX + i;
				id_from  = idz + idy + idx ; 			
				target[id_to] = tmp0[id_from];
			}//x
		}//y
	}//z	
}


//----- type conversion
template<class T>
void double2int(std::valarray<T>& out, std::valarray<double>& in, DoseGeometry& geo ){
	const double in_min = in.min();
	const double in_max = in.max();
	check_overflow(out, in_min, in_max);
    if (in_max <= 1.0){
        const double slope = (in_max-in_min)/10000.0; 
        //or maximum step size for given data type
	    //const double slope = (in_max-in_min)/(double( std::numeric_limits<T>::highest() - std::numeric_limits<T>::lowest()));
        std::cout<<"#Double 2 integer digitizing step: "<<slope<<std::endl;
	    geo.contents.push_back(string("#!!!! all values are less than 1.0, output is rescaled. (max-min)/1e4") );
	    geo.contents.push_back(string("# double2int digitizing step: "+ to_string(slope) )  );
        for(std::size_t i =0; i< in.size(); ++i ) out[i] =  (T) ((in[i] - in_min)/slope) ; 
	}else{
        //integer values great enough to direct convert
        for(std::size_t i =0; i< in.size(); ++i ) out[i] =  floor( (T) (in[i]+0.5) );
    }//if
}
	
void double2float(std::valarray<float>& out, std::valarray<double>& in, DoseGeometry& geo ){
	const double in_min = in.min();
	const double in_max = in.max();
	check_overflow(out, in_min, in_max);
	for(std::size_t i =0; i< in.size(); ++i ){
		out[i] = (float) (in[i]) ;
	}
}

//----- check overflow 
template<class T>
void check_overflow(std::valarray<T>& out, double min, double max){
	std::cout<< "\tChecking your output ranges..."<<std::endl;

	if( min < std::numeric_limits<T>::lowest()  ){
		std::cout<< "\t lower bound "<< std::numeric_limits<T>::lowest() <<" is greather than your min " << min << " exit(1)"<<std::endl;
		exit(1);
	}
	if( max > std::numeric_limits<T>::max()){
		std::cout<< "\t upper bound "<< std::numeric_limits<T>::max() <<" is smaller than your max" << max << " exit(1)"<<std::endl;
		exit(1);
	}
	std::cout<<"\t---no overflow is detected. \n your output values are in the range ["<< std::numeric_limits<T>::lowest() <<", "
			 << std::numeric_limits<T>::max()<<"]" <<std::endl;

}


//----- output file
template<class T>
void write_out(std::valarray<T>& output, const string filename){
	ofstream file1( filename, ios::out | ofstream::binary);
    file1.write(reinterpret_cast<const char *>(& output[0]), output.size() * sizeof(T));
    file1.close();
}

//----- output file
template<class T>
void write_out(std::valarray<T>& output, const string filename, DoseGeometry header){
	write_out(output, filename);	
	
	header.contents.push_back( string("#\t min:" + to_string(output.min())) ); 
	header.contents.push_back( string("#\t max:" + to_string(output.max())) ); 
	header.contents.push_back( string("#\t sum:" + to_string(output.sum())) ); 

	write_header(filename, header);
}

//----- calculate variance
//n starts 0 thus 
/*
//from TOPAS doc
	for x in data:
        n = n+1
        delta = x - mean
        mean = mean + delta/n
        M2 = M2 + delta*(x - mean)
    sum = n * mean
    variance = M2/(n - 1)
    standard deviation = sqrt(variance)
    population_variance = M2/n
*/
void calculate_std(int n,
                   std::valarray<double>& data, 
                   std::valarray<double>& mean,
                   std::valarray<double>& M2
                   ){
    int j = n + 1 ;
    std::valarray<double> delta  = data - mean;
    mean     += (delta/double(j));
    M2  += (delta*(data - mean));
}


void finalize_std(int n,
                  std::valarray<double>& mu, 
                  std::valarray<double>& s 
                  ){
    //percent standard deviation, s is M2 in calculate_std
    s = s/(n - 1) ;
    s =  100.0*sqrt(s)/mu ;
    s[s != s ] = 0.0;
	std::cout<<"Final std max,min:"<< s.max()<< ", " <<s.min() <<std::endl;
}
