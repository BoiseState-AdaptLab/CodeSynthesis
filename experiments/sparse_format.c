#include "sparse_format.h"

using namespace std;

void read_sparse_coo (const std::string& filename, 
		coo_d& coo_data){
	bool past_comments = false;
	int i = 0;
	char * elems [10];
	std::string buffer;
	ifstream mat_d;
	mat_d.open(filename);
	if (mat_d){
	  while (getline(mat_d,buffer)){
             if (buffer[0]!='%'){
	       char s[buffer.length()];
               strcpy(s,buffer.c_str());
	       split(s," ",elems);
               if(!past_comments){
		 coo_data.nr = atoi(elems[0]);
		 coo_data.nc = atoi(elems[1]);
		 coo_data.nz = atoi(elems[2]);
		 coo_data.nnz = atoi(elems[3]);
		 coo_data.rows = new int [coo_data.nnz];
		 coo_data.cols = new int [coo_data.nnz];
		 coo_data.zs = new int [coo_data.nnz];
		 coo_data.vals = new float [coo_data.nnz];
		 past_comments = true;
	       }else{
                 coo_data.rows[i] = atoi(elems[0]);
		 coo_data.cols[i] = atoi(elems[1]);
		 coo_data.zs[i] = atoi(elems[2]);
		 coo_data.vals[i] = atof(elems[3]);
		 i++;
               }
	     }
	  }
	}
 	
}


void split(char *s, char* delim, char* result [10]) {
   char * tok = strtok(s,delim);
   int i = 0;
   while (tok!=NULL){
      result[i] = tok;
      tok = strtok(NULL,delim);
      i++;
   }
}
