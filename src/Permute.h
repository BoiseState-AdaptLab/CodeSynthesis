#ifndef SYNTH_PERMUTE_HEADER
#define SYNTH_PERMUTE_HEADER
#include <functional>
#include <algorithm>
#include <set>
#include <vector>
#include <map>
#include <math.h>
#include <assert.h>
#include <iostream>
#include <string>
#include <sstream>
#include <math.h>
#include <time.h>
#include <stdlib.h>
// Permute Using Grouped Vectors
template < typename T, typename L, typename C>
class Permute{
   // Linearization
   L lin;
   // Comparator 
   C comp;
   int bitCap = 0b0;
   // Buckets
   std::map<uint32_t,std::vector<std::vector<T>>> buckets;
   int binarySearch(const std::vector<std::vector<T>>& v, std::vector<T>& val){
      int l = 0;
      int r = v.size();
      while (  r >= l){
         int m  = l + (r - l) / 2;
	 if (v[m] == val) return m;
	 if (comp(val,v[m])) r = m - 1;
	 else l = m + 1;
      }
      return -1;
   }
   uint32_t mask = 0b0;
public:
   Permute(L lin,C comp, int bitCap): lin(lin),
	comp(comp),bitCap(bitCap) {
      mask = ~((int)pow(2.0,(double)bitCap) - 1);
   }
   void insert (std::vector<T> tup){
      uint32_t linValue = lin(tup);
      uint32_t bucketValue = linValue & mask;
      //Shift binary to the right so the 
      uint32_t bucketPos = (bucketValue >> bitCap);
      auto& bucket = buckets[bucketPos];
      int insertPos = 0;
      while( insertPos < bucket.size() &&
		comp(bucket[insertPos],tup)) insertPos++;
      // Cant have duplicate tuple as this DS
      // must be a function.
      if (bucket.size() > 0 && bucket[insertPos] == tup) return;
      bucket.insert(bucket.begin() + insertPos,tup); 
   }
  
   uint32_t get (std::vector<T> tup){
      uint32_t linValue = lin(tup);
      uint32_t bucketValue = linValue & mask;
      uint32_t bucketPos = bucketValue >> bitCap;
      auto& bucket = buckets[bucketPos];
      auto localPos = binarySearch(bucket,tup);
      assert(localPos !=-1 && "Tuple does not exist");
      // Count up bucket sizes with bucketPositions
      //
      for( auto it = buckets.begin(); it != buckets.end() &&
		     it->first < bucketPos ; it++){
         localPos+=it->second.size();
      }

      return localPos;
   }

   std::vector<T> getInv(uint32_t pos){
      uint32_t offset = 0;
      for( auto  it = buckets.begin() ;
		      it != buckets.end() ; it++){
	 if ((offset + it->second.size()) > pos ) { 
 	    auto localPos = pos - offset;
	    return it->second[localPos];
       	 }
         offset+= it->second.size();
      }  
      assert(0 && "Position not found");
      return {};
   }

   
};

#define MAXLEVEL 32
// Permute using IndexableSkipList
template <typename T> struct SkipNode{
   SkipNode<T>** forward;
   int * width;
   T value;
   int level;
   SkipNode(int level,T value):level(level),
	value(value){
      forward = new SkipNode <T>* [level+1];
      width = new int [level+1];
      memset(forward,0,sizeof(SkipNode <T>*) * (level+1));
      memset(width,0,sizeof(int) * (level+1));
   }
   ~SkipNode(){
      delete [] forward;
      delete [] width;
   }
};


template < typename T, typename L, typename C>
class PermuteSL{
   uint32_t size = 0;
   // Linearization
   L lin;
   // Comparator 
   C comp;
   typedef SkipNode<std::vector<int>> Node;
   Node * header;
   uint32_t maxLevel;
   uint32_t currentLevel = 0;
   double p = 0.5f;
   inline double random(){
      return (double)rand() / RAND_MAX;
   }
   int randomLevel(){
      int level = 0;
      double ranNum;
      while ((ranNum = random()) < p && level < maxLevel-1){
         level++;
      }
      std::cerr << "level: "<<level << "\n";
      return level;
   }
public:
   PermuteSL(L lin,C comp): lin(lin),
	comp(comp) {
	   maxLevel = MAXLEVEL;
	   header = new Node(maxLevel-1,std::vector<int>());
	}

   uint32_t get (std::vector<T> tup){
      Node* temp = header;
      uint32_t pos = 0;
      for(int i = currentLevel; i >= 0; i--){
            while( temp->forward[i]!= NULL && 
			    comp(temp->forward[i]->value,tup)){
	       //Increment position by width 
	       // before skipping.
	       pos += temp->width[i];
	       temp = temp->forward[i];
	    }
      }
      // The value should be the next value on 
      // the same level
      // -- x→key < searchKey ≤ x→forward[1]→key
      temp = temp->forward[0];
      if (temp!= NULL && temp->value == tup){
         return pos;
      }
      assert( 0 && "tuple does not exist");
      return pos;
   }
   std::vector<T> getInv(uint32_t index){
      // Dont count the head
      uint32_t pos = index+1;
      Node* temp = header;
      for(int i = currentLevel; i >= 0; i--){
           while(temp->forward[i] != NULL && pos >= temp->width[i]){
	       //Increment position by width 
	       // before skipping.
	       pos -= temp->width[i];
	       temp = temp->forward[i];
	    }
      }
      assert(pos >= 0 && temp!=header
		      && "Position does not exit");
      return temp->value;
      
   }
   void insert (std::vector<T> tup){
      int newLevel = randomLevel();
      
      Node*    update[currentLevel+1] = {0};
      uint32_t width[currentLevel+1] = {0};
      Node* temp = header;
      
      for(int i = currentLevel; i >= 0; i--){
	  uint32_t off = 0;

          while(temp->forward[i] != NULL && comp(temp->forward[i]->value,tup)){
	     off += temp->width[i];  
	     temp = temp->forward[i];
	  }
	  update[i] = temp;
	  width[i]  = off;
      }
      
      // Cant have duplicates Permute must be a function
      if (temp!= NULL && temp!=header && temp->forward[0]
		      ->value == tup){
         return;
      }
      temp = new Node(newLevel,tup);
      
      if (newLevel > currentLevel){
         
	 for(int i = currentLevel + 1; i <= newLevel; i++){
	    header->width[i] = size + 1;
	    update[i]        = header; 
	 }
	 currentLevel = newLevel;
      }
     
      uint32_t cumWidth = 0;
      for (int i =0; i <= currentLevel; i++){
	 if ( i > 0){
	    cumWidth += width[i-1];
	    update[i]->width[i]-=cumWidth;
	 }else{
	    update[i]->width[i] = 1;
	 }

         temp->forward[i] = update[i]->forward[i];
	 update[i]->forward[i] = temp;
      } 
      size++;
   }
   std::string toString(){
      std::stringstream ss;
      ss << "SkipList Size: "<< size << "\n"; 
      for(int i = currentLevel; i >= 0; i--){
          Node* temp = header;
	  ss << "Level: " << i << ": ";
          while(temp!= NULL){
	     ss << "{[";
	     for( auto v : temp->value){
	        ss << "," << v;
	     }
	     ss << "],"<<
		   " W: " <<temp->width[i] <<
		  " L:" << temp->level << " } --->";
	     temp = temp->forward[i];
	  }
	  ss << "\n";
      }
      return ss.str();
   }
};

#endif   
