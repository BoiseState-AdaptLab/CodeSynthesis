#ifndef SYNTH_HEADER
#define SYNTH_HEADER
#include <functional>
#include <algorithm>
#include <set>
#include <vector>
#include <assert.h>
#include <iostream>
#include <string>
#include <sstream>
#include <math.h>

template < typename T, typename L, typename C>
class Permute{
   // Linearization
   L lin;
   // Comparator 
   C comp;
   int bitCap = 0b0;
   int msb = 0b0;
   // Buckets
   std::vector<std::vector<std::vector<T>>> buckets;
   uint32_t binarySearch(const std::vector<std::vector<T>>& v, std::vector<T>& val){
      int l = 0;
      int r = v.size() - 1;
      while ( l <  r){
         int m  = l + (r - l) / 2;
	 if (v[m] == val) return m;
	 if (comp(val,v[m])) r = m - 1;
	 if (comp(v[m],val)) l = m + 1;
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
      if (bucketPos > buckets.size() - 1)
	    buckets.resize(bucketPos+1);
      auto& bucket = buckets[bucketPos];
      int insertPos = 0;
      while( insertPos < bucket.size() &&
		comp(bucket[insertPos],tup)) insertPos++;
      
      bucket.resize(bucket.size() + 1);
      int swapIndex = bucket.size() - 1;
      while(swapIndex != insertPos){
         // SHift forward
         bucket[swapIndex] = bucket[swapIndex -1];
	 swapIndex--;
      }
      bucket[insertPos] = tup;
   }
  
   uint32_t get (std::vector<T> tup){
      uint32_t linValue = lin(tup);
      uint32_t bucketValue = linValue & mask;
      uint32_t bucketPos = bucketValue >> bitCap;
      auto& bucket = buckets[bucketPos];
      auto localPos = binarySearch(bucket,linValue);
      assert(localPos!=-1 && "Tuple does not exist");
      for( int i = 0; i < bucketValue; i++){
         localPos+=buckets[i].size();
      }
      return localPos;
   }

   std::vector<T> getInv(uint32_t pos){
      uint32_t offset = 0;
      for( int i = 0 ; i <buckets.size() ; i++){
	 if ((offset + buckets[i].size()) > pos ) { 
 	    auto localPos = pos - offset;
	    return buckets[i][localPos];
	 }
         offset+= buckets[i].size();
      }  
      assert(0 && "Position not found");
   }

   
};
#endif   
