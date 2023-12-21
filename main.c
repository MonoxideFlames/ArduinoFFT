//
//  main.cpp
//  CPPProject
//
//  Created by Adith Talupuru on 11/30/22.
//

#include <stdio.h>
#include "NTT.h"
int main(int argc, const char** argv){
	U32 a[] = {10000,10000,10000,10000,10000,10000,10000,10000};
	U32 m  = 10001;
	FFT_evaluate(a, 8, m);
	for(int i = 0; i < 8; i++){
		a[i] = (a[i])% m;
	}
	FFT_interpolate(a, 8, m);

	for(int i = 0; i < 8; i++){
		printf("%d, ", (a[i] * 8751) % m);
	}
	return 0;
}
