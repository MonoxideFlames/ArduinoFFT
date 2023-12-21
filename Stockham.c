//
//  Stockham.c
//  CProjects
//
//  Created by Adith Talupuru on 11/30/22.
//

#include "Stockham.h"
#include <stdlib.h>
typedef unsigned int uint;// easier to write and read

static double principalRoots[24] = {
	-1, 0,
	0, 1,
	0.7071067812, 0.7071067812,
	0.9238795325, 0.3826834324,
	0.9807852804, 0.195090322 ,
	0.9951847267, 0.09801714033,
	0.9987954562, 0.04906767433,
	0.9996988187, 0.02454122852,
	0.9999247018, 0.01227153829,
	0.9999811753, 0.006135884649,
	0.9999952938, 0.003067956763,
	0.9999988235, 0.001533980186,
};

Complex* fft_internal_stockham(Complex* in, Complex* out, uint length, int t){
	uint x = length / 2;
	uint k = 1;
	uint count = 0;
	uint end;
	
	Complex w;
	Complex z;
	Complex odd;
	
	double temp;
	Complex* pOut;
	Complex* pIn;
	
	while(k < length) {
		w.real = principalRoots[count];
		w.imag = principalRoots[count + 1] * t;
		for(uint i = 0; i < x; i += k){
			z.real = 1;
			z.imag = 0;
			end = i + k;
			for(uint j = i; j < end; j++){
				//odd = mul(z, in[j + x]);
				pIn = in + j + x;
				odd.real = pIn->real * z.real - pIn->imag * z.imag;
				odd.imag = pIn->imag * z.real + pIn->real * z.imag;
				
				//out[2i + j] = add(in[j], odd);
				pIn -= x;
				pOut = out + i + j;
				pOut->real = pIn->real + odd.real;
				pOut->imag = pIn->imag + odd.imag;
				
				//out[2i + j + k] = sub(in[j], odd);
				pOut += k;
				pOut->real = pIn->real - odd.real;
				pOut->imag = pIn->imag - odd.imag;
				
				//z = mul(z, w);
				temp = z.real * w.real - w.imag * z.imag;
				z.imag = z.real * w.imag + w.real * z.imag;
				z.real = temp;
			}
		}
		k = k << 1;
		
		pIn = in;
		in = out;
		out = pIn;
		
		count += 2;
	}
	return in;//since we switch on that last iteration
}
void fft_evaluate(Complex* input, uint length){
	Complex* temp = alloca(length * sizeof(Complex));
	Complex* out = fft_internal_stockham(input, temp, length, 1);
	
	if(out == input){
		return;
	}
	for(uint i = 0; i < length; i++)
		input[i] = temp[i];
}

void fft_DIF(Complex* input, uint length){
	Complex* temp = alloca(length * sizeof(Complex));
	Complex* out = fft_internal_stockham(input, temp, length, -1);
	
	if(out == input){
		return;
	}
	for(uint i = 0; i < length; i++)
		input[i] = temp[i];
}

void fft_interpolate(Complex* input, uint length){
	Complex* temp = alloca(length * sizeof(Complex));
	Complex* out = fft_internal_stockham(input, temp, length, -1);
	
	double norm = 1.0 / length;
	if(out == input){
		for(uint i = 0; i < length; i++) {
			input[i].real *= norm;
			input[i].imag *= norm;
		}
		return;
	}
	for(uint i = 0; i < length; i++){
		input[i].real = temp[i].real * norm;
		input[i].imag = temp[i].imag * norm;
	}
}

void fft_DIT(Complex* input, uint length){
	Complex* temp = alloca(length * sizeof(Complex));
	Complex* out = fft_internal_stockham(input, temp, length, 1);
	
	double norm = 1.0 / length;
	if(out == input){
		for(uint i = 0; i < length; i++) {
			input[i].real *= norm;
			input[i].imag *= norm;
		}
		return;
	}
	for(uint i = 0; i < length; i++){
		input[i].real = temp[i].real * norm;
		input[i].imag = temp[i].imag * norm;
	}
}


