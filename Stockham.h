//
//  Stockham.h
//  CProjects
//
//  Created by Adith Talupuru on 11/30/22.
//

#ifndef Stockham_h
#define Stockham_h

typedef struct{
	double real;
	double imag;
}Complex;


//Complex* fft_internal_stockham(Complex* input, Complex* output, uint length, int traversalDirection);
void fft_evaluate(Complex* input, unsigned int length);
void fft_interpolate(Complex* input, unsigned int length);
void fft_DIT(Complex* input, unsigned int length);
void fft_DIF(Complex* input, unsigned int length);

#endif /* Stockham_h */
