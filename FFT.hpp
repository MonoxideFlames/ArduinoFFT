//
//  FFT.hpp
//  CS Assignment workspace
//
//  Created by Adith Talupuru on 9/27/22.
//

#ifndef FFT_hpp
#define FFT_hpp
#include <cmath>
#include <thread>
struct Complex{
   double real;
   double imag;
   Complex(double, double);
   Complex(double, bool);
   Complex();
   Complex operator*(Complex);
   Complex operator*(double);
   Complex operator+(Complex);
   Complex operator-(Complex);
   Complex operator/(Complex);
   double abs();
   double arg();
	Complex(double);
};
namespace FFT {
   namespace internal {
      void bitReversalAndSwap(Complex*, uint_fast32_t, uint_fast32_t);
      void transform(Complex*, uint_fast32_t, uint_fast32_t, uint_fast32_t, int);
      void fft(Complex*, uint_fast32_t, uint_fast32_t, int);
   }
   void evaluate(Complex*, uint_fast32_t);//actually difft
   void interpolate(Complex*, uint_fast32_t);//actually ditft
   void DITFT(Complex*, uint_fast32_t);//actually interpolation
   void DIFFT(Complex*, uint_fast32_t);//actually evaluation
   void storeAmplitudes(Complex*, uint_fast32_t, double*);
   void storePhaseShifts(Complex*, uint_fast32_t, double*);
}

#endif /* FFT_hpp */
