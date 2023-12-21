//
//  FFT.cpp
//  CS Assignment workspace
//
//  Created by Adith Talupuru on 9/27/22.
//

#include "FFT.hpp"
static const Complex principalRoots[24] = {
   Complex(M_PI/1, true),Complex(M_PI/2, true),Complex(M_PI/4, true),Complex(M_PI/8, true),
   Complex(M_PI/16, true),Complex(M_PI/32, true),Complex(M_PI/64, true),Complex(M_PI/128, true),
   Complex(M_PI/256, true),Complex(M_PI/512, true),Complex(M_PI/1024, true),Complex(M_PI/2048, true),
   Complex(M_PI/4096, true),Complex(M_PI/8192, true),Complex(M_PI/16384, true),Complex(M_PI/32768, true),
   Complex(M_PI/65536, true),Complex(M_PI/131072, true),Complex(M_PI/262144, true),Complex(M_PI/524288, true),
   Complex(M_PI/1048576, true),Complex(M_PI/2097152, true),Complex(M_PI/4194304, true),Complex(M_PI/8388608, true)
};

Complex::Complex(double re, double im){
   real = re;
   imag = im;
}
Complex::Complex(double x){
	real = x;
	imag = 0;
}
Complex::Complex(double arg, bool t){
   real = cos(arg);
   imag = sin(arg);
}

Complex::Complex(){
   
}

Complex Complex::operator*(Complex other){
   Complex m;
   m.real = real * other.real - imag * other.imag;
   m.imag = real * other.imag + imag * other.real;
   return m;
}

Complex Complex::operator*(double scale){
   Complex m;
   m.real = real * scale;
   m.imag = imag * scale;
   return m;
}

Complex Complex::operator+(Complex other){
   Complex m;
   m.real = real + other.real;
   m.imag = imag + other.imag;
   return m;
}

Complex Complex::operator-(Complex other){
   Complex m;
   m.real = real - other.real;
   m.imag = imag - other.imag;
   return m;
}

Complex Complex::operator/(Complex other){
   Complex m;
   double s = 1.0 / (other.real * other.real + other.imag + other.imag);
   m.real = (real * other.real + imag * other.imag) * s;
   m.imag = (imag * other.real - real * other.imag) * s;
   return m;
}

double Complex::abs(){
   return sqrt(real * real + imag * imag);
}

double Complex::arg(){
   if(real == 0){
      if(imag > 0){
         return M_PI_2;
      }
      else{
         return -M_PI_2;
      }
   }
   return atan(imag / real) + ((real < 0)? M_PI : 0);
}

void FFT::internal::bitReversalAndSwap(Complex* p,  unsigned int offset,  unsigned int n){
   double swap1 = 0;
   double swap2 = 0;
	unsigned int k = __builtin_clz(n - 1);

   for(int i = 0; i < n; i++){
		//reversing the bits of i:
		unsigned int a = i;
		a = ((a & 0xaaaaaaaa) >> 1) | ((a & 0x55555555) << 1);
		a = ((a & 0xcccccccc) >> 2) | ((a & 0x33333333) << 2);
		a = ((a & 0xf0f0f0f0) >> 4) | ((a & 0x0f0f0f0f) << 4);
		a = ((a & 0xff00ff00) >> 8) | ((a & 0x00ff00ff) << 8);
		a = ((a & 0xffff0000) >> 16) | ((a & 0x0000ffff) << 16);
		//a = __builtin_bitreverse32(i);
		a = a >> k;
      if(a > i){
         a += offset;
         i += offset;
         
         swap1 = p[i].real;
         p[i].real = p[a].real;
         p[a].real = swap1;
         swap2 = p[i].imag;
         p[i].imag = p[a].imag;
         p[a].imag = swap2;
         
         i -= offset;
      }
   }
}

void FFT::internal::transform(Complex* p,  unsigned int offset,  unsigned int n,  unsigned int run, int traverseDirection){
   int i = 30 - __builtin_clz(run);
   Complex root;
   Complex acc;
   Complex pOdd;
   Complex pEven;
    unsigned int halfRun = run >> 1;
    unsigned int endpoint = n + offset;
   for(; run <= n; run = run << 1, i++){
      for(int j = offset; j < endpoint; j += run) {
         root = principalRoots[i];
         root.imag *= traverseDirection;
         acc.real = 1.0;
         acc.imag = 0.0;
         
          unsigned int oddShift = j + halfRun;
         for( unsigned int k = 0; k < halfRun; k++){
            pOdd = acc * p[k + oddShift];
            pEven = p[k + j];
            p[k + j] = pEven + pOdd;
            p[k + oddShift] = pEven - pOdd;
            acc = acc * root;
         }
      }
      halfRun = halfRun << 1;
   }
}

void FFT::internal::fft(Complex* p,  unsigned int offset,  unsigned int n, int traverseDirection){
   bitReversalAndSwap(p, offset, n);
   if(n > 32){
      unsigned int m = n >> 2;//make a smaller size value 1/8th of the size
      std::thread thread0(transform, p, offset        , m, 2, traverseDirection);//run fft on each 1/8 chunk of the data
      std::thread thread1(transform, p, offset +     m, m, 2, traverseDirection);
      std::thread thread2(transform, p, offset + 2 * m, m, 2, traverseDirection);
      transform(p, offset + 3 * m, m, 2, traverseDirection);
      thread0.join();
      thread1.join();
      thread2.join();
      transform( p, offset, n, (n >> 1), traverseDirection);
   }
   else{
      transform(p, offset, n, 2, traverseDirection);
   }
} // fft core operation - very fiddly, so it's obfuscated
// fft core operation - very fiddly, so it's obfuscated
void FFT::DIFFT(Complex* p,  unsigned int length){//easier function calls containing very little reference to confusing parameters
   FFT::internal::fft(p, 0, length, -1);
}

void FFT::DITFT(Complex* p,  unsigned int length){
   FFT::internal::fft(p, 0, length, 1);
   double s = 1.0 / length;
   for(int i = 0; i < length; i++){
      p[i] = p[i] * s;
   }
}

void FFT::interpolate(Complex* signal,  unsigned int numSamples){
   FFT::internal::fft(signal, 0, numSamples, -1);
   double s = 1.0 / numSamples;
   for(int i = 0; i < numSamples; i++){
      signal[i] = signal[i] * s;
   }
}

void FFT::evaluate(Complex* frequencies,  unsigned int numSamples){
   FFT::internal::fft(frequencies, 0, numSamples, 1);
}
void FFT::storeAmplitudes(Complex* frequencies,  unsigned int numSamples, double* amplitudes){
   for(int i = 0; i < numSamples; i++)
      amplitudes[i] = frequencies[i].abs();
}
void FFT::storePhaseShifts(Complex* frequencies,  unsigned int numSamples, double* phaseShifts){
   for(int i = 0; i< numSamples; i++)
      phaseShifts[i] = frequencies[i].arg();
}
