#ifndef SPHC_H_
#define SPHC_H_

double poly6kernel(int, int);
void spikeykernel(int, int, double, double, double);
double viscositykernel(int, int);

void step();
void run();
void init();

int main(int, char);

#endif // SPHC_H_
