#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <time.h>

#define PI 3.14159265358979323846

//Calcula la DFT para una señal
//Ojo, no es fft
double complex* DFT(double complex signal[],int size){
    double complex* ft = (double complex*)calloc(size,sizeof(double complex));
    for(int k=0;k<size;k++){
        double complex s = 0;
        for(int i=0;i<size;i++){
            s += signal[i]*cexp(-I*PI*k*i);

        }
        ft[k]=s;
    }
    return ft;
}
double** mSignal(double signal[],int size){
    double* mdata = (double*)calloc(size/2,sizeof(double));
    double* ndata = (double*)calloc(size/2,sizeof(double));
    for(int i=0;i<size/2;i++){
        mdata[i] = signal[2*i];
        ndata[i] = signal[2*i+1];
    }
    double** dt = (double**)calloc(2,sizeof(double*));
    dt[0] = mdata;dt[1]=ndata;
    return dt;
}
double complex** cmSignal(double complex signal[],int size){
    double complex* mdata = (double complex*)calloc(size/2,sizeof(double complex));
    double complex* ndata = (double complex*)calloc(size/2,sizeof(double complex));
    for(int i=0;i<size/2;i++){
        mdata[i] = signal[2*i];
        ndata[i] = signal[2*i+1];
    }
    double complex** dt = (double complex**)calloc(2,sizeof(double complex*));
    dt[0] = mdata;dt[1]=ndata;
    return dt;
}

double complex* FFTi(double complex signal[],double complex m,int k,int size,int r){
    double complex* ift = (double complex*)calloc(2,sizeof(double complex));
    if(size>2){
        double complex** dcd = cmSignal(signal,size);
        double complex m1;
        if(k==0){
            m1 = 1;
        }
        else{
            double complex m1 = cexp(-4*PI*I*k/size);
        }
        double complex* Eks = FFTi(dcd[0],m1,k,size/2,1);
        double complex* Oks = FFTi(dcd[1],m1,k,size/2,1);

        double complex Ek = Eks[0];
        double complex Ok = Oks[0];
        ift[0] = Ek+m*Ok;
        if(r==0){
            ift[1] = Ek-m*Ok;
        }

        free(dcd);
        free(Eks);
        free(Oks);
    }
    else{
        ift[0] += signal[0];
        if(k%2==0){
            ift[1] += signal[1];
        }
        else{
            ift[1] -= signal[1];
        }   
    }
    return ift;
}


double complex* FFT(double complex signal[],int size){
    double complex* ft = (double complex*)calloc(size,sizeof(double complex));
    double complex mq = cexp(-2*PI*I/size);
    for(int k=0;k<size/2;k++){
        double complex m = pow(mq,k);
        double complex* ftq = FFTi(signal,m,k,size,0);
        ft[k] = ftq[0];
        ft[k+size/2] = ftq[1];

        free(ftq);
    }
    return ft;
}


double complex* FFT1(double complex signal[],int size){
    double complex* ft = (double complex*)calloc(size,sizeof(double complex));
    if(size>2){
        double complex** dcd = cmSignal(signal,size);
        double complex* dft1 = FFT1(dcd[0],size/2);
        double complex* dft2 = FFT1(dcd[1],size/2);

        for(int i=0;i<size/2;i++){
            if(i==0){
                ft[i] = dft1[i] + dft2[i];    
                ft[i+size/2] = dft1[i] - dft2[i];
            }
            else{
                double complex M = cexp(-2*PI*i*I/size);
            
                ft[i] = dft1[i] + M*dft2[i];
                ft[i+size/2] = dft1[i] - M*dft2[i];
            }   
        }
        free(dcd);
        free(dft1);
        free(dft2);
    }
    else{
        ft = DFT(signal,size);
    }
    return ft;
}

double complex* FFT2(double complex signal[], int size) {
    double complex* ft = (double complex*)calloc(size, sizeof(double complex));

    if (size > 2) {
        double complex* dft1 = FFT2(signal, size / 2);
        double complex* dft2 = FFT2(signal + size / 2, size / 2);

        for (int i = 0; i < size / 2; i++) {
            double complex M = cexp(-2 * PI * I * i / size);

            ft[i] = dft1[i] + M * dft2[i];
            ft[i + size / 2] = dft1[i] - M * dft2[i];
        }

        free(dft1);
        free(dft2);
    } else {
        ft = DFT(signal, size); // Considera cambiar esto por la implementación real de DFT
    }

    return ft;
}


void printcn(double complex number){
    if(cimag(number)>=0){
        printf("%.16f+%.16fi",creal(number),cimag(number));
    }
    else{
        printf("%.16f%.16fi",creal(number),cimag(number));
    }
}

void show_array(double array[],int size){
    for(int i=0;i<size;i++){
        printf("%f\n",array[i]);
    }
}

void show_complex_array(double complex array[],int size){
    for(int i=0;i<size;i++){
        printcn(array[i]);
        printf("\n");
    }
}

int main(){
    FILE *archivo;
    char filename[] = "data.txt";
    archivo = fopen(filename,"r");
        
    int size = 5;
    double complex* signal = (double complex*)calloc(size,sizeof(double complex));

    double complex numero;
    int cantidad = 0;
    while(fscanf(archivo,"%lf",&numero)==1){
        if(cantidad>=size){
            size +=1;
            signal = (double complex*)realloc(signal,size*sizeof(double complex));
        }
        signal[cantidad] = numero;
        cantidad++;
    }
    FFT1(signal,size);

    clock_t c1,c2,c3;
    c1 = clock();
    double complex* ftf = FFT1(signal,size);
    c2 = clock();
    double complex* ftf1 = FFT2(signal,size);
    c3 = clock();
    show_complex_array(ftf,size);
    printf("\n-------------------\n");
    printf("Tiempo 1:%f",(double)(c2-c1)/CLOCKS_PER_SEC);
    printf("\n-------------------\n");
    //show_complex_array(ftf1,size);
    printf("Tiempo 2:%f",(double)(c3-c2)/CLOCKS_PER_SEC);
    
}