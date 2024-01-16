#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define CL_TARGET_OPENCL_VERSION 300
#include <CL/cl.h>
#include <CL/opencl.h>
#include <math.h>
#include <complex.h>
#include <time.h>

#define PI 3.14159265358979323846265358979323846

float** read_file(char filename[],int *asize){
    FILE* file;
    file = fopen(filename,"r");
    int size = 4; 
    float* rdata = (float*)calloc(size,sizeof(float));
    float* idata = (float*)calloc(size,sizeof(float));
    float n;
    float in;
    int num_rfiles = 0;
    while(fscanf(file,"%f",&n)==1){
	float qn = 0;
	if(fscanf(file,"+%fi",&in)==1){
	     qn = in;
	}
	else if(fscanf(file,"-%fi",&in)==1){
	    qn = -in;
	}
	else{
	     qn=0;
	}
	
        if(num_rfiles>=size){
            size+=1;
            rdata = (float *)realloc(rdata,size*sizeof(float));
	        idata = (float *)realloc(idata,size*sizeof(float));

        }
            rdata[num_rfiles] = n;
	    idata[num_rfiles] = qn;
            num_rfiles++;
    }
    *asize = size;
    fclose(file);
    float** cdata = (float**)calloc(2,sizeof(float*));
    cdata[0] = rdata;cdata[1] = idata;
    return cdata;
}
void printcn(float complex number){
    if(cimag(number)>=0){
        printf("%.16f+%.16fi",creal(number),cimag(number));
    }
    else{
        printf("%.16f%.16fi",creal(number),cimag(number));
    }
}

void show_array(float array[],int size){
    for(int i=0;i<size;i++){
        printf("%f\n",array[i]);
    }
}

float** getFFTs(cl_program program,cl_command_queue queue,cl_context context,float rsignal[],float isignal[],int size){
    cl_int csize = size;
    cl_int qcsize = size/2;
    float* rfreq = (float*)calloc(size,sizeof(float));
    float* ifreq = (float*)calloc(size,sizeof(float));

    cl_mem buffer1 = clCreateBuffer(context,CL_MEM_READ_ONLY|CL_MEM_COPY_HOST_PTR,sizeof(float)*size,rsignal,NULL);
    cl_mem buffer2 = clCreateBuffer(context,CL_MEM_READ_ONLY|CL_MEM_COPY_HOST_PTR,sizeof(float)*size,isignal,NULL);

    cl_mem buffer3 = clCreateBuffer(context,CL_MEM_READ_WRITE,sizeof(float)*size,NULL,NULL);
    cl_mem buffer4 = clCreateBuffer(context,CL_MEM_READ_WRITE,sizeof(float)*size,NULL,NULL);

    cl_int err;
    cl_kernel kernel = clCreateKernel(program,"DFT",&err);
    if(err != CL_SUCCESS){
        printf("Error al crear el Kernel: %d \n",err);
    }
    clSetKernelArg(kernel,0,sizeof(cl_mem),(void*)&buffer1);
    clSetKernelArg(kernel,1,sizeof(cl_mem),(void*)&buffer2);
    clSetKernelArg(kernel,2,sizeof(cl_mem),(void*)&buffer3);
    clSetKernelArg(kernel,3,sizeof(cl_mem),(void*)&buffer4);
    clSetKernelArg(kernel,4,sizeof(cl_int),&qcsize);
    size_t gsize = size/2;
    clEnqueueNDRangeKernel(queue,kernel,1,NULL,&gsize,NULL,0,NULL,NULL);
    clEnqueueReadBuffer(queue,buffer3,CL_TRUE,0,sizeof(float)*size,rfreq,0,NULL,NULL);
    clEnqueueReadBuffer(queue,buffer4,CL_TRUE,0,sizeof(float)*size,ifreq,0,NULL,NULL);

    float** arr = (float **)calloc(size,sizeof(float* ));
    arr[0] = rfreq;
    arr[1] = ifreq;

    return arr;
}
float complex** cmSignal(float complex signal[],int size){
    float complex* mdata = (float complex*)calloc(size/2,sizeof(float complex));
    float complex* ndata = (float complex*)calloc(size/2,sizeof(float complex));
    for(int i=0;i<size/2;i++){
        mdata[i] = signal[2*i];
        ndata[i] = signal[2*i+1];
    }
    float complex** dt = (float complex**)calloc(2,sizeof(float complex*));
    dt[0] = mdata;dt[1]=ndata;
    return dt;
}
float** rmSignal(float rsignal[],float isignal[],int size){
    float * rmdata = (float *)calloc(size/2,sizeof(float));
    float * rndata = (float *)calloc(size/2,sizeof(float));
    float * imdata = (float *)calloc(size/2,sizeof(float));
    float * indata = (float *)calloc(size/2,sizeof(float));
    for(int i=0;i<size/2;i++){
        rmdata[i] = rsignal[2*i];
        rndata[i] = rsignal[2*i+1];
        imdata[i] = isignal[2*i];
        indata[i] = isignal[2*i+1];
    }
    float** dt = (float**)calloc(4,sizeof(float*));
    dt[0] = rmdata;dt[1]=rndata;
    dt[2] = imdata;dt[3]=indata;
    return dt;
}

float complex* iFFT(float* rffts,float* iffts,int size){
    float complex* ft = (float complex*)calloc(size,sizeof(float complex*));
    if(size>2){
        float** dcd = rmSignal(rffts,iffts,size);
        float complex* dft1 = iFFT(dcd[0],dcd[2],size/2);
        float complex* dft2 = iFFT(dcd[1],dcd[3],size/2);

        for(int i=0;i<size/2;i++){
            if(i==0){
                ft[i] = dft1[i] + dft2[i];    
                ft[i+size/2] = dft1[i] - dft2[i];
            }
            else{
                float complex M = cexp(-2*PI*i*I/size);
            
                ft[i] = dft1[i] + M*dft2[i];
                ft[i+size/2] = dft1[i] - M*dft2[i];
            }   
        }
        free(dcd);
        free(dft1);
        free(dft2);
    }
    else{
        ft[0] = rffts[0]+I*iffts[0];ft[1]=rffts[1]+I*iffts[1];
    }
    return ft;
}
float complex* FFT(cl_program program,cl_command_queue queue,cl_context context,float rsignal[],float isignal[],int size){
    float ** ffts = getFFTs(program,queue,context,rsignal,isignal,size);
    float complex* ft = iFFT(ffts[0],ffts[1],size);
    return ft;
}

void show_complex_array(float complex array[],int size){
    for(int i=0;i<size;i++){
        printcn(array[i]);
        printf("\n");
    }
}
void write_data(char filename[],float complex* array,int size){
	int max = 100;
	char* datafile = (char*)malloc(max*size);
	datafile[0]='\0';
	for(int i=0;i<size;i++){
		char rsnumber[max];
        sprintf(rsnumber,"%.5f",creal(array[i]));
		//gcvt(creal(array[i]),5,rsnumber);
		char isnumber[max];
		sprintf(isnumber,"%.5f",cimag(array[i]));
        //gcvt(cimag(array[i]),5,isnumber);
		
		strcat(datafile,rsnumber);
		char carac = '-';
		char* res = strchr(isnumber,carac);
		if(res==NULL){
			strcat(datafile,"+");
		}

		strcat(datafile,isnumber);
		strcat(datafile,"i\n");

    }
	printf("%s",datafile);
	
	FILE* file;
	file = fopen(filename,"w+");
    if(file == NULL){
        printf("No se pudo escribir \n");
    }
    rewind(file);
	fprintf(file,"%s",datafile);
    fclose(file);
}
float complex* IFFT(cl_program program, cl_command_queue queue, cl_context context,float rfreqs[],float ifreqs[],int size){
	for(int i=0;i<size;i++){
		ifreqs[i] = -ifreqs[i];
	}
	float complex* recsignal = FFT(program,queue,context,rfreqs,ifreqs,size);
	for(int i=0;i<size;i++){
		recsignal[i] = conjf(recsignal[i])/(float)size;
	}
	return recsignal;
}
//Archivo para leer el kernel
char* readTextFile(char filename[]){
    FILE* file = fopen(filename,"r");
    if(file == NULL){
        perror("Error al abrir el archivo");
    }
    fseek(file,0,SEEK_END);
    long file_size = ftell(file);
    fseek(file,0,SEEK_SET);
    char* KernelS = (char*)malloc((file_size+1)*sizeof(char));
    fread(KernelS,1,file_size,file);
    fclose(file);
    KernelS[file_size]='\0';
    return KernelS;
}


cl_int qerror = CL_SUCCESS;
cl_int cerror = CL_SUCCESS;
int main(){
    const char* KernelSource = readTextFile("kernel.cl");
    //Obtenemos el id del dispositivo
    const cl_uint num = 1;
    cl_device_type devt = CL_DEVICE_TYPE_CPU;
    clGetDeviceIDs(NULL,devt,0,NULL,(cl_uint*)&num);
    cl_device_id devices[1];
    clGetDeviceIDs(NULL,devt,num,devices,NULL);
    //Creamos el context del dispositivo

    cl_context context = clCreateContextFromType(NULL,devt,NULL,NULL,&cerror);
    if(cerror != CL_SUCCESS){
        printf("Error --- \n");
    }
    //Creamos la linea de comunicacion con el dipositivo
    clGetDeviceIDs(NULL,devt,1,devices,NULL);
    cl_command_queue queue = clCreateCommandQueueWithProperties(context,devices[0],NULL,&qerror);
    //Creamos el programa
    cl_program program = clCreateProgramWithSource(context,1,(const char**)&KernelSource,NULL,NULL);
    //No olvidar esto jaja
    clBuildProgram(program,0,NULL,NULL,NULL,NULL);

    cl_build_status buildStatus;
    clGetProgramBuildInfo(program, devices[0], CL_PROGRAM_BUILD_STATUS, sizeof(cl_build_status), &buildStatus, NULL);
    if (buildStatus != CL_SUCCESS) {
        size_t logSize;
        clGetProgramBuildInfo(program, devices[0], CL_PROGRAM_BUILD_LOG, 0, NULL, &logSize);
        char *log = (char *)malloc(logSize);
        clGetProgramBuildInfo(program, devices[0], CL_PROGRAM_BUILD_LOG, logSize, log, NULL);
        printf("Error de compilación:\n%s\n", log);
        free(log);
        return -1;
    }
    /*
    Nota-
    el 1 en clCreateProgramWithSource nos indica el numero de funciones __kernel que tenemos en el
    string, en caso de añadir mas modificar el numero
    */
    //Creamos el buffer de los datos
    int size;
	int max_name_size = 50;
	char in_filename[max_name_size];
	printf("Escribe el nombre del archivo con los datos: ");
	scanf("%s",in_filename);
    float** data = read_file(in_filename,&size);
    clock_t t0,t1;
	printf("¿Que quieres calcular?\n1:Transformada rapida de Fourier\n2:Transformada inversa\n:"); 
    int opt1;
	scanf("%d",&opt1);
	float complex* fts;
	if(opt1==1){
		t0 = clock();
    	fts = FFT(program,queue,context,data[0],data[1],size);
    	t1 = clock();
	}
	else if(opt1 ==2){
		t0 = clock();
		fts = IFFT(program,queue,context,data[0],data[1],size);
		t1 = clock();
	}
	show_complex_array(fts,size);
	printf("Tiempo de calculo:\n");
    printf("%f  segundos \n",(float)(t1-t0)/CLOCKS_PER_SEC);
	char opt[1];
	printf("¿Deseas guardar los datos transformados?[y/n]: ");
	scanf("%s",opt);
	if(strcmp(opt,"y")==0){
		char out_filename[max_name_size];
		printf("Escribe el nombre del archivo: ");
		scanf("%s",out_filename);
		write_data(out_filename,fts,size);
	}
}
