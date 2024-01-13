R"(

__kernel void DFT(
        __global float* rsignal,
        __global float* isignal,
        __global float* rfreq,
        __global float* ifreq,
        __const int size){
    int gid = get_global_id(0);
    if(gid<size){
    	int id = gid;
    	float PI = 3.14159265358979323846265358979323846;

        for(int k=0;k<2;k++){
            float rn = 0;
            float in = 0;
            float marg = -PI*k;
            for(int n=0;n<2;n++){
                int idx = id+size*n;
                float arg = marg*n;
                float carg = cos(arg);
                float sarg = sin(arg);
                rn += rsignal[idx]*carg-isignal[idx]*sarg; //Real
                in += rsignal[idx]*sarg+isignal[idx]*carg; //Imaginario
            }
            int qidx = id+size*k;
            rfreq[qidx] = rn;
            ifreq[qidx] = in;
        }
    }   
}
)"
