
  ///////////////////////////////////////////////////////////////////////
  // y = Ax                                                            //
  // A : m-by-n matrix, x : n elements vector, y : m elements vector   // 
  // m and n are arbitrary positive integers                           //
  ///////////////////////////////////////////////////////////////////////

  texture<float4, 2, cudaReadModeElementType> texRefA;

__global__ void mv_kernel( float* y, cudaArray* A, float* x, int m, int n);

//********************************************
//********************************************
//********************************************
//********************************************
//********************************************
//********************************************
//********************************************
//********************************************
//********************************************
//********************************************
//********************************************
//********************************************

extern "C" void matrix_vector_real_optimized(float *y, float *A, float *x, int m, int n)
{
    int blkNum = (m >> 4) + ((m & 15) ? 1 : 0); 
    int height = blkNum << 4;
    int width = (n & 255) ? (((n >> 8) + 1) << 8) : n;

    dim3 threads(16, 16);
    dim3 grid(blkNum, 1);

    cudaArray *d_A;
    float *d_x, *d_y;

    cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float4>();
    cudaMallocArray(&d_A, &channelDesc, width >> 2, height);
    cudaMemcpy2DToArray(d_A, 0, 0, A, n * sizeof(float), n * sizeof(float), m, cudaMemcpyHostToDevice);
    cudaBindTextureToArray(texRefA, d_A);
    cudaMalloc((void **) &d_x, n * sizeof(float));
    cudaMalloc((void **) &d_y, m * sizeof(float));

    cudaMemcpy(d_x, x, n * sizeof(float), cudaMemcpyHostToDevice);
    mv_kernel<<< grid, threads >>>(d_y, d_A, d_x, m, n);
    cudaMemcpy(y, d_y, m * sizeof(float), cudaMemcpyDeviceToHost);

    cudaFree(d_y); cudaFree(d_x); cudaUnbindTexture(texRefA); cudaFreeArray(d_A);
}

//********************************************
//********************************************
//********************************************
//********************************************
//********************************************
//********************************************
//********************************************
//********************************************
//********************************************
//********************************************
//********************************************
//********************************************

#define bx blockIdx.x
#define tx threadIdx.x
#define ty threadIdx.y

__global__ void mv_kernel( float* y, cudaArray* A, float* x, int m, int n)
{
    __shared__ float xs[16][16];
    __shared__ float Ps[16][16];

    float4 a;
    float *Psptr = (float *) Ps + (ty << 4) + tx;
    int ay = (bx << 4) + ty;
    float *xptr = x + (ty << 4) + tx;
    float *xsptr = (float *) xs + (tx << 2);

    *Psptr = 0.0f;
    int i;
    for (i = 0; i < (n & ~255); i += 256, xptr += 256) {
        xs[ty][tx] = *xptr;
        __syncthreads();
        int ax = tx + (i >> 2);
        a = tex2D(texRefA, ax     , ay);
        *Psptr += a.x * *xsptr         + a.y * *(xsptr +   1) + a.z * *(xsptr +   2) + a.w * *(xsptr +   3); 
        a = tex2D(texRefA, ax + 16, ay);
        *Psptr += a.x * *(xsptr +  64) + a.y * *(xsptr +  65) + a.z * *(xsptr +  66) + a.w * *(xsptr +  67); 
        a = tex2D(texRefA, ax + 32, ay);
        *Psptr += a.x * *(xsptr + 128) + a.y * *(xsptr + 129) + a.z * *(xsptr + 130) + a.w * *(xsptr + 131); 
        a = tex2D(texRefA, ax + 48, ay);
        *Psptr += a.x * *(xsptr + 192) + a.y * *(xsptr + 193) + a.z * *(xsptr + 194) + a.w * *(xsptr + 195); 
        __syncthreads();
    }
    
    if (i + (ty << 4) + tx < n) {
        xs[ty][tx] = *xptr;
    }
    __syncthreads();
    int j;
    for (j = 0; j < ((n - i) >> 6); j++, xsptr += 61) {
        a = tex2D(texRefA, tx + (i >> 2) + (j << 4), ay);
        *Psptr += a.x * *xsptr++ + a.y * *xsptr++ + a.z * *xsptr++ + a.w * *xsptr; 
    }
    __syncthreads();
    int remain = (n - i) & 63;
    if ((tx << 2) < remain) {
        a = tex2D(texRefA, tx + (i >> 2) + (j << 4), ay);
        *Psptr += a.x * *xsptr++;
    }
    if ((tx << 2) + 1 < remain) *Psptr += a.y * *xsptr++;
    if ((tx << 2) + 2 < remain) *Psptr += a.z * *xsptr++;
    if ((tx << 2) + 3 < remain) *Psptr += a.w * *xsptr;
    __syncthreads();

    if (tx < 8) *Psptr += *(Psptr + 8);
    if (tx < 4) *Psptr += *(Psptr + 4);
    if (tx < 2) *Psptr += *(Psptr + 2);
    if (tx < 1) *Psptr += *(Psptr + 1);

    __syncthreads();
    if (ty == 0 && (bx << 4) + tx < m) y[(bx << 4) + tx] = Ps[tx][0];
}


//********************************************
//********************************************
//********************************************
//********************************************
//********************************************
//********************************************
//********************************************
//********************************************
//********************************************
//********************************************
//********************************************
//********************************************

