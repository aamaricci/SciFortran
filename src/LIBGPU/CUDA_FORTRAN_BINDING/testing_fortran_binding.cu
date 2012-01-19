/*
 * A simple example that scales values in data1, with the corresponding value
 * from data2.
 */

extern "C" __global__ void scale(int *data1, int *data2)
{
	int idx = blockDim.x * blockIdx.y + threadIdx.x;
	data1[idx] *= data2[blockIdx.y];
}

