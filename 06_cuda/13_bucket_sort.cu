#include <cstdio>
#include <cstdlib>

__global__ void sort(int *key, int *bucket, int range)
{
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  atomicAdd(&bucket[key[i]], 1);
  extern __shared__ int count[];
  for (int j = 1; j < range; j <<= 1)
  {
    count[i] = bucket[i];
    __syncthreads();
    if (i >= j)
      bucket[i] += count[i - j];
    __syncthreads();
  }
  for (int j = 0; j < range; j++)
  {
    if (i < bucket[j])
    {
      key[i] = j;
      return;
    }
  }
}

int main()
{
  int n = 50;
  int range = 5;
  int *key, *bucket;
  cudaMallocManaged(&key, n * sizeof(int));
  cudaMallocManaged(&bucket, range * sizeof(int));
  for (int i = 0; i < n; i++)
  {
    key[i] = rand() % range;
    printf("%d ", key[i]);
  }
  printf("\n");

  for (int i = 0; i < range; i++)
    bucket[i] = 0;

  sort<<<1, n, range> > >(key, bucket, range);
  cudaDeviceSynchronize();

  for (int i = 0; i < n; i++)
  {
    printf("%d ", key[i]);
  }
  printf("\n");
  cudaFree(key);
  cudaFree(bucket);
}
