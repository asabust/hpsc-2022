#include <cstdio>
#include <cmath>
#include <cstring>
#include <mpi.h>

const int nx = 40;
const int ny = 40;
const int nt = 30;
const int nit = 30;

void print(float u[][nx], float v[]][nx], float p[][nx])
{
    for (int j = 1; j < ny - 1; j++)
    {
        for (int i = 1; i < nx - 1; i++)
        {
            // printf("(%.3f) ", p[j][i]);
            printf("(%.2f, %.2f, %.2f) ", u[j][i] * 1000, v[j][i] * 1000, p[j][i] * 1000);
        }
        printf("\n");
    }
    printf("\n");
}

struct Point
{
    float u, v, p;
};

int test_main(int argc, char **argv)
{
    float dx = 2.0 / (nx - 1);
    float dy = 2.0 / (ny - 1);
    float dt = 0.01;
    float rho = 1;
    float nu = .02;

    MPI_Init(&argc, &argv);
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // float u[ny][nx], v[ny][nx], p[ny][nx];
    // float pn[ny][nx], un[ny][nx], vn[ny][nx];
    Point point[nx / size][ny / size], rev_point[nx / size][ny / size];

    for (int i = 0; i < ny / size; ++i)
    {
        for (int j = 0; j < nx / size; ++j)
        {
            point[i][j].u = rev_point[i][j].u = 0;
            point[i][j].v = rev_point[i][j].v = 0;
            point[i][j].p = rev_point[i][j].p = 0;
        }
    }

    int recv_from = (rank + 1) % size;
    int send_to = (rank - 1 + size) % size;
    MPI_Datatype MPI_POINT;
    MPI_Type_contiguous(3, MPI_FLOAT, &MPI_POINT);
    MPI_Type_commit(&MPI_POINT);

    int bx = rank * (nx / size) + 1;
    int ex = (rank + 1) * (nx / size);
    int by = rank * (ny / size) + 1;
    int ey = (rank + 1) * (ny / size);

    MPI_Win win;
    MPI_Win_create(rev_point, (nx / size) * (ny / size) * sizeof(MPI_POINT), sizeof(MPI_POINT), MPI_INFO_NULL, MPI_COMM_WORLD, &win);

    for (int n = 0; n < nt; n++)
    {
        for (int it = 0; it < nit; it++)
        {
            MPI_Win_fence(0, win);
            MPI_Put(rev_point, (nx / size) * (ny / size), MPI_POINT, send_to, 0, (nx / size) * (ny / size), MPI_POINT, win);
            MPI_Win_fence(0, win);
            for (int j = 0; j < ny - 1; j++)
            {
                for (int i = 1; i < nx - 1; i++)
                {
                    float b = 1 / dt * ((point[j][i + 1].u - point[j][i - 1].u) / (2 * dx) + (point[j + 1][i].v - point[j - 1][i].v) / (2 * dy)) - 2 * ((point[j + 1][i].u - point[j - 1][i].u) / (2 * dy) * (point[j][i + 1].v - point[j][i - 1].v) / (2 * dx)) - pow(((point[j][i + 1].u - point[j][i - 1].u) / (2 * dx)), 2) - pow(((point[j + 1][i].v - point[j - 1][i].v) / (2 * dy)), 2);
                    point[j][i].p = (dy * dy * (rev_point[j][i + 1].p + rev_point[j][i - 1].p) +
                                     dx * dx * (rev_point[j + 1][i].p + rev_point[j - 1][i].p) -
                                     b * rho * dx * dx * dy * dy) /
                                    (2 * (dx * dy + dy * dy));
                }
            }

            //
            for (int i = 0; i < nx; i++)
            {
                point[ny - 1][i].p = 0;
                point[0][i].p = point[1][i].p;
            }
            for (int i = 0; i < ny; i++)
            {
                point[i][nx - 1].p = point[i][nx - 2].p;
                point[i][0].p = point[i][1].p;
            }
        }
        MPI_Win_fence(0, win);
        MPI_Put(rev_point, (nx / size) * (ny / size), MPI_POINT, send_to, 0, (nx / size) * (ny / size), MPI_POINT, win);
        MPI_Win_fence(0, win);
        for (int j = 1; j < ny - 1; j++)
        {
            for (int i = 1; i < nx - 1; i++)
            {
                point[j][i].u = rev_point[j][i].u - rev_point[j][i].u * dt / dx * (rev_point[j][i].u - rev_point[j][i - 1].u) -
                                rev_point[j][i].u * dt / dy * (rev_point[j][i].u - rev_point[j - 1][i].u) -
                                dt / (2 * rho * dx) * (rev_point[j][i + 1].p - rev_point[j][i - 1].p) +
                                nu * dt / dx * dx * (rev_point[j][i + 1].u - 2 * rev_point[j][i].u + rev_point[j][i - 1].u) +
                                nu * dt / dy * dy * (rev_point[j + 1][i].u - 2 * rev_point[j][i].u + rev_point[j - 1][i].u);
                point[j][i].v = rev_point[j][i].v - rev_point[j][i].v * dt / dx * (rev_point[j][i].v - rev_point[j][i - 1].v) -
                                rev_point[j][i].v * dt / dy * (rev_point[j][i].v - rev_point[j - 1][i].v) -
                                dt / (2 * rho * dx) * (rev_point[j + 1][i].p - rev_point[j - 1][i].p) +
                                nu * dt / dx * dx * (rev_point[j][i + 1].v - 2 * rev_point[j][i].v + rev_point[j][i - 1].v) +
                                nu * dt / dy * dy * (rev_point[j + 1][i].v - 2 * rev_point[j][i].v + rev_point[j - 1][i].v);
            }
            // printf("%d / u:%f / v:%.2f \n ", j, point[j][10], v.u[j][10]);
        }
        for (int j = 0; j < ny; j++)
        {
            point[j][0].u = 0;
            point[j][ny - 1].u = 0;
            point[j][0].v = 0;
            point[j][ny - 1].v = 0;
        }
        for (int i = 0; i < nx; i++)
        {
            point[0][i].u = 0;
            point[ny - 1][i].u = 1;
            point[0][i].v = 0;
            point[ny - 1][i].v = 0;
        }
        // print(u, v, p);
    }

    MPI_Win_free(&win);
    MPI_Finalize();
    return 0;
}

// mpirun -np 4 ./main-mpi