#include <stdio.h>
#include <math.h>
#include <string.h>
#include <omp.h>

int main()
{
    const int nx = 21;
    const int ny = 21;
    const int nt = 30;
    const int nit = 30;
    float dx = 2.0 / (nx - 1);
    float dy = 2.0 / (ny - 1);
    float dt = 0.01;
    float rho = 1;
    float nu = .02;

    float u[ny][nx], v[ny][nx], p[ny][nx];
    float pn[ny][nx], un[ny][nx], vn[ny][nx];
    for (int j = 1; j < ny - 1; j++)
    {
        for (int i = 1; i < nx - 1; i++)
        {
            pn[j][i] = p[j][i] = 0;
            un[j][i] = u[j][i] = 0;
            vn[j][i] = v[j][i] = 0;
        }
    }

#pragma omp parallel shared(p, u, v)
    for (int n = 0; n < nt; n++)
    {
        for (int it = 0; it < nit; it++)
        {
            for (int j = 1; j < ny - 1; j++)
            {
                for (int i = 1; i < nx - 1; i++)
                {
                    pn[j][i] = p[j][i];
                }
            }
#pragma omp parallel for
            for (int j = 1; j < ny - 1; j++)
            {
                for (int i = 1; i < nx - 1; i++)
                {
                    float b = 1 / dt * ((u[j][i + 1] - u[j][i - 1]) / (2 * dx) + (v[j + 1][i] - v[j - 1][i]) / (2 * dy)) - 2 * ((u[j + 1][i] - u[j - 1][i]) / (2 * dy) * (v[j][i + 1] - v[j][i - 1]) / (2 * dx)) - ((u[j][i + 1] - u[j][i - 1]) / (2 * dx)) * ((u[j][i + 1] - u[j][i - 1]) / (2 * dx)) - (((v[j + 1][i] - v[j - 1][i]) / (2 * dy)) * ((v[j + 1][i] - v[j - 1][i]) / (2 * dy)));
                    p[j][i] = (dy * dy * (pn[j][i + 1] + pn[j][i - 1]) +
                               dx * dx * (pn[j + 1][i] + pn[j - 1][i]) -
                               b * rho * dx * dx * dy * dy) /
                              (2 * (dx * dy + dy * dy));
                }
            }

            for (int i = 0; i < nx; i++)
            {
                p[ny - 1][i] = 0;
                p[0][i] = p[1][i];
            }

            for (int i = 0; i < ny; i++)
            {
                p[i][nx - 1] = p[i][nx - 2];
                p[i][0] = p[i][1];
            }
        }

        for (int j = 1; j < ny - 1; j++)
        {
            for (int i = 1; i < nx - 1; i++)
            {
                un[j][i] = u[j][i];
                vn[j][i] = v[j][i];
            }
        }

#pragma omp parallel for
        for (int j = 1; j < ny - 1; j++)
        {
            for (int i = 1; i < nx - 1; i++)
            {
                u[j][i] = un[j][i] - un[j][i] * dt / dx * (un[j][i] - un[j][i - 1]) -
                          un[j][i] * dt / dy * (un[j][i] - un[j - 1][i]) -
                          dt / (2 * rho * dx) * (p[j][i + 1] - p[j][i - 1]) +
                          nu * dt / dx * dx * (un[j][i + 1] - 2 * un[j][i] + un[j][i - 1]) +
                          nu * dt / dy * dy * (un[j + 1][i] - 2 * un[j][i] + un[j - 1][i]);
                v[j][i] = vn[j][i] - vn[j][i] * dt / dx * (vn[j][i] - vn[j][i - 1]) -
                          vn[j][i] * dt / dy * (vn[j][i] - vn[j - 1][i]) -
                          dt / (2 * rho * dx) * (p[j + 1][i] - p[j - 1][i]) +
                          nu * dt / dx * dx * (vn[j][i + 1] - 2 * vn[j][i] + vn[j][i - 1]) +
                          nu * dt / dy * dy * (vn[j + 1][i] - 2 * vn[j][i] + vn[j - 1][i]);
            }
            // printf("%d / u:%f / v:%.2f \n ", j, u[j][10], v[j][10]);
        }

        for (int j = 0; j < ny; j++)
        {
            u[j][0] = 0;
            u[j][ny - 1] = 0;
            v[j][0] = 0;
            v[j][ny - 1] = 0;
        }

        for (int i = 0; i < nx; i++)
        {
            u[0][i] = 0;
            u[ny - 1][i] = 1;
            v[0][i] = 0;
            v[ny - 1][i] = 0;
        }
    }

    //     for (int j = 1; j < ny - 1; j++)
    //     {
    //         for (int i = 1; i < nx - 1; i++)
    //         {
    //             //printf("(%.3f) ", p[j][i]);
    //             printf("(%.2f, %.2f, %.2f) ", u[j][i] * 1000, v[j][i]* 1000, p[j][i]* 1000);
    //         }
    //         printf("\n");
    //     }
    return 0;
}