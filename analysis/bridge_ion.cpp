#include<iostream>
#include<fstream>
#include<cmath>
#define N_tot 440
#define N_PEP 80
#define N_NA 180
#define N_CL 180
#define R_NA 3
#define R_CL 4
using namespace std;

float sc1_x[60000][50], sc1_y[60000][50], sc1_z[60000][50], sc2_x[60000][50], sc2_y[60000][50], sc2_z[60000][50];
float density_NA[50][50][60000], density_NA_avg[50][50], density_NA_std[50][50];
float density_CL[50][50][60000], density_CL_avg[50][50], density_CL_std[50][50];

struct example
{
    char name[5];
    float x, y, z;
};

example atom_info[60000][500];

int main()
{
    ifstream infile;
    infile.open("../ion_env.pdb");

    if (!infile)
    {
        cerr << "Open infile failure!" << endl;
        return -1;
    }

    // 将原子坐标格式化存入结构数组
    int i = 1, j, k, l, t;
    float aa;
    while (!infile.eof())
    {
        aa = (i + 0.0) / (N_tot + 0.0); 
        k = fmod(i, N_tot);
        if (k == 0) 
        {
            k = N_tot;
            j = floor(aa);
        }
        else 
            j = floor(aa) + 1;
        infile >> atom_info[j][k].name >> atom_info[j][k].x >> atom_info[j][k].y >> atom_info[j][k].z;
        i++;
    }
    t = j - 1;
    cout << t << endl;
    infile.close();

    ofstream outfile_1, outfile_2;
    outfile_1.open("bridge_NA.dat", ios::out | ios::ate);
    outfile_2.open("bridge_CL.dat", ios::out | ios::ate);

    // 计算侧链电荷中心1随时间的变化
    for (i = 1; i <= t; i++)
    {
        for (j = 1; j <= N_PEP; j = j + 2)
        {
            sc1_x[i][(j + 1) / 2] = atom_info[i][j].x;
            sc1_y[i][(j + 1) / 2] = atom_info[i][j].y;
            sc1_z[i][(j + 1) / 2] = atom_info[i][j].z;
        }
    }

    // 计算侧链电荷中心2随时间的变化
    for (i = 1; i <= t; i++)
    {
        for (j = 2; j <= N_PEP; j = j + 2)
        {
            sc2_x[i][(j + 0) / 2] = atom_info[i][j].x;
            sc2_y[i][(j + 0) / 2] = atom_info[i][j].y;
            sc2_z[i][(j + 0) / 2] = atom_info[i][j].z;
        }
    }
    
    // 计算两两电荷中心共享NA离子个数随时间的变化
    float d_i1, d_i2, d_j1, d_j2;
    for (i = 1; i <= 40; i++)
    {
        for (j = 1; j <= 40; j++)
        {
            for (k = 1; k <= t; k++)
            {
                density_NA[i][j][k] = 0;
                for (l = N_PEP + 1; l <= N_PEP + N_NA; l++)
                {
                    d_i1 = sqrt((atom_info[k][l].x - sc1_x[k][i]) * (atom_info[k][l].x - sc1_x[k][i]) + (atom_info[k][l].y - sc1_y[k][i]) * (atom_info[k][l].y - sc1_y[k][i]) + (atom_info[k][l].z - sc1_z[k][i]) * (atom_info[k][l].z - sc1_z[k][i]));
                    d_i2 = sqrt((atom_info[k][l].x - sc2_x[k][i]) * (atom_info[k][l].x - sc2_x[k][i]) + (atom_info[k][l].y - sc2_y[k][i]) * (atom_info[k][l].y - sc2_y[k][i]) + (atom_info[k][l].z - sc2_z[k][i]) * (atom_info[k][l].z - sc2_z[k][i]));
                    d_j1 = sqrt((atom_info[k][l].x - sc1_x[k][j]) * (atom_info[k][l].x - sc1_x[k][j]) + (atom_info[k][l].y - sc1_y[k][j]) * (atom_info[k][l].y - sc1_y[k][j]) + (atom_info[k][l].z - sc1_z[k][j]) * (atom_info[k][l].z - sc1_z[k][j]));
                    d_j2 = sqrt((atom_info[k][l].x - sc2_x[k][j]) * (atom_info[k][l].x - sc2_x[k][j]) + (atom_info[k][l].y - sc2_y[k][j]) * (atom_info[k][l].y - sc2_y[k][j]) + (atom_info[k][l].z - sc2_z[k][j]) * (atom_info[k][l].z - sc2_z[k][j]));
                    if (d_i1 < R_NA || d_i2 < R_NA)
                    {
                        if (d_j1 < R_NA || d_j2 < R_NA)
                            density_NA[i][j][k]++;
                    }
                }
            }
        }
    }

    // 计算两两电荷中心共享NA离子个数的平均值
    for (i = 1; i <= 40; i++)
    {
        for (j = 1; j <= 40; j++)
        {
            density_NA_avg[i][j] = 0;
            for (k = 1; k <= t; k++)
            {
                density_NA_avg[i][j] += density_NA[i][j][k];
            }
            density_NA_avg[i][j] = density_NA_avg[i][j] / t;
        }
    }

    // 计算两两电荷中心共享NA离子个数的标准差
    for (i = 1; i <= 40; i++)
    {
        for (j = 1; j <= 40; j++)
        {        
            density_NA_std[i][j] = 0;
            for (k = 1; k <= t; k++)
            {
                density_NA_std[i][j] += (density_NA[i][j][k] - density_NA_avg[i][j]) * (density_NA[i][j][k] - density_NA_avg[i][j]);
            }
            density_NA_std[i][j] = sqrt(density_NA_std[i][j] / t);
            outfile_1 << i << " " << j << " " << density_NA_avg[i][j] << " " << density_NA_std[i][j] << endl;
        }
        outfile_1 << endl;
    }
    outfile_1.close();

    // 计算两两电荷中心共享CL离子个数随时间的变化
    for (i = 1; i <= 40; i++)
    {
        for (j = 1; j <= 40; j++)
        {
            for (k = 1; k <= t; k++)
            {
                density_CL[i][j][k] = 0;
                for (l = N_PEP + N_NA + 1; l <= N_PEP + N_NA + N_CL; l++)
                {
                    d_i1 = sqrt((atom_info[k][l].x - sc1_x[k][i]) * (atom_info[k][l].x - sc1_x[k][i]) + (atom_info[k][l].y - sc1_y[k][i]) * (atom_info[k][l].y - sc1_y[k][i]) + (atom_info[k][l].z - sc1_z[k][i]) * (atom_info[k][l].z - sc1_z[k][i]));
                    d_i2 = sqrt((atom_info[k][l].x - sc2_x[k][i]) * (atom_info[k][l].x - sc2_x[k][i]) + (atom_info[k][l].y - sc2_y[k][i]) * (atom_info[k][l].y - sc2_y[k][i]) + (atom_info[k][l].z - sc2_z[k][i]) * (atom_info[k][l].z - sc2_z[k][i]));
                    d_j1 = sqrt((atom_info[k][l].x - sc1_x[k][j]) * (atom_info[k][l].x - sc1_x[k][j]) + (atom_info[k][l].y - sc1_y[k][j]) * (atom_info[k][l].y - sc1_y[k][j]) + (atom_info[k][l].z - sc1_z[k][j]) * (atom_info[k][l].z - sc1_z[k][j]));
                    d_j2 = sqrt((atom_info[k][l].x - sc2_x[k][j]) * (atom_info[k][l].x - sc2_x[k][j]) + (atom_info[k][l].y - sc2_y[k][j]) * (atom_info[k][l].y - sc2_y[k][j]) + (atom_info[k][l].z - sc2_z[k][j]) * (atom_info[k][l].z - sc2_z[k][j]));
                    if (d_i1 < R_CL || d_i2 < R_CL)
                    {
                        if (d_j1 < R_CL || d_j2 < R_CL)
                            density_CL[i][j][k]++;
                    }
                }
            }
        }
    }

    // 计算两两电荷中心共享CL离子个数的平均值
    for (i = 1; i <= 40; i++)
    {
        for (j = 1; j <= 40; j++)
        {
            density_CL_avg[i][j] = 0;
            for (k = 1; k <= t; k++)
            {
                density_CL_avg[i][j] += density_CL[i][j][k];
            }
            density_CL_avg[i][j] = density_CL_avg[i][j] / t;
        }
    }

    // 计算两两电荷中心共享CL离子个数的标准差
    for (i = 1; i <= 40; i++)
    {
        for (j = 1; j <= 40; j++)
        {        
            density_CL_std[i][j] = 0;
            for (k = 1; k <= t; k++)
            {
                density_CL_std[i][j] += (density_CL[i][j][k] - density_CL_avg[i][j]) * (density_CL[i][j][k] - density_CL_avg[i][j]);
            }
            density_CL_std[i][j] = sqrt(density_CL_std[i][j] / t);
            outfile_2 << i << " " << j << " " << density_CL_avg[i][j] << " " << density_CL_std[i][j] << endl;
        }
        outfile_2 << endl;
    }
    outfile_2.close();

    return 0;
}
