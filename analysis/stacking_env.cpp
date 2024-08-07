#include<iostream>
#include<fstream>
#include<cmath>
#define N_ele_env 80
#define N_CZ 20
#define interval 10
#define stacking_dist 4.2
#define cut 10
using namespace std;

float p_center_x[60000][30], p_center_y[60000][30], p_center_z[60000][30];
float n_center_x[60000][30], n_center_y[60000][30], n_center_z[60000][30];
float CZ_ARG_env[60000][30], CZ_ASP_env[60000][30];
float CZ_pair_env[60000][30][30], CZ_pair_dist[60000][30][30];
float CZ_pair_dist_interval[60000][30][30], CZ_pair_env_interval[60000][30][30];
float p_dist_env[100][100];
float CZ_pair_dist_avg[30][30], CZ_pair_env_avg[30][30];

struct example
{
    char name[5];
    float x, y, z;
};

example CZ_info[60000][30];
example ele_env_info[60000][100];

int main()
{
    ifstream infile;
    infile.open("CZ.pdb");

    if (!infile)
    {
        cerr << "Open infile failure!" << endl;
        return -1;
    }

    int i = 1, j, k, l, m, t;
    float aa;
    while (!infile.eof())
    {
        aa = (i + 0.0) / (N_CZ + 0.0); 
        k = fmod(i, N_CZ);
        if (k == 0) 
        {
            k = N_CZ;
            j = floor(aa);
        }
        else 
            j = floor(aa) + 1;
        infile >> CZ_info[j][k].name >> CZ_info[j][k].x >> CZ_info[j][k].y >> CZ_info[j][k].z;
        i++;
    }
    t = j - 1;
    cout << t << endl;
    infile.close();

    infile.open("../res_env.pdb");

    if (!infile)
    {
        cerr << "Open infile failure!" << endl;
        return -1;
    }

    i = 1;
    while (!infile.eof())
    {
        aa = (i + 0.0) / (N_ele_env + 0.0); 
        k = fmod(i, N_ele_env);
        if (k == 0) 
        {
            k = N_ele_env;
            j = floor(aa);
        }
        else 
            j = floor(aa) + 1;
        infile >> ele_env_info[j][k].name >> ele_env_info[j][k].x >> ele_env_info[j][k].y >> ele_env_info[j][k].z;
        i++;
    }
    t = j - 1;
    cout << t << endl;
    infile.close();

    ofstream outfile_1, outfile_2, outfile_3, outfile_4;
    outfile_1.open("stacking_env.dat", ios::out | ios::ate);
    outfile_2.open("p_dist_env.dat", ios::out | ios::ate);
    //outfile_3.open("CZ_pair_dist.dat", ios::out | ios::ate);
    //outfile_4.open("CZ_pair_env.dat", ios::out | ios::ate);

    for (i = 1; i <= t; i++)
    {
        for (j = 1; j <= 20; j = j + 2)
        {
            p_center_x[i][(j + 1) / 2] = (ele_env_info[i][j].x + ele_env_info[i][j + 1].x) / 2;
            p_center_y[i][(j + 1) / 2] = (ele_env_info[i][j].y + ele_env_info[i][j + 1].y) / 2;
            p_center_z[i][(j + 1) / 2] = (ele_env_info[i][j].z + ele_env_info[i][j + 1].z) / 2;
        }
        for (j = 41; j <= 60; j = j + 2)
        {
            p_center_x[i][(j - 19) / 2] = (ele_env_info[i][j].x + ele_env_info[i][j + 1].x) / 2;
            p_center_y[i][(j - 19) / 2] = (ele_env_info[i][j].y + ele_env_info[i][j + 1].y) / 2;
            p_center_z[i][(j - 19) / 2] = (ele_env_info[i][j].z + ele_env_info[i][j + 1].z) / 2;
        }
    }
    
    for (i = 1; i <= t; i++)
    {
        for (j = 21; j <= 40; j = j + 2)
        {
            n_center_x[i][(j - 19) / 2] = (ele_env_info[i][j].x + ele_env_info[i][j + 1].x) / 2;
            n_center_y[i][(j - 19) / 2] = (ele_env_info[i][j].y + ele_env_info[i][j + 1].y) / 2;
            n_center_z[i][(j - 19) / 2] = (ele_env_info[i][j].z + ele_env_info[i][j + 1].z) / 2;
        }
        for (j = 61; j <= 80; j = j + 2)
        {
            n_center_x[i][(j - 39) / 2] = (ele_env_info[i][j].x + ele_env_info[i][j + 1].x) / 2;
            n_center_y[i][(j - 39) / 2] = (ele_env_info[i][j].y + ele_env_info[i][j + 1].y) / 2;
            n_center_z[i][(j - 39) / 2] = (ele_env_info[i][j].z + ele_env_info[i][j + 1].z) / 2;
        }
    }

    for (i = 1; i <= t; i++)
    {
        for (j = 1; j <= 20; j++)
        {
            CZ_ARG_env[i][j] = 0;
            for (k = 1; k <= 20; k++)
            {
                CZ_ARG_env[i][j] += 1 / sqrt((p_center_x[i][k] - CZ_info[i][j].x) * (p_center_x[i][k] - CZ_info[i][j].x) + (p_center_y[i][k] - CZ_info[i][j].y) * (p_center_y[i][k] - CZ_info[i][j].y) + (p_center_z[i][k] - CZ_info[i][j].z) * (p_center_z[i][k] - CZ_info[i][j].z));
            }
        }
    }

    for (i = 1; i <= t; i++)
    {
        for (j = 1; j <= 20; j++)
        {
            CZ_ASP_env[i][j] = 0;
            for (k = 1; k <= 20; k++)
            {
                CZ_ASP_env[i][j] += -1 / sqrt((n_center_x[i][k] - CZ_info[i][j].x) * (n_center_x[i][k] - CZ_info[i][j].x) + (n_center_y[i][k] - CZ_info[i][j].y) * (n_center_y[i][k] - CZ_info[i][j].y) + (n_center_z[i][k] - CZ_info[i][j].z) * (n_center_z[i][k] - CZ_info[i][j].z));
            }
        }
    }

    for (i = 1; i <= t; i++)
    {
        for (j = 1; j <= 20; j++)
        {
            for (k = 1; k <= 20; k++)
            {
                CZ_pair_dist[i][j][k] = sqrt((CZ_info[i][j].x - CZ_info[i][k].x) * (CZ_info[i][j].x - CZ_info[i][k].x) + (CZ_info[i][j].y - CZ_info[i][k].y) * (CZ_info[i][j].y - CZ_info[i][k].y) + (CZ_info[i][j].z - CZ_info[i][k].z) * (CZ_info[i][j].z - CZ_info[i][k].z));
                CZ_pair_env[i][j][k] = CZ_ASP_env[i][j] + CZ_ASP_env[i][k];
            }
        }
    }

    float CZ_pair_dist_tmp, CZ_pair_env_tmp;
    for (i = 1; i <= t; i = i + interval)
    {
        for (j = 1; j <= 20; j++)
        {
            for (k = 1; k <= 20; k++)
            {
                if (j > k)
                {
                    CZ_pair_dist_tmp = 0, CZ_pair_env_tmp = 0;
                    for (l = 1; l <= interval; l++)
                    {
                        CZ_pair_dist_tmp += CZ_pair_dist[i + l - 1][j][k];
                        CZ_pair_env_tmp += CZ_pair_env[i + l - 1][j][k];
                    }
                    CZ_pair_dist_interval[(i + interval - 1) / interval][j][k] = CZ_pair_dist_tmp / interval;
                    CZ_pair_env_interval[(i + interval - 1) / interval][j][k] = CZ_pair_env_tmp / interval;
                    outfile_1 << CZ_pair_dist_interval[(i + interval - 1) / interval][j][k] << " " << CZ_pair_env_interval[(i + interval - 1) / interval][j][k] << endl;
                }
            }
        }
    }
    outfile_1.close();

    int count;
    for (i = 0; i <= 60; i++)
    {
        for (j = -5 * cut; j <= -1 * cut; j++)
        {
            count = 0;
            for (k = 1; k <= t / interval - 1; k++)
            {
                for (l = 1; l <= 20; l++)
                {
                    for (m = 1; m <= 20; m++)
                    {
                        if (l > m)
                        {
                            if (CZ_pair_dist_interval[k][l][m] >= i && CZ_pair_dist_interval[k][l][m] < i + 1)
                            {
                                if (CZ_pair_env_interval[k][l][m] >= (float)j / cut && CZ_pair_env_interval[k][l][m] < (float)(j + 1) / cut)
                                count++;
                            }
                        }
                    }
                }
            }
            p_dist_env[i][j] = (float)count / (t / interval) / 190;
            outfile_2 << i << " " << (float)j / cut << " " << p_dist_env[i][j] << endl;
        }
        outfile_2 << endl;
    }
    outfile_2.close();

/*    for (i = 1; i <= 20; i++)
    {
        for (j = 1; j <= 20; j++)
        {
            for (k = 1; k <= t / interval - 1; k++)
            {
                if (CZ_pair_dist_interval[k][i][j] < stacking_dist)
                    CZ_pair_dist_avg[i][j]++;
                CZ_pair_env_avg[i][j] += CZ_pair_env_interval[k][i][j];
            }
            CZ_pair_dist_avg[i][j] /= t / interval;
            CZ_pair_env_avg[i][j] /= t / interval;
            outfile_3 << i << " " << j << " " << CZ_pair_dist_avg[i][j] << endl;
            outfile_4 << i << " " << j << " " << CZ_pair_env_avg[i][j] << endl;
        }
        outfile_3 << endl;
        outfile_4 << endl;
    }
    outfile_3.close();
    outfile_4.close();*/

    return 0;
}
