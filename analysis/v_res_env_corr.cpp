#include<iostream>
#include<fstream>
#include<cmath>
#define N_tot 80
#define r_cut 4
#define cut 10
using namespace std;

float p1_x[60000][30], p1_y[60000][30], p1_z[60000][30];
float p2_x[60000][30], p2_y[60000][30], p2_z[60000][30];
float n1_x[60000][30], n1_y[60000][30], n1_z[60000][30];
float n2_x[60000][30], n2_y[60000][30], n2_z[60000][30];
float v_p1[60000][30], v_p2[60000][30], v_n1[60000][30], v_n2[60000][30];
float p_center_n_env[60000][30], n_center_p_env[60000][30];
float P_p_center_n_env[100], P_n_center_p_env[100];

struct example
{
    char name[5];
    float x, y, z;
};

example atom_info[60000][200];

int main()
{
    ifstream infile;
    infile.open("../res_env.pdb");

    if (!infile)
    {
        cerr << "Open infile failure!" << endl;
        return -1;
    }

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

    ofstream outfile_1, outfile_2, outfile_3, outfile_4;
    outfile_1.open("p_p_center_n_env.dat", ios::out | ios::ate);
    outfile_2.open("p_n_center_p_env.dat", ios::out | ios::ate);
    outfile_3.open("p_env_series.dat", ios::out | ios::ate);
    outfile_4.open("n_env_series.dat", ios::out | ios::ate);

    for (i = 1; i <= t; i++)
    {
        for (j = 1; j <= 20; j = j + 2)
        {
            p1_x[i][(j + 1) / 2] = atom_info[i][j].x;
            p1_y[i][(j + 1) / 2] = atom_info[i][j].y;
            p1_z[i][(j + 1) / 2] = atom_info[i][j].z;
        }
        for (j = 41; j <= 60; j = j + 2)
        {
            p1_x[i][(j - 19) / 2] = atom_info[i][j].x;
            p1_y[i][(j - 19) / 2] = atom_info[i][j].y;
            p1_z[i][(j - 19) / 2] = atom_info[i][j].z;
        }
    }

    for (i = 1; i <= t; i++)
    {
        for (j = 2; j <= 20; j = j + 2)
        {
            p2_x[i][(j + 0) / 2] = atom_info[i][j].x;
            p2_y[i][(j + 0) / 2] = atom_info[i][j].y;
            p2_z[i][(j + 0) / 2] = atom_info[i][j].z;
        }
        for (j = 42; j <= 60; j = j + 2)
        {
            p2_x[i][(j - 20) / 2] = atom_info[i][j].x;
            p2_y[i][(j - 20) / 2] = atom_info[i][j].y;
            p2_z[i][(j - 20) / 2] = atom_info[i][j].z;
        }
    }
    
    for (i = 1; i <= t; i++)
    {
        for (j = 21; j <= 40; j = j + 2)
        {
            n1_x[i][(j - 19) / 2] = atom_info[i][j].x;
            n1_y[i][(j - 19) / 2] = atom_info[i][j].y;
            n1_z[i][(j - 19) / 2] = atom_info[i][j].z;
        }
        for (j = 61; j <= 80; j = j + 2)
        {
            n1_x[i][(j - 39) / 2] = atom_info[i][j].x;
            n1_y[i][(j - 39) / 2] = atom_info[i][j].y;
            n1_z[i][(j - 39) / 2] = atom_info[i][j].z;
        }
    }

    for (i = 1; i <= t; i++)
    {
        for (j = 22; j <= 40; j = j + 2)
        {
            n2_x[i][(j - 20) / 2] = atom_info[i][j].x;
            n2_y[i][(j - 20) / 2] = atom_info[i][j].y;
            n2_z[i][(j - 20) / 2] = atom_info[i][j].z;
        }
        for (j = 62; j <= 80; j = j + 2)
        {
            n2_x[i][(j - 40) / 2] = atom_info[i][j].x;
            n2_y[i][(j - 40) / 2] = atom_info[i][j].y;
            n2_z[i][(j - 40) / 2] = atom_info[i][j].z;
        }
    }

    for (i = 1; i <= t - 1; i++)
    {
        for (j = 1; j <= 20; j++)
        {
            v_p1[i][j] = sqrt((p1_x[i + 1][j] - p1_x[i][j]) * (p1_x[i + 1][j] - p1_x[i][j]) + (p1_y[i + 1][j] - p1_y[i][j]) * (p1_y[i + 1][j] - p1_y[i][j]) + (p1_z[i + 1][j] - p1_z[i][j]) * (p1_z[i + 1][j] - p1_z[i][j]));
        }
    }

    for (i = 1; i <= t - 1; i++)
    {
        for (j = 1; j <= 20; j++)
        {
            v_p2[i][j] = sqrt((p2_x[i + 1][j] - p2_x[i][j]) * (p2_x[i + 1][j] - p2_x[i][j]) + (p2_y[i + 1][j] - p2_y[i][j]) * (p2_y[i + 1][j] - p2_y[i][j]) + (p2_z[i + 1][j] - p2_z[i][j]) * (p2_z[i + 1][j] - p2_z[i][j]));
        }
    }

    for (i = 1; i <= t - 1; i++)
    {
        for (j = 1; j <= 20; j++)
        {
            v_n1[i][j] = sqrt((n1_x[i + 1][j] - n1_x[i][j]) * (n1_x[i + 1][j] - n1_x[i][j]) + (n1_y[i + 1][j] - n1_y[i][j]) * (n1_y[i + 1][j] - n1_y[i][j]) + (n1_z[i + 1][j] - n1_z[i][j]) * (n1_z[i + 1][j] - n1_z[i][j]));
        }
    }

    for (i = 1; i <= t - 1; i++)
    {
        for (j = 1; j <= 20; j++)
        {
            v_n2[i][j] = sqrt((n2_x[i + 1][j] - n2_x[i][j]) * (n2_x[i + 1][j] - n2_x[i][j]) + (n2_y[i + 1][j] - n2_y[i][j]) * (n2_y[i + 1][j] - n2_y[i][j]) + (n2_z[i + 1][j] - n2_z[i][j]) * (n2_z[i + 1][j] - n2_z[i][j]));
        }
    }
    
    int count;
    float d_p1n1, d_p1n2, d_p2n1, d_p2n2;
    for (i = 1; i <= t; i++)
    {
        for (j = 1; j <= 20; j++)
        {
            count = 0;
            for (k = 1; k <= 20; k++)
            {
                d_p1n1 = sqrt((n1_x[i][k] - p1_x[i][j]) * (n1_x[i][k] - p1_x[i][j]) + (n1_y[i][k] - p1_y[i][j]) * (n1_y[i][k] - p1_y[i][j]) + (n1_z[i][k] - p1_z[i][j]) * (n1_z[i][k] - p1_z[i][j]));
                d_p1n2 = sqrt((n2_x[i][k] - p1_x[i][j]) * (n2_x[i][k] - p1_x[i][j]) + (n2_y[i][k] - p1_y[i][j]) * (n2_y[i][k] - p1_y[i][j]) + (n2_z[i][k] - p1_z[i][j]) * (n2_z[i][k] - p1_z[i][j]));
                d_p2n1 = sqrt((n1_x[i][k] - p2_x[i][j]) * (n1_x[i][k] - p2_x[i][j]) + (n1_y[i][k] - p2_y[i][j]) * (n1_y[i][k] - p2_y[i][j]) + (n1_z[i][k] - p2_z[i][j]) * (n1_z[i][k] - p2_z[i][j]));
                d_p2n2 = sqrt((n2_x[i][k] - p2_x[i][j]) * (n2_x[i][k] - p2_x[i][j]) + (n2_y[i][k] - p2_y[i][j]) * (n2_y[i][k] - p2_y[i][j]) + (n2_z[i][k] - p2_z[i][j]) * (n2_z[i][k] - p2_z[i][j]));
                if (d_p1n1 < r_cut || d_p1n2 < r_cut || d_p2n1 < r_cut || d_p2n2 < r_cut)
                    count++;
            }
            p_center_n_env[i][j] = count;
        }
    }
    
    float d_n1p1, d_n1p2, d_n2p1, d_n2p2;
    for (i = 1; i <= t; i++)
    {
        for (j = 1; j <= 20; j++)
        {
            count = 0;
            for (k = 1; k <= 20; k++)
            {
                d_n1p1 = sqrt((p1_x[i][k] - n1_x[i][j]) * (p1_x[i][k] - n1_x[i][j]) + (p1_y[i][k] - n1_y[i][j]) * (p1_y[i][k] - n1_y[i][j]) + (p1_z[i][k] - n1_z[i][j]) * (p1_z[i][k] - n1_z[i][j]));
                d_n1p2 = sqrt((p2_x[i][k] - n1_x[i][j]) * (p2_x[i][k] - n1_x[i][j]) + (p2_y[i][k] - n1_y[i][j]) * (p2_y[i][k] - n1_y[i][j]) + (p2_z[i][k] - n1_z[i][j]) * (p2_z[i][k] - n1_z[i][j]));
                d_n2p1 = sqrt((p1_x[i][k] - n2_x[i][j]) * (p1_x[i][k] - n2_x[i][j]) + (p1_y[i][k] - n2_y[i][j]) * (p1_y[i][k] - n2_y[i][j]) + (p1_z[i][k] - n2_z[i][j]) * (p1_z[i][k] - n2_z[i][j]));
                d_n2p2 = sqrt((p2_x[i][k] - n2_x[i][j]) * (p2_x[i][k] - n2_x[i][j]) + (p2_y[i][k] - n2_y[i][j]) * (p2_y[i][k] - n2_y[i][j]) + (p2_z[i][k] - n2_z[i][j]) * (p2_z[i][k] - n2_z[i][j]));
                if (d_n1p1 < r_cut || d_n1p2 < r_cut || d_n2p1 < r_cut || d_n2p2 < r_cut)
                    count++;
            }
            n_center_p_env[i][j] = count;
        }
    }

    for (i = 0; i <= 4; i++)
    {
        count = 0;
        for (j = 10001; j <= t; j++)
        {
            for (k = 1; k <= 20; k++)
            {
                if (p_center_n_env[j][k] >= (float)i && p_center_n_env[j][k] < (float)(i + 1))
                    count++;
            }
        }
        P_p_center_n_env[i] = (float)count / (t - 10000) / 20;
        outfile_1 << (float)i << " " << P_p_center_n_env[i] << endl;
    }
    outfile_1.close();

    for (i = 0; i <= 4; i++)
    {
        count = 0;
        for (j = 10001; j <= t; j++)
        {
            for (k = 1; k <= 20; k++)
            {
                if (n_center_p_env[j][k] >= (float)i && n_center_p_env[j][k] < (float)(i + 1))
                    count++;
            }
        }
        P_n_center_p_env[i] = (float)count / (t - 10000) / 20;
        outfile_2 << (float)i << " " << P_n_center_p_env[i] << endl;
    }
    outfile_2.close();

    for (i = 1; i <= t; i = i + 10)
    {
        for (j = 1; j <= 20; j++)
        {
            outfile_3 << i << " " << j << " " << p_center_n_env[i][j] << endl;
            outfile_4 << i << " " << j << " " << n_center_p_env[i][j] << endl;
        }
        outfile_3 << endl;
        outfile_4 << endl;
    }
    outfile_3.close();
    outfile_4.close();

    return 0;
}
