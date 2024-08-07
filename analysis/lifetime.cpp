#include<iostream>
#include<fstream>
#include<cmath>
#define N_tot 100
#define r_cut 2.8
using namespace std;

float p1_x[60000][40], p1_y[60000][40], p1_z[60000][40];
float p2_x[60000][40], p2_y[60000][40], p2_z[60000][40];
float p3_x[60000][40], p3_y[60000][40], p3_z[60000][40];
float n1_x[60000][30], n1_y[60000][30], n1_z[60000][30];
float n2_x[60000][30], n2_y[60000][30], n2_z[60000][30];
float s_b[60000][30][30];
float lifetime[1000000];
int frq[50000];

struct example
{
    char name[5];
    float x, y, z;
};

example atom_info[60000][200];

int main()
{
    ifstream infile;
    infile.open("../side_chain.pdb");

    if (!infile)
    {
        cerr << "Open infile failure!" << endl;
        return -1;
    }

    // Format the coordinates of P atoms and store them in a structural array
    int i = 1, j, k, l, m, t, count;
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

    ofstream outfile_1, outfile_2, outfile_3, outfile_4, outfile_5;
    outfile_1.open("lifetime.dat", ios::out | ios::ate);
    outfile_2.open("lifetime_dis.dat", ios::out | ios::ate);
    outfile_3.open("lifetime_2.dat", ios::out | ios::ate);
    outfile_4.open("lifetime_dis_2.dat", ios::out | ios::ate);
    outfile_5.open("s_b_series.dat", ios::out | ios::ate);

    for (i = 1; i <= t; i++)
    {
        for (j = 1; j <= 30; j = j + 3)
        {
            p1_x[i][(j + 2) / 3] = atom_info[i][j].x;
            p1_y[i][(j + 2) / 3] = atom_info[i][j].y;
            p1_z[i][(j + 2) / 3] = atom_info[i][j].z;
        }
        for (j = 51; j <= 80; j = j + 3)
        {
            p1_x[i][(j - 18) / 3] = atom_info[i][j].x;
            p1_y[i][(j - 18) / 3] = atom_info[i][j].y;
            p1_z[i][(j - 18) / 3] = atom_info[i][j].z;
        }
    }

    for (i = 1; i <= t; i++)
    {
        for (j = 2; j <= 30; j = j + 3)
        {
            p2_x[i][(j + 1) / 3] = atom_info[i][j].x;
            p2_y[i][(j + 1) / 3] = atom_info[i][j].y;
            p2_z[i][(j + 1) / 3] = atom_info[i][j].z;
        }
        for (j = 52; j <= 80; j = j + 3)
        {
            p2_x[i][(j - 19) / 3] = atom_info[i][j].x;
            p2_y[i][(j - 19) / 3] = atom_info[i][j].y;
            p2_z[i][(j - 19) / 3] = atom_info[i][j].z;
        }
    }

    for (i = 1; i <= t; i++)
    {
        for (j = 3; j <= 30; j = j + 3)
        {
            p3_x[i][(j + 0) / 3] = atom_info[i][j].x;
            p3_y[i][(j + 0) / 3] = atom_info[i][j].y;
            p3_z[i][(j + 0) / 3] = atom_info[i][j].z;
        }
        for (j = 53; j <= 80; j = j + 3)
        {
            p3_x[i][(j - 20) / 3] = atom_info[i][j].x;
            p3_y[i][(j - 20) / 3] = atom_info[i][j].y;
            p3_z[i][(j - 20) / 3] = atom_info[i][j].z;
        }
    }
    
    for (i = 1; i <= t; i++)
    {
        for (j = 31; j <= 50; j = j + 2)
        {
            n1_x[i][(j - 29) / 2] = atom_info[i][j].x;
            n1_y[i][(j - 29) / 2] = atom_info[i][j].y;
            n1_z[i][(j - 29) / 2] = atom_info[i][j].z;
        }
        for (j = 81; j <= 100; j = j + 2)
        {
            n1_x[i][(j - 59) / 2] = atom_info[i][j].x;
            n1_y[i][(j - 59) / 2] = atom_info[i][j].y;
            n1_z[i][(j - 59) / 2] = atom_info[i][j].z;
        }
    }

    for (i = 1; i <= t; i++)
    {
        for (j = 32; j <= 50; j = j + 2)
        {
            n2_x[i][(j - 30) / 2] = atom_info[i][j].x;
            n2_y[i][(j - 30) / 2] = atom_info[i][j].y;
            n2_z[i][(j - 30) / 2] = atom_info[i][j].z;
        }
        for (j = 82; j <= 100; j = j + 2)
        {
            n2_x[i][(j - 60) / 2] = atom_info[i][j].x;
            n2_y[i][(j - 60) / 2] = atom_info[i][j].y;
            n2_z[i][(j - 60) / 2] = atom_info[i][j].z;
        }
    }
    
    float d[5];
    for (i = 1; i <= t; i++)
    {
        for (j = 1; j <= 20; j++)
        {
            for (k = 1; k <= 20; k++)
            {
                d[1] = sqrt((n1_x[i][k] - p1_x[i][j]) * (n1_x[i][k] - p1_x[i][j]) + (n1_y[i][k] - p1_y[i][j]) * (n1_y[i][k] - p1_y[i][j]) + (n1_z[i][k] - p1_z[i][j]) * (n1_z[i][k] - p1_z[i][j]));
                d[2] = sqrt((n1_x[i][k] - p2_x[i][j]) * (n1_x[i][k] - p2_x[i][j]) + (n1_y[i][k] - p2_y[i][j]) * (n1_y[i][k] - p2_y[i][j]) + (n1_z[i][k] - p2_z[i][j]) * (n1_z[i][k] - p2_z[i][j]));
                d[3] = sqrt((n1_x[i][k] - p3_x[i][j]) * (n1_x[i][k] - p3_x[i][j]) + (n1_y[i][k] - p3_y[i][j]) * (n1_y[i][k] - p3_y[i][j]) + (n1_z[i][k] - p3_z[i][j]) * (n1_z[i][k] - p3_z[i][j]));
                d[4] = sqrt((n2_x[i][k] - p1_x[i][j]) * (n2_x[i][k] - p1_x[i][j]) + (n2_y[i][k] - p1_y[i][j]) * (n2_y[i][k] - p1_y[i][j]) + (n2_z[i][k] - p1_z[i][j]) * (n2_z[i][k] - p1_z[i][j]));
                d[5] = sqrt((n2_x[i][k] - p2_x[i][j]) * (n2_x[i][k] - p2_x[i][j]) + (n2_y[i][k] - p2_y[i][j]) * (n2_y[i][k] - p2_y[i][j]) + (n2_z[i][k] - p2_z[i][j]) * (n2_z[i][k] - p2_z[i][j]));
                d[6] = sqrt((n2_x[i][k] - p3_x[i][j]) * (n2_x[i][k] - p3_x[i][j]) + (n2_y[i][k] - p3_y[i][j]) * (n2_y[i][k] - p3_y[i][j]) + (n2_z[i][k] - p3_z[i][j]) * (n2_z[i][k] - p3_z[i][j]));
                count = 0;
                for (l = 1; l <= 6; l++)
                {
                    if (d[l] < r_cut)
                        count++;
                }
                s_b[i][j][k] = count;
            }
        }
    }

    m = 1;
    for (i = 1; i <= 20; i++)
    {
        for (j = 1; j <= 20; j++)
        {
            for (k = 1; k <= t; k = k + count)
            {
                count = 0;
                for (l = k; l <= t; l++)
                {
                    if (s_b[l][i][j] == 0 && s_b[l + 1][i][j] == 0)
                    {
                        count++;
                        continue;
                    }
                    if (s_b[l][i][j] == 0 && s_b[l + 1][i][j] > 0)
                    {
                        count++;
                        break;
                    }
                    if (s_b[l][i][j] > 0 && s_b[l + 1][i][j] > 0)
                    {
                        count++;
                        continue;
                    }
                    if (s_b[l][i][j] > 0 && s_b[l + 1][i][j] == 0)
                    {
                        count++;
                        lifetime[m] = count;
                        if (count > 40000)
                            cout << i << " " << j << endl;
                        m++;
                        break;
                    }
                }   
            }
        }
    }

    for (i = 1; i <= m - 1; i++)
    {
        outfile_1 << lifetime[i] << endl;
    }
    outfile_1.close();

    for (i = 1; i <= 50000; i++)
    {
        frq[i] = 0;
        for (j = 1; j <= m - 1; j++)
        {
            if (lifetime[j] >= i && lifetime[j] < i + 1)
            frq[i]++;
        }
        outfile_2 << i * 0.04 << " " << frq[i] << endl;
    }
    outfile_2.close();

    m = 1;
    for (i = 1; i <= 20; i++)
    {
        for (j = 1; j <= 20; j++)
        {
            for (k = 1; k <= t; k = k + count)
            {
                count = 0;
                for (l = k; l <= t; l++)
                {
                    if (s_b[l][i][j] < 2 && s_b[l + 1][i][j] < 2)
                    {
                        count++;
                        continue;
                    }
                    if (s_b[l][i][j] < 2 && s_b[l + 1][i][j] >= 2)
                    {
                        count++;
                        break;
                    }
                    if (s_b[l][i][j] >= 2 && s_b[l + 1][i][j] >= 2)
                    {
                        count++;
                        continue;
                    }
                    if (s_b[l][i][j] >= 2 && s_b[l + 1][i][j] < 2)
                    {
                        count++;
                        lifetime[m] = count;
                        if (count > 40000)
                            cout << i << " " << j << endl;
                        m++;
                        break;
                    }
                }   
            }
        }
    }

    for (i = 1; i <= m - 1; i++)
    {
        outfile_3 << lifetime[i] << endl;
    }
    outfile_3.close();

    for (i = 1; i <= 50000; i++)
    {
        frq[i] = 0;
        for (j = 1; j <= m - 1; j++)
        {
            if (lifetime[j] >= i && lifetime[j] < i + 1)
            frq[i]++;
        }
        outfile_4 << i * 0.04 << " " << frq[i] << endl;
    }
    outfile_4.close();

    for (i = 1; i <= t; i++)
    {
        outfile_5 << s_b[i][10][2] << " " << s_b[i][11][2] << " " << s_b[i][14][2] << endl;
    }
    outfile_5.close();

    return 0;
}
