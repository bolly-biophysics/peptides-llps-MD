#include<iostream>
#include<fstream>
#include<cmath>
#define N_tot 480
#define N_PEP 120
#define N_NA 180
#define N_CL 180
#define r_NA 3
#define r_CL 4
using namespace std;

float p1_x[60000][30], p1_y[60000][30], p1_z[60000][30];
float p2_x[60000][30], p2_y[60000][30], p2_z[60000][30];
float po_x[60000][30], po_y[60000][30], po_z[60000][30];
float n1_x[60000][30], n1_y[60000][30], n1_z[60000][30];
float n2_x[60000][30], n2_y[60000][30], n2_z[60000][30];
float no_x[60000][30], no_y[60000][30], no_z[60000][30];
float v_NA[60000][200], v_CL[60000][200];
float NA_n_env[60000][200], CL_p_env[60000][200];
float P_v_env_NA[200][100], P_v_env_CL[200][100];
float P_v_NA[200], P_v_CL[200];

struct example
{
    char name[5];
    float x, y, z;
};

example atom_info[60000][600];

int main()
{
    ifstream infile;
    infile.open("../ion_env_all.pdb");

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

    ofstream outfile_1, outfile_2, outfile_3, outfile_4;
    outfile_1.open("p_v_env_NA.dat", ios::out | ios::ate);
    outfile_2.open("p_v_env_CL.dat", ios::out | ios::ate);
    outfile_3.open("v_NA_dis.dat", ios::out | ios::ate);
    outfile_4.open("v_CL_dis.dat", ios::out | ios::ate);

    // 计算NA离子速度随时间的变化
    for (i = 1; i <= t - 1; i++)
    {
        for (j = N_PEP + 1; j <= N_PEP + N_NA; j++)
        {
            v_NA[i][j - N_PEP] = sqrt((atom_info[i + 1][j].x - atom_info[i][j].x) * (atom_info[i + 1][j].x - atom_info[i][j].x) + (atom_info[i + 1][j].y - atom_info[i][j].y) * (atom_info[i + 1][j].y - atom_info[i][j].y) + (atom_info[i + 1][j].z - atom_info[i][j].z) * (atom_info[i + 1][j].z - atom_info[i][j].z));
        }
    }

    // 计算CL离子速度随时间的变化
    for (i = 1; i <= t - 1; i++)
    {
        for (j = N_PEP + N_NA + 1; j <= N_PEP + N_NA + N_CL; j++)
        {
            v_CL[i][j - N_PEP - N_NA] = sqrt((atom_info[i + 1][j].x - atom_info[i][j].x) * (atom_info[i + 1][j].x - atom_info[i][j].x) + (atom_info[i + 1][j].y - atom_info[i][j].y) * (atom_info[i + 1][j].y - atom_info[i][j].y) + (atom_info[i + 1][j].z - atom_info[i][j].z) * (atom_info[i + 1][j].z - atom_info[i][j].z));
        }
    }

    // 计算ARG电正中心1随时间的变化
    for (i = 1; i <= t; i++)
    {
        for (j = 1; j <= 30; j = j + 3)
        {
            p1_x[i][(j + 2) / 3] = atom_info[i][j].x;
            p1_y[i][(j + 2) / 3] = atom_info[i][j].y;
            p1_z[i][(j + 2) / 3] = atom_info[i][j].z;
        }
        for (j = 61; j <= 90; j = j + 3)
        {
            p1_x[i][(j - 28) / 3] = atom_info[i][j].x;
            p1_y[i][(j - 28) / 3] = atom_info[i][j].y;
            p1_z[i][(j - 28) / 3] = atom_info[i][j].z;
        }
    }

    // 计算ARG电正中心2随时间的变化
    for (i = 1; i <= t; i++)
    {
        for (j = 2; j <= 30; j = j + 3)
        {
            p2_x[i][(j + 1) / 3] = atom_info[i][j].x;
            p2_y[i][(j + 1) / 3] = atom_info[i][j].y;
            p2_z[i][(j + 1) / 3] = atom_info[i][j].z;
        }
        for (j = 62; j <= 90; j = j + 3)
        {
            p2_x[i][(j - 29) / 3] = atom_info[i][j].x;
            p2_y[i][(j - 29) / 3] = atom_info[i][j].y;
            p2_z[i][(j - 29) / 3] = atom_info[i][j].z;
        }
    }

    // 计算ARG中O原子随时间的变化
    for (i = 1; i <= t; i++)
    {
        for (j = 3; j <= 30; j = j + 3)
        {
            po_x[i][(j + 0) / 3] = atom_info[i][j].x;
            po_y[i][(j + 0) / 3] = atom_info[i][j].y;
            po_z[i][(j + 0) / 3] = atom_info[i][j].z;
        }
        for (j = 63; j <= 90; j = j + 3)
        {
            po_x[i][(j - 30) / 3] = atom_info[i][j].x;
            po_y[i][(j - 30) / 3] = atom_info[i][j].y;
            po_z[i][(j - 30) / 3] = atom_info[i][j].z;
        }
    }
    
    // 计算ASP电负中心1随时间的变化
    for (i = 1; i <= t; i++)
    {
        for (j = 31; j <= 60; j = j + 3)
        {
            n1_x[i][(j - 28) / 3] = atom_info[i][j].x;
            n1_y[i][(j - 28) / 3] = atom_info[i][j].y;
            n1_z[i][(j - 28) / 3] = atom_info[i][j].z;
        }
        for (j = 91; j <= 120; j = j + 3)
        {
            n1_x[i][(j - 58) / 3] = atom_info[i][j].x;
            n1_y[i][(j - 58) / 3] = atom_info[i][j].y;
            n1_z[i][(j - 58) / 3] = atom_info[i][j].z;
        }
    }

    // 计算ASP电负中心2随时间的变化
    for (i = 1; i <= t; i++)
    {
        for (j = 32; j <= 60; j = j + 3)
        {
            n2_x[i][(j - 29) / 3] = atom_info[i][j].x;
            n2_y[i][(j - 29) / 3] = atom_info[i][j].y;
            n2_z[i][(j - 29) / 3] = atom_info[i][j].z;
        }
        for (j = 92; j <= 120; j = j + 3)
        {
            n2_x[i][(j - 59) / 3] = atom_info[i][j].x;
            n2_y[i][(j - 59) / 3] = atom_info[i][j].y;
            n2_z[i][(j - 59) / 3] = atom_info[i][j].z;
        }
    }

    // 计算ASP中O原子随时间的变化
    for (i = 1; i <= t; i++)
    {
        for (j = 33; j <= 60; j = j + 3)
        {
            no_x[i][(j - 30) / 3] = atom_info[i][j].x;
            no_y[i][(j - 30) / 3] = atom_info[i][j].y;
            no_z[i][(j - 30) / 3] = atom_info[i][j].z;
        }
        for (j = 93; j <= 120; j = j + 3)
        {
            no_x[i][(j - 60) / 3] = atom_info[i][j].x;
            no_y[i][(j - 60) / 3] = atom_info[i][j].y;
            no_z[i][(j - 60) / 3] = atom_info[i][j].z;
        }
    }
    
    // 计算NA离子周围电负环境随时间的变化
    int count;
    float dist_n1, dist_n2, dist_no, dist_po;
    for (i = 1; i <= t; i++)
    {
        for (j = N_PEP + 1; j <= N_PEP + N_NA; j++)
        {
            count = 0;
            for (k = 1; k <= 20; k++)
            {
                dist_n1 = sqrt((n1_x[i][k] - atom_info[i][j].x) * (n1_x[i][k] - atom_info[i][j].x) + (n1_y[i][k] - atom_info[i][j].y) * (n1_y[i][k] - atom_info[i][j].y) + (n1_z[i][k] - atom_info[i][j].z) * (n1_z[i][k] - atom_info[i][j].z));
                dist_n2 = sqrt((n2_x[i][k] - atom_info[i][j].x) * (n2_x[i][k] - atom_info[i][j].x) + (n2_y[i][k] - atom_info[i][j].y) * (n2_y[i][k] - atom_info[i][j].y) + (n2_z[i][k] - atom_info[i][j].z) * (n2_z[i][k] - atom_info[i][j].z));
                dist_no = sqrt((no_x[i][k] - atom_info[i][j].x) * (no_x[i][k] - atom_info[i][j].x) + (no_y[i][k] - atom_info[i][j].y) * (no_y[i][k] - atom_info[i][j].y) + (no_z[i][k] - atom_info[i][j].z) * (no_z[i][k] - atom_info[i][j].z));
                if (dist_n1 < r_NA || dist_n2 < r_NA || dist_no < r_NA)
                    count++;
            }
            for (k = 1; k <= 20; k++)
            {
                dist_po = sqrt((po_x[i][k] - atom_info[i][j].x) * (po_x[i][k] - atom_info[i][j].x) + (po_y[i][k] - atom_info[i][j].y) * (po_y[i][k] - atom_info[i][j].y) + (po_z[i][k] - atom_info[i][j].z) * (po_z[i][k] - atom_info[i][j].z));
                if (dist_po < r_NA)
                    count++;
            }
            NA_n_env[i][j - N_PEP] = count;
        }
    }
    
    // 计算CL离子周围电正环境随时间的变化
    float dist_p1, dist_p2;
    for (i = 1; i <= t; i++)
    {
        for (j = N_PEP + N_NA + 1; j <= N_PEP + N_NA + N_CL; j++)
        {
            count = 0;
            for (k = 1; k <= 20; k++)
            {
                dist_p1 = sqrt((p1_x[i][k] - atom_info[i][j].x) * (p1_x[i][k] - atom_info[i][j].x) + (p1_y[i][k] - atom_info[i][j].y) * (p1_y[i][k] - atom_info[i][j].y) + (p1_z[i][k] - atom_info[i][j].z) * (p1_z[i][k] - atom_info[i][j].z));
                dist_p2 = sqrt((p2_x[i][k] - atom_info[i][j].x) * (p2_x[i][k] - atom_info[i][j].x) + (p2_y[i][k] - atom_info[i][j].y) * (p2_y[i][k] - atom_info[i][j].y) + (p2_z[i][k] - atom_info[i][j].z) * (p2_z[i][k] - atom_info[i][j].z));
                if (dist_p1 < r_CL || dist_p2 < r_CL)
                    count++;
            }
            CL_p_env[i][j - N_PEP - N_NA] = count;
        }
    }
    
    // 计算NA离子速度与电负环境关联的二维分布
    for (i = 0; i <= 30; i++)
    {
        for (j = 0; j <= 6; j++)
        {
            count = 0;
            for (k = 1; k <= t - 1; k++)
            {
                for (l = 1; l <= N_NA; l++)
                {
                    if (v_NA[k][l] >= i && v_NA[k][l] < i + 1)
                    {
                        if (NA_n_env[k][l] >= (float)j && NA_n_env[k][l] < (float)(j + 1))
                            count++;
                    }
                }
            }
            P_v_env_NA[i][j] = (float)count / (t - 1) / N_NA;
            outfile_1 << i << " " << j << " " << P_v_env_NA[i][j] << endl;
        }
        outfile_1 << endl;
    }
    outfile_1.close();
    
    // 计算CL离子速度与电正环境关联的二维分布
    for (i = 0; i <= 30; i++)
    {
        for (j = 0; j <= 6; j++)
        {
            count = 0;
            for (k = 1; k <= t - 1; k++)
            {
                for (l = 1; l <= N_CL; l++)
                {
                    if (v_CL[k][l] >= i && v_CL[k][l] < i + 1)
                    {
                        if (CL_p_env[k][l] >= (float)j && CL_p_env[k][l] < (float)(j + 1))
                            count++;
                    }
                }
            }
            P_v_env_CL[i][j] = (float)count / (t - 1) / N_CL;
            outfile_2 << i << " " << j << " " << P_v_env_CL[i][j] << endl;
        }
        outfile_2 << endl;
    }
    outfile_2.close();

    // 统计NA离子的速度分布
    for (i = 0; i <= 120; i++)
    {
        count = 0;
        for (j = 1; j <= t - 1; j++)
        {
            for (k = 1; k <= N_NA; k++)
            {
                if (v_NA[j][k] >= (float)i && v_NA[j][k] < (float)(i + 1))
                    count++;
            }
        }
        P_v_NA[i] = (float)count / (t - 1) / N_NA;
        outfile_3 << i << " " << P_v_NA[i] << endl;
    }
    outfile_3.close();

    // 统计CL离子的速度分布
    for (i = 0; i <= 120; i++)
    {
        count = 0;
        for (j = 1; j <= t - 1; j++)
        {
            for (k = 1; k <= N_CL; k++)
            {
                if (v_CL[j][k] >= (float)i && v_CL[j][k] < (float)(i + 1))
                    count++;
            }
        }
        P_v_CL[i] = (float)count / (t - 1) / N_CL;
        outfile_4 << i << " " << P_v_CL[i] << endl;
    }
    outfile_4.close();

    return 0;
}
