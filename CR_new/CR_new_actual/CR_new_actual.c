#include<stdio.h>
#include<malloc.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>

int NETWORK_SIZE = 50;
double PROBABILITY_OF_EAGE = 0.8;
double** adjacentMatrix_w;
double* sum_w;
double pd;
double lambda;
int k = 15;//每次遇到的人
int theta = 4;//免疫力较强集群的受感染力
double alpha = 0.4;//,免疫力薄弱集群个数占比
double mu = 0.2;
int tmax = 501;
int t_test = 1;


FILE * fp3;
FILE * fp4;
FILE *fp1;
FILE *fp2;
FILE * fp5;
int *city_people_num;

double ** peopleTodes_num;
double *I_pro_ini;
double *P;
double *PI;
double *rho;
double *rho_t;
double ** adjacentMatrix_w_alpha;
double** contact_M;
double *sum;
struct people
{
	int sta;
	char ini_state;
	int des;
	char fin_state;
}** city;

void load_city_people_num();
void initial();
void load_Matrix_w_fenzi();
void load_Matrix_w();//将邻接矩阵变为带权的
void load_struct_people();
void load_peopleTodes_num();
void load_I_pro_ini();
void load_P();
void load_PI();
void load_rho(int t);//100个步长后的每个集群的感染占比
void load_rho_t(int t);
void load_limit_rho();
double load_I_pro_zong();
void load_write_file();
//void load_Matrix_w_alpha();
//void load_contact_M();
//void load_write_M();

int main()
{
	int i, j;
	load_city_people_num();
	initial();
	load_struct_people();//y
	load_I_pro_ini();//y
	load_Matrix_w_fenzi();
	load_Matrix_w();
	
	load_write_file();
	//load_write_M();
	return 0;
}
void load_city_people_num()
{
	city_people_num = (int *)malloc(sizeof(int)*NETWORK_SIZE);
	int i, a;
	errno_t err;
	err = fopen_s(&fp1, "peoplenum_a.txt", "r");
	if (err != 0)
	{
		puts("不能打开文件");
	}
	for (i = 0; i < NETWORK_SIZE; i++)
	{
		fscanf_s(fp1, "%d ", &a);
		city_people_num[i] = a;
	}
	fclose(fp1);
}
void initial()
{
	adjacentMatrix_w = (double**)malloc(sizeof(double *) * NETWORK_SIZE); //分配指针数组
	int i;
	for (i = 0; i < NETWORK_SIZE; i++)
	{
		adjacentMatrix_w[i] = (double *)malloc(sizeof(double) * NETWORK_SIZE);//分配每个指针指向的数组
	}
	city = (struct people**)malloc(sizeof(struct people *) * NETWORK_SIZE);//有多少个城市,city为指向指针的指针
	for (i = 0; i < NETWORK_SIZE; i++)
	{
		city[i] = (struct people*)malloc(sizeof(struct people)* city_people_num[i]);//叶子城市里有多少个人
	}
	sum_w = (double *)malloc(sizeof(double)*NETWORK_SIZE);
	peopleTodes_num = (double**)malloc(sizeof(double*) * NETWORK_SIZE); //分配指针数组
	for (i = 0; i < NETWORK_SIZE; i++)
	{
		peopleTodes_num[i] = (double*)malloc(sizeof(double) * NETWORK_SIZE);//分配每个指针指向的数组
	}
	I_pro_ini = (double *)malloc(sizeof(double)*NETWORK_SIZE);
	P = (double *)malloc(sizeof(double)*NETWORK_SIZE);
	PI = (double *)malloc(sizeof(double)*NETWORK_SIZE);
	rho = (double *)malloc(sizeof(double)*NETWORK_SIZE);
	rho_t = (double *)malloc(sizeof(double)*t_test);
	contact_M = (double**)malloc(sizeof(double *) * NETWORK_SIZE); //分配指针数组
	for (i = 0; i < NETWORK_SIZE; i++)
	{
		contact_M[i] = (double *)malloc(sizeof(double) * NETWORK_SIZE);//分配每个指针指向的数组
	}
	adjacentMatrix_w_alpha = (double**)malloc(sizeof(double *) * NETWORK_SIZE); //分配指针数组
	for (i = 0; i < NETWORK_SIZE; i++)
	{
		adjacentMatrix_w_alpha[i] = (double *)malloc(sizeof(double) * NETWORK_SIZE);//分配每个指针指向的数组
	}

	sum = (double *)malloc(sizeof(double)*NETWORK_SIZE);
}
void load_struct_people()
{
	int i, j, m;
	for (i = 0; i < NETWORK_SIZE; i++)
	{
		for (j = 0; j < city_people_num[i]; j++)
		{
			city[i][j].sta = i;
			if (i == 0 && j < 10)
			{
				city[i][j].ini_state = 'I';
			}
			else
			{
				city[i][j].ini_state = 'S';
			}
			city[i][j].des = 100;
			city[i][j].fin_state = 'N';
		}
	}
}
void load_I_pro_ini()
{
	int i, j, m, S, I;
	for (m = 0; m < NETWORK_SIZE; m++)
	{
		S = 0;//记录每个种群的S人数
		I = 0;//记录每个种群的I人数
		for (i = 0; i < NETWORK_SIZE; i++)
		{
			for (j = 0; j < city_people_num[i]; j++)
			{
				if (city[i][j].sta == m && city[i][j].ini_state == 'S')
					S++;
				if (city[i][j].sta == m && city[i][j].ini_state == 'I')
					I++;
			}
		}
		I_pro_ini[m] = (double)I / (I + S);
	}
}
void load_Matrix_w_fenzi()
{
	int i, j;
	int a;
	errno_t err;
	err = fopen_s(&fp2, "w_a.txt", "r");
	if (err != 0)
	{
		puts("不能打开文件");
	}
	for (i = 0; i < NETWORK_SIZE; i++)
	{
		for (j = 0; j < NETWORK_SIZE; j++)
		{
			fscanf_s(fp2, "%d ", &a);
			adjacentMatrix_w[i][j] = a;
		}
	}
	fclose(fp2);
}
void load_Matrix_w()
{
	int i, j;
	double s;
	for (i = 0; i < NETWORK_SIZE; i++)
	{
		s = 0;
		for (j = 0; j < NETWORK_SIZE; j++)
		{
			s = adjacentMatrix_w[i][j] + s;
		}
		sum_w[i] = s;
	}
	for (i = 0; i < NETWORK_SIZE; i++)
	{
		for (j = 0; j < NETWORK_SIZE; j++)
		{
			adjacentMatrix_w[i][j] = adjacentMatrix_w[i][j] / sum_w[i];
		}
	}
}
void load_peopleTodes_num()
{
	int i, j;
	for (i = 0; i < NETWORK_SIZE; i++)
	{
		for (j = 0; j < NETWORK_SIZE; j++)
		{
			if (i == j)
			{
				peopleTodes_num[j][i] = (1 - pd)*city_people_num[i];
			}
			if (i != j)
			{
				peopleTodes_num[j][i] = pd * adjacentMatrix_w[j][i] * city_people_num[j];
			}
		}
	}
}
double load_C(int m, int n)//排列组合 C(rho_[j]*k,theta)
{
	double a;
	int i, j, k,b=1,c=1,d=1;
	for (i = m; i > 0; i--)
	{
		b = b * i;
	}
	for (j = n; j > 0; j--)
	{
		c = c * j;
	}
	for (k = (m-n); k > 0; k--)
	{
		d = d * k;
	}
	a = (double)b / (c * d);
	return a;
}
double load_binomial(double lambda ,int m,int n )//lambda 不取1  二项分布
{
	double a, b;
	a = load_C(m, n);
	b = pow((1 - lambda), (m - n))*pow(lambda, n)*a;
	return b;
}
double load_binomial_double(double lambda, double m, int n)//lambda 不取1  二项分布
{
	double a, b;
	a =0.6*load_C(floor(m), n)+0.4*load_C(ceil(m),n);
	b = pow((1 - lambda), (m - n))*pow(lambda, n)*a;
	return b;
} 
double load_rho_nji(int i)
{
	int j;
	double a=0,b=0,c;
	for (j = 0; j < NETWORK_SIZE; j++)
	{
		a = peopleTodes_num[j][i] * rho[j]+a;
	}
	for (j = 0; j < NETWORK_SIZE; j++)
	{
		b = peopleTodes_num[j][i] + b;
	}
	return c = a / b;
}
void load_P()
{
	int i, j,m;
	double a,b1,b2;
	for (i = 0; i < NETWORK_SIZE; i++)
	{
		b1 = 0;
		a = load_rho_nji(i)*k;
		if ((int)a == a)
		{
			for (m = 0; m < theta - 1; m++)
			{
				b1 = load_binomial(lambda, a, m) + b1;
			}
			b2 = load_binomial(lambda, a, 0);
		}
		else
		{
			for (m = 0; m < theta - 1; m++)
			{
				b1 = load_binomial_double(lambda, a, m) + b1;
			}
			b2 = load_binomial_double(lambda, a, 0);
		}
		P[i] = alpha*(1-b2)+(1-alpha)*(1-b1);
	}
}
void load_PI()
{
	int i, j;
	double sum;
	for (i = 0; i < NETWORK_SIZE; i++)
	{
		sum = 0;
		for (j = 0; j < NETWORK_SIZE; j++)
		{
			sum = (adjacentMatrix_w[i][j] * P[j]) + sum;
		}
		PI[i] = (1 - pd)*P[i] + pd * sum;
	}
}
void load_rho(int t)
{
	int i;
	for (i = 0; i < NETWORK_SIZE; i++)
	{
		if (t == 0)
			rho[i] = I_pro_ini[i];
		else
		{
			rho[i] = (1 - mu)*rho[i] + (1 - rho[i])*PI[i];
		}
	}
}
void load_rho_t(int t)
{
	int i;
	double sum = 0;
	for (i = 0; i < NETWORK_SIZE; i++)
	{
		sum = sum + rho[i];
	}
	rho_t[t - (tmax - t_test)] = (double)sum / NETWORK_SIZE;
}
void load_limit_rho()
{
	int t, i, j;
	for (t = 0; t < tmax; t++)
	{
		//printf(" pd:%If lambda:%If t:%d \n", pd,lambda,t);
		load_rho(t);
		/*printf("rho\n");
		for (i = 0; i < NETWORK_SIZE; i++)
		{
			printf("%.20lf ", rho[i]);

		}
		printf("\n");*/
		//printf("\n");
		if (t >= (tmax - t_test))
		{
			load_rho_t(t);
		}
		/*printf("W:\n");
		for (j = 0; j < NETWORK_SIZE; j++)
		{
			for (i = 0; i < NETWORK_SIZE; i++)
			{
				printf("%f ", adjacentMatrix_w[j][i]);
			}
			printf("\n");
		}
		printf("\n");*/
		
		/*printf("R_fenzi:\n");
		for (j = 0; j < NETWORK_SIZE; j++)
		{
			for (i = 0; i < NETWORK_SIZE; i++)
			{
				printf("%f ", adjacentMatrix_R[j][i]);
				if (adjacentMatrix_R[j][i] < 0)
				{
					printf("%f ", adjacentMatrix_R[j][i]);
					printf("%d,%d\n", j, i);
				}

			}
			printf("\n");
		}
		printf("\n");*/
		load_peopleTodes_num();
		/*printf("peopletodes:\n");
		for (j = 0; j < NETWORK_SIZE; j++)
		{
			for (i = 0; i < NETWORK_SIZE; i++)
			{
				printf("%f ", peopleTodes_num[j][i]);
			}
			printf("\n");
		}
		printf("\n");*/
		load_P();
		load_PI();
	}
}

double load_I_pro_zong()
{
	int i;
	double pro_zong, zong = 0;
	/*for (i = 0; i < t_test; i++)
	{
		printf("%d:%f\n",i, rho_t[i]);
	}*/
	for (i = 0; i < t_test; i++)
	{
		zong = rho_t[i] + zong;
	}
	pro_zong = (double)zong / t_test;
	return pro_zong;
}
void load_write_file()
{
	double pro;
	int i, j;
	errno_t err;
	err = fopen_s(&fp3, "I_L_pd.txt", "w");
	if (err != 0)
	{
		puts("不能打开文件");
	}
	for (pd = 0; pd <= 1.1; pd = pd + 0.5)//50
	{
		fprintf_s(fp3, "pd=%f\n", pd);
		for (lambda = 0; lambda <= 0.101; lambda = lambda + 0.002)//50
		{
			load_limit_rho();
			pro = load_I_pro_zong();
			fprintf_s(fp3, "%f ", pro);
			//printf("%f %f %f %.2f", pd, la, lambda, pro);
		}
		fprintf_s(fp3, "\n");
		fprintf_s(fp3, "\n");
		//printf("\n");
	}
	fclose(fp3);
}

		
/*void load_write_file()
{
	double pro;
	int i, j;
	errno_t err;
	err = fopen_s(&fp3, "help1.txt", "w");
	if (err != 0)
	{
		puts("不能打开文件");
	}
	for (pd = 0.3; pd <= 0.3; pd = pd + 0.2)//11
	{
		fprintf_s(fp3, "pd=%f\n", pd);
		for (la = -6; la <= 6; la = la + 3)//50
		{
			for (lambda = 0; lambda <= 0.01; lambda = lambda + 0.0002)//50
			{
				load_limit_rho();
				pro = load_I_pro_zong();
				fprintf_s(fp3, "%f ", pro);
				//printf("%f %f %f %.2f", pd, la, lambda, pro);
			}
			fprintf_s(fp3, "\n");
			//printf("\n");
		}
		fprintf_s(fp3, "\n");
		fprintf_s(fp3, "\n");
		fprintf_s(fp3, "\n");
	}
	fclose(fp3);
}

void load_Matrix_w_alpha()
{
	int i, j;
	double s;
	for (i = 0; i < NETWORK_SIZE; i++)
	{
		s = 0;
		for (j = 0; j < NETWORK_SIZE; j++)
		{
			if (adjacentMatrix_w[i][j] != 0)
			{
				s = s + pow(adjacentMatrix_w[i][j], la);
			}
		}
		sum[i] = s;//存储每一行的和，共i行
	}
	for (i = 0; i < NETWORK_SIZE; i++)
	{
		for (j = 0; j < NETWORK_SIZE; j++)
		{
			if (adjacentMatrix_w[i][j] == 0)
			{
				adjacentMatrix_w_alpha[i][j] = 0;
			}
			else
				adjacentMatrix_w_alpha[i][j] = pow(adjacentMatrix_w[i][j], la) / sum[i];//出现误差
		}
	}
}
void load_contact_M()
{
	double sum;
	int i, j, l, m, n;
	for (i = 0; i < NETWORK_SIZE; i++)
	{
		for (j = 0; j < NETWORK_SIZE; j++)
		{
			sum = 0;
			for (l = 0; l < NETWORK_SIZE; l++)
			{
				sum = sum + (adjacentMatrix_w_alpha[i][l] * adjacentMatrix_w_alpha[j][l]);
			}
			if (i == j)
			{
				contact_M[i][j] = (1 - pd)*(1 - pd) *city_people_num[j] + pd * (1 - pd)*city_people_num[j] * (adjacentMatrix_w_alpha[i][j] + adjacentMatrix_w_alpha[j][i]) + pd * pd*city_people_num[j] * sum;
			}
			if (i != j)
			{
				contact_M[i][j] = pd * (1 - pd)*city_people_num[j] * (adjacentMatrix_w_alpha[i][j] + adjacentMatrix_w_alpha[j][i]) + pd * pd*city_people_num[j] * sum;
			}
		}
	}
}
void load_write_M()
{
	double pro;
	int i, j;
	errno_t err;
	err = fopen_s(&fp4, "M_a_la.txt", "w");
	if (err != 0)
	{
		puts("不能打开文件");
	}
	for (pd = 0.3; pd <= 0.8; pd = pd + 0.2)
	{
		fprintf_s(fp4, "pd=%f\n", pd);
		for (la = -6; la <= 6; la = la + 3)
		{
			load_Matrix_w_alpha();
			load_contact_M();
			//fprintf_s(fp3, "pd=%f\n", pd);
			for (i = 0; i < NETWORK_SIZE; i++)
			{
				for (j = 0; j < NETWORK_SIZE; j++)
				{
					fprintf_s(fp4, "%f ", contact_M[i][j]);
				}
				fprintf_s(fp4, "\n");
			}
		}
		fprintf_s(fp4, "\n");
		fprintf_s(fp4, "\n");
		fprintf_s(fp4, "\n");
	}
	fclose(fp4);
}*/

/*for (j = 0; j < NETWORK_SIZE; j++)
		{
			for (i = 0; i < NETWORK_SIZE; i++)
			{
				printf("%f ", adjacentMatrix[j][i]);
			}
			printf("\n");
		}
		printf("\n");*/

		/*
		printf("rho:\n");
				for (i = 0; i < NETWORK_SIZE; i++)
				{
					printf("%f ", rho[i]);
				}
				printf("\n");
				printf("\n");
		*/

		/*
		void load_city_people_num()
		{
			city_people_num = (int *)malloc(sizeof(int)*NETWORK_SIZE);
			int i, a, num;
			errno_t err;
			err = fopen_s(&fp, "people3.txt", "r");
			if (err != 0)
			{
				puts("不能打开文件");
			}
			for (i = 0; i < NETWORK_SIZE; i++)
			{
				num = 0;
				while (1)
				{
					fscanf_s(fp, "%d", &a);
					if (feof(fp))
						break;
					if (a == i)
					{
						num++;
						fseek(fp, 3, 1);
					}
					else
						break;
				}
				city_people_num[i] = num;
			}
			fclose(fp);
		}
		*/

		/*
		void load_struct_people()
		{
			int i, j;
			errno_t err;
			err = fopen_s(&fp3, "people3.txt", "r");
			if (err != 0)
			{
				puts("不能打开文件");
			}
			for (i = 0; i < NETWORK_SIZE; i++)
			{
				for (j = 0; j < city_people_num[i]; j++)
				{
					fscanf_s(fp3, "%d", &city[i][j].sta);
					fscanf_s(fp3, "%c", &city[i][j].ini_state, sizeof(char));
					city[i][j].des = 22;
					city[i][j].fin_state = 'N';
				}
			}
			fclose(fp3);
		}
		*/