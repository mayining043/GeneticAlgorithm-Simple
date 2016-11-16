#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#define POPSIZE 500
#define maximization 1
#define minimization 2
#define cmax 100
#define cmin 0
#define length1 10
#define length2 10
#define chromlength length1+length2  //染色体长度
int functionmode=maximization;
int popsize;        //种群大小
int maxgeneration;  //最大世代数
double pc;          //交叉率
double pm;          //变异率
struct individual
{
	char  chrom[chromlength+1];
	double value;         
	double fitness;      //适应度
};
int generation;      //世代数
int best_index;
int worst_index;
struct individual bestindividual;  //最佳个体
struct individual worstindividual; //最差个体
struct individual currentbest;
struct individual population[POPSIZE];
//函数声明                                       
void generateinitialpopulation();                 
void generatenextpopulation();
void evaluatepopulation();
long decodechromosome(char *,int,int);
void calculateobjectvalue();
void calculatefitnessvalue();
void findbestandworstindividual();
void performevolution();
void selectoperator();
void crossoveroperator();
void mutationoperator();
void input();
void outputtextreport();

void generateinitialpopulation( )  //种群初始化
{
	int i,j;

	for (i=0;i<popsize; i++)
	{
		for(j=0;j<chromlength;j++)
		{
			population[i].chrom[j]=(rand()%10<5)?'0':'1';
		}
		population[i].chrom[chromlength]='\0';
	}
}
void generatenextpopulation()  //生成下一代
{
	selectoperator();
	crossoveroperator();
	mutationoperator();
}
void evaluatepopulation()   //评价个体，求最佳个体
{
	calculateobjectvalue();
	calculatefitnessvalue();
	findbestandworstindividual();
}
long decodechromosome(char *string ,int point,int length) //给染色体解码
{
	int i;
	long decimal=0;
	char*pointer;
	for(i=0,pointer=string+point;i<length;i++,pointer++)
		if(*pointer-'0')

		{decimal +=long(2<<i);
	}
	return (decimal);

}
void calculateobjectvalue()  //计算函数值
{
	int i;
	long temp1,temp2;
	double x1,x2;

	for (i=0; i<popsize; i++)
	{
		temp1=decodechromosome(population[i].chrom,0,length1);
		temp2=decodechromosome(population[i].chrom,length1,length2);
		x1=4.096*temp1/1023.0-2.048;
		x2=4.096*temp2/1023.0-2.048;
		population[i].value=100*(x1*x1-x2)* (x1*x1-x2)+(1-x1)*(1-x1);
	}
}
void calculatefitnessvalue()//计算适应度
{
	int  i;
	double temp;
	for(i=0;i<popsize;i++)
	{
		if(functionmode==maximization)
		{if((population[i].value+cmin)>0.0)
		{temp=cmin+population[i].value;}
		else
		{temp=0.0;
		}
		}
		else if (functionmode==minimization)
		{
			if(population[i].value<cmax)
			{temp=cmax-population[i].value;}
			else{ temp=0.0;}
		}
		population[i].fitness=temp;
	}
}
void findbestandworstindividual( ) //求最佳个体和最差个体
{
	int i;
	double sum=0.0;

	bestindividual=population[0];
	worstindividual=population[0];

	for (i=1;i<popsize; i++){
		if (population[i].fitness>bestindividual.fitness){
			bestindividual=population[i];
			best_index=i;
		}
		else if (population[i].fitness<worstindividual.fitness)
		{
			worstindividual=population[i];
			worst_index=i;

		}
		sum+=population[i].fitness;
	}
	if (generation==0){
		currentbest=bestindividual;
	}
	else{
		if(bestindividual.fitness>=currentbest.fitness){
			currentbest=bestindividual;
		}
	}
}
void performevolution() //演示评价结果
{
	if (bestindividual.fitness>currentbest.fitness){
		currentbest=population[best_index];
	}
	else{
		population[worst_index]=currentbest;
	}
}
void selectoperator() //比例选择算法
{
	int i,index;
	double p,sum=0.0;
	double cfitness[POPSIZE];

	struct individual newpopulation[POPSIZE];
	for(i=0;i<popsize;i++)
	{sum+=population[i].fitness;}

	for(i=0;i<popsize; i++){
		cfitness[i]=population[i].fitness/sum;
	}

	for(i=1;i<popsize; i++){
		cfitness[i]=cfitness[i-1]+cfitness[i];
	}

	for (i=0;i<popsize;i++)
	{
		p=rand()%1000/1000.0;
		index=0;
		while (p>cfitness[index])
		{
			index++;
		}
		newpopulation[i]=population[index];
	}
	for(i=0;i<popsize; i++){
		population[i]=newpopulation[i];
	}
}
void crossoveroperator() //交叉算法
{
	int i,j;
	int index[POPSIZE];
	int point,temp;
	double p;
	char ch;


	for (i=0;i<popsize;i++){
		index[i]=i;
	}
	for (i=0;i<popsize;i++){
		point=rand()%(popsize-i);
		temp=index[i];
		index[i]=index[point+i];
		index[point+i]=temp;
	}

	for (i=0;i<popsize-1;i+=2){
		p=rand()%1000/1000.0;
		if (p<pc){
			ppoint=rand()%(chromlength-1)+1;
			for (j=point; j<chromlength;j++){
				ch=population[index[i]].chrom[j];
				population[index[i]].chrom[j]=population[index[i+1]].chrom[j];
				population[index[i+1]].chrom[j]=ch;
			}
		}
	}

}
void mutationoperator() //变异操作
{
	int i,j;
	double p;

	for (i=0;i<popsize;i++){
		for(j=0;j<chromlength;j++){
			p=rand()%1000/1000.0;
			if (p<pm){
				population[i].chrom[j]=(population[i].chrom[j]=='0')?'1':'0';
			}
		}
	}
}
void input() //数据输入
{    
}
void outputtextreport()//数据输出
{
	int i;
	double sum;
	double average;
	sum=0.0;
	for(i=0;i<popsize;i++)
	{sum+=population[i].value;}
	average=sum/popsize;

	printf("当前世代=%d\n当前世代平均函数值=%f\n当前世代最高函数值=%f\n",generation,average,population[best_index].value);

}
void Main()    //主函数
{   int i;
printf("本程序为求函数y=100*(x1*x1-x2)*(x1*x2-x2)+(1-x1)*(1-x1)的最大值 \n其中-2.048<=x1,x2<=2.048\n");
generation=0;
input();
generateinitialpopulation();
evaluatepopulation();
while(generation<maxgeneration)
{
	generation++;
	generatenextpopulation();
	evaluatepopulation();
	performevolution();
	outputtextreport();
}
printf("\n");
printf("         统计结果:        ");
printf("\n");
printf("最大函数值等于：%f\n",currentbest.fitness);
printf("其染色体编码为：");
for (i=0;i<chromlength;i++)
{
	printf("%c",currentbest.chrom[i]);

}
printf("\n"); 
system("pause");
}
