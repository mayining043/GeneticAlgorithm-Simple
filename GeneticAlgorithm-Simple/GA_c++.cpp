#include<iostream>
#include<cstdlib>
#include<cstdio>
#include<string>
#include<ctime>
#include<cmath>
#include<vector>
#include<list>
#include<algorithm>
#include<fstream>
using namespace std;

ofstream out;
long long  _sum;
long long  _max_sum;

//个体模板类，使用浮点数编码(FPR)方式
template <class T>class individual{
public:
	vector<T> chromosome;//染色体
	double fitness;//适应度
	double val;//函数值
	void print(){
		vector<T>::iterator it=this->chromosome.begin();	
		cout<<"个体有"<<chromosome.size()<<"个基因：";
		for(;it!=this->chromosome.end();it++)
			cout<<*it<<" ";
		cout<<endl;
	}
	void write(){
		vector<T>::iterator it=this->chromosome.begin();	
		cout<<"个体有"<<chromosome.size()<<"个基因：\n";
		for(;it!=this->chromosome.end();it++)
			out<<*it<<" ";
		out<<endl;
		out<<val<<endl;
	}
};

//GA遗传算法类模板
template <class T>class GAalg{

public:
	//算法初始化
	GAalg(int sg,int sc,int mg,double pc,double pm):
	  sizeof_generation(sg),sizeof_chrom(sc),max_generation(mg),p_crossover(pc),p_mutation(pm){
		//文件初始化
		out.open("data.txt");
		cout<<"Start GA algorithm...\n";
		generationbox.clear();
		//更新代数
		cur_generation=0;
		_sum=0;
		//随机数撒种
		srand(time(NULL));
	}
	//算法结束
	virtual ~GAalg(){
		cout<<"End GA algorithm...\n";
		cout<<endl<<endl<<"历史最优个体为:\n";
		bestone.print();
		generationbox.clear();
		system("pause");
	}
	//打印所有个体
	virtual void print(){
		cout<<"共有"<<sizeof_generation<<"个个体，每个个体基因如下：\n";
		vector<individual<T> >::iterator it=generationbox.begin();
		for(;it!=generationbox.end();it++)
			(*it).print();
	}
	//初始化种群中的个体
	virtual void InitIndividual(individual<T>&my){
		do{
			my.chromosome.clear();
			for(int i=0;i<sizeof_chrom;i++){
				cout<<"初始化第"<<i+1<<"个基因"<<endl;
				T temp;
				InitChromosome(temp);
				my.chromosome.push_back(temp);
			}
		}while(!IsOk(my));
		return;
	}
	//种群初始化，初始化基因
	virtual void GenerateInit(){
		for(int i=0;i<sizeof_generation;i++){
			cout<<"正在初始化第"<<i+1<<"个个体"<<endl;
			individual<T> obj;
			InitIndividual(obj);
			generationbox.push_back(obj);
		}
		cur_generation=0;
		return;
	}
	//生成下一代
	virtual void NextPopulation(){
		cur_generation++;
		cout<<endl<<endl<<"正在处理第"<<cur_generation<<"代...\n";
		SelectOperator();//选择算子
		CrossoverOperator();//交配算子
		MutationOperator();//变异算子
	}
	//评价个体
	virtual void ValuatePopulation(){
		CalculateFitnessValue();//计算适应度
		FindBestWorstIndividual();//计算最佳个体和最差个体
	}
	//计算最佳个体和最差个体
	virtual void FindBestWorstIndividual(){
		//初始化当前种群极值
		cur_bestone=*generationbox.begin();
		cur_worstone=*generationbox.begin();
		//遍历所用iterator
		vector<individual<T> >::iterator it=generationbox.begin();
		int i=0;
		//寻找当前最佳最差
		for (it=generationbox.begin();it!=generationbox.end();it++){  
			if(it->fitness>=cur_bestone.fitness){
				cur_bestone=*it;
				best_index=i;
			}
			else if(it->fitness<=cur_worstone.fitness){
				cur_worstone=*it;
				worst_index=i;
			}
			i++;
		}  
		//更新历史最好最差
		if (cur_generation==0){
			bestone=cur_bestone;
			worstone=cur_worstone;
		}
		else{
			if(cur_bestone.fitness>=bestone.fitness)
				bestone=cur_bestone;
			if(cur_worstone.fitness<=worstone.fitness)
				worstone=cur_worstone;
		}
			return;
	}
	//交配方法，交换自从位置pos之后的x与y的染色体
	void cross(int x,int y,int pos,int end){
		//临时变量
		vector<individual<T> >::iterator itx=generationbox.begin();
		vector<individual<T> >::iterator ity=generationbox.begin();
		individual<T> valx;
		individual<T> valy;
		//寻找暂存
		for(int i=0;i!=x;i++)
			itx++;
		for(int i=0;i!=y;i++)
			ity++;
		valx=*itx;
		valy=*ity;
		//cout<<"个体"<<x+1<<"和"<<y+1<<"在基因位置"<<pos+1<<"开始交配直到位置"<<end+1<<endl;
		//开始交配x
		vector<T>::iterator from,targe;
		from=valy.chromosome.begin();
		targe=(*itx).chromosome.begin();
		for(int i=0;i!=pos;i++){
			from++;
			targe++;
		}
		for(int p=pos;p<=end;p++){
			*targe=*from;
			from++;
			targe++;
		}
		//开始交配y
		from=valx.chromosome.begin();
		targe=(*ity).chromosome.begin();
		for(int i=0;i!=pos;i++){
			from++;
			targe++;
		}
		for(int p=pos;p<=end;p++){
			*targe=*from;
			from++;
			targe++;
		}
		//判断是否为有效解
		if(!IsOk(*itx)){
			*itx=valx;
			*ity=valy;
			cout<<"fail\n";
			return;
		}
		if(!IsOk(*ity)){
			*itx=valx;
			*ity=valy;
			cout<<"fail\n";
			return;
		}
		return;
	}
	//选择算子*轮盘赌选算法ROUTE WHEEL SELECTION*
	virtual void SelectOperator(){
		//cout<<"选择算子*轮盘赌选算法，正在更新种群...\n";
		//cout<<"更新前：\n";
		//this->print();
		int index,i=0;
		//随机概率以及总适应度
		double p,sum=0.0;
		//新种群temp
		vector<individual<T> >newgenerationbox;
		//前缀和适应度
		double *addfit=new double[sizeof_generation];
		//遍历所用iterator
		vector<individual<T> >::iterator it=generationbox.begin();
		//求总适应度
		for (it=generationbox.begin();it!=generationbox.end();it++)
			sum+=it->fitness;
		//求轮盘区域概率
		for (it=generationbox.begin();it!=generationbox.end();it++)
			addfit[i++]=it->fitness/sum;
		
		//求轮盘区域概率数组前缀和
		i=1;
		for (it=++generationbox.begin();it!=generationbox.end();it++,i++)
			addfit[i]=addfit[i-1]+addfit[i];
		/*//输出前缀和数组
		for(int j=0;j<sizeof_generation;j++)
			cout<<addfit[j]<<" ";
		cout<<"\n"<<p<<"\n";*/
		//模拟轮盘赌选算法
		for(i=0;i<sizeof_generation;i++){
			p=rand()%1000/1000.0;
			index=0;
			while (p>=addfit[index])
				index++;
			it=generationbox.begin();
			for(int j=0;j<index;j++)
				it++;
			newgenerationbox.push_back(*it);
		}
		//更新当代种群
		generationbox.clear();
		it=newgenerationbox.begin();
		for(i=0;i<sizeof_generation;i++){
			generationbox.push_back(*it);
			it++;
		}
		//cout<<"更新后：\n";
		//this->print();
		return;
	}
	//交配算子
	virtual void CrossoverOperator(){
		//cout<<"交配算子，正在更新种群...\n";
		//产生交配概率，判断是否交配
		bool *pi=new bool[sizeof_generation];
		for(int i=0;i<sizeof_generation;i++){
			double p=rand()%1000/1000.0;
			pi[i]=(p<p_crossover);
		}
		//产生一个交配队列
		list<int>cl;
		for(int i=0;i<sizeof_generation;i++)
			if(pi[i])
				cl.push_back(i);
		//判断是否为偶数
		if(cl.size()%2!=0)
			cl.pop_back();
		//开始交配
		while(!cl.empty()){
			int x=cl.front();
			cl.pop_front();
			int y=cl.front();
			cl.pop_front();
			double p1=rand()%sizeof_chrom,p2=rand()%sizeof_chrom;
			if(p2<p1)
				swap(p1,p2);
			cross(x,y,p1,p2);
		}
		return;
	}
	//变异算子
	virtual void MutationOperator(){
		//cout<<"变异算子，正在更新种群...\n";
		int i,j;
		double p;
		//获取每一个个体
		vector<individual<T> >::iterator iti;
		for (iti=generationbox.begin();iti!=generationbox.end();iti++){
			//获取个体的每一个基因
			individual<T> temp=*iti;
			vector<T>::iterator itj=temp.chromosome.begin();
			for(;itj!=temp.chromosome.end();itj++){
				//判断是否变异
				p=rand()%1000/1000.0;
				if(p<p_mutation){
					//开始变异
					//cout<<"发生变异"<<endl;
					InitChromosome(*itj);
					//cout<<"有数据变异为"<<*itj<<endl;
				} 
				//else
					//cout<<"不发生变异"<<endl;
			}
			//是否保存变异
			if(IsOk(temp))
				*iti=temp;
		}
		return;
	}
	//数据输出展示函数
	virtual void OutputPrint(){
		out<<bestone.fitness<<endl;
/*
		cout<<"正在展示第"<<cur_generation<<"代...\n";
		cout<<"本次最优个体为:\n";
		cur_bestone.print();
		out<<bestone.val<<endl;
		cout<<"函数值为："<<cur_bestone.val<<endl;
		cout<<"适应度为："<<cur_bestone.fitness<<endl;
		cout<<"当代平均适应度：";
		double sum=0;
		for(vector<individual<T> >::iterator it=generationbox.begin();it!=generationbox.end();it++)
			sum+=it->fitness;
		cout<<sum/(double)sizeof_generation<<endl;
		cout<<"历史最优个体为:\n";
		bestone.print();
		cout<<"函数值为："<<bestone.val<<endl;
		cout<<"适应度为："<<bestone.fitness<<endl;*/
	}
	//输出最终答案
	virtual void ans(){
		out<<endl<<endl<<"历史最优个体为:\n";
		bestone.write();
		out<<"所用代数："<<cur_generation<<endl;
		out<<"适应值："<<bestone.fitness<<endl;
		out<<"函数值："<<bestone.val<<endl;
		out.close();
	}

	//初始化个体中的基因
	virtual void InitChromosome(T &)=0;
	//计算函数值
	virtual void CalculateVal(individual<T>&)=0;
	//计算适应度
	virtual void CalculateFitnessValue()=0;
	//终止条件
	virtual bool IsEnd()=0;
	//是否为有效解
	virtual bool IsOk(individual<T>cur){return 1;}
public:
	//染色体容器
	vector<individual<T> >generationbox;
	//最佳个体
	individual<T> cur_bestone;
	//最差个体
	individual<T> cur_worstone;
	//历史最佳个体
	individual<T> bestone;
	//历史最差个体
	individual<T> worstone;

	//数据成员
	int sizeof_generation;//种群大小
	int sizeof_chrom;//基因个数
	int max_generation;//最大遗传世代数
	int cur_generation;//当前世代数
	double p_crossover;//交配概率
	double p_mutation;//变异概率
	int best_index;//最佳个体位置
	int worst_index;//最差个体位置
};


class CalcMaxVal:public GAalg<double>{
public:
	CalcMaxVal(int sg,int sc,int mg,double pc,double pm):GAalg(sg,sc,mg,pc,pm){
		//innit the class
	}
	//初始化个体中的基因
	void InitChromosome(double &my){
		my=rand()%20100/100.0-100; 
		//cout<<"基因初始值为"<<my<<endl;
	}
	//计算适应度
	void CalculateFitnessValue(){
		//遍历所用iterator
		vector<individual<double> >::iterator it=generationbox.begin();
		//求函数值
		for (it=generationbox.begin();it!=generationbox.end();it++){
			CalculateVal(*it);
				it->fitness=pow(1/(it->val),8);
			_sum++;
		}
	}
	//计算函数值
	void CalculateVal(individual<double>&my){
		vector<double>::iterator it=my.chromosome.begin();
		double ans=0;
		for(;it!=my.chromosome.end();it++){
			ans+=(*it)*(*it);
		}
		my.val=ans;
		return;
	}
	//终止条件
	bool IsEnd(){
		if(cur_generation>=max_generation || _sum>_max_sum)
			return 1;
		else
			return 0;
	}
	//是否为有效解
	bool IsOk(individual<double> cur)
	{
		return 1;
	}
};

int main(){
	int sg,sc,mg;//种群个数，基因个数，最大代数 
	double pc,pm;//交配概率，变异概率
	/*printf("初始化全局变量:\n");
	printf("\t种群大小:");
	cin>>sg;  
	if((sg%2) != 0){
		printf("种群大小已设置为偶数\n");
		sg++;
	};
	printf("\t最大世代数：");
	cin>>mg;
	printf("\t交配率：");
	cin>>pc;
	printf("\t变异率：");
	cin>>pm;
	printf("\t最大评估次数：");
	cin>>_max_sum;*/
	_max_sum=10000000000;
	CalcMaxVal a(30,30,800000,0.60,0.005);//种群个数，基因个数，最大代数，交配概率，变异概率 
	a.GenerateInit();
	//a.print();
	cout<<"开始遗传算法："<<endl;
	a.ValuatePopulation();
	while(!a.IsEnd()){
		a.NextPopulation();
		a.ValuatePopulation();
		a.OutputPrint();
	}
	a.ans();
}