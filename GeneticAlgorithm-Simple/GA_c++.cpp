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

//����ģ���࣬ʹ�ø���������(FPR)��ʽ
template <class T>class individual{
public:
	vector<T> chromosome;//Ⱦɫ��
	double fitness;//��Ӧ��
	double val;//����ֵ
	void print(){
		vector<T>::iterator it=this->chromosome.begin();	
		cout<<"������"<<chromosome.size()<<"������";
		for(;it!=this->chromosome.end();it++)
			cout<<*it<<" ";
		cout<<endl;
	}
	void write(){
		vector<T>::iterator it=this->chromosome.begin();	
		cout<<"������"<<chromosome.size()<<"������\n";
		for(;it!=this->chromosome.end();it++)
			out<<*it<<" ";
		out<<endl;
		out<<val<<endl;
	}
};

//GA�Ŵ��㷨��ģ��
template <class T>class GAalg{

public:
	//�㷨��ʼ��
	GAalg(int sg,int sc,int mg,double pc,double pm):
	  sizeof_generation(sg),sizeof_chrom(sc),max_generation(mg),p_crossover(pc),p_mutation(pm){
		//�ļ���ʼ��
		out.open("data.txt");
		cout<<"Start GA algorithm...\n";
		generationbox.clear();
		//���´���
		cur_generation=0;
		_sum=0;
		//���������
		srand(time(NULL));
	}
	//�㷨����
	virtual ~GAalg(){
		cout<<"End GA algorithm...\n";
		cout<<endl<<endl<<"��ʷ���Ÿ���Ϊ:\n";
		bestone.print();
		generationbox.clear();
		system("pause");
	}
	//��ӡ���и���
	virtual void print(){
		cout<<"����"<<sizeof_generation<<"�����壬ÿ������������£�\n";
		vector<individual<T> >::iterator it=generationbox.begin();
		for(;it!=generationbox.end();it++)
			(*it).print();
	}
	//��ʼ����Ⱥ�еĸ���
	virtual void InitIndividual(individual<T>&my){
		do{
			my.chromosome.clear();
			for(int i=0;i<sizeof_chrom;i++){
				cout<<"��ʼ����"<<i+1<<"������"<<endl;
				T temp;
				InitChromosome(temp);
				my.chromosome.push_back(temp);
			}
		}while(!IsOk(my));
		return;
	}
	//��Ⱥ��ʼ������ʼ������
	virtual void GenerateInit(){
		for(int i=0;i<sizeof_generation;i++){
			cout<<"���ڳ�ʼ����"<<i+1<<"������"<<endl;
			individual<T> obj;
			InitIndividual(obj);
			generationbox.push_back(obj);
		}
		cur_generation=0;
		return;
	}
	//������һ��
	virtual void NextPopulation(){
		cur_generation++;
		cout<<endl<<endl<<"���ڴ����"<<cur_generation<<"��...\n";
		SelectOperator();//ѡ������
		CrossoverOperator();//��������
		MutationOperator();//��������
	}
	//���۸���
	virtual void ValuatePopulation(){
		CalculateFitnessValue();//������Ӧ��
		FindBestWorstIndividual();//������Ѹ����������
	}
	//������Ѹ����������
	virtual void FindBestWorstIndividual(){
		//��ʼ����ǰ��Ⱥ��ֵ
		cur_bestone=*generationbox.begin();
		cur_worstone=*generationbox.begin();
		//��������iterator
		vector<individual<T> >::iterator it=generationbox.begin();
		int i=0;
		//Ѱ�ҵ�ǰ������
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
		//������ʷ������
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
	//���䷽���������Դ�λ��pos֮���x��y��Ⱦɫ��
	void cross(int x,int y,int pos,int end){
		//��ʱ����
		vector<individual<T> >::iterator itx=generationbox.begin();
		vector<individual<T> >::iterator ity=generationbox.begin();
		individual<T> valx;
		individual<T> valy;
		//Ѱ���ݴ�
		for(int i=0;i!=x;i++)
			itx++;
		for(int i=0;i!=y;i++)
			ity++;
		valx=*itx;
		valy=*ity;
		//cout<<"����"<<x+1<<"��"<<y+1<<"�ڻ���λ��"<<pos+1<<"��ʼ����ֱ��λ��"<<end+1<<endl;
		//��ʼ����x
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
		//��ʼ����y
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
		//�ж��Ƿ�Ϊ��Ч��
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
	//ѡ������*���̶�ѡ�㷨ROUTE WHEEL SELECTION*
	virtual void SelectOperator(){
		//cout<<"ѡ������*���̶�ѡ�㷨�����ڸ�����Ⱥ...\n";
		//cout<<"����ǰ��\n";
		//this->print();
		int index,i=0;
		//��������Լ�����Ӧ��
		double p,sum=0.0;
		//����Ⱥtemp
		vector<individual<T> >newgenerationbox;
		//ǰ׺����Ӧ��
		double *addfit=new double[sizeof_generation];
		//��������iterator
		vector<individual<T> >::iterator it=generationbox.begin();
		//������Ӧ��
		for (it=generationbox.begin();it!=generationbox.end();it++)
			sum+=it->fitness;
		//�������������
		for (it=generationbox.begin();it!=generationbox.end();it++)
			addfit[i++]=it->fitness/sum;
		
		//�����������������ǰ׺��
		i=1;
		for (it=++generationbox.begin();it!=generationbox.end();it++,i++)
			addfit[i]=addfit[i-1]+addfit[i];
		/*//���ǰ׺������
		for(int j=0;j<sizeof_generation;j++)
			cout<<addfit[j]<<" ";
		cout<<"\n"<<p<<"\n";*/
		//ģ�����̶�ѡ�㷨
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
		//���µ�����Ⱥ
		generationbox.clear();
		it=newgenerationbox.begin();
		for(i=0;i<sizeof_generation;i++){
			generationbox.push_back(*it);
			it++;
		}
		//cout<<"���º�\n";
		//this->print();
		return;
	}
	//��������
	virtual void CrossoverOperator(){
		//cout<<"�������ӣ����ڸ�����Ⱥ...\n";
		//����������ʣ��ж��Ƿ���
		bool *pi=new bool[sizeof_generation];
		for(int i=0;i<sizeof_generation;i++){
			double p=rand()%1000/1000.0;
			pi[i]=(p<p_crossover);
		}
		//����һ���������
		list<int>cl;
		for(int i=0;i<sizeof_generation;i++)
			if(pi[i])
				cl.push_back(i);
		//�ж��Ƿ�Ϊż��
		if(cl.size()%2!=0)
			cl.pop_back();
		//��ʼ����
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
	//��������
	virtual void MutationOperator(){
		//cout<<"�������ӣ����ڸ�����Ⱥ...\n";
		int i,j;
		double p;
		//��ȡÿһ������
		vector<individual<T> >::iterator iti;
		for (iti=generationbox.begin();iti!=generationbox.end();iti++){
			//��ȡ�����ÿһ������
			individual<T> temp=*iti;
			vector<T>::iterator itj=temp.chromosome.begin();
			for(;itj!=temp.chromosome.end();itj++){
				//�ж��Ƿ����
				p=rand()%1000/1000.0;
				if(p<p_mutation){
					//��ʼ����
					//cout<<"��������"<<endl;
					InitChromosome(*itj);
					//cout<<"�����ݱ���Ϊ"<<*itj<<endl;
				} 
				//else
					//cout<<"����������"<<endl;
			}
			//�Ƿ񱣴����
			if(IsOk(temp))
				*iti=temp;
		}
		return;
	}
	//�������չʾ����
	virtual void OutputPrint(){
		out<<bestone.fitness<<endl;
/*
		cout<<"����չʾ��"<<cur_generation<<"��...\n";
		cout<<"�������Ÿ���Ϊ:\n";
		cur_bestone.print();
		out<<bestone.val<<endl;
		cout<<"����ֵΪ��"<<cur_bestone.val<<endl;
		cout<<"��Ӧ��Ϊ��"<<cur_bestone.fitness<<endl;
		cout<<"����ƽ����Ӧ�ȣ�";
		double sum=0;
		for(vector<individual<T> >::iterator it=generationbox.begin();it!=generationbox.end();it++)
			sum+=it->fitness;
		cout<<sum/(double)sizeof_generation<<endl;
		cout<<"��ʷ���Ÿ���Ϊ:\n";
		bestone.print();
		cout<<"����ֵΪ��"<<bestone.val<<endl;
		cout<<"��Ӧ��Ϊ��"<<bestone.fitness<<endl;*/
	}
	//������մ�
	virtual void ans(){
		out<<endl<<endl<<"��ʷ���Ÿ���Ϊ:\n";
		bestone.write();
		out<<"���ô�����"<<cur_generation<<endl;
		out<<"��Ӧֵ��"<<bestone.fitness<<endl;
		out<<"����ֵ��"<<bestone.val<<endl;
		out.close();
	}

	//��ʼ�������еĻ���
	virtual void InitChromosome(T &)=0;
	//���㺯��ֵ
	virtual void CalculateVal(individual<T>&)=0;
	//������Ӧ��
	virtual void CalculateFitnessValue()=0;
	//��ֹ����
	virtual bool IsEnd()=0;
	//�Ƿ�Ϊ��Ч��
	virtual bool IsOk(individual<T>cur){return 1;}
public:
	//Ⱦɫ������
	vector<individual<T> >generationbox;
	//��Ѹ���
	individual<T> cur_bestone;
	//������
	individual<T> cur_worstone;
	//��ʷ��Ѹ���
	individual<T> bestone;
	//��ʷ������
	individual<T> worstone;

	//���ݳ�Ա
	int sizeof_generation;//��Ⱥ��С
	int sizeof_chrom;//�������
	int max_generation;//����Ŵ�������
	int cur_generation;//��ǰ������
	double p_crossover;//�������
	double p_mutation;//�������
	int best_index;//��Ѹ���λ��
	int worst_index;//������λ��
};


class CalcMaxVal:public GAalg<double>{
public:
	CalcMaxVal(int sg,int sc,int mg,double pc,double pm):GAalg(sg,sc,mg,pc,pm){
		//innit the class
	}
	//��ʼ�������еĻ���
	void InitChromosome(double &my){
		my=rand()%20100/100.0-100; 
		//cout<<"�����ʼֵΪ"<<my<<endl;
	}
	//������Ӧ��
	void CalculateFitnessValue(){
		//��������iterator
		vector<individual<double> >::iterator it=generationbox.begin();
		//����ֵ
		for (it=generationbox.begin();it!=generationbox.end();it++){
			CalculateVal(*it);
				it->fitness=pow(1/(it->val),8);
			_sum++;
		}
	}
	//���㺯��ֵ
	void CalculateVal(individual<double>&my){
		vector<double>::iterator it=my.chromosome.begin();
		double ans=0;
		for(;it!=my.chromosome.end();it++){
			ans+=(*it)*(*it);
		}
		my.val=ans;
		return;
	}
	//��ֹ����
	bool IsEnd(){
		if(cur_generation>=max_generation || _sum>_max_sum)
			return 1;
		else
			return 0;
	}
	//�Ƿ�Ϊ��Ч��
	bool IsOk(individual<double> cur)
	{
		return 1;
	}
};

int main(){
	int sg,sc,mg;//��Ⱥ��������������������� 
	double pc,pm;//������ʣ��������
	/*printf("��ʼ��ȫ�ֱ���:\n");
	printf("\t��Ⱥ��С:");
	cin>>sg;  
	if((sg%2) != 0){
		printf("��Ⱥ��С������Ϊż��\n");
		sg++;
	};
	printf("\t�����������");
	cin>>mg;
	printf("\t�����ʣ�");
	cin>>pc;
	printf("\t�����ʣ�");
	cin>>pm;
	printf("\t�������������");
	cin>>_max_sum;*/
	_max_sum=10000000000;
	CalcMaxVal a(30,30,800000,0.60,0.005);//��Ⱥ�����������������������������ʣ�������� 
	a.GenerateInit();
	//a.print();
	cout<<"��ʼ�Ŵ��㷨��"<<endl;
	a.ValuatePopulation();
	while(!a.IsEnd()){
		a.NextPopulation();
		a.ValuatePopulation();
		a.OutputPrint();
	}
	a.ans();
}