#include <iostream>
#include <cstdlib>
#include <random>
#include <string>
#include <vector>
#include <fstream>
#include <math.h>
using namespace std;

vector<vector<bool> > population;
vector<double> norm_fitness;
vector<double> running_tot;
vector<double> fitness;

vector<double> gen_fit_avg;
vector<double> gen_fit_max;
vector<int> gen_count_max;

void Initialize_Population(int l,int N) {
	vector<bool> individual(l);
	
	random_device rd;
	mt19937 gen(rd());
	bernoulli_distribution d(0.7);

	for(int h=0; h<N; h++) {
		for(int i=0; i<l; i++)
			individual[i] = d(gen);
		
		population.push_back(individual);
	}
}

void Calculate_Fitness(int l, int N) {
	double tot2 = 0.0;
	double tot = 0.0;
	vector<int> tmp;
	double power;
	double norm;
	double fit;
	int val;
	
	fitness.clear();
	//calculate integer value for each individual
	for(int i=0; i<N; i++) {
		val = 0;
		for(int g=0; g<l;g++) {
			if(population[i][g] == 1){
				power = (l - (g+1.0));
				val += pow(2.0,power);
			}
		}
		tmp.push_back(val);
	}

	//calculate fitness for each individual
	for(int i=0; i<N; i++) {
		fit = pow((tmp[i] / (pow(2.0,(double)l))),10.0);
		fitness.push_back(fit);
		tot += fit;
	}
	
	//normalize fitness values for each individual
	for(int i=0; i<N; i++) {
		norm = (fitness[i] / tot);
		tot2 += norm;

		norm_fitness.push_back(norm);
		running_tot.push_back(tot2);
	}
}

double Bounded_Random(double x0, double x1) {
	return (x0 + (x1 - x0) * rand() / ((double)RAND_MAX));
}

bool Is_Between(double low, double high, double x) {
	return ((low <= x)&&(x <= high));
}

void Find_Parents(int *p1, int *p2, int N){
	double r1,r2;
	
	//two random numbers
	r1 = Bounded_Random(0,1);
	r2 = Bounded_Random(0,1);

	//find parent 1
	for(int i=0; i<N; i++) {
		if(Is_Between(running_tot[i],running_tot[i+1],r1)) {
			(*p1) = (i+1);
			break;
		}
	}

	//find parent 2
	for(int i=0; i<N; i++) {
		if(Is_Between(running_tot[i],running_tot[i+1],r2)) {
			if((*p1) == (i+1)) {
				r2 = Bounded_Random(0,1);
				continue;
			}else{
				(*p2) = i+1;
				break;
			}
		}
	}
}

void Process_Offspring(int l, int N, double Pc, double Pm){
	cout<<"\n	starting..";
	random_device rd;
	mt19937 gen(rd());
	uniform_int_distribution<> dis(0,(l-1));

	vector<vector<bool> >tmp;
	vector<bool> c1;
	vector<bool> c2;
	double point;
	double r1;
	int p1,p2;

	//	cout<<"\n	loop O' DOOM";
	for(int i=0; i<(N/2); i++) {

		//get parents
		//		cout<<"\n		finding parents ("<<i<<")...";
		Find_Parents(&p1, &p2, N);
		//		cout<<"done p1: "<<p1<<" p2: "<<p2;
		r1 = Bounded_Random(0,1);
		//		cout<<"\nr1: "<<r1<<" Pc: "<<Pc;

		//		cout<<"\n		crossover ("<<i<<")...";
		if(r1 < Pc){
			//			cout<<"we crossing!\n";
			//cross
			point = dis(gen);
			for(int i=0; i<l; i++){
				if(i < point){
					c1.push_back(population[p1][i]);
					c2.push_back(population[p2][i]);
				}else{
					c1.push_back(population[p2][i]);
					c2.push_back(population[p1][i]);
				}
			}

		}else{
			//			cout<<"we not crossing!\n";
			//dont cross - copy over bits
			for(int i=0; i<l; i++) {
				c1.push_back(population[p1][i]);
				c2.push_back(population[p2][i]);
			}
		}

		//		cout<<"done\n		mutate 1 ("<<i<<")...";
		//mutate child 1
		for(int i=0; i<l; i++) {
			r1 = Bounded_Random(0,1);
			if(r1 < Pm){
				if(c1[i] == 0)
					c1[i] = 1;
				else
					c1[i] = 0;
			}
		}

		//		cout<<"done\n		mutate 2 ("<<i<<")...";
		//mutate child 2
		for(int i=0; i<l; i++) {
			r1 = Bounded_Random(0,1);
			if(r1 < Pm){
				if(c2[i] == 0)
					c2[i] = 1;
				else
					c2[i] = 0;
			}
		}

		tmp.push_back(c1);
		tmp.push_back(c2);

		c1.clear();
		c2.clear();
	}

	//update population
	population.clear();
	cout<<"pop cleared: "<<population.size();
	for(int i=0; i<N; i++)
		population.push_back(tmp[i]);
	cout<<"pop nsize: "<<population.size();
}

void Process_Stats(int l, int N) {
	double sum = 0.0;
	double max = 0.0;
	double avg;
	int count;
	int maxC = 0;

	for(int i=0; i<N; i++) {
		sum += fitness[i];
		if(fitness[i] > max)
			max = fitness[i];
	}
	//pushback avg and max
	gen_fit_avg.push_back(sum/N);
	gen_fit_max.push_back(max);

	for(int i=0; i<N; i++){
		count = 0;
		for(int g=0; g<l; g++){
			if(population[i][g] == 1)
				count++;
		}
		if(count > maxC)
			maxC = count;
	}
	//pushback maxC	
	gen_count_max.push_back(maxC);
}

//print generational stats to external file
void Print_GenStats(int G) {
	string file = "data.csv";
	ofstream data;
	data.open(file,ios_base::app);

	for(int g = 0; g < G; g++)
		data<<g<<","<<gen_fit_avg[g]<<","<<gen_fit_max[g]<<","<<gen_count_max[g]<<endl;
}

void Print_Population(int l, int N, int c) {
	cout<<"generation "<<c<<endl;
	cout<<"pop: "<<population.size()<<endl;
	for(int i=0; i<N; i++){
		cout<<i<<"   ";
		for(int g=0; g<l; g++){
			cout<<population[i][g];
		}
		cout<<"   "<<fitness[i]<<endl;
	}
	cout<<endl;
}

int main(int argc, char *argv[]) {

	//argument check real quick
	if(argc != 6){
		cout<<"error: a.out [# genes] [pop size] [# gen] [mut p] [cros p]\n";
		return 1;
	}
	
	srand(time(0));
	int l = atoi(argv[1]);
	int N = atoi(argv[2]);
	int G = atoi(argv[3]);
	double Pm = atof(argv[4]);
	double Pc = atof(argv[5]);

//	cout<<"N: "<<N<<" l: "<<l<<endl;

//	cout<<"initializing...";
	Initialize_Population(l,N);
//	cout<<"done\nCalculating fitness...";
	Calculate_Fitness(l,N);
//	cout<<"done\nProcessing Stats...";
	Process_Stats(l,N);
//	cout<<"done\n";
	Print_Population(l,N,0);
	
	for(int gen=1; gen < G; gen++) {
		
//		cout<<"processing offspring (gen "<<gen<<")...";
		Process_Offspring(l,N,Pc,Pm);
//		cout<<"done\nCalculating fitness...";
		Calculate_Fitness(l,N);
//		cout<<"done\nProcessing Stats...";
		Process_Stats(l,N);
		cout<<"done\n";
//		Print_Population(l,N,gen);
	}
	cout<<"Printing data...";
	Print_GenStats(G);
	cout<<"done\nexiting..\n";

	return 0;
}
