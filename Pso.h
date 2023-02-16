#pragma once
#include<iostream>
#include<map>
#include<vector>
#include<fstream>
#include<algorithm>
#include<functional>
#include<random>
#include<cmath>
#include<math.h>

using namespace std;

class Pso
{
private:
	int _number_of_particles;//������Ŀ
	int _number_of_dimensions;//ά����Ŀ
	int _maximum_of_iteration;//�����������ֵ

	double _upper_bound;//�����Ͻ�
	double _lower_bound;//�����½�

	double _velocity_max;//�ٶ��Ͻ�
	double _velocity_min;//�ٶ��½�

	double _w;//inertia coefficient
	double _c1;//local coefficient
	double _c2;//global coefficient

	vector<double>_particles;//����
	vector<double>_particles_velocity;//�����ٶȼ�¼

	//ȫ������������key���� value.first����ǰȫ��������Ӧ��ֵ��value.second ������ǰȫ�����ɶȵ�����
	map<int, pair<double,double>,less<int>>_global_best;


	//�ֲ���������
	//key ���ӱ��
	//value ������ʷ��¼
	//value.key ������Ӧ��ֵ
	//value.value.first ����ֵ
	//value.value.second ��
	map<int,multimap<double,pair<double,int>,less<double>>,less<int>>_individual_best;

public:
	Pso(int NOP = 1000, int NOD = 1, int MOI = 1000, double UB = 50.0, double LB = 0.0,
		double VMAX = 1, double VMIN = -1, double W = 0.4, double C1 = 0.5, double C2 = 0.2) :
		_number_of_particles(NOP),
		_number_of_dimensions(NOD),
		_maximum_of_iteration(MOI),
		_upper_bound(UB),
		_lower_bound(LB),
		_velocity_min(VMIN),
		_velocity_max(VMAX),
		_w(W),
		_c1(C1),
		_c2(C2)
	{
		this->_particles.resize(NOP, 0);
		this->_particles_velocity.resize(NOP, this->_velocity_min);
	}

	void set_nop(const int& nop);
	void set_nod(const int& nod);
	void set_moi(const int& moi);

	void set_ub(const double& ub);
	void set_lb(const double& lb);

	void set_vmin(const double& vmin);
	void set_vmax(const double& vmax);

	void set_w(const double& w);
	void set_c1(const double& c1);
	void set_c2(const double& c2);

	//��ʼ������
	void init();

	void calc_v();

	double calc_fit_value(double num);

	void particle_init();

	void individual_best_init();

	void global_best_init();

	void evolve();

	void optimize();

	void update_individual_best();

	void update_global_best();

	void print_particles();

	void print_individual_best();

	void print_global_best();

	void print_individual_results();
	
	void print_global_result();

	void save_global_best();

	void save_individual_best();
};