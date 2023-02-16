#include"Pso.h"

void Pso::set_nop(const int& nop)
{
	this->_number_of_particles = nop;
}

void Pso::set_nod(const int& nod)
{
	this->_number_of_dimensions = nod;
}

void Pso::set_moi(const int& moi)
{
	this->_maximum_of_iteration = moi;
}

void Pso::set_ub(const double& ub)
{
	this->_upper_bound = ub;
}

void Pso::set_lb(const double& lb)
{
	this->_lower_bound = lb;
}

void Pso::set_vmin(const double& vmin)
{
	this->_velocity_min = vmin;
}

void Pso::set_vmax(const double& vmax)
{
	this->_velocity_max = vmax;
}

void Pso::set_w(const double& w)
{
	this->_w = w;
}

void Pso::set_c1(const double& c1)
{
	this->_c1 = c1;
}

void Pso::set_c2(const double& c2)
{
	this->_c2 = c2;
}

void Pso::init()
{
	particle_init();
	individual_best_init();
	global_best_init();
}

double Pso::calc_fit_value(double num) {
	return num * sin(num) * cos(2 * num) - 2 * (num)*sin(3 * num) + 3 * num * sin(4 * num);
}

void Pso::particle_init()
{
	cout << "Particles are initing..." << endl;
	mt19937 mersenne(static_cast<unsigned int>(time(nullptr)));
	for (vector<double>::iterator iter = this->_particles.begin();iter != this->_particles.end();iter++) {
			uniform_real_distribution<>rnd_particle_component(0,1);//�������ȷֲ�
			*iter = this->_lower_bound + (this->_upper_bound - this->_lower_bound) * (rnd_particle_component(mersenne));
		}
}

void Pso::individual_best_init()
{
	cout << "Individual best is initing..." << endl;
	int i = 1;
	this->_individual_best;
	for (vector<double>::const_iterator citer = this->_particles.begin();citer != this->_particles.end();citer++, i++) {
			multimap<double, pair<double,int>, less<double>>m;

			double fit_val = calc_fit_value(*citer);

			m.insert(make_pair(fit_val, pair<double, int>(*citer, 1)));

			this->_individual_best.insert(make_pair(i, m));
		}
}

void Pso::global_best_init()
{
	cout << "Global best is initing..." << endl;
	mt19937 mersenne(static_cast<unsigned int>(time(nullptr)));
	uniform_real_distribution<double>rnd_global(0,1);

	double g = this->_lower_bound + (this->_upper_bound - this->_lower_bound) * (rnd_global(mersenne));

	double test_value_g = g * sin(g) * cos(2 * g) - 2 * (g)*sin(3 * g) + 3 * g * sin(4 * g);//ȫ������ֵ
		//map<int,map<double,double,less<double>>,less<int>>global_best;//ȫ�����ż�¼ int ������
	pair<double, double>p(make_pair(test_value_g,g ));

	this->_global_best.insert(make_pair(1,p));
}

void Pso::calc_v()
{
	map<int, multimap<double, pair<double, int>, less<double>>, less<int>>::iterator iiter = this->_individual_best.begin();

	auto giter = this->_global_best.rbegin();
	
	vector<double>::iterator piter = this->_particles.begin();

	for (vector<double>::iterator iter = this->_particles_velocity.begin(); 
		
		iter != this->_particles_velocity.end()&&
		piter!=this->_particles.end()&&
		iiter!=this->_individual_best.end();
		
		iter++,
		piter++,
		iiter++) {

		mt19937 mersenne(static_cast<unsigned int>(time(nullptr)));
		uniform_real_distribution<double>rnd_global(0, 1);

		double c_arr[2] = { 0.0, 0.0 };

		for (int i = 0; i < 2; i++)
		{
			c_arr[i] = rnd_global(mersenne);
		}

		*iter = (*iter) * this->_w +
			this->_c1 * c_arr[0] * (iiter->second.begin()->first - *(piter)) +
			this->_c2 * c_arr[1] * ((giter->second.first)-(*piter));

		if (*iter>this->_velocity_max)
		{
			*iter = this->_velocity_max;
		}
		else if (*iter<this->_velocity_min)
		{
			*iter = this->_velocity_min;
		}
	}

}




void Pso::evolve()
{
	this->calc_v();
	vector<double>::iterator viter = this->_particles_velocity.begin();

	for (vector<double>::iterator piter = this->_particles.begin(); piter != this->_particles.end()&&viter!=this->_particles_velocity.end();piter++,viter++ ) {
		*piter += *viter;

		if (*piter>this->_upper_bound)
		{
			*piter = this->_upper_bound;
		}
		else if(*piter<this->_lower_bound)
		{
			*piter = this->_lower_bound;
		}
	}
}

void Pso::optimize()
{
	for (int i = 1; i <=this->_maximum_of_iteration; i++)
	{
		evolve();
		update_individual_best();
		update_global_best();
	}
}

void Pso::update_individual_best()//�˺�����update_global_best֮ǰ
{
	map<int, multimap<double, pair<double, int>, less<double>>, less<int>>::iterator miter = this->_individual_best.begin();

	for (vector<double>::iterator iter = this->_particles.begin();
		iter != this->_particles.end() && miter != this->_individual_best.end();
		iter++,miter++) {
		double fit_val = calc_fit_value(*iter);//�������ӱ䶯���Ӧ����Ӧ��

		pair<double, int>p = miter->second.begin()->second;//p�洢ֵ������
		p.second = miter->second.size() + 1;//�����������Ӽ�¼������++1 ÿ�ֶ������

		double current_best_fit_val = miter->second.begin()->first;//��ǰ����ֵ

		if (fit_val<miter->second.begin()->first)//����µ����ӵ���Ӧ��ֵ��С
		{
			p.first = *iter;//��������ֵ

			miter->second.insert(make_pair(fit_val, p));//����map
		}
		else
		{
			miter->second.insert(make_pair(current_best_fit_val, p));//������Ǽ�1����������ֵ����Ӧ��ֵ�������仯
		}
	}
}

void Pso::update_global_best()
{	


	//����individualbest������Ѱ����Ӧ����Сֵ

	pair<double, double>current_global_best = this->_global_best.rbegin()->second;
	
	map<int, multimap<double, pair<double, int>, less<double>>, less<int>>::iterator iter;

	for (iter=this->_individual_best.begin();
		iter!=this->_individual_best.end(); 
		iter++)
	{
		//cout << current_global_best.first << " " << iter->second.begin()->first << endl;
		if (current_global_best.first>iter->second.begin()->first)
		{
			current_global_best.first = iter->second.begin()->first;//������Ӧ��
			current_global_best.second = iter->second.begin()->second.first;//����ֵ
		}

	}
	this->_global_best.insert(make_pair(this->_global_best.size() + 1, current_global_best));
	
}

void Pso::print_particles()
{
	for (auto& p : this->_particles) {
		cout << p << " ";
	}
	cout << endl;
}

void Pso::print_individual_best()
{
	cout << "���Ӹ�������" << endl;
	cout << "���ӱ��" << ' ' << "������Ӧ��ֵ" << ' ' << "����ֵ" << ' ' << "��" << endl;
	
	for (map<int, multimap<double, pair<double, int>, less<double>>, less<int>>::iterator particle_it = this->_individual_best.begin(); particle_it != this->_individual_best.end();particle_it++) {
		for (auto& deatil : particle_it->second) {
			cout << particle_it->first <<' ';
			cout << deatil.first << ' ';
			cout << deatil.second.first << ' ';
			cout << deatil.second.second <<endl;
		}
	}
}

void Pso::print_global_best()
{
	cout << "ȫ������" << endl;
	cout << "��" << ' ' << "ȫ��������Ӧ��" << ' ' << "����ֵ" << endl;
	for (auto& rec : this->_global_best) {
		cout << rec.first << ' ';
		cout << rec.second.first << ' ';
		cout << rec.second.second << endl;
	}
}

void Pso::print_individual_results()
{
	cout << "���ӱ�� " << "������Ӧ�� " << "����ֵ" << endl;
	for (auto& c : this->_individual_best) {
		cout << c.first << " ";
		cout << c.second.begin()->first << " ";
		cout << c.second.begin()->second.first << endl;
	}
}

void Pso::print_global_result()
{
	cout <<"������Ӧ�� " << "����ֵ" << endl;
	cout << this->_global_best.rbegin()->second.first << " " << this->_global_best.rbegin()->second.second << endl;
}

void Pso::save_global_best()
{
	ofstream ofs("global_best_recs.txt", ios::out|ios::trunc);

	for (auto& rec : this->_global_best) {
		ofs<< rec.first << ' ';
		ofs << rec.second.first << ' ';
		ofs << rec.second.second << endl;
	}

	ofs.close();
}

void Pso::save_individual_best()
{
	ofstream ofs("individual_best_recs.txt", ios::out|ios::trunc);
	for (map<int, multimap<double, pair<double, int>, less<double>>, less<int>>::iterator particle_it = this->_individual_best.begin(); particle_it != this->_individual_best.end(); particle_it++) {
		for (auto& deatil : particle_it->second) {
			ofs << particle_it->first << ' ';
			ofs << deatil.first << ' ';
			ofs << deatil.second.first << ' ';
			ofs << deatil.second.second << endl;
		}
	}
	ofs.close();
}


