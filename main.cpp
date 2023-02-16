#include"Pso.h"

void test() {

	mt19937 mersenne(static_cast<unsigned int>(time(nullptr)));
	uniform_real_distribution<double>rnd_global(0, 1);

	double c_arr[2] = { 0.0, 0.0 };

	for (int i = 0; i < 2; i++)
	{
		c_arr[i] = rnd_global(mersenne);
	}
	cout << c_arr[0] << " " << c_arr[1] << endl;
}

int main() {
	Pso pso;
	pso.init();
	pso.print_particles();
	pso.print_individual_best();
	pso.print_global_best();
	pso.optimize();
	//pso.print_particles();
	//pso.print_individual_best();
	//pso.print_global_best();
	pso.print_individual_results();
	pso.print_global_result();
	pso.save_global_best();
	pso.save_individual_best();


	//test();


}