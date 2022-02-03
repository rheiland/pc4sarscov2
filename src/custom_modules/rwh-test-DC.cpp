// self-contained test of DC_history
//
// Using whatever your OpenMP-enabled compiler is, e.g.:
// g++-10 -O3 -fomit-frame-pointer -mfpmath=both -fopenmp -m64 -std=c++11 test_DC_history.cpp -o test_DC_history
//
// clang++ -g -fomit-frame-pointer -Xpreprocessor -fopenmp -m64 -std=c++11 test_DC_history.cpp -o test_DC_history -L/usr/local/opt/libomp/lib -lomp

// Usage: <max_cells> <max_dt_iters> <# threads>
//
// ~/git/dev-MG/PhysiCell/custom_modules$ time test_DC_history 100 1000 2
// real	0m3.669s
// user	0m7.234s
// sys	0m0.055s
//

// #include "./DC_history.h" 
#include <omp.h>
#include <sstream>

#include <iostream> 
#include <vector> 
#include <random> 
#include <algorithm> 


class Phenotype
{
 public:
	bool dead; 
};

class Cell
{
 public:
    int type;
    int index;
	// Phenotype phenotype; 
    bool dead;

    Cell();
};

Cell::Cell()
{
	type = 13; 
    dead = false;
}


int DC_type = 13;
std::vector<int> history(144000);
double DCAMOUNT = 0; //counter
std::vector<Cell*> all_cells[10000];


//rwh
std::random_device rd;
std::mt19937 gen(rd());

double UniformRandom()
{
	return std::generate_canonical<double, 10>(gen);
}

// void DC_history_model( Cell* pCell, Phenotype& phenotype, double dt )  //rwh: dt never used?
void DC_history_model( Cell* pCell )
{
	static double DCprob = 0.0000033;

	// if not DC, do nothing 
	if( pCell->type != DC_type )  //rwh: do in calling function
	{ return; } 
	
    //rwh: may be faster to reverse these 2 tests, if the 2nd happens less often
	// if( pCell->custom_data["activated_immune_cell"] >  0.5 && UniformRandom() < DCprob)

	if( UniformRandom() < DCprob)
	{
		// (Adrianne) DC leaves the tissue and so we lyse that DC
		std::cout<<"DC leaves tissue"<<std::endl;
		// pCell->lyse_cell(); 

		//rwh: could avoid 'critical' and do a thread-dependent (array) counter, then sum at end.
		#pragma omp critical   
		{ DCAMOUNT++; } // add one
		return;
	}
	return; 
}

// void DC_history_main_model( double dt )
int main(int argc, char *argv[])
{
    int max_dt_iters = 100;
    int max_cells = 10;
    int num_threads = 1;

    std::cout << "argc = " << argc << std::endl;
    if (argc < 4)
    {
        std::cout << "usage: provide <max_cells> <max_dt_iters> <num_threads>" << std::endl;
        std::exit(1);
    }
    else  // parse the args
    {
        std::istringstream iss1( argv[1] );
        if (iss1 >> max_cells)
        {
            std::cout << "max_cells = " << max_cells << std::endl;  // total agents: 2871
        }
        std::istringstream iss2( argv[2] );
        if (iss2 >> max_dt_iters)
        {
            std::cout << "max_dt_iters= " << max_dt_iters<< std::endl;
        }
        std::istringstream iss3( argv[3] );
        if (iss3 >> num_threads)
        {
            std::cout << "num_threads= " << num_threads<< std::endl;
        }
    }

	gen.seed(0.0);

    omp_set_num_threads(num_threads);

    Cell* pNew;
    for (int idx=0; idx < max_cells; idx++)  //rwh: create cells
    {
        // std::cout << "--- create cell # " << idx << std::endl;
        pNew = new Cell;
	    (*all_cells).push_back( pNew ); 
	    pNew->index = (*all_cells).size()-1;
    }

    for (int itime=0; itime < max_dt_iters; itime++)  //rwh: simulate being called every dt_diffusion
    {
        #pragma omp parallel for 
        for( int n=0; n < (*all_cells).size() ; n++ )
        {
            Cell* pC = (*all_cells)[n]; 
            // if( pC->phenotype.death.dead == false )
            if( pC->dead == false )  //rwh
            // { DC_history_model( pC, pC->phenotype , dt ); }  //rwh: faster to not incur function call overhead; do in-place
            { DC_history_model( pC ); }  //rwh: faster to not incur function call overhead; do in-place
        }
        std::rotate(history.rbegin(),history.rbegin()+1,history.rend());
        history.front() = DCAMOUNT;
        
        // std::copy(history.begin(), history.end(), std::ostream_iterator<int>(std::cout, " "));
        // std::cout << std::endl; 
    }
    std::cout << "DCAMOUNT = " << DCAMOUNT<< std::endl; 
	return 0; 
}
