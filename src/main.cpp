#include "helper_functions.h"
#include "TEST_GATES.h"
#include "TEST_CIRC_COMB.h"
#include "TEST_CIRC_SEQ.h"
#include "TEST_CIRC_ARITHM.h"
#include "QRS_DETECTION.h"

int main(int argc, char **argv)
{
	cout << "======================================\n";
	bool verbose = true;

	Timing t_all("Overall");
	t_all.start();
	Errors e("test");
	
    vector<double> samples = get_samples_from_file("MIT_BIH_Records/100_1s.txt", 2, verbose); 

	QRS_DETECTION q_detection(samples, 360, verbose);
	e = q_detection.test();
	e.display();
	
    t_all.end();
	cout << "======================================\n";
	return 0;
}
