#include "helper_functions.h"
#include "qrs_detection.h"

int main(int argc, char **argv)
{
    print_banner("Beginning QRS detection");
	bool verbose = false;

	Timing t_all("Overall");
	t_all.start();
	Errors e("test");
	
    vector<double> raw_samples = get_samples_from_file("MIT_BIH_Records/100_10s.txt", 1, verbose); // read files from wfdb generated sample file, read values from 1st column. See helper_functions for details 
	vector<long> samples = scale_samples(raw_samples, 1000); // scale samples to be integers. If precision != 3, adjust "1000" to be 10^precision
   
    QRS_Detection qrs_det(samples, 360, verbose);

	e = qrs_det.test_all();
	e.display();

    t_all.end();
    print_banner("Finished QRS detection");

    return 0;
}
