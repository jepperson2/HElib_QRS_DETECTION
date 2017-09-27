#include "qrs_detection.h"

int main(int argc, char **argv)
{
    print_banner("Beginning QRS detection");
	bool verbose = false;

	Timing t_all("Overall");
	t_all.start();
	Errors e("test");
    
    // read files from wfdb generated sample file, read values from 1st column. See helper_functions for details
    string filename = "207_10s";
    string sample_path = "MIT_BIH_Records/" + filename + ".txt";
    string annotation_path = "MIT_BIH_Records/ANNOTATIONS/" + filename + "_ann.txt";

    vector<double> raw_samples = get_samples_from_file(sample_path, 1, verbose); 
    vector<int> annotations = get_annotations_from_file(annotation_path, verbose); 
    QRS_Detection qrs_det(raw_samples, annotations, 360, verbose);

	e = qrs_det.test_all();
	e.display();

    t_all.end();
    print_banner("Finished QRS detection");

    return 0;
}
