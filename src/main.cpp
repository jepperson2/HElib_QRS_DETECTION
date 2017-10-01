#include "qrs_detection.h"

int main(int argc, char **argv)
{
    print_banner("Beginning QRS detection");
	bool verbose = false;

	Timing t_all("Overall");
	t_all.start();
	Errors e("test");
    
    string filename = "207_10s";

    string mV_sample_path = "MIT_BIH_Records/" + filename + ".txt";
    string digital_sample_path = "MIT_BIH_Records/DIGITAL/" + filename + ".txt";
    string annotation_path = "MIT_BIH_Records/ANNOTATIONS/" + filename + ".txt";

    vector<double> mV_samples = get_mV_samples(mV_sample_path, 1, 2, verbose);
    // scale samples to be integers. If precision != 3, adjust "1000" to be 10^precision
    vector<long> scaled_mV_samples = scale_samples(mV_samples, 1000); 

    vector<long> digital_samples = get_digital_samples(digital_sample_path, 1, 1, verbose);
    vector<int> annotations = get_annotations(annotation_path, verbose); 
/*
    boost::circular_buffer<int> cb(8);
    
    for (int i = 0; i < 12; i++){
        cb.push_back(i);
        int sum = std::accumulate(cb.begin(), cb.end(), 0);
        cout << "sum = " << sum << endl;
    }
*/
    QRS_Detection qrs_det(digital_samples, annotations, 360, verbose);
	e = qrs_det.test_all();
	e.display();

    t_all.end();
    print_banner("Finished QRS detection");

    return 0;
}
