#include "qrs_detection.h"

int main(int argc, char **argv)
{
    print_banner("Beginning QRS detection");
	bool verbose = false;
    bool file_io_verbose = false;

	Timing t_all("Overall");
	t_all.start();
	Errors e("test");
    
    vector<string> files = { "100","101","102","103","104","105","106","107","108","109","111","112","113","114","115","116","117","118","119","121","122","123","124","200","201","202","203","205","207","208","209","210","212","213","214","215","217","219","220","221","222","223","228","230","231","232","233","234" };
    for (int i = 0; i < 1; i++){
        int file_index = i;
        string filename = files[file_index];// + "_10s";

        string digital_sample_path = "MIT_BIH_Records/DIGITAL/" + filename + ".txt";
        string annotation_path = "MIT_BIH_Records/ANNOTATIONS/" + filename + ".txt";

        vector<long> digital_samples = get_digital_samples(digital_sample_path, 1, 1, file_io_verbose);
        vector<int> annotations = get_annotations(annotation_path, file_io_verbose); 

        QRS_Detection qrs_det(digital_samples, annotations, 360, verbose);
        e = qrs_det.test_all(file_index);
        e.display();
    }

    t_all.end();
    print_banner("Finished QRS detection");

    return 0;
}
