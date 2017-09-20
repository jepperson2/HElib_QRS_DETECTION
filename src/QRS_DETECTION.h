#ifndef QRS_DETECTION_H
#define QRS_DETECTION_H
#include "helper_functions.h"
#include "he.h"
#include <bitset>
#include <stdlib.h> //rand

class QRS_DETECTION {
	public:
		QRS_DETECTION(vector<double> ecgs, int sampling_frequency, bool dbg);
		Errors test();

	private:
		const string className();
		bool debug;
		HE he;
		Timing t;
		Conversion conv;

        int fs;                 // Sampling Frequency
        vector<double> samples;
        vector<double> sample_difference_widths;

        int a;                  // minimum distance away from considered sample = 0.027*fs
        int b;                  // maximum distance away from considered sample = 0.063*fs
        int n_samples;          // number of samples processed at one time
        int lr_size;            // size of "window" of values being compared to considered sample

        double diff_threshold;  // Theta_diff
        double min_threshold;   // Theta_min
        double avg_height;      // H_ave 

        unsigned bits; 
		unsigned N_numbers;
		long nslots;
		key_params params;
		vector< vector < long > > inputs;
		vector< vector < vector<long> > >  v_in;
		vector< vector < mkt > > k_constant, k;
		void t_start();
		void t_end(string name);
		void set_params();
		void initialize();
		void make_copies();

		bool test_RCADDER();
        bool test_RBSUBER();
        bool test_RCMP();
		bool test_DUALSLOPE();
        bool test_REQ();
        bool test_SHIFTR();
        bool test_SHIFTL();
        bool test_NMUX();
};
#endif
