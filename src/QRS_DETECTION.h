#ifndef QRS_DETECTION_H
#define QRS_DETECTION_H
#include "helper_functions.h"
#include "he.h"
#include <bitset>
#include <stdlib.h> //rand

class QRS_Detection {
	public:
		QRS_Detection(vector<long> ecgs, int sampling_frequency, bool dbg);
		Errors test_all();

	private:
		const string className();
		bool debug;
		HE he;                  // instance of HE, our API for homomorphic binary circuits 
		Timing t;
		Conversion conv;

        int fs;                 // Sampling Frequency
        vector<long> samples;
        vector<double> sample_difference_widths;

        int a;                  // minimum distance away from considered sample = 0.027*fs
        int b;                  // maximum distance away from considered sample = 0.063*fs
        int n_considered;       // number of samples processed at one time
        int lr_size;            // size of "window" of values being compared to considered sample

        double diff_threshold;  // Theta_diff
        double min_threshold;   // Theta_min
        double avg_height;      // H_ave 

        vector<mkt> encrypted_thresholds; // Encrypted thresholds

        unsigned bits;          // number of bits in each sample
		unsigned n_samples;     // total number of samples
		long nslots;            // number of slots used in each Ciphertext (computed internally)
		key_params params;

        vector < long > scaling_factors;   // how much to scale samples by so that calculations always done on whole numbers. 

        vector < long > samples_x;         // samples scaled by x (which is set in initialize()) 
        vector < long > samples_x_10;      // samples scaled by x / 10
        vector < long > samples_x_11;      // samples scaled by x / 11
        vector < long > samples_x_12;      // ...
        vector < long > samples_x_13;
        vector < long > samples_x_14;
        vector < long > samples_x_15;
        vector < long > samples_x_16;
        vector < long > samples_x_17;
        vector < long > samples_x_18;
        vector < long > samples_x_19;
        vector < long > samples_x_20;
        vector < long > samples_x_21;
        vector < long > samples_x_22;
        vector < long > samples_x_23;      // samples scaled by x / 23

		vector < vector<long> > samples_bits_x;        // samples scaled by x converted to binary 
		vector < vector<long> > samples_bits_x_10;     // samples scaled by x / 10 converted to binary
		vector < vector<long> > samples_bits_x_11;     // ...
		vector < vector<long> > samples_bits_x_12;
		vector < vector<long> > samples_bits_x_13;
		vector < vector<long> > samples_bits_x_14;
		vector < vector<long> > samples_bits_x_15;
		vector < vector<long> > samples_bits_x_16;
		vector < vector<long> > samples_bits_x_17;
		vector < vector<long> > samples_bits_x_18;
		vector < vector<long> > samples_bits_x_19;
		vector < vector<long> > samples_bits_x_20;
		vector < vector<long> > samples_bits_x_21;
		vector < vector<long> > samples_bits_x_22;
		vector < vector<long> > samples_bits_x_23;     // samples scaled by x / 23 converted to binary 

		vector < mkt > k_constant_x;        // samples_x encrypted - a vector of keys (see helper_functions.h)
		vector < mkt > k_constant_x_10;     // samples_x_10 encrypted
		vector < mkt > k_constant_x_11;     // ...
		vector < mkt > k_constant_x_12;
		vector < mkt > k_constant_x_13;
		vector < mkt > k_constant_x_14;
		vector < mkt > k_constant_x_15;
		vector < mkt > k_constant_x_16;
		vector < mkt > k_constant_x_17;
		vector < mkt > k_constant_x_18;
		vector < mkt > k_constant_x_19;
		vector < mkt > k_constant_x_20;
		vector < mkt > k_constant_x_21;
		vector < mkt > k_constant_x_22;
		vector < mkt > k_constant_x_23;     // samples_x_23 encrypted

		vector< vector < long > > inputs;          //variable passed to binary circuits
		vector< vector < vector<long> > >  v_in;   //inputs converted to bitsets to be passed to circuits
		vector< vector < mkt > > k, k_constant;    //encrypted values passed to circuits

		void t_start();                     // start the timer
		void t_end(string name);            // end the timer
		void set_params();                  // set parameters for BGV cryptosystem
		void initialize();                  // prepare for calculations
		void make_copies(vector< vector<mkt> > input, vector< vector<mkt> > destination); // copy contents of one vector<vector<mkt>> to another one

        void prepare_data(vector<long> sample_subset); // prepare data for homomorphic operations. 

        // Components of Dualslope algorithm broken into pieces (operating on FHE encrypted data)
        void compute_lr_slopes(vector< vector<mkt> > encrypted_pairs, bool simd);
        vector<mkt> compute_min_max(vector<mkt> encrypted_list);
        vector<mkt> compute_diff_max(vector<mkt> encrypted_mins_maxs);
        vector<mkt> compare_to_thresholds(vector<mkt> diff_maxs);
        mkt check_peak_closeness(vector<mkt> peaks);
        vector<mkt> update_thresholds(vector<mkt> diff_maxs); 

        // Components of Dualslope algorithm broken into pieces (operating on plaintext data)
        void compute_lr_slopes(vector< vector<long> > plain_pairs);
        vector<long> compute_min_max(vector<long> plain_list);
        vector<long> compute_diff_max(vector<long> plain_mins_maxs);
        vector<bool> compare_to_thresholds(vector<long> diff_maxs);
        long check_peak_closeness(vector<long> peaks);
        void update_thresholds(long diff_max); 

        void ds_fhe();
        void ds_fhe(int iterations);
        void ds_unpacked_fhe();
        void ds_unpacked_fhe(int iterations);
        void ds_plain();
        void ds_plain(int iterations);

        bool test_ds_fhe();
        bool test_ds_fhe(int iterations);
        bool test_ds_unpacked_fhe();
        bool test_ds_unpacked_fhe(int iterations);
        bool test_ds_plain();
        bool test_ds_plain(int iterations);        
};
#endif
