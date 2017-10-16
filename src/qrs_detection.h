#ifndef QRS_DETECTION_H
#define QRS_DETECTION_H
#include "helper_functions.h"
#include "he.h"
#include <bitset>
#include <boost/circular_buffer.hpp>

class QRS_Detection {
	public:
        QRS_Detection(vector<long> digital_ecgs, vector<int> anns, int sampling_frequency, bool dbg);
		Errors test_all(int file_index);

	private:
        
/***************************************************************************/
/************************ Classes / Helper Methods *************************/
/***************************************************************************/

        // start the timer
		void t_start();                     
        // end the timer
		void t_end(string name);            
        // set parameters for BGV cryptosystem
		void set_params();                  
        // prepare for calculations
		void initialize();                  
        // copy contents of input vector<vector<mkt>> to destination
		void make_copies(vector< vector<mkt> > input, vector< vector<mkt> > destination); 
        // return the name of the current class
		const string className();


        // instance of HE, our API for homomorphic binary circuits 
		HE he;                      
        // instance of Timing, a class used for...timing
		Timing t;
        // instance of Conversion, a class used for converting between data types
		Conversion conv;
        
/***************************************************************************/
/*************************** Plaintext Variables ***************************/
/***************************************************************************/

        // samples for use in computations (ADC values)
        vector<long> samples;       
        // annotations of qrs locations by doctors - used to test accuracy of algorithms
        vector<int> annotations;    
        // sampling frequency
        double fs;                     
        // set whether debugging information should be printed 
		bool debug;


        // minimum distance away from considered sample = 0.027*fs
        int a;                      
        // maximum distance away from considered sample = 0.063*fs
        int b;                      
        // number of samples processed at one time
        int n_considered;           
        // size of "window" of values being compared to considered sample
        int lr_size;                
        // total number of samples
		unsigned n_samples;         


        // possible "widths" between two considered samples. From 1/a to 1/b
        vector<double> sample_difference_widths;
        // false negatives reported by Wang in paper
        vector<int> false_negatives;
        // false positives reported by Wang in paper
        vector<int> false_positives;

        // Theta_diff
        double diff_threshold;      
        // Theta_min
        double min_threshold;       


        // s_ave = s_diff_max for peak. s_aves holds last 8 s_ave values
        boost::circular_buffer<double> s_aves;  
        // h_cur values for last 8 peaks  
        boost::circular_buffer<double> h_aves;  
        // h_cur = height of current sample in mV (converted from ADU) 
        double h_cur;                           
        // s_diff_max = max of slope differences from left and right sides of considered sample
        double s_diff_max;


        // locations of all potential QRS locations by sample number
        vector<int> qrs_locations;

/***************************************************************************/
/****************************** FHE Variables ******************************/
/***************************************************************************/

        // number of bits needed to represent each scaled sample
        unsigned bits;
        // parameters used by HE to run BGV
		key_params params;
        // number of slots used in each Ciphertext (computed by HE's keyGen)
		long nslots;


        // low threshold for encrypted_diff_threshold. hard coded so we don't have to recompute every update
        vector<mkt> encrypted_low_threshold;
        // mid threshold for encrypted_diff_threshold. hard coded so we don't have to recompute every update
        vector<mkt> encrypted_mid_threshold;
        // high threshold for encrypted_diff_threshold. hard coded so we don't have to recompute every update
        vector<mkt> encrypted_high_threshold;

        // low S_ave value for comparing against current S_ave to see if we need to update diff_threshold
        mkt encrypted_low_ave;
        // low S_ave value for comparing against current S_ave to see if we need to update diff_threshold
        mkt encrypted_high_ave;

        // Theta_diff encrypted
        vector<mkt> encrypted_diff_threshold;
        // Theta_min encrypted
        mkt encrypted_min_threshold;


        // s_ave = s_diff_max for peak. s_aves holds last 8 encrypted s_ave values
        boost::circular_buffer<mkt> encrypted_s_aves;
        // encrypted h_cur values for last 8 peaks
        boost::circular_buffer<mkt> encrypted_h_aves;
        // encrypted h_cur
        mkt encrypted_h_cur;
        // encrypted s_diff_max
        mkt encrypted_s_diff_max;


        // least common multiple of a, a+1, a+2, ..., b
        // Note: If fs != 360, then likely a != 10 and b != 23
        // This means lcm should be recomputed to be the minimum value such that 
        // lcm / a, lcm / (a+1), ... , lcm / b  all yield whole numbers 
        long lcm = 5354228880;
        // how much to scale samples by so that calculations always done on whole numbers. 
        vector < long > scaling_factors;   


        // samples scaled by lcm 
        vector < long > samples_lcm;         
        // samples scaled by lcm / 10
        vector < long > samples_lcm_10;      
        // samples scaled by lcm / 11
        vector < long > samples_lcm_11;      
        // ...
        vector < long > samples_lcm_12;      
        vector < long > samples_lcm_13;
        vector < long > samples_lcm_14;
        vector < long > samples_lcm_15;
        vector < long > samples_lcm_16;
        vector < long > samples_lcm_17;
        vector < long > samples_lcm_18;
        vector < long > samples_lcm_19;
        vector < long > samples_lcm_20;
        vector < long > samples_lcm_21;
        vector < long > samples_lcm_22;
        // samples scaled by lcm / 23
        vector < long > samples_lcm_23;      
        // Note: if fs != 360, compute a, b and change above to samples_lcm_a, ..., samples_lcm_b


        // samples scaled by lcm converted to binary 
		vector < vector<long> > samples_lcm_bits;        
        // samples scaled by lcm / 10 converted to binary
		vector < vector<long> > samples_lcm_10_bits;    
        // ...
		vector < vector<long> > samples_lcm_11_bits;     
		vector < vector<long> > samples_lcm_12_bits;
		vector < vector<long> > samples_lcm_13_bits;
		vector < vector<long> > samples_lcm_14_bits;
		vector < vector<long> > samples_lcm_15_bits;
		vector < vector<long> > samples_lcm_16_bits;
		vector < vector<long> > samples_lcm_17_bits;
		vector < vector<long> > samples_lcm_18_bits;
		vector < vector<long> > samples_lcm_19_bits;
		vector < vector<long> > samples_lcm_20_bits;
		vector < vector<long> > samples_lcm_21_bits;
		vector < vector<long> > samples_lcm_22_bits;
        // samples scaled by lcm / 23 converted to binary 
		vector < vector<long> > samples_lcm_23_bits;     


        // samples_lcm encrypted - a vector of keys (see helper_functions.h)
		vector < mkt > k_constant_lcm;        
        // samples_lcm_10 encrypted
		vector < mkt > k_constant_lcm_10;     
        // ...
		vector < mkt > k_constant_lcm_11;     
		vector < mkt > k_constant_lcm_12;
		vector < mkt > k_constant_lcm_13;
		vector < mkt > k_constant_lcm_14;
		vector < mkt > k_constant_lcm_15;
		vector < mkt > k_constant_lcm_16;
		vector < mkt > k_constant_lcm_17;
		vector < mkt > k_constant_lcm_18;
		vector < mkt > k_constant_lcm_19;
		vector < mkt > k_constant_lcm_20;
		vector < mkt > k_constant_lcm_21;
		vector < mkt > k_constant_lcm_22;
        // samples_lcm_23 encrypted
		vector < mkt > k_constant_lcm_23;     
        

        //variable passed to binary circuits
		vector< vector < long > > inputs;          
        //inputs converted to bitsets to be passed to circuits
		vector< vector < vector<long> > >  v_in;   
        //encrypted values passed to circuits
		vector< vector < mkt > > k, k_constant;    


        vector<mkt> encrypted_qrs_locations;


/***************************************************************************/
/************************** FHE Dualslope Methods **************************/
/***************************************************************************/

        void prepare_data(int iteration); 
        vector< vector<mkt> > compute_lr_slopes(vector< vector<mkt> > encrypted_pairs, bool simd);
        vector<mkt> compute_mins_maxs(vector<mkt> encrypted_list);
        vector<mkt> compute_diff_maxs(vector<mkt> encrypted_mins_maxs);
        vector<mkt> compare_to_thresholds(vector<mkt> diff_maxs);
        mkt check_peak_closeness(vector<mkt> peaks);
        vector<mkt> update_thresholds(vector<mkt> diff_maxs); 

/***************************************************************************/
/*********************** Plaintext Dualslope Methods ***********************/
/***************************************************************************/

        vector< vector<double> > compute_lr_slopes(int index);
        vector< vector<double> > compute_mins_maxs(vector< vector<double> > plain_list);
        bool compare_to_thresholds(vector< vector<double> > mins_maxs, int index);
        bool check_local_extremes(int index);
        void update_thresholds(); 

/***************************************************************************/
/************************** Dualslope Algorithms ***************************/
/***************************************************************************/

        void ds_fhe();
        void ds_fhe(int iteration);
        void ds_unpacked_fhe();
        void ds_unpacked_fhe(int iteration);
        void ds_plain();
        void ds_plain(int index);

/***************************************************************************/
/******************************** Testing **********************************/
/***************************************************************************/

        bool test_ds_fhe(int file_index);
        bool test_ds_unpacked_fhe(int file_index);
        bool test_ds_plain(int file_index);
};
#endif
