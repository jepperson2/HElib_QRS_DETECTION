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
/************************ Helper Methods / Classes *************************/
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
        // annotations of qrs locations - used to test accuracy of algorithms
        vector<int> annotations;    
        // Sampling Frequency
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
        // H_cur values for last 8 peaks  
        boost::circular_buffer<double> h_aves;  
        // h_cur 
        double h_cur;                           
        // s_diff_max
        double s_diff_max;


        // locations of all potential QRS locations by sample number
        vector<int> qrs_locations;

/***************************************************************************/
/****************************** FHE Variables ******************************/
/***************************************************************************/

        // Encrypted thresholds
        vector<mkt> encrypted_thresholds; 


        // number of bits in each sample
        unsigned bits;              
        // number of slots used in each Ciphertext (computed internally)
		long nslots;                
		key_params params;


        // how much to scale samples by so that calculations always done on whole numbers. 
        vector < long > scaling_factors;   


        // samples scaled by x (which is set in initialize()) 
        vector < long > samples_x;         
        // samples scaled by x / 10
        vector < long > samples_x_10;      
        // samples scaled by x / 11
        vector < long > samples_x_11;      
        // ...
        vector < long > samples_x_12;      
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
        // samples scaled by x / 23
        vector < long > samples_x_23;      


        // samples scaled by x converted to binary 
		vector < vector<long> > samples_bits_x;        
        // samples scaled by x / 10 converted to binary
		vector < vector<long> > samples_bits_x_10;    
        // ...
		vector < vector<long> > samples_bits_x_11;     
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
        // samples scaled by x / 23 converted to binary 
		vector < vector<long> > samples_bits_x_23;     


        // samples_x encrypted - a vector of keys (see helper_functions.h)
		vector < mkt > k_constant_x;        
        // samples_x_10 encrypted
		vector < mkt > k_constant_x_10;     
        // ...
		vector < mkt > k_constant_x_11;     
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
        // samples_x_23 encrypted
		vector < mkt > k_constant_x_23;     
        

        //variable passed to binary circuits
		vector< vector < long > > inputs;          
        //inputs converted to bitsets to be passed to circuits
		vector< vector < vector<long> > >  v_in;   
        //encrypted values passed to circuits
		vector< vector < mkt > > k, k_constant;    


/***************************************************************************/
/************************** FHE Dualslope Methods **************************/
/***************************************************************************/

        void prepare_data(int iteration, int leftovers); 
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
        void ds_fhe(int iterations, int leftovers);
        void ds_unpacked_fhe();
        void ds_unpacked_fhe(int iterations, int leftovers);
        void ds_plain();
        void ds_plain(int index);

/***************************************************************************/
/******************************** Testing **********************************/
/***************************************************************************/

        bool test_ds_fhe();
        bool test_ds_fhe(int iterations, int leftovers);
        bool test_ds_unpacked_fhe();
        bool test_ds_unpacked_fhe(int iterations, int leftovers);
        bool test_ds_plain(int file_index);
        bool test_ds_plain(int file_index, int sample_index);        
};
#endif
