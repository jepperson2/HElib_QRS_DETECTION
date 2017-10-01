#include "qrs_detection.h"

/*TODO: 
    move stuff in initialize that prepares/encrypts data
    create methods to match .h file
    finish making methods in .h file (calls to actual algorithms)
    write calls to algorithsm to call functions defined in h files (constructor() -> loop through specific (int iterations = samples.size/(nslots-n_considered)
    rethink through parameters/outputs in .h file new methods
    implement baby example in plain
    think through security of methods in .h -> reveal info by return values? can hide. 
    implement min_max based on that one paper's method 3
    change parameters (params.L calculation) based on expected levels needed for subtration/comparison. 
    double check that comparison works on neg numbers 
    if not, scale everything...? need to make sure that works... 
    

  */
QRS_Detection::QRS_Detection(vector<long> digital_ecgs, vector<int> anns, int sampling_frequency, bool dbg){
	samples = digital_ecgs;
    annotations = anns;
    fs = sampling_frequency;
    debug = dbg;

    // set here because set_params depends on them and initialize depends on set_params
    bits = 64; 

	//set_params();
	initialize();
}

Errors QRS_Detection::test_all(){
	Errors e("QRS_Detection");
//	e.add("ds_fhe()", test_ds_fhe());
//	e.add("ds_fhe(5,200)", test_ds_fhe(5,200));
//	e.add("ds_fhe_unpacked()", test_ds_unpacked_fhe());
//	e.add("ds_fhe_unpacked(5,200)", test_ds_unpacked_fhe(5,200));
//	e.add("ds_plain()", test_ds_plain());
	e.add("ds_plain(46)", test_ds_plain(46));
	return e;
}

const string QRS_Detection::className(){
	return "QRS_Detection";
}

void QRS_Detection::t_start(){
	if (debug){
		t.start();
	}
}

void QRS_Detection::t_end(string name){
	if (debug){
		double duration = t.end("silent");
		cout << className() << ": " << name << " - " << nslots << " operations done in " << duration << "s" << endl;
		cout << className() << ": " << name << " - Time for each operation: " << 1000*duration/nslots << "ms" << endl;
	}
}

void QRS_Detection::set_params(){
	params.p = 2;
	params.r = 1;
    
	params.d = 0;   // field p^d 
	params.k = 128;
	params.slb = 800;
	params.L = 40;
	
	if ((params.L > 42) && (params.slb > 600)){
		params.L *= 1.2;
	}

	params.c = 3;
	params.w = 64;
	
	params.m = FindM(params.k,params.L,params.c,params.p,params.d,params.slb,0);
}

void QRS_Detection::initialize(){
    if (debug){
        cout << className() << ": Initializing..." << endl;
    }   

/***************************************************************************/
/******************** Initializing Plaintext Variables *********************/
/***************************************************************************/

    // set as in Wang's paper, an estimate of the min width of QRS complex
    double double_a = 0.027 * fs;   
    a = round(double_a);
    // set as in Wang's paper, an estimate of the max width of QRS complex
    double double_b = 0.063 * fs;   
    b = round(double_b);
    // number of samples considered in one iteration (goes from 0 to -2b)
    n_considered = (2 * b) + 1;     
    // convenience for looping through each side 
    lr_size = (b - a) + 1;          


    // Constant values of 1/a, 1/(a+1), ... , 1/b. These are the "run" in the rise over run slope calculation 
    double width;
    for (double i = a; i < (b + 1); i++){
        width = 1/i;
        sample_difference_widths.push_back(width);
    }
    //Total number of samples given to process 
	n_samples = samples.size(); 


    // initial diff_threshold value. to be updated later based on S_ave
    diff_threshold = 3840 / (double)fs; 
    // threshold value given in paper
    min_threshold = 1536 / (double)fs;  


    // initialize s_aves to size 8.
    s_aves = boost::circular_buffer<double>(8); 
    // initialize h_aves to size 8.
    h_aves = boost::circular_buffer<double>(8); 
    // set average height to 0. to be updated later based on sample heights
    h_cur = 0.0;                   


    if (debug){
        cout << "a,b = " << a << "," << b << endl;
        cout << "n_considered = " << n_considered << endl;
        cout << "lr_size = " << lr_size << endl;

        for (int i = 0; i < lr_size; i++){
            cout << "sample_difference_widths[" << i << "]: " << sample_difference_widths[i] << endl;
        }
    }

/***************************************************************************/
/*********************** Initializing FHE Variables ************************/
/***************************************************************************/

// TODO - initialize encrypted thresholds
/*	he.debug_on(debug);
	cout << className() << ": Number of bits n was set to " << bits << endl; 
	nslots = he.keyGen(params);
	mkt k_ones = he.setOnes(nslots);
	he.set01(k_ones);
*/
    //TEMP! PUT BACK THE ABOVE WHEN RUNNING NON_FHE TESTS
    nslots = 1024; 

    /* Compute scaling_factors based on a and b values. */
    // x = 5,354,228,880 = 2^4*3^2*5*7*11*13*17*19*23 = Minimum number divisible by 10, 11, ... , 23
    // Note: If fs != 360, then likely a != 10 and b != 23
    // This means x should be recomputed to be the minimum value such that 
    // x / a, x / (a+1), ... , x / b  all yield whole numbers 
    long x = 5354228880;

    scaling_factors.push_back(x);
    for (int i = a; i <= b; i++){
        scaling_factors.push_back(x/i);
    }

    samples_x.resize(nslots,0);
    samples_x_10.resize(nslots,0);
    samples_x_11.resize(nslots,0);
    samples_x_12.resize(nslots,0);
    samples_x_13.resize(nslots,0);
    samples_x_14.resize(nslots,0);
    samples_x_15.resize(nslots,0);
    samples_x_16.resize(nslots,0);
    samples_x_17.resize(nslots,0);
    samples_x_18.resize(nslots,0);
    samples_x_19.resize(nslots,0);
    samples_x_20.resize(nslots,0);
    samples_x_21.resize(nslots,0);
    samples_x_22.resize(nslots,0);
    samples_x_23.resize(nslots,0);

    samples_bits_x.resize(bits, vector<long> (nslots,0));
    samples_bits_x_10.resize(bits, vector<long> (nslots,0));
    samples_bits_x_11.resize(bits, vector<long> (nslots,0));
    samples_bits_x_12.resize(bits, vector<long> (nslots,0));
    samples_bits_x_13.resize(bits, vector<long> (nslots,0));
    samples_bits_x_14.resize(bits, vector<long> (nslots,0));
    samples_bits_x_15.resize(bits, vector<long> (nslots,0));
    samples_bits_x_16.resize(bits, vector<long> (nslots,0));
    samples_bits_x_17.resize(bits, vector<long> (nslots,0));
    samples_bits_x_18.resize(bits, vector<long> (nslots,0));
    samples_bits_x_19.resize(bits, vector<long> (nslots,0));
    samples_bits_x_20.resize(bits, vector<long> (nslots,0));
    samples_bits_x_21.resize(bits, vector<long> (nslots,0));
    samples_bits_x_22.resize(bits, vector<long> (nslots,0));
    samples_bits_x_23.resize(bits, vector<long> (nslots,0));
}

void QRS_Detection::make_copies(vector< vector<mkt> > input, vector< vector<mkt> > destination){
	for(unsigned n = 0; n < destination.size(); n++){
		for (unsigned b = 0; b < destination[n].size(); b++){
			he.erase(destination[n][b]);
		}
	}

	destination = vector< vector<mkt> >(input.size(), vector<mkt>(input[0].size())); 
	for(unsigned n = 0; n < input.size(); n++){
		for (unsigned b = 0; b < input[n].size(); b++){
			destination[n][b] = he.copy(input[n][b]);
		}
	}
}

/****************************************************************************/
/************************* FHE Dualslope Methods ****************************/
/****************************************************************************/

void QRS_Detection::prepare_data(int iteration, int leftovers){
    // where to start pulling samples from 
    int index = ((nslots - n_considered + 1) * iteration); 
    
    if (index + nslots - 1 > samples.size()){
        // read leftovers in full batch of nslots data (meaning we re-read some data)
        index = samples.size() - nslots;    
    }

    // Compute samples*x, samples*x/a, ... , samples*x/b
	for (int i = 0; i < nslots; i++){
        samples_x[i] = samples[index + i]*scaling_factors[0];
        samples_x_10[i] = samples[index + i]*scaling_factors[1];
        samples_x_11[i] = samples[index + i]*scaling_factors[2];
        samples_x_12[i] = samples[index + i]*scaling_factors[3];
        samples_x_13[i] = samples[index + i]*scaling_factors[4];
        samples_x_14[i] = samples[index + i]*scaling_factors[5];
        samples_x_15[i] = samples[index + i]*scaling_factors[6];
        samples_x_16[i] = samples[index + i]*scaling_factors[7];
        samples_x_17[i] = samples[index + i]*scaling_factors[8];
        samples_x_18[i] = samples[index + i]*scaling_factors[9];
        samples_x_19[i] = samples[index + i]*scaling_factors[10];
        samples_x_20[i] = samples[index + i]*scaling_factors[11];
        samples_x_21[i] = samples[index + i]*scaling_factors[12];
        samples_x_22[i] = samples[index + i]*scaling_factors[13];
        samples_x_23[i] = samples[index + i]*scaling_factors[14];
    }

    for (int i = 0; i < nslots; i++){
        bitset<64> b_x(samples_x[i]);
        bitset<64> b_x_10(samples_x_10[i]);
        bitset<64> b_x_11(samples_x_11[i]);
        bitset<64> b_x_12(samples_x_12[i]);
        bitset<64> b_x_13(samples_x_13[i]);
        bitset<64> b_x_14(samples_x_14[i]);
        bitset<64> b_x_15(samples_x_15[i]);
        bitset<64> b_x_16(samples_x_16[i]);
        bitset<64> b_x_17(samples_x_17[i]);
        bitset<64> b_x_18(samples_x_18[i]);
        bitset<64> b_x_19(samples_x_19[i]);
        bitset<64> b_x_20(samples_x_20[i]);
        bitset<64> b_x_21(samples_x_21[i]);
        bitset<64> b_x_22(samples_x_22[i]);
        bitset<64> b_x_23(samples_x_23[i]);
        for (int b = 0; b < bits; b++){
            samples_bits_x[b][i] = b_x[b];
            samples_bits_x_10[b][i] = b_x_10[b];
            samples_bits_x_11[b][i] = b_x_11[b];
            samples_bits_x_12[b][i] = b_x_12[b];
            samples_bits_x_13[b][i] = b_x_13[b];
            samples_bits_x_14[b][i] = b_x_14[b];
            samples_bits_x_15[b][i] = b_x_15[b];
            samples_bits_x_16[b][i] = b_x_16[b];
            samples_bits_x_17[b][i] = b_x_17[b];
            samples_bits_x_18[b][i] = b_x_18[b];
            samples_bits_x_19[b][i] = b_x_19[b];
            samples_bits_x_20[b][i] = b_x_20[b];
            samples_bits_x_21[b][i] = b_x_21[b];
            samples_bits_x_22[b][i] = b_x_22[b];
            samples_bits_x_23[b][i] = b_x_23[b];
        }
    }
//TODO - implement any more data prep? set first two things to be comp'ed into k/k_constants? 
    k_constant_x.resize(bits);
/*
    for (int i = 0; i < bits; i++){
        k_constant_x[i] = he.encrypt(samples_bits_x[i]);
        cout << "encrypted as: " << k_constant_x[i] << endl;
    }
*/

    inputs.resize(n_considered, vector < long > (nslots,0));
    v_in.resize(n_considered,vector< vector<long> >(bits,vector<long>(nslots,0)));
    k_constant.resize(n_considered, vector < mkt>(bits));

    //inputs to N bit circuits
    for(unsigned j = 0; j < nslots; j++){
            inputs[0][j] = samples_x[j];
            inputs[1][j] = samples_x_10[j];
    }

    if(debug){
        for(unsigned j = 0; j < nslots; j++){
                cout << "inputs[0][" << j << "]: " << inputs[0][j] << endl;
                cout << "inputs[1][" << j << "]: " << inputs[1][j] << endl;
        }
    }
}

vector< vector<mkt> > QRS_Detection::compute_lr_slopes(vector< vector<mkt> > encrypted_pairs, bool simd){
    vector< vector<mkt> > lr_slopes;

    return lr_slopes;
}

vector<mkt> QRS_Detection::compute_mins_maxs(vector<mkt> encrypted_list){
    vector<mkt> result;


    return result;
}

vector<mkt> QRS_Detection::compute_diff_maxs(vector<mkt> encrypted_mins_maxs){
    vector<mkt> result;

    
    return result;
}

vector<mkt> QRS_Detection::compare_to_thresholds(vector<mkt> diff_maxs){
    vector<mkt> result;

    
    return result;
}

mkt QRS_Detection::check_peak_closeness(vector<mkt> peaks){
    mkt true_peak;


    return true_peak;
}

vector<mkt> QRS_Detection::update_thresholds(vector<mkt> diff_maxs){
    vector<mkt> results;


    return results;
}

/***************************************************************************/
/********************* Cleartext Dualslope Methods *************************/
/***************************************************************************/

vector< vector<double> > QRS_Detection::compute_lr_slopes(int index){
    vector< vector<double> > lr_slopes(2, vector<double>(lr_size,0));
    
    if (index < n_considered - 1){
        cout << "Error! Must choose index >= n_considered - 1. Returning empty vector." << endl;
        return lr_slopes; 
    }

    // compute slopes for values a to b samples away from index-b, our center sample. 
    for (int i = 0; i < lr_size; i++){
        double l_diff = samples[index-b] - samples[index-(b-a)+i];
        double r_diff = samples[index-b] - samples[index-(b+a)-i];

        l_diff *= sample_difference_widths[i];
        r_diff *= sample_difference_widths[i];
      /* 
        cout << "samples[" << index-b << "]: " << samples[index-b] << endl;
        cout << "samples[" << index-(b-a)+i << "]: " << samples[index-(b-a)+i] << endl; 
        cout << "samples[" << index-(b+a)-i << "]: " << samples[index-(b+a)-i] << endl;  
        cout << "l_diff, r_diff: " << l_diff << ", " << r_diff << endl;
      */  

        lr_slopes[0][i] = l_diff;
        lr_slopes[1][i] = r_diff;
    }
    
    return lr_slopes;
}

vector< vector<double> > QRS_Detection::compute_mins_maxs(vector< vector<double> > plain_list){
    vector< vector<double> > result(2, vector<double>(2,0)); 

    double l_min = *min_element(begin(plain_list[0]),end(plain_list[0]));
    double l_max = *max_element(begin(plain_list[0]),end(plain_list[0]));

    double r_min = *min_element(begin(plain_list[1]),end(plain_list[1]));
    double r_max = *max_element(begin(plain_list[1]),end(plain_list[1]));

    result[0][0] = l_min;
    result[0][1] = l_max;
    result[1][0] = r_min;
    result[1][1] = r_max;

    if (debug){
        cout << "l_min, l_max = " << l_min << ", " << l_max << endl;
        cout << "r_min, r_max = " << r_min << ", " << r_max << endl;
    }

    return result;
}

bool QRS_Detection::compare_to_thresholds(vector< vector<double> > mins_maxs){
    double r_max_l_min = mins_maxs[1][1] - mins_maxs[0][0];
    double l_max_r_min = mins_maxs[0][1] - mins_maxs[1][0];
    
    double diff_max;
    double s_min;
    bool signs_different;

    if (r_max_l_min > l_max_r_min){
        diff_max = r_max_l_min;
        s_min = min(mins_maxs[1][1], mins_maxs[0][0]);
        signs_different = (signbit(mins_maxs[1][1]) != signbit(mins_maxs[0][0]));
    } else {
        diff_max = l_max_r_min;
        s_min = min(mins_maxs[0][1], mins_maxs[1][0]);
        signs_different = (signbit(mins_maxs[0][1]) != signbit(mins_maxs[1][0]));
    }
    
    if (debug){
        cout << "diff_max = " << diff_max << endl;
        cout << "s_min = " << s_min << endl;
        cout << "signs_different = " << signs_different << endl;
        cout << "diff_threshold = " << diff_threshold << endl; 
        cout << "min_threshold = " << min_threshold << endl;
        cout << "h_cur = " << h_cur << endl; 
    }

    if (diff_max <= diff_threshold){
        return false;
    }

    if (s_min <= min_threshold){
        return false;
    }
    
    if (!signs_different){
        return false;
    }

    double H_ave = std::accumulate(h_aves.begin(), h_aves.end(), 0);
    if (h_cur <= 0.4 * H_ave){
        return false;
    }

    return true;
}

long QRS_Detection::check_peak_closeness(vector<long> peaks){
    long true_peak;


    return true_peak;
}

void QRS_Detection::update_thresholds(long diff_max){

}

/***************************************************************************/
/************************ FHE Dualslope Algorithms *************************/
/***************************************************************************/

// Process whole sample file using the dualslope (ds) algorithm run using fhe (HElib + hbc API) 
void QRS_Detection::ds_fhe(){
    
    // number of ciphertext vectors with nslots of encrypted samples needed to run algorithm on all samples.
 	// note: subtracting (n_considered - 1) because of overlap needed. Algorithm looks at (n_considered - 1) 
    // samples previous to the current sample, thus we must "repeat" these values 
    // in the new packed ciphertext in order to process the next sample. 
    int iterations_needed = n_samples / (nslots - n_considered + 1);
    int leftovers =  n_samples % (nslots - n_considered + 1); 

    if (leftovers != 0){
        // if leftovers exist, add an iteration to process the leftover samples
        iterations_needed++; 
    }
    if (debug){
        cout << "iterations_needed = " << iterations_needed << endl;
        cout << "leftovers = " << leftovers << endl;
    }
    
    ds_fhe(iterations_needed, leftovers);
}

// Process sample file from beginning x = iterations number of times (1 iteration = nslots samples processed) 
void QRS_Detection::ds_fhe(int iterations, int leftovers){    
    for (int i = 0; i < iterations; i++){ 
        prepare_data(i, leftovers); 

        compute_lr_slopes(k, true);
    }

}

// Process whole sample file using dualslope (ds) algortithm but not using packed ciphertexts
void QRS_Detection::ds_unpacked_fhe(){

}

// Process sample file from beginning x = iterations number of times (1 iteration = nslots samples processed) but not using packed ciphertexts  
void QRS_Detection::ds_unpacked_fhe(int iterations, int leftovers){

}

/***************************************************************************/
/******************** Cleartext Dualslope Algorithms ***********************/
/***************************************************************************/

// Process whole sample file using the dualslope (ds) algorithm unencrypted 
void QRS_Detection::ds_plain(){
// We start at n_considered because the algorithm looks "back" at the previous n_considered - 1 samples 
    for (int i = n_considered; i < samples.size(); i++){
        ds_plain(i);
    }
}

// Process sample file from beginning x = iterations number of times (1 iteration = nslots samples processed) unencrypted 
void QRS_Detection::ds_plain(int index){
    h_cur = samples[index - b];
    vector< vector<double> > slopes = compute_lr_slopes(index);
    vector< vector<double> > mins_maxs = compute_mins_maxs(slopes);
    bool need_update = compare_to_thresholds(mins_maxs);
    if (need_update){
        cout << "Yo! We need an update up in hurr!" << endl;
    } else {
        cout << "No update needed..." << endl;
    }
}





































/***************************************************************************/
/******************************* Testing ***********************************/
/***************************************************************************/

bool QRS_Detection::test_ds_fhe(){
    ds_fhe();
    //TODO implement check of values, return true if error occured
    // if (calculated_qrs_locations[everywhere] != annotations[everywhere]) return true;
    return false;
}

bool QRS_Detection::test_ds_fhe(int iterations, int leftovers){
    ds_fhe(iterations, leftovers);
    //TODO implement check of values, return true if error occured
    // if (calculated_qrs_locations[everywhere] != annotations[everywhere]) return true;
    return false;
}

bool QRS_Detection::test_ds_unpacked_fhe(){
    ds_unpacked_fhe();
    //TODO implement check of values, return true if error occured
    // if (calculated_qrs_locations[everywhere] != annotations[everywhere]) return true;
    return false;
}

bool QRS_Detection::test_ds_unpacked_fhe(int iterations, int leftovers){
    ds_unpacked_fhe(iterations, leftovers);
    //TODO implement check of values, return true if error occured
    // if (calculated_qrs_locations[everywhere] != annotations[everywhere]) return true;
    return false;
}

bool QRS_Detection::test_ds_plain(){
    ds_plain();
    //TODO implement check of values, return true if error occured
    // if (calculated_qrs_locations[everywhere] != annotations[everywhere]) return true;
    return false;
}

bool QRS_Detection::test_ds_plain(int index){
    ds_plain(index);
    //TODO implement check of values, return true if error occured
    // if (calculated_qrs_locations[everywhere] != annotations[everywhere]) return true;
    return false;
}
