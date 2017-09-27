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
QRS_Detection::QRS_Detection(vector<double> raw_ecgs, vector<int> anns, int sampling_frequency, bool dbg){
	raw_samples = raw_ecgs;
    annotations = anns;
    fs = sampling_frequency;
    debug = dbg;

    bits = 64; // set here because set_params depends on them and initialize depends on set_params

	//set_params();
	initialize();
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

Errors QRS_Detection::test_all(){
	Errors e("QRS_Detection");
//	e.add("ds_fhe()", test_ds_fhe());
//	e.add("ds_fhe(5,200)", test_ds_fhe(5,200));
//	e.add("ds_fhe_unpacked()", test_ds_unpacked_fhe());
//	e.add("ds_fhe_unpacked(5,200)", test_ds_unpacked_fhe(5,200));
//	e.add("ds_plain()", test_ds_plain());
	e.add("ds_plain(393)", test_ds_plain(393));
	return e;
}

void QRS_Detection::set_params(){
	params.p = 2;
	params.r = 1;
	params.d = 0; // field p^d
	params.k = 128;
	params.slb = 800;
	params.L = 40;
	
	if ((params.L > 42)&&(params.slb > 600)){
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

    samples = scale_samples(raw_samples, 1000); // scale samples to be integers. If precision != 3, adjust "1000" to be 10^precision

    double double_a = 0.027 * fs;   // set as in Wang's paper, an estimate of the min width of QRS complex
    a = round(double_a);
    double double_b = 0.063 * fs;   // set as in Wang's paper, an estimate of the max width of QRS complex
    b = round(double_b);
    
    n_considered = (2 * b) + 1;     // number of samples considered in one iteration (goes from 0 to -2b)
    lr_size = (b - a) + 1;          // convenience for looping through each side 

    diff_threshold = 3840 / fs;     // initial diff_threshold value. to be updated later based on S_ave
    min_threshold = 1536 / fs;      // threshold value given in paper
    avg_height = 0.0;               // set average height to 0. to be updated later based on sample heights

    // TODO - initialize encrypted thresholds

    // Constant values of 1/a, 1/(a+1), ... , 1/b. These are the "run" in the rise over run calculation of slope
    double width;
    double rounded_width;

    for (double i = a; i < (b + 1); i++){
        width = 1/i;
        rounded_width = round(100000 * width)/100000;

        sample_difference_widths.push_back(rounded_width);
    }

    if (debug){
        cout << "a,b = " << a << "," << b << endl;
        cout << "n_considered = " << n_considered << endl;
        cout << "lr_size = " << lr_size << endl;

        for (int i = 0; i < lr_size; i++){
            cout << "sample_difference_widths[" << i << "]: " << sample_difference_widths[i] << endl;
        }
    }
	n_samples = samples.size(); //Total number of samples given to process 
/*	he.debug_on(debug);
	cout << className() << ": Number of bits n was set to " << bits << endl; 
	nslots = he.keyGen(params);
	mkt k_ones = he.setOnes(nslots);
	he.set01(k_ones);
*/
    nslots = 1024; //TEMP! PUT BACK THE ABOVE WHEN RUNNING NON_FHE TESTS

	// Compute scaling_factors based on a and b values. 
    long x = 5354228880;
    // x = 5,354,228,880 = 2^4*3^2*5*7*11*13*17*19*23 = Minimum number divisible by 10, 11, ... , 23
    // Note: If fs != 360, then likely a != 10 and b != 23, which means x should be recomputed to be the minimum value such that x / a, x / (a+1), ... , x / b  all yield whole numbers 

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
	destination = vector< vector<mkt> >(input.size(), vector<mkt>(input[0].size())); // used to be (n_considered, vector<mkt>(bits))
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
    int index = ((nslots - n_considered + 1) * iteration); // where to start pulling samples from 
    
    if (index + nslots - 1 > samples.size()){
        index = samples.size() - nslots;    // read leftovers in full batch of nslots data (meaning we re-read some data)
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

vector<mkt> QRS_Detection::compute_min_max(vector<mkt> encrypted_list){
    vector<mkt> result;


    return result;
}

vector<mkt> QRS_Detection::compute_diff_max(vector<mkt> encrypted_mins_maxs){
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
    
    if (index < n_considered){
        cout << "Error! Must choose index >= n_considered" << endl;
        return lr_slopes; 
    }

    // compute slopes for values a to b samples away from index-b, our center sample. 
    for (int i = 0; i < lr_size; i++){
        // index offset by 1 because of raw_samples' 0-based indexing 
        double l_diff = raw_samples[index-b] - raw_samples[index-(b-a)+i];
        double r_diff = raw_samples[index-b] - raw_samples[index-(b+a)-i];

        l_diff *= sample_difference_widths[i];
        r_diff *= sample_difference_widths[i];
      /* 
        cout << "raw_samples[" << index-b << "]: " << raw_samples[index-b] << endl;
        cout << "raw_samples[" << index-(b-a)+i << "]: " << raw_samples[index-(b-a)+i] << endl; 
        cout << "raw_samples[" << index-(b+a)-i << "]: " << raw_samples[index-(b+a)-i] << endl;  
        cout << "l_diff, r_diff: " << l_diff << ", " << r_diff << endl;
        */

        lr_slopes[0][i] = l_diff;
        lr_slopes[1][i] = r_diff;
    }
    
  /*  for (int i = index - n_considered + 1; i < index + 1; i++){
        cout << "raw_samples[" << i << "]: " << raw_samples[i] << endl;
    }
    */

    return lr_slopes;
}

vector<long> QRS_Detection::compute_min_max(vector<long> plain_list){
    vector<long> result; 


    return result;
}

vector<long> QRS_Detection::compute_diff_max(vector<long> plain_mins_maxs){
    vector<long> result;


    return result;
}

vector<bool> QRS_Detection::compare_to_thresholds(vector<long> diff_maxs){
    vector<bool> result;


    return result;
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
    
    // number of ciphertext vectors with nslots of encrypted samples needed to run algorithm on whole sample set.
 	// note: subtracting n_considered and 1 because of overlap needed. Algorithm looks at (n_considered - 1) 
    // samples previous to the currently considered sample, thus we must "repeat" these values in the new packed 
    // ciphertext in order to process the next sample. 
    int iterations_needed = n_samples / (nslots - n_considered + 1);
    int leftovers =  n_samples % (nslots - n_considered + 1); 

    if (leftovers != 0){
        iterations_needed++; // if leftovers exist, add an iteration to process the leftover samples
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
    for (int i = n_considered; i < raw_samples.size(); i++){
        ds_plain(i);
    }
}

// Process sample file from beginning x = iterations number of times (1 iteration = nslots samples processed) unencrypted 
void QRS_Detection::ds_plain(int index){
    vector< vector<double> > slopes = compute_lr_slopes(index);

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
