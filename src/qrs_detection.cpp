#include "qrs_detection.h"

/*TODO: 
    rethink through parameters/outputs in .h file new methods
    think through security of methods in .h -> reveal info by return values? can hide. 
    implement min_max based on that one paper's method 3
    change parameters (params.L calculation) based on expected levels needed for subtration/comparison. 
    double check that comparison works on neg numbers 
*/

QRS_Detection::QRS_Detection(vector<long> digital_ecgs, vector<int> anns, int sampling_frequency, bool dbg){
	samples = digital_ecgs;
    annotations = anns;
    fs = (double)sampling_frequency;
    debug = dbg;

    // number of bits needed based on the largest sample scaled by lcm and multiplied by 8 = h_aves/s_aves size
    // this is because we may need to average these 8 values, which would involve adding up to 8 together
    int max_sample = *max_element(samples.begin(), samples.end());
    long highest_value = (long)max_sample * lcm * 8;
    long double bits_needed = log2((long double)highest_value);

    // set here because set_params depends on bits and initialize depends on set_params
    bits = ceil(bits_needed); 
	set_params();
	initialize();
}

Errors QRS_Detection::test_all(int file_index){
	Errors e("QRS_Detection");
    //string plain = "ds_plain(" + to_string(file_index) + ")";
	//e.add(plain, test_ds_plain(file_index));

    string fhe = "ds_fhe(" + to_string(file_index) + ")";
    e.add(fhe, test_ds_fhe(file_index));
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
	params.k = 128; // could set to 256
	params.slb = 800;
	params.L = 20; // TODO - change params.L to equal a level just high enough to run calculations without error
	
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
    //Total number of samples given to process 
	n_samples = samples.size(); 


    // Constant values of 1/a, 1/(a+1), ... , 1/b. These are the "run" in the rise over run slope calculation 
    double width;
    for (double i = a; i < (b + 1); i++){
        width = 1/i;
        sample_difference_widths.push_back(width);
    }
    // false negatives reported by Wang in paper. values = # of FN in sample 100, 101, ... 
    false_negatives = { 1,0,0,1,0,1,3,2,23,0,1,0,1,1,0,24,0,0,0,2,0,0,2,0,51,4,3,0,15,39,0,4,0,3,5,0,1,0,0,4,0,1,3,0,0,0,2,2 };
    // false positives reported by Wang in paper. values = # of FP in sample 100, 101, ... 
    false_positives = { 0,1,0,0,58,72,3,0,21,0,2,0,0,5,0,2,1,6,0,0,0,1,1,30,0,1,79,0,15,5,4,18,4,0,3,0,17,0,0,2,5,1,44,2,0,2,0,0 };


    // lowest diff_threshold value in paper. we check s_diff_max against diff_threshold to determine if peak. 
    diff_threshold = 3840 / fs; 
    // min_threshold a constant as given in paper. check smallest slope in s_diff_max against this. 
    min_threshold = 1536 / fs; 


    // 8 = the max number used in computing S_ave, the average of the previous 8 peaks' s_diff_max.
    s_aves = boost::circular_buffer<double>(8); 
    // 8 = the max number used in computing H_ave, the average of the previous 8 peaks' h_cur.
    h_aves = boost::circular_buffer<double>(8); 
    // updated later based on currently considered sample's value in mV
    h_cur = 0.0;
    // updated later based on currently considered sample's max slope difference between left and right side
    s_diff_max = 0.0;


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
	he.debug_on(true); //TODO - change to debug when done testing
	cout << className() << ": Number of bits was set to " << bits << endl; 
	nslots = he.keyGen(params);
	mkt k_ones = he.setOnes(nslots);
	he.set01(k_ones);

  //  /*//TEMP! USE WHEN RUNNING NON_FHE TESTS
  //   nslots = 1024; 
  //  */
  
    // create nslots with the encrypted, scaled low diff_threshold value found in paper 
    //long scaled_low_threshold = 3840 * lcm / (long)fs;
    long scaled_low_threshold = 3840 / (long)fs;
    long long sd = 3840 * lcm / (long long)fs;
    vector<long> vector_low_thresholds(nslots, scaled_low_threshold);
    vector< vector<long> > low_thresholds_bits = conv.longVec2Matrix(vector_low_thresholds);
    unsigned bits_low = low_thresholds_bits.size();

    encrypted_low_threshold.resize(bits_low);
    for (unsigned b = 0; b < bits_low; b++){
        encrypted_low_threshold[b] = he.encrypt(low_thresholds_bits[b]);
    }

    // create nslots with the encrypted, scaled mid diff_threshold value found in paper TODO CHANGE vector<mkt>  
    long scaled_mid_threshold = 4352 * lcm / (long)fs;
    vector<long> vector_mid_thresholds(nslots, scaled_mid_threshold);
    vector< vector<long> > mid_thresholds_bits = conv.longVec2Matrix(vector_mid_thresholds);
    unsigned bits_mid = mid_thresholds_bits.size();

    encrypted_mid_threshold.resize(bits_mid);
    for (unsigned b = 0; b < bits_mid; b++){
        encrypted_mid_threshold[b] = he.encrypt(mid_thresholds_bits[b]);
    }

    // create nslots with the encrypted, scaled high diff_threshold value found in paper 
    //long scaled_high_threshold = 7680 * lcm / (long)fs;
    long scaled_high_threshold = 7680 / (long)fs;
    vector<long> vector_high_thresholds(nslots, scaled_high_threshold);
    vector< vector<long> > high_thresholds_bits = conv.longVec2Matrix(vector_high_thresholds);
    unsigned bits_high = high_thresholds_bits.size();

    encrypted_high_threshold.resize(bits_high);
    for (unsigned b = 0; b < bits_high; b++){
        encrypted_high_threshold[b] = he.encrypt(high_thresholds_bits[b]);
    }

    // create nslots with the encrypted, scaled low S_ave value found in paper
    long scaled_low_ave = 12800 * lcm / (long)fs;
    vector<long> low_aves(nslots, scaled_low_ave);
    mkt encrypted_low_ave = he.encrypt(low_aves);
    // create nslots with the encrypted, scaled high S_ave value found in paper
    long scaled_high_ave = 20480 * lcm / (long)fs;
    vector<long> high_aves(nslots, scaled_high_ave);
    mkt encrypted_high_ave = he.encrypt(high_aves);

    // set encrypted_diff_threshold to the low_threshold value to start. will change later in update_thresholds()
    encrypted_diff_threshold.resize(bits_low);
    for (unsigned b = 0; b < bits_low; b++){
        encrypted_diff_threshold[b] = he.copy(encrypted_low_threshold[b]);
    }

    // set encrypted_min_threshold to encrypted, scaled version of value found in paper. Doesn't change.
    long scaled_min_threshold = 1536 * lcm / (long)fs;
    vector<long> v_min_thresholds(nslots, scaled_min_threshold);
    encrypted_min_threshold = he.encrypt(v_min_thresholds);


    // 8 = the max number used in computing S_ave, the average of the previous 8 peaks' s_diff_max.
    encrypted_s_aves = boost::circular_buffer<mkt>(8); 
    // 8 = the max number used in computing H_ave, the average of the previous 8 peaks' h_cur.
    encrypted_h_aves = boost::circular_buffer<mkt>(8); 
    // create nslots of encrypted zeros
    vector<long> zeros(nslots,0);
    mkt encrypted_zeros = he.encrypt(zeros);
    // initialize encrypted_h_cur to nslots 0s. to be updated later based on sample value in mV
    encrypted_h_cur = encrypted_zeros;
    // initialize encrypted_s_diff_max to nslots 0s. to be updated later based on current sample's s_diff_max
    encrypted_s_diff_max = encrypted_zeros;

    /* Compute scaling_factors based on a and b values. */
    scaling_factors.push_back(lcm);
    for (long i = a; i <= b; i++){
        scaling_factors.push_back(lcm/i);
    }

    // resize scaled sample vectors. _lcm = scaled by lcm, _lcm_xx = scaled by lcm / xx
    samples_lcm.resize(nslots,0);
    samples_lcm_10.resize(nslots,0);
    samples_lcm_11.resize(nslots,0);
    samples_lcm_12.resize(nslots,0);
    samples_lcm_13.resize(nslots,0);
    samples_lcm_14.resize(nslots,0);
    samples_lcm_15.resize(nslots,0);
    samples_lcm_16.resize(nslots,0);
    samples_lcm_17.resize(nslots,0);
    samples_lcm_18.resize(nslots,0);
    samples_lcm_19.resize(nslots,0);
    samples_lcm_20.resize(nslots,0);
    samples_lcm_21.resize(nslots,0);
    samples_lcm_22.resize(nslots,0);
    samples_lcm_23.resize(nslots,0);

    // resize the scaled sample bit vectors. these hold the binary representations of the above vectors
    samples_lcm_bits.resize(bits, vector<long> (nslots,0));
    samples_lcm_10_bits.resize(bits, vector<long> (nslots,0));
    samples_lcm_11_bits.resize(bits, vector<long> (nslots,0));
    samples_lcm_12_bits.resize(bits, vector<long> (nslots,0));
    samples_lcm_13_bits.resize(bits, vector<long> (nslots,0));
    samples_lcm_14_bits.resize(bits, vector<long> (nslots,0));
    samples_lcm_15_bits.resize(bits, vector<long> (nslots,0));
    samples_lcm_16_bits.resize(bits, vector<long> (nslots,0));
    samples_lcm_17_bits.resize(bits, vector<long> (nslots,0));
    samples_lcm_18_bits.resize(bits, vector<long> (nslots,0));
    samples_lcm_19_bits.resize(bits, vector<long> (nslots,0));
    samples_lcm_20_bits.resize(bits, vector<long> (nslots,0));
    samples_lcm_21_bits.resize(bits, vector<long> (nslots,0));
    samples_lcm_22_bits.resize(bits, vector<long> (nslots,0));
    samples_lcm_23_bits.resize(bits, vector<long> (nslots,0));

    // resize the k_constant vectors. these hold the permanent encrypted values of different scaled samples
    k_constant_lcm.resize(bits);
    k_constant_lcm_10.resize(bits);
    k_constant_lcm_11.resize(bits);
    k_constant_lcm_12.resize(bits);
    k_constant_lcm_13.resize(bits);
    k_constant_lcm_14.resize(bits);
    k_constant_lcm_15.resize(bits);
    k_constant_lcm_16.resize(bits);
    k_constant_lcm_17.resize(bits);
    k_constant_lcm_18.resize(bits);
    k_constant_lcm_19.resize(bits);
    k_constant_lcm_20.resize(bits);
    k_constant_lcm_21.resize(bits);
    k_constant_lcm_22.resize(bits);
    k_constant_lcm_23.resize(bits);

    //TODO - initialize/resize inputs, v_in, k, k_constant?
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

void QRS_Detection::prepare_data(int iteration){
    // where to start pulling samples from 
    int index = ((nslots - n_considered + 1) * iteration); 
    
    if (index + nslots - 1 > samples.size()){
        // read leftovers in full batch of nslots data (meaning we re-read some data)
        index = samples.size() - nslots;    
    }

    // Compute samples*lcm, samples*lcm/a, ... , samples*lcm/b
	for (int i = 0; i < nslots; i++){
        samples_lcm[i] = samples[index + i]*scaling_factors[0];
        samples_lcm_10[i] = samples[index + i]*scaling_factors[1];
        samples_lcm_11[i] = samples[index + i]*scaling_factors[2];
        samples_lcm_12[i] = samples[index + i]*scaling_factors[3];
        samples_lcm_13[i] = samples[index + i]*scaling_factors[4];
        samples_lcm_14[i] = samples[index + i]*scaling_factors[5];
        samples_lcm_15[i] = samples[index + i]*scaling_factors[6];
        samples_lcm_16[i] = samples[index + i]*scaling_factors[7];
        samples_lcm_17[i] = samples[index + i]*scaling_factors[8];
        samples_lcm_18[i] = samples[index + i]*scaling_factors[9];
        samples_lcm_19[i] = samples[index + i]*scaling_factors[10];
        samples_lcm_20[i] = samples[index + i]*scaling_factors[11];
        samples_lcm_21[i] = samples[index + i]*scaling_factors[12];
        samples_lcm_22[i] = samples[index + i]*scaling_factors[13];
        samples_lcm_23[i] = samples[index + i]*scaling_factors[14];
    }

    for (int i = 0; i < nslots; i++){
        bitset<64> b_lcm(samples_lcm[i]);
        bitset<64> b_lcm_10(samples_lcm_10[i]);
        bitset<64> b_lcm_11(samples_lcm_11[i]);
        bitset<64> b_lcm_12(samples_lcm_12[i]);
        bitset<64> b_lcm_13(samples_lcm_13[i]);
        bitset<64> b_lcm_14(samples_lcm_14[i]);
        bitset<64> b_lcm_15(samples_lcm_15[i]);
        bitset<64> b_lcm_16(samples_lcm_16[i]);
        bitset<64> b_lcm_17(samples_lcm_17[i]);
        bitset<64> b_lcm_18(samples_lcm_18[i]);
        bitset<64> b_lcm_19(samples_lcm_19[i]);
        bitset<64> b_lcm_20(samples_lcm_20[i]);
        bitset<64> b_lcm_21(samples_lcm_21[i]);
        bitset<64> b_lcm_22(samples_lcm_22[i]);
        bitset<64> b_lcm_23(samples_lcm_23[i]);
        for (unsigned b = 0; b < bits; b++){
            samples_lcm_bits[b][i] = b_lcm[b];
            samples_lcm_10_bits[b][i] = b_lcm_10[b];
            samples_lcm_11_bits[b][i] = b_lcm_11[b];
            samples_lcm_12_bits[b][i] = b_lcm_12[b];
            samples_lcm_13_bits[b][i] = b_lcm_13[b];
            samples_lcm_14_bits[b][i] = b_lcm_14[b];
            samples_lcm_15_bits[b][i] = b_lcm_15[b];
            samples_lcm_16_bits[b][i] = b_lcm_16[b];
            samples_lcm_17_bits[b][i] = b_lcm_17[b];
            samples_lcm_18_bits[b][i] = b_lcm_18[b];
            samples_lcm_19_bits[b][i] = b_lcm_19[b];
            samples_lcm_20_bits[b][i] = b_lcm_20[b];
            samples_lcm_21_bits[b][i] = b_lcm_21[b];
            samples_lcm_22_bits[b][i] = b_lcm_22[b];
            samples_lcm_23_bits[b][i] = b_lcm_23[b];
        }
    }
//TODO - implement any more data prep? set first two things to be comp'ed into k/k_constants? 
    k_constant_lcm.resize(bits);
/*
    for (int i = 0; i < bits; i++){
        k_constant_lcm[i] = he.encrypt(samples_lcm_bits[i]);
        cout << "encrypted as: " << k_constant_lcm[i] << endl;
    }
*/

    inputs.resize(n_considered, vector<long>(nslots,0));
    v_in.resize(n_considered, vector< vector<long> >(bits,vector<long>(nslots,0)));
    k_constant.resize(n_considered, vector<mkt>(bits));

    //inputs to N bit circuits
    for(unsigned j = 0; j < nslots; j++){
            inputs[0][j] = samples_lcm[j];
            inputs[1][j] = samples_lcm_10[j];
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

/*
    double S_ave = std::accumulate(s_aves.begin(), s_aves.end(), 0.0);
    if (s_aves.size() != 0){
        S_ave /= s_aves.size();
    }
    
    double tmp = diff_threshold;
    // high_ave and low_ave are given in Wang's paper
    double high_ave = 20480 / fs; 
    double low_ave = 12800 / fs;
    
    if (S_ave > high_ave){
        diff_threshold = 7680 / fs;
    } else if (S_ave > low_ave){
        diff_threshold = 4352 / fs;
    } else {
        diff_threshold = 3840 / fs;
    }

    if (debug){
        if (tmp != diff_threshold){
            cout << "Diff_threshold changing! Old = " << tmp << ", New = " << diff_threshold << endl;
        }
    }*/


    vector<long> tmp(nslots); 
    tmp = conv.matrix2LongVec(he.decryptNbits(encrypted_diff_threshold));
    cout << "tmp[0] = " << tmp[0] << endl;

    unsigned bits_high = encrypted_high_threshold.size();
    encrypted_diff_threshold.resize(bits_high);
    for (unsigned b = 0; b < bits_high; b++){
        encrypted_diff_threshold[b] = he.copy(encrypted_high_threshold[b]);
    }

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

    return result;
}

bool QRS_Detection::compare_to_thresholds(vector< vector<double> > mins_maxs, int index){
    double r_max_l_min = mins_maxs[1][1] - mins_maxs[0][0];
    double l_max_r_min = mins_maxs[0][1] - mins_maxs[1][0];
    
    double s_min;
    bool signs_different;

    if (r_max_l_min > l_max_r_min){
        s_diff_max = r_max_l_min;
        s_min = min(mins_maxs[1][1], mins_maxs[0][0]);
        signs_different = (signbit(mins_maxs[1][1]) != signbit(mins_maxs[0][0]));
    } else {
        s_diff_max = l_max_r_min;
        s_min = min(mins_maxs[0][1], mins_maxs[1][0]);
        signs_different = (signbit(mins_maxs[0][1]) != signbit(mins_maxs[1][0]));
    }
    
    if (debug){
        cout << "s_diff_max = " << s_diff_max << endl;
        cout << "diff_threshold = " << diff_threshold << endl; 
        cout << "s_min = " << s_min << endl;
        cout << "min_threshold = " << min_threshold << endl;
        cout << "signs_different = " << signs_different << endl;
        cout << "h_cur = " << h_cur << endl; 
    }

    if (s_diff_max <= diff_threshold){
        return false;
    }
    if (debug){
        cout << "s_diff_max > diff_threshold!" << "\t\t\t" << "index - b = " << index - b << endl;
    }

    if (s_min <= min_threshold){
        return false;
    }
    if (debug){
        cout << "s_min > min_threshold!" << "\t\t\t\t\t" << "index - b = " << index - b << endl;
    }
   /* 
    if (signs_different){
        return false;
    }
    if (debug){
        cout << "!signs_different, so the same?" << "\t\t\t" << "index - b = " << index - b << endl;
    }*/

    double H_ave = std::accumulate(h_aves.begin(), h_aves.end(), 0.0);
    if (h_aves.size() != 0){
        // H_ave = average sample height in mV for the previous 8 peaks  
        H_ave /= h_aves.size(); 
    }

    if (h_cur <= 0.4 * H_ave){
        return false;
    }

    if (debug) {
        cout << "h_cur > 0.4 * H_ave (" << H_ave << ")" << "\t\t\t\t" << "index - b = " << index - b << endl;
    }
    return true;
}

bool QRS_Detection::check_local_extremes(int index){
    // if peak index is within b samples of previous peak, update peak to be the one with the greater s_diff_max 
    // (current one stored in global s_diff_max, the other is the most recently added value to s_aves)
    if (qrs_locations.size() != 0){
        int last_peak = qrs_locations.back();
        // if last_peak is within the considered sample range, eliminate the peak with lower s_diff_max
        if ((index - n_considered) <= last_peak){ 
            double last_peak_diff_max = s_aves.back();
            if (s_diff_max > last_peak_diff_max){
                qrs_locations[qrs_locations.size() - 1] = (index - b);
                s_aves[s_aves.size() - 1] = s_diff_max;
                h_aves[h_aves.size() - 1] = h_cur;
            }
            // return true because there was a local_extreme in range and appropriate action has been taken
            return true;
        }
    }
    // return false because there was no local_extreme in range and no action has been taken 
    return false;
}

void QRS_Detection::update_thresholds(){
    double S_ave = std::accumulate(s_aves.begin(), s_aves.end(), 0.0);
    if (s_aves.size() != 0){
        S_ave /= s_aves.size();
    }
    
    double tmp = diff_threshold;
    // high_ave and low_ave are given in Wang's paper
    double high_ave = 20480 / fs; 
    double low_ave = 12800 / fs;
    
    if (S_ave > high_ave){
        diff_threshold = 7680 / fs;
    } else if (S_ave > low_ave){
        diff_threshold = 4352 / fs;
    } else {
        diff_threshold = 3840 / fs;
    }

    if (debug){
        if (tmp != diff_threshold){
            cout << "Diff_threshold changing! Old = " << tmp << ", New = " << diff_threshold << endl;
        }
    }
}

/***************************************************************************/
/************************ FHE Dualslope Algorithms *************************/
/***************************************************************************/

/* Process whole sample file using the dualslope (ds) algorithm run using fhe (HElib + hbc API) */
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
/*
    for (int i = 0; i < iterations_needed; i++){
        ds_fhe(i);
    }
*/
    ds_fhe(0);
}

/* Process sample file from beginning x = iterations number of times (1 iteration = nslots samples processed) */
void QRS_Detection::ds_fhe(int iteration){    
    prepare_data(iteration); 

    vector< vector<mkt> > slopes = compute_lr_slopes(k, true);
/*
    // Convert samples from ADU counts to mV and store in h_cur. Database is 11-bit over ~10mV. 1024 ADU = 0mV
    h_cur = (samples[index - b] - 1024) * 0.005; 
    vector< vector<double> > slopes = compute_lr_slopes(index);
    vector< vector<double> > mins_maxs = compute_mins_maxs(slopes);
    bool is_peak = compare_to_thresholds(mins_maxs, index);
    if (is_peak){
        if (debug){
            cout << "we found a peak! ------------------------------" << endl;
            cout << "index - b, s_diff_max, h_cur = " << index-b << ", " << s_diff_max << ", " << h_cur << endl;
        }

        bool local_extreme_found = check_local_extremes(index); 
        if (!local_extreme_found){
            qrs_locations.push_back(index-b);
            s_aves.push_back(s_diff_max);
            h_aves.push_back(h_cur);
        }
        update_thresholds();
    }
*/
}

/* Process whole sample file using dualslope (ds) algortithm but not using packed ciphertexts */
void QRS_Detection::ds_unpacked_fhe(){

}

/* Process sample file from beginning x = iterations number of times (1 iteration = nslots samples processed) but not using packed ciphertexts */ 
void QRS_Detection::ds_unpacked_fhe(int iteration){

}

/***************************************************************************/
/******************** Cleartext Dualslope Algorithms ***********************/
/***************************************************************************/

/* Process whole sample file using the dualslope (ds) algorithm unencrypted */ 
void QRS_Detection::ds_plain(){
// We start at n_considered because the algorithm looks "back" at the previous n_considered - 1 samples 
    for (int i = n_considered; i < samples.size(); i++){
        ds_plain(i);
    }
}

/* Process sample file from beginning x = iterations number of times (1 iteration = nslots samples processed) unencrypted */ 
void QRS_Detection::ds_plain(int index){
    /*
       int start = 33777;
       int end = 33825;
       if ((index > start) && (index < end)){
        debug = true;
        cout << "WE'RE ON INDEX = " << index << ", ok? ==================, index - b = " << index - b << endl;
    } else {
        debug = false;
    }
   */ 
    // Convert samples from ADU counts to mV and store in h_cur. Database is 11-bit over ~10mV. 1024 ADU = 0mV
    h_cur = (samples[index - b] - 1024) * 0.005; 
    vector< vector<double> > slopes = compute_lr_slopes(index);
    vector< vector<double> > mins_maxs = compute_mins_maxs(slopes);
    bool is_peak = compare_to_thresholds(mins_maxs, index);
    if (is_peak){
        if (debug){
            cout << "we found a peak! ------------------------------" << endl;
            cout << "index - b, s_diff_max, h_cur = " << index-b << ", " << s_diff_max << ", " << h_cur << endl;
        }

        bool local_extreme_found = check_local_extremes(index); 
        if (!local_extreme_found){
            qrs_locations.push_back(index-b);
            s_aves.push_back(s_diff_max);
            h_aves.push_back(h_cur);
        }
        update_thresholds();
    }
    
    return;
}




















/***************************************************************************/
/******************************* Testing ***********************************/
/***************************************************************************/

bool QRS_Detection::test_ds_fhe(int file_index){
    ds_fhe();

    int false_neg = 2011; // TODO set to zero once implemented 
    int false_pos = 2011; // TODO set to zero once implemented 

    vector<mkt> throwaway;
    vector<mkt> resulting = update_thresholds(throwaway);

    // decrypt, check
    vector<long> decrypted_diff_thresholds(nslots); 
    decrypted_diff_thresholds = conv.matrix2LongVec(he.decryptNbits(encrypted_diff_threshold));
    cout << "decrypted_diff_thresholds[0]: " << decrypted_diff_thresholds[0] << endl;

    bool hasError = false;

    if (false_neg != false_negatives[file_index]){
        cout << "ERROR! Expected " << false_negatives[file_index] << " false negative(s), but found " << false_neg << endl;
        hasError = true; // found a different number of false negatives than paper
    }

    if (false_pos != false_positives[file_index]){
        cout << "ERROR! Expected " << false_positives[file_index] << " false positive(s), but found " << false_pos << endl;
        hasError = true; // found a different number of false positives than paper
    }

    return hasError;
}
/****************************************************************************************************************/
/****************************************************************************************************************/
bool QRS_Detection::test_ds_unpacked_fhe(int file_index){
    ds_unpacked_fhe();

    int false_neg = 2011; // TODO set to zero once implemented 
    int false_pos = 2011; // TODO set to zero once implemented 

    // decrypt, check

    bool hasError = false;

    if (false_neg != false_negatives[file_index]){
        cout << "ERROR! Expected " << false_negatives[file_index] << " false negative(s), but found " << false_neg << endl;
        hasError = true; // found a different number of false negatives than paper
    }

    if (false_pos != false_positives[file_index]){
        cout << "ERROR! Expected " << false_positives[file_index] << " false positive(s), but found " << false_pos << endl;
        hasError = true; // found a different number of false positives than paper
    }

    return hasError;
}

/****************************************************************************************************************/
/****************************************************************************************************************/
bool QRS_Detection::test_ds_plain(int file_index){
    ds_plain(); // run the qrs_detection algorithm
    
    int min_size = min(qrs_locations.size(), annotations.size());

    int false_neg = 0;
    int false_pos = 0;
    
    // we say we have either a false_pos or false_neg location when diff is > 5 or < -5. Threshold can be changed
    for (int i = 0; i < min_size; i++){
        int diff = annotations[i] - qrs_locations[i];
        if (diff > 5){ // thought there was a beat here. oops!
            annotations.insert(annotations.begin() + i, -1); // insert dummy value to make vector sizes equal
            false_pos++;
        } 
        if (diff < -5){ // missed an annotated beat!
            qrs_locations.insert(qrs_locations.begin() + i, -1); // insert dummy value to make vector sizes equal
            false_neg++;
        }
    }

    if (qrs_locations.size() < annotations.size()){ // missed annotated beats at the end of the file
        for (int i = qrs_locations.size(); i < annotations.size(); i++){
            qrs_locations.insert(qrs_locations.begin() + i, -1); // insert dummy value to make vector sizes equal
            false_neg++;
        }
    } else { // found too many "beats", or... sizes are equal, in which case for loop won't run 
        for (int i = annotations.size(); i < qrs_locations.size(); i++){
            annotations.insert(annotations.begin() + i, -1); // insert dummy value to make vector sizes equal
            false_pos++;
        }
    }
    
    if (debug){
        for (int i = 0; i < qrs_locations.size(); i++){ // print pairs. qrs_locations.size() = annotations.size() 
            cout << qrs_locations[i] << ", " << annotations[i] << endl;
        }
        cout << "false_neg, false_pos = " << false_neg << ", " << false_pos << endl;
    }

    bool hasError = false;

    if (false_neg != false_negatives[file_index]){
        cout << "ERROR! Expected " << false_negatives[file_index] << " false negative(s), but found " << false_neg << endl;
        hasError = true; // found a different number of false negatives than paper
    }

    if (false_pos != false_positives[file_index]){
        cout << "ERROR! Expected " << false_positives[file_index] << " false positive(s), but found " << false_pos << endl;
        hasError = true; // found a different number of false positives than paper
    }

    return hasError;
}
