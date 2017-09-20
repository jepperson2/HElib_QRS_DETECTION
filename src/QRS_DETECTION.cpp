#include "QRS_DETECTION.h"

QRS_DETECTION::QRS_DETECTION(vector<double> ecgs, int sampling_frequency, bool dbg){
	samples = ecgs;
    fs = sampling_frequency;
    debug = dbg;

    bits = 16;

	set_params();
	initialize();
}
const string QRS_DETECTION::className(){
	return "QRS_DETECTION";
}
void QRS_DETECTION::t_start(){
	if (debug){
		t.start();
	}
}
void QRS_DETECTION::t_end(string name){
	if (debug){
		double duration = t.end("silent");
		cout << className() << ": " << name << " - " << nslots << " operations done in " << duration << "s" << endl;
		cout << className() << ": " << name << " - Time for each operation: " << 1000*duration/nslots << "ms" << endl;
	}
}

Errors QRS_DETECTION::test(){
	Errors e("QRS_DETECTION");
//	e.add("RIPPLE CARRY ADDER sequential circuit", test_RCADDER());
//	e.add("RIPPLE BORROW SUBTRACTOR sequential circuit", test_RBSUBER());
//	e.add("RIPPLE COMPARATOR sequential circuit", test_RCMP());
	e.add("DUAL SLOPE algorithm test", test_DUALSLOPE());
//	e.add("RIPPLE EQUALITY COMPARATOR sequential circuit", test_REQ());
//	e.add("SHIFT RIGHT sequential circuit", test_SHIFTR());
//	e.add("SHIFT LEFT sequential circuit", test_SHIFTL());
//	e.add("MUX sequential circuit", test_NMUX());
	return e;
}

void QRS_DETECTION::set_params(){
	params.p = 2;
	params.r = 1;
	params.d = 0; // field p^d
	params.k = 128;
	params.slb = 800;
	params.L = 0;
	if(params.L == 0){ //L not set 
	//for the ripple carry adder so for circuits with complexity up to 3n+1
	//not valid for ripple comparator !
		switch(bits){
			case 1:
				params.L = 5;
				break;
			case 2:
				params.L = 7;
				break;
			case 3:
				params.L = 9;
				break;
			case 4:
				params.L = 12;
				break;
			case 5:
				params.L = 13;
				break;
			case 6:
				params.L = 15;
				break;
			case 7:
				params.L = 17;
				break;
			case 8:
				params.L = 19;
				break;
			case 9:
				params.L = 21;
				break;
			case 10:
				params.L = 23;
				break;
			case 11:
				params.L = 27;
				break;
			case 12:
				params.L = 29;
				break;
			case 13:
				params.L = 31;
				break;
			case 14:
				params.L = 34;
				break;
			case 15:
				params.L = 35;
				break;
			case 16:
				params.L = 37;
				break;
			if(params.L == 0){ // bits not in 1 .. 16
				params.L = 44; //should work with everything
			}
		}
	}
	
	
	
	if ((params.L > 42)&&(params.slb > 600)){
		params.L *= 1.2;
	}
	params.c = 3;
	params.w = 64;
	
	params.m = FindM(params.k,params.L,params.c,params.p,params.d,params.slb,0);
}
void QRS_DETECTION::initialize(){
    if (debug){
        cout << className() << ": Initializing..." << endl;
    }   

    double double_a = 0.027 * fs;   // set as in Wang's paper, an estimate of the min width of QRS complex
    a = round(double_a);
    double double_b = 0.063 * fs;   // set as in Wang's paper, an estimate of the max width of QRS complex
    b = round(double_b);
    
    n_samples = (2 * b) + 1;        // number of samples considered in one iteration (goes from 0 to -2b)
    lr_size = (b - a) + 1;          // convenience for looping through each side 

    diff_threshold = 3840 / fs;     // initial diff_threshold value. to be updated later based on S_ave
    min_threshold = 1536 / fs;      // threshold value given in paper
    avg_height = 0.0;               // set average height to 0. to be updated later based on sample heights

    // Constant values of 1/a, 1/(a+1), ... , 1/b. These are the values of the "run" in the rise over run calculation of slope
    double width;
    double rounded_width;

    for (double i = a; i < (b + 1); i++){
        width = 1/i;
        rounded_width = round(100000 * width)/100000;

        sample_difference_widths.push_back(rounded_width);
    }

    if (debug){
        cout << "a,b = " << a << "," << b << endl;
        cout << "n_samples = " << n_samples << endl;
        cout << "lr_size = " << lr_size << endl;

        for (int i = 0; i < lr_size; i++){
            cout << "sample_difference_widths[" << i << "]: " << sample_difference_widths[i] << endl;
        }
    }
	N_numbers = n_samples;
	he.debug_on(debug);
	cout << className() << ": Number of bits n was set to " << bits << endl; 
	nslots = he.keyGen(params);
	mkt k_ones = he.setOnes(nslots);
	he.set01(k_ones);
	inputs.resize(N_numbers, vector < long > (nslots,0));
	v_in.resize(N_numbers,vector< vector<long> >(bits,vector<long>(nslots,0)));
    cout << "v_in size: " << v_in.size() << endl;
    cout << "v_in[0] size: " << v_in[0].size() << endl;
    cout << "v_in[0][0] size: " << v_in[0][0].size() << endl;
	k_constant.resize(N_numbers, vector < mkt>(bits));

	//inputs to N bit circuits
	for(unsigned i = 0; i < nslots; i++){
		inputs[0][i] = rand() % (unsigned)pow(2,bits);
		inputs[1][i] = rand() % (unsigned)pow(2,bits);

		if(debug){
			cout << "inputs[0][" << i << "], inputs[1][" << i << "]: " << inputs[0][i] << ", " << inputs[1][i] << endl;
		}
	}
	
	//Converts inputs to bits into v_in for parallel ciphertexts
	for(unsigned n = 0; n < N_numbers; n++){
		for(unsigned j = 0; j < nslots; j++){
			bitset<64> bin(inputs[n][j]); //max is 2^64 so max nbits = 64
			for(unsigned b = 0; b < bits; b++){
				v_in[n][b][j] = bin[b]; //first ctxt (b = 0) is LSB
			}
		}
	}
	
	//Encrypts all the vectors into ciphertexts
	if(debug){
		cout << className() << ": Encrypting input vectors (" << N_numbers * bits << " vectors)" << endl;
	}
	for(unsigned n = 0; n < N_numbers; n++){
		for (unsigned b = 0; b < bits; b++){
			k_constant[n][b] = he.encrypt(v_in[n][b]);
			if(debug){
				cout << "k_constant[" << n << "][" << b << "]: " << k_constant[n][b] << endl; 
			}
		}
	}
}
void QRS_DETECTION::make_copies(){
	for(unsigned n = 0; n < k.size(); n++){
		for (unsigned b = 0; b < k[n].size(); b++){
			he.erase(k[n][b]);
		}
	}
	k = vector< vector<mkt> >(N_numbers, vector<mkt>(bits));
	for(unsigned n = 0; n < N_numbers; n++){
		for (unsigned b = 0; b < bits; b++){
			k[n][b] = he.copy(k_constant[n][b]);
		}
	}
}
bool QRS_DETECTION::test_RCADDER(){
	//NO 2's complement, MSB is used for "overflow" (last carry out)
	make_copies();
	t_start();
	he.RCADDER(k[0],k[1]);
	t_end(__FUNCTION__);
	vector<long> results(nslots);
	results = conv.matrix2LongVec(he.decryptNbits(k[0]));
	for (unsigned i = 0; i < nslots; i++){ //bit level
		if (results[i] != (inputs[0][i] + inputs[1][i])){
			return true;
		}
	}
	return false;
}

bool QRS_DETECTION::test_RBSUBER(){
	make_copies();
	string bin_str;
	t_start();
	he.RBSUBER(k[0],k[1]);
	t_end(__FUNCTION__);
	vector<long> results(nslots);
	results = conv.matrix2SignedLongVec(he.decryptNbits(k[0]));
	for (unsigned i = 0; i < nslots; i++){ //bit level
		if (results[i] != (inputs[0][i] - inputs[1][i])){
			return true;
		}
	}
	return false;
}
bool QRS_DETECTION::test_RCMP(){
	make_copies();
	t_start();
	he.RCMP(k[0],k[1]);
	t_end(__FUNCTION__);
	
	vector< vector<long> > v_out1, v_out2;
	v_out1 = he.decryptNbits(k[0]);
	v_out2 = he.decryptNbits(k[1]);

	if(debug){
		for(unsigned i = 0; i < nslots; i++){
			cout << "v_out1[0][" << i << "], v_out2[0][" << i << "]: " << v_out1[0][i] << ", " << v_out2[0][i] <<  endl;  
		}
	}

	for (unsigned i = 0; i < nslots; i++){
		if (v_out1[0][i] != (inputs[0][i] == inputs[1][i])){
			return true;
		}
		if (v_out2[0][i] != (inputs[0][i] > inputs[1][i])){
			return true;
		}
		if ((!v_out1[0][i] && !v_out2[0][i]) != (inputs[0][i] < inputs[1][i])){
			return true;
		}
	}
	return false;
}

bool QRS_DETECTION::test_DUALSLOPE(){
	make_copies();
	t_start();
//	he.EFF_RCMP(k[0],k[1]);
	he.RCMP(k[0],k[1]);
	t_end(__FUNCTION__);	

	vector< vector<long> > v_out1, v_out2;
	v_out1 = he.decryptNbits(k[0]);
	v_out2 = he.decryptNbits(k[1]);

	if(debug){
                for(unsigned i = 0; i < nslots; i++){
                        cout << "v_out1[0][" << i << "], v_out2[0][" << i << "]: " << v_out1[0][i] << ", " << v_out2[0][i] <<  endl;
                }
        }

	for (unsigned i = 0; i < nslots; i++){
                if (v_out1[0][i] != (inputs[0][i] == inputs[1][i])){
                        return true;
                }
                if (v_out2[0][i] != (inputs[0][i] > inputs[1][i])){
                        return true;
                }
                if ((!v_out1[0][i] && !v_out2[0][i]) != (inputs[0][i] < inputs[1][i])){
                        return true;
                }
        }

	return false;
}

bool QRS_DETECTION::test_REQ(){
	make_copies();
	t_start();
	he.REQ(k[0],k[1]);
	t_end(__FUNCTION__);
	vector<long> results(nslots);
	results = conv.matrix2LongVec(he.decryptNbits(k[0]));
	for (unsigned i = 0; i < nslots; i++){
		if (results[i] != (inputs[0][i] == inputs[1][i])){
			return true;
		}
	}
	return false;
}
bool QRS_DETECTION::test_SHIFTR(){
	if(debug){
		cout << className() << ": Running " << __FUNCTION__ << "..." << endl;
	}
	//complexity negligible
	make_copies();
	vector<long> results(nslots);
	unsigned nbits = k[0].size();
	for (unsigned shift = 0; shift <= nbits; shift++){
		make_copies();
		he.SHIFTR(k[0],shift);
		results = conv.matrix2LongVec(he.decryptNbits(k[0]));
		for (unsigned i = 0; i < nslots; i++){
			if ((inputs[0][i] >> shift) != results[i]){
				return true;
			}
		}
	}
	return false;
}
bool QRS_DETECTION::test_SHIFTL(){
	if(debug){
		cout << className() << ": Running " << __FUNCTION__ << "..." << endl;
	}
	//complexity negligible
	make_copies();
	vector<long> results(nslots);
	for (unsigned shift = 0; shift <= bits; shift++){
		make_copies();
		he.SHIFTL(k[0], shift);
		results = conv.matrix2LongVec(he.decryptNbits(k[0]));
		for (unsigned i = 0; i < nslots; i++){
			if (results[i] != (inputs[0][i] << shift)){	
				return true;
			}
		}
	}
	return false;
}
bool QRS_DETECTION::test_NMUX(){
	vector<long> v_zeros(nslots, 0), v_ones(nslots, 1);
	mkt k_zeros = he.encrypt(v_zeros);
	mkt k_ones = he.encrypt(v_ones);
	vector<long> results(nslots);
	
	make_copies();
	t_start();
	he.NMUX(k[0], k[1], k_zeros);
	t_end(__FUNCTION__);
	results = conv.matrix2LongVec(he.decryptNbits(k[0]));
	for (unsigned i = 0; i < nslots; i++){
		if (results[i] != inputs[1][i]){	
			return true;
		}
	}
	make_copies();
	he.NMUX(k[0], k[1], k_ones);
	results = conv.matrix2LongVec(he.decryptNbits(k[0]));
	for (unsigned i = 0; i < nslots; i++){
		if (results[i] != inputs[0][i]){	
			return true;
		}
	}
	return false;
}

