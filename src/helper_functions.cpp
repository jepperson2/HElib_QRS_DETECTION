#include "helper_functions.h"

string generate_string(int length){
	const string alphanum = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
	//"0123456789!@#$%^&*ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";
	string s;
	for(unsigned int i = 0; i < length; ++i){
		s += alphanum[rand() % (alphanum.size() - 1)];
	}
	return s;
}

vector<double> get_mV_samples(string filename, int channel, int num_header_lines, bool debug){
    vector<int> sample_numbers;
    vector<double> channel_1_values;
    vector<double> channel_2_values;

    ifstream infile(filename);

    if (!infile.is_open()){
        cout << "Error! Could not open file: " << filename << ". Returning empty vector" << endl;
        return channel_1_values;
    }

    vector<string> header(num_header_lines);

    for (int i = 0; i < num_header_lines; i++){
        getline(infile, header[i]);
    }

    string line; 
    int tmp_samp_num;
    double tmp_chan_1_val;
    double tmp_chan_2_val;
    
    while (getline(infile,line)){
        istringstream iss(line);
        if (!(iss >> tmp_samp_num >> tmp_chan_1_val >> tmp_chan_2_val)){
            cout << "Error! We expected a different file format" << endl;
            cout << "Something more like: sample_number\t channel_1 value (mV) \t channel_2 value (mV)" << endl;
            break;
        }

        sample_numbers.push_back(tmp_samp_num);
        channel_1_values.push_back(tmp_chan_1_val);
        channel_2_values.push_back(tmp_chan_2_val);
    }

    if (debug){
        for (int i = 0; i < num_header_lines; i++){
            cout << header[i] << endl;
        }

        for (int i = 0; i < sample_numbers.size(); i++){
            cout << sample_numbers[i] << "\t\t" << channel_1_values[i] << "\t" << channel_2_values[i] << endl;
        }   
    }   

    infile.close();

    if (channel == 1){ 
        return channel_1_values;
    } else {
        return channel_2_values;
    }   
}

vector<long> get_digital_samples(string filename, int channel, int num_header_lines, bool debug){
    vector<int> sample_numbers;
    vector<long> channel_1_values;
    vector<long> channel_2_values;

    ifstream infile(filename);

    if (!infile.is_open()){
        cout << "Error! Could not open file: " << filename << ". Returning empty vector" << endl;
        return channel_1_values;
    }

    vector<string> header(num_header_lines);

    for (int i = 0; i < num_header_lines; i++){
        getline(infile, header[i]);
    }

    string line; 
    int tmp_samp_num;
    long tmp_chan_1_val;
    long tmp_chan_2_val;
    
    while (getline(infile,line)){
        istringstream iss(line);
        if (!(iss >> tmp_samp_num >> tmp_chan_1_val >> tmp_chan_2_val)){
            cout << "Error! We expected a different file format" << endl;
            cout << "Something more like: sample_number\t channel_1 value (ADC unit) \t channel_2 value (ADC unit)" << endl;
            break;
        }

        sample_numbers.push_back(tmp_samp_num);
        channel_1_values.push_back(tmp_chan_1_val);
        channel_2_values.push_back(tmp_chan_2_val);
    }

    if (debug){
        for (int i = 0; i < num_header_lines; i++){
            cout << header[i] << endl;
        }

        for (int i = 0; i < sample_numbers.size(); i++){
            cout << sample_numbers[i] << "\t\t" << channel_1_values[i] << "\t" << channel_2_values[i] << endl;
        }   
    }   

    infile.close();

    if (channel == 1){ 
        return channel_1_values;
    } else {
        return channel_2_values;
    }   
}

vector<int> get_annotations(string filename, bool debug){
    vector<string> times;
    vector<int> sample_numbers;
    vector<char> types;
    vector<int> subs;
    vector<int> chans;
    vector<int> nums;
    vector<string> auxs;

    ifstream infile(filename);
    
    if (!infile.is_open()){
        cout << "Error! Could not open file: " << filename << ". Returning empty vector" << endl;
        return sample_numbers;
    }

    string file_heading;
    getline(infile, file_heading);

    string line;
    string tmp_time;
    int tmp_samp_num;
    char tmp_char;
    int tmp_sub;
    int tmp_chan;
    int tmp_num;
    string tmp_aux;

    while (getline(infile,line)){
        istringstream iss(line);
        if(!(iss >> tmp_time >> tmp_samp_num >> tmp_char >> tmp_sub >> tmp_chan >> tmp_num)){
            cout << "Error! Expecting different file format" << endl;
            cout << "Something more like: 0:00.214 \t 77 \t N \t 0 \t 0 \t 0 \t (VT" << endl;
            break;
        }
        if(iss >> tmp_aux){
            times.push_back(tmp_time);
            sample_numbers.push_back(tmp_samp_num);
            types.push_back(tmp_char);
            subs.push_back(tmp_sub);
            chans.push_back(tmp_chan);
            nums.push_back(tmp_num);
            auxs.push_back(tmp_aux);
        } else {
            times.push_back(tmp_time);
            sample_numbers.push_back(tmp_samp_num);
            types.push_back(tmp_char);
            subs.push_back(tmp_sub);
            chans.push_back(tmp_chan);
            nums.push_back(tmp_num);
            auxs.push_back("");
        }
    }   
    
    if (debug){
        cout << file_heading << endl;

        for (int i = 0; i < sample_numbers.size(); i++){
            cout << times[i] << "\t" << sample_numbers[i] << "\t" << types[i] << "\t" << subs[i] << "\t" << chans[i] << "\t" << nums[i] << "\t" << auxs[i] << endl;
        }   
    }   
    
    return sample_numbers;
}

vector<long> scale_samples(vector<double> samples, long scaling_factor){
    vector<long> scaled;

    for (int i = 0; i < samples.size(); i++){
        scaled.push_back((long)(scaling_factor*samples[i]));
    } 

    return scaled;
}

void print_banner(string title){
    if (!title.empty()){
        size_t title_length = title.length();
        size_t banner_length = title_length + 2 + 2 * 10;
        string banner_top(banner_length, '*');
        string banner_middle = string(10, '*') + " " + title + " " + string(10, '*');

        cout << endl
            << banner_top << endl
            << banner_middle << endl
            << banner_top << endl
            << endl;
    }
}

Errors::Errors(string t){
	title = t;
}
void Errors::add(string name, bool error){
	names.push_back(name);
	errors.push_back(error);
}
void Errors::display(){
	bool no_error = true;
	for (int i = 0; i < errors.size(); i++){
		if (errors[i]){
			no_error = false;
			cout << title << ": Error occured for test " << names[i] << endl;
		}
	}
	if (no_error){
		cout << title << ": ALL TESTS PASSED" << endl;
	}
}

Timing::Timing(){
	title = "";
	measure_id = 0;
}
Timing::Timing(string t){
	title = t + ": ";
	measure_id = 0;
}
void Timing::start(){
	measure_id++;
	a = clock();
}
void Timing::end(){
	b = clock();
	duration = double(b - a)/CLOCKS_PER_SEC;
	cout << title << "Duration of measurement " << measure_id << ": ";
	cout << duration << " seconds" << endl;
}
double Timing::end(string silent){
	b = clock();
	duration = double(b - a)/CLOCKS_PER_SEC;
	return duration;
}

bool Conversion::str2Bool(string bit){
	istringstream is(bit);
	bool ret;
	is >> ret;
	return ret;
}
string Conversion::bool2Str(bool value){
	ostringstream os ;
	os << value;
	return os.str();
}
string Conversion::long2Str(long value){
	ostringstream os;
	os << value;
	return os.str();
}
string Conversion::long2bitStr(long value){
	//Only for positive values (no 2's complement)
	unsigned bits = 0;
	while (value > 0) { 
		bits++;
		value = value >> 1;
	}
	bitset<32> bin(value);
	string bin_str = "";
	for (int i = bits - 1; i >= 0; i--){
		bin_str += long2Str(bin[i]);
	}
	return bin_str;
}
long Conversion::bitStr2Long(string bin_str){ //max 32 bits
	bitset<32> bin(bin_str);
	return bin.to_ulong();
}
string Conversion::bitStr2LongStr(string bin_str){ //max 32 bits
	//Only for positive values (no 2's complement)
	return long2Str(bitStr2Long(bin_str));
}
long Conversion::signedBitStr2Long(string bin_str){
	string substring = bin_str.substr(0,1);
	if (str2Bool(substring)){ //2s complement MSB
		int inv;
		string bin_inv_str = "";
		for (int i = 0; i < bin_str.length(); i++){
			substring = bin_str.substr(i,1);
			inv = !str2Bool(substring);
			bin_inv_str += bool2Str(inv);
		}
		bitset<32> binary(bin_inv_str);
		return - (binary.to_ulong() + 1);
	}
	else{
		bitset<32> binary(bin_str);
		return binary.to_ulong();
	}
}
string Conversion::signedBitStr2LongStr(string bin_str){
	return long2Str(signedBitStr2Long(bin_str));
}
vector<string> Conversion::matrix2bitStrVec(vector< vector<long> > v_out){
	unsigned nbits = v_out.size();
	unsigned nslots = v_out[0].size();
	string bin_str;
	vector<string> v_out_binstr(nslots);
	for (unsigned j = 0; j < nslots; j++){ //bit level
		bin_str = "";
		for (int i = nbits - 1; i >= 0; i--){
			//array of independent bits (ciphertexts)
			bin_str += long2Str(v_out[i][j]);
		}
		v_out_binstr[j] = bin_str;
	}
	return v_out_binstr;
}
vector<long> Conversion::matrix2LongVec(vector< vector<long> > v_out){
	unsigned nbits = v_out.size();
	unsigned nslots = v_out[0].size();
	string bin_str;
	vector<long> v_out_long(nslots);
	for (unsigned j = 0; j < nslots; j++){ //bit level
		bin_str = "";
		for (int i = nbits - 1; i >= 0; i--){
			//array of independent bits (ciphertexts)
			bin_str += long2Str(v_out[i][j]);
		}
		long x = bitStr2Long(bin_str);
		v_out_long[j] = x;
	}
	return v_out_long;
}
vector<long> Conversion::matrix2SignedLongVec(vector< vector<long> > v_out){
	unsigned nbits = v_out.size();
	unsigned nslots = v_out[0].size();
	string bin_str;
	vector<long> v_out_long(nslots);
	for (unsigned j = 0; j < nslots; j++){ //bit level
		bin_str = "";
		for (int i = nbits - 1; i >= 0; i--){
			//array of independent bits (ciphertexts)
			bin_str += long2Str(v_out[i][j]);
		}
		long x = signedBitStr2Long(bin_str);
		v_out_long[j] = x;
	}
	return v_out_long;
}
vector< vector<long> > Conversion::longVec2Matrix(vector< long > v_in){
	//Finds the number of bits necessary to represent the integers
	unsigned temp, bits = 0;
	vector< long > temp_vec = v_in;
	for(unsigned i = 0; i < temp_vec.size(); i++){
		temp = 0;
		while (temp_vec[i] > 0) { 
			temp++;
			temp_vec[i] = temp_vec[i] >> 1;
		}
		if (temp > bits){
			bits = temp;
		}
	}
	
	unsigned nslots = v_in.size();
	vector< vector<long> > v_in_bits;
	v_in_bits.resize(bits, vector<long>(nslots,0));
	
	
	for(int j = 0; j < nslots; j++){
		bitset<64> bin(v_in[j]); //max is 2^64
		for(int b = 0; b < bits; b++){
			v_in_bits[b][j] = bin[b]; //first ctxt (bits = 0) is LSB
		}
	}
	
	
	return v_in_bits;
}


