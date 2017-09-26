// Just code excerpts from what used to be TEST_CIRC_SEQ (I think), testing the different ripple circuits

/************************.h****************************

	bool test_RCADDER();
    bool test_RBSUBER();
    bool test_RCMP();
    bool test_REQ();
    bool test_SHIFTR();
    bool test_SHIFTL();
    bool test_NMUX();
	bool test_DUALSLOPE();

*/

/***********************.cpp***************************

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
bool QRS_DETECTION::test_DUALSLOPE(){
	make_copies();
    t_start();
//  he.EFF_RCMP(k[0],k[1]);
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
*/
