void Decoder::ADMMDecoder_float()
{
    int maxIter          = maxIteration;
    const float mu       = para_mu;
    const float rho      = para_rho;
    const float un_m_rho = 1.0 - rho;
    const float _amu_    = alpha/mu;
    const float _2_amu_  = _amu_+ _amu_;

    for (int i = 0;i < mPCheckMapSize; i++)
    {
        Lambda  [i] = 0;
        zReplica[i] = 0.5;
    }

    for(int i = 0; i < maxIter; i++)
    {
        int ptr    = 0;
        mIteration = i + 1;
        //
        // VN processing kernel
        //
        for (int j = 0; j < mBlocklength; j++)
        {
            float temp = 0;
            for(int k = 0; k < VariableDegree[j]; k++)
            {
                temp += (zReplica[ t_row[ptr] ] + Lambda[ t_row[ptr] ]);
                ptr  +=1;
            }

            float llr  = _LogLikelihoodRatio[j];
            float t    = temp - llr / mu;
            float deg  = (float)VariableDegree[j];
            float xx   = (t  -  _amu_)/(deg - _2_amu_);
            float vMax = std::min(xx,   1.0f);
            float vMin = std::max(vMax, 0.0f);
            OutputFromDecoder[j] = vMin;
        }



        //
        // CN processing kernel
        //
        int CumSumCheckDegree = 0; // cumulative position of currect edge in factor graph
        int allVerified       = 0;
        float vector_before_proj[32] __attribute__((aligned(64)));
        for(int j = 0; j < mNChecks; j++)
        {
            int syndrom = 0;
#pragma unroll
            for(int k = 0; k < CheckDegree[j]; k++) // fill out the vector
            {
                const auto ind         = CumSumCheckDegree + k;

                //
                // CALCUL TON VN...
                //

                const auto xpred       = OutputFromDecoder[ t_col1[ind] ];
                syndrom               += (xpred - 0.5) > 0.0f;
                vector_before_proj[k]  = rho * xpred + un_m_rho * zReplica[ind] - Lambda[ind];
            }

            allVerified   += ( syndrom & 0x01 );
            _type_ *ztemp  = mp.mProjectionPPd(vector_before_proj, CheckDegree[j]);

            for(int k = 0; k < CheckDegree[j]; k++)
            {
                const auto l   = CumSumCheckDegree + k;
                Lambda[l]      = Lambda[l] + (rho * (ztemp[k] - OutputFromDecoder[ t_col1[l] ]) + un_m_rho * (ztemp[k] - zReplica[l]));
            }

#pragma unroll
            for(int k = 0; k < CheckDegree[j]; k++)
            {
                zReplica[CumSumCheckDegree + k] = ztemp[k];
            }
            CumSumCheckDegree += CheckDegree[j];

            //
            // MISE A JOUR DE VN
            //
        }

        if(allVerified == 0)
        {
            mAlgorithmConverge = true;
            mValidCodeword     = true;
            break;
        }
    }
}


void Decoder::ADMMDecoder_float_coeffs()
{
    int maxIter          = maxIteration;
    const float mu       = para_mu;
    const float rho      = para_rho;
    const float un_m_rho = 1.0 - rho;
    const float _amu_    = alpha/mu;
    const float _2_amu_  = _amu_+ _amu_;

    for (int i = 0;i < mPCheckMapSize; i++)
    {
        Lambda  [i] = 0;
        zReplica[i] = 0.5;
    }
    
    for(int i = 0; i < maxIter; i++)
    {
        int ptr    = 0;
        mIteration = i + 1;
        //
        // VN processing kernel
        //
        for (int j = 0; j < mBlocklength; j++)
        {
            float temp = 0;
            const int degVn = VariableDegree[j];
            for(int k = 0; k < degVn; k++)
            {
                temp += (zReplica[ t_row[ptr] ] + Lambda[ t_row[ptr] ]);
                ptr  +=1;
            }

            const float _amu_    = o_value[ degVn ]; // o_value = tableau des coeffs (attribut de classe)
            const float _2_amu_  = _amu_+ _amu_;
            float llr  = _LogLikelihoodRatio[j];
    		float t    = temp - llr / mu;
            float xx   = (t  -  _amu_)/(degVn - _2_amu_);
    		float vMax = std::min(xx,   1.0f);
    		float vMin = std::max(vMax, 0.0f);
    		OutputFromDecoder[j] = vMin;
/*
            float llr  = _LogLikelihoodRatio[j];
            float t    = temp - llr / mu;
            float deg  = (float)VariableDegree[j];
            float xx   = (t  -  _amu_)/(deg - _2_amu_);
            float vMax = std::min(xx,   1.0f);
            float vMin = std::max(vMax, 0.0f);
            OutputFromDecoder[j] = vMin;
*/
        }



        //
        // CN processing kernel
        //
        int CumSumCheckDegree = 0; // cumulative position of currect edge in factor graph
        int allVerified       = 0;
        float vector_before_proj[32] __attribute__((aligned(64)));
        for(int j = 0; j < mNChecks; j++)
        {
            int syndrom = 0;
#pragma unroll
            for(int k = 0; k < CheckDegree[j]; k++) // fill out the vector
            {
                const auto ind         = CumSumCheckDegree + k;
                
                //
                // CALCUL TON VN...
                //
                
                const auto xpred       = OutputFromDecoder[ t_col1[ind] ];
                syndrom               += (xpred - 0.5) > 0.0f;
                vector_before_proj[k]  = rho * xpred + un_m_rho * zReplica[ind] - Lambda[ind];
            }

            allVerified   += ( syndrom & 0x01 );
            _type_ *ztemp  = mp.mProjectionPPd(vector_before_proj, CheckDegree[j]);

            for(int k = 0; k < CheckDegree[j]; k++)
            {
                const auto l   = CumSumCheckDegree + k;
                Lambda[l]      = Lambda[l] + (rho * (ztemp[k] - OutputFromDecoder[ t_col1[l] ]) + un_m_rho * (ztemp[k] - zReplica[l]));
            }

#pragma unroll
            for(int k = 0; k < CheckDegree[j]; k++)
            {
                zReplica[CumSumCheckDegree + k] = ztemp[k];
            }
            CumSumCheckDegree += CheckDegree[j];
            
            //
            // MISE A JOUR DE VN
            //
        }

        if(allVerified == 0)
        {
            mAlgorithmConverge = true;
            mValidCodeword     = true;
            break;
        }
    }
}
