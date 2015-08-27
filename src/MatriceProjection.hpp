/*
 * Decoder.hpp
 *
 *  Created on: 17 avr. 2015
 *      Author: legal
 *  Modified: june 2015
 */

//#define PERFORMANCE_ANALYSIS



#ifndef MatriceProjection_HPP_
#define MatriceProjection_HPP_

#include "./sorting/SortDeg6.hpp"
#include "./sorting/SortDeg7.hpp"
#include "./sorting/SortDeg8.hpp"
#include "./sorting/SortDegN.hpp"



inline bool the_compare_fx(const NODE & a, const NODE & b){return (a.value > b.value);}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////


template <class T = float, int N=32>
class MatriceProjection{
private:

	float results [32]; // DECLARATION STATIQUE POUR GAGNER DU TEMPS !
	float results_0 [32] = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f}; // DECLARATION STATIQUE POUR GAGNER DU TEMPS !
	float results_1 [32] = { 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f}; // DECLARATION STATIQUE POUR GAGNER DU TEMPS !

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

private:
    long long int counter_1;
    long long int counter_2;
    long long int counter_3;
    long long int counter_4;
    long long int temps_1;
    long long int temps_2;
    long long int temps_3;
    long long int temps_4;

#ifdef __AVX2__
    __m256 invTab  [8];
    __m256 iClipTab[8];
    __m256 minusOne[8];
#endif

public:
	MatriceProjection()
	{
        counter_1 = 0;
        counter_2 = 0;
        counter_3 = 0;
        counter_4 = 0;
        temps_1   = 0;
        temps_2   = 0;
        temps_3   = 0;
        temps_4   = 0;

				#ifdef __AVX2__
        minusOne[0] =  _mm256_set_ps (+0.0f, +0.0f, +0.0f, +0.0f, +0.0f, +0.0f, +0.0f, -1.0f);
        minusOne[1] =  _mm256_set_ps (+0.0f, +0.0f, +0.0f, +0.0f, +0.0f, +0.0f, -1.0f, -1.0f);
        minusOne[2] =  _mm256_set_ps (+0.0f, +0.0f, +0.0f, +0.0f, +0.0f, -1.0f, -1.0f, -1.0f);
        minusOne[3] =  _mm256_set_ps (+0.0f, +0.0f, +0.0f, +0.0f, -1.0f, -1.0f, -1.0f, -1.0f);
        minusOne[4] =  _mm256_set_ps (+0.0f, +0.0f, +0.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f);
        minusOne[5] =  _mm256_set_ps (+0.0f, +0.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f);
        minusOne[6] =  _mm256_set_ps (+0.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f);
        minusOne[7] =  _mm256_set_ps (-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f);
				#endif
				#ifdef __AVX2__
        iClipTab[0] =  _mm256_castsi256_ps( _mm256_set_epi32 (0x80000000, 0x80000000, 0x80000000, 0x80000000, 0x80000000, 0x80000000, 0x80000000, 0x00000000) );
        iClipTab[1] =  _mm256_castsi256_ps( _mm256_set_epi32 (0x80000000, 0x80000000, 0x80000000, 0x80000000, 0x80000000, 0x80000000, 0x00000000, 0x00000000) );
        iClipTab[2] =  _mm256_castsi256_ps( _mm256_set_epi32 (0x80000000, 0x80000000, 0x80000000, 0x80000000, 0x80000000, 0x00000000, 0x00000000, 0x00000000) );
        iClipTab[3] =  _mm256_castsi256_ps( _mm256_set_epi32 (0x80000000, 0x80000000, 0x80000000, 0x80000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000) );
        iClipTab[4] =  _mm256_castsi256_ps( _mm256_set_epi32 (0x80000000, 0x80000000, 0x80000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000) );
        iClipTab[5] =  _mm256_castsi256_ps( _mm256_set_epi32 (0x80000000, 0x80000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000) );
        iClipTab[6] =  _mm256_castsi256_ps( _mm256_set_epi32 (0x80000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000) );
        iClipTab[7] =  _mm256_castsi256_ps( _mm256_set_epi32 (0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000) );
				#endif
				#ifdef __AVX2__
        invTab[0] =  _mm256_castsi256_ps( _mm256_set_epi32 (0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x80000000) );
        invTab[1] =  _mm256_castsi256_ps( _mm256_set_epi32 (0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x80000000, 0x80000000) );
        invTab[2] =  _mm256_castsi256_ps( _mm256_set_epi32 (0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x80000000, 0x80000000, 0x80000000) );
        invTab[3] =  _mm256_castsi256_ps( _mm256_set_epi32 (0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x80000000, 0x80000000, 0x80000000, 0x80000000) );
        invTab[4] =  _mm256_castsi256_ps( _mm256_set_epi32 (0x00000000, 0x00000000, 0x00000000, 0x80000000, 0x80000000, 0x80000000, 0x80000000, 0x80000000) );
        invTab[5] =  _mm256_castsi256_ps( _mm256_set_epi32 (0x00000000, 0x00000000, 0x80000000, 0x80000000, 0x80000000, 0x80000000, 0x80000000, 0x80000000) );
        invTab[6] =  _mm256_castsi256_ps( _mm256_set_epi32 (0x00000000, 0x80000000, 0x80000000, 0x80000000, 0x80000000, 0x80000000, 0x80000000, 0x80000000) );
        invTab[7] =  _mm256_castsi256_ps( _mm256_set_epi32 (0x80000000, 0x80000000, 0x80000000, 0x80000000, 0x80000000, 0x80000000, 0x80000000, 0x80000000) );
				#endif
	}

	~MatriceProjection()
	{
#ifdef PERFORMANCE_ANALYSIS
        counter_1 ++;
        counter_2 ++;
        counter_3 ++;
        counter_4 ++;
		printf("(PP) Counter 1 value = %12ld, time %12ld (%ld)\n", counter_1, temps_1, temps_1/counter_1);
		printf("(PP) Counter 2 value = %12ld, time %12ld (%ld)\n", counter_2, temps_2, temps_2/counter_2);
		printf("(PP) Counter 3 value = %12ld, time %12ld (%ld)\n", counter_3, temps_3, temps_3/counter_3);
		printf("(PP) Counter 4 value = %12ld, time %12ld (%ld)\n", counter_4, temps_4, temps_4/counter_4);
#endif
	}



    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#define OPTIMISATION true

	T* mProjectionPPd(float v[], int length)
	{
		if( length == 6 )
		{
			#ifdef __AVX2__
      __m256 rIn = _mm256_loadu_ps(v);
      __m256 rOu = projection_deg6( rIn );
      _mm256_storeu_ps(results, rOu);
      return results;
			#else
      return mProjectionPPd<6>(v);
			#endif
    } else if( length == 7 ) {
					#ifdef __AVX2__
          __m256 rIn = _mm256_loadu_ps(v);
          __m256 rOu = projection_deg7( rIn );
          _mm256_storeu_ps(results, rOu);
          return results;
					#else
          return mProjectionPPd<7>(v);
					#endif
    } else if( length == 8 ) {
					#ifdef __AVX2__
          __m256 rIn = _mm256_loadu_ps(v);
          __m256 rOu = projection_deg8( rIn );
          _mm256_storeu_ps(results, rOu);
          return results;
					#else
          return mProjectionPPd<8>(v);
					#endif
        } else if( length == 9 ) {
        	return mProjectionPPd<9>(v);

        } else if( length == 10 ) {
        	return mProjectionPPd<10>(v);

        } else if( length == 11 ) {
        	return mProjectionPPd<11>(v);

        } else if( length == 12 ) {
        	return mProjectionPPd<12>(v);

        } else if( length == 13 ) {
        	return mProjectionPPd<13>(v);

        } else if( length == 14 ) {
        	return mProjectionPPd<14>(v);

        } else if( length == 20 ) {
        	return mProjectionPPd<20>(v);

        } else {
                cout<<"ADMM not supported"<<endl;
                exit ( 0 );
		}
		return NULL;
	}



    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



	template <int length=6> T* mProjectionPPd(float llr[])
	{
	int AllZero = 0;
	int AllOne  = 0;
	#pragma unroll
	for(int i = 0; i < length; i++)
	{
		AllZero = (llr[i] <= 0) ? AllZero + 1 : AllZero;
		AllOne  = AllOne  + (llr[i] >  1);
	}

        ///////////// VERSION ACCELEREEE DU CALCUL DES ZEROS ET UN //////////////
		if(AllZero == length) // exit if the vector has all negative values, which means that the projection is all zero vector
		{
			return results_0;
		}

		/////////////////////////////////////////////////////////////////////////////////////////////////////////

		if ( (AllOne==length) && ((length&0x01) == 0) ) // exit if the vector should be all one vector.
		{
			return results_1;
		}

		/////////////////////////////////////////////////////////////////////////////////////////////////////////
#if 1
		int zSorti[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
		if( length == 6 )
		{
			//function_sort666(llr, zSorti);
			sort6_rank_order_reg(llr, zSorti);
		}
        else if( length == 7 )
        {
            sort7_rank_order_reg(llr, zSorti);
        }
        else if( length == 8 )
        {
            sort8_rank_order_reg(llr, zSorti);
        }
        else
        {
            function_sort<length>(llr, zSorti/*, length*/);
        }

#else

        //	NODE *zSort = new NODE[length]; // vector that save the sorted input vector
		#pragma unroll
		for(int i = 0; i < length; i++) // keep the original indices while sorting,
		{
			zSort[i].index = i;
			zSort[i].value = llr[i];
		}

		if( length == 6 )
		{
			function_sort_6VNs(zSort);
		}
		else if( length == 7 )
		{
			function_sort_7VNs(zSort);
		}
		else
		{
			sort(zSort, zSort + length, the_compare_fx);//sort input vector, decreasing order
		}
#endif


		int clip_idx = 0, zero_idx= 0;
		T constituent = 0;
	//	NODE *zClip = new NODE[length];

		float llrClip[N];
		#pragma unroll
		for(int i = 0;i < length; i++)// project on the [0,1]^d cube
		{
			//zClip[i].value = min( max( zSort[i].value, 0.0f ), 1.0f);
			const auto tmp  = llr[i];
			const auto vMin = std::max(tmp,  0.0f);
			const auto vMax = std::min(vMin, 1.0f);
			llrClip[i]      = vMax;
			constituent    += vMax;//zClip[i].value;
		}

//		int r = (int)floor(constituent);
		int r = (int)constituent;
        r = r - (r & 0x01);
//		if (r & 1)
//			r--;

		// calculate constituent parity
		// calculate sum_Clip = $f_r^T z$
		T sum_Clip = 0;
		for(int i = 0; i < r+1; i++)
			sum_Clip += llrClip[i];

		for(int i = r + 1; i < length; i++)
			sum_Clip -= llrClip[i];


		if (sum_Clip <= r) // then done, return projection, beta = 0
		{
			#pragma unroll
			for(int i = 0; i < 8; i++)
				results[zSorti[i]] = llrClip[i];
			#ifdef PERFORMANCE_ANALYSIS
					counter_3 += 1;
					temps_3   += (timer() - start);
			#endif

			return results;
		}

		/////////////////////////////////////////////////////////////////////////////////////////////////////////

		T beta     = 0;
		T beta_max = 0;
		if (r + 2 <= length)
			beta_max = (llr[r] - llr[r+1])/2; // assign beta_max
		else
			beta_max = llr[r];

        // sorting zBetaRep
#if 0
		// merge beta, save sort
		int left_idx  = r, right_idx = r + 1;
		int count_idx = 0;
		while(count_idx < length)
		{
			if(left_idx < 0)
			{
				while(count_idx < length)
				{
					zBetaRep[count_idx].index = right_idx;
					zBetaRep[count_idx].value = -llr[right_idx];
					right_idx++;count_idx++;
				}
				break;
			}
			if(right_idx >= length)
			{
				while(count_idx < length)
				{
					zBetaRep[count_idx].index = left_idx;
					zBetaRep[count_idx].value = llr[left_idx] - 1;
					left_idx--;count_idx++;
				}
				break;
			}

			T temp_a =  llr[ left_idx] - 1;
			T temp_b = -llr[right_idx];
			if( (temp_a > temp_b) || (left_idx < 0) )
			{
				zBetaRep[count_idx].index = right_idx;
				zBetaRep[count_idx].value = temp_b;
				right_idx++;
				count_idx++;
			}
			else
			{
				zBetaRep[count_idx].index = left_idx;
				zBetaRep[count_idx].value = temp_a;
				left_idx --;
				count_idx++;
			}
		}
#else
		int   zSorti_m[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
		double T_in[N];
        double T_out[N];
        int   order_out[N];

		for(int i = 0; i < r+1; i++)
            T_in[i] = llr[i] - 1.0f;

		for(int i = r+1; i < length; i++)
            T_in[i] = -llr[i];

		if( length == 6 )
		{
			sort6_rank_order_reg_modif(T_in, T_out, zSorti_m, order_out);
		}
        else if( length == 7 )
        {
            sort7_rank_order_reg_modif(T_in, T_out, zSorti_m, order_out);
        }
        else if( length == 8 )
        {
            sort8_rank_order_reg_modif(T_in, T_out, zSorti_m, order_out);
        }
        else
        {
            function_sort<length>(T_in, T_out, zSorti_m, order_out);
        }
#endif


		#pragma unroll
		for(int i = 0; i < length; i++)
		{
			if (llr[i] > 1)
				clip_idx++;
			if (llr[i] >= -1e-10)
				zero_idx++;
		}

		clip_idx--;

		int idx_start = 0;
		int idx_end   = 0;

        // les beta sont trié dans zBetaRep, il faut considérer uniquement ceux entre 0 et beta_max
		#pragma unroll
		for(int i = 0; i < length; i++)
		{
			if(T_out[i] < 1e-10 )
			{
				idx_start++;
			}
			if(T_out[i] < beta_max )
			{
				idx_end++;
			}
		}
		idx_end--;
		T active_sum = 0;

		#pragma unroll
		for(int i = 0;i < length; i++)
		{
			if(i > clip_idx && i <= r)
				active_sum += llr[i];
			if(i > r && i < zero_idx)
				active_sum -= llr[i];
		}

		T total_sum = 0;
		total_sum = active_sum + clip_idx + 1;

		int previous_clip_idx, previous_zero_idx;
		T previous_active_sum;
		bool change_pre     = true;

//		previous_clip_idx   = clip_idx;
//		previous_zero_idx   = zero_idx;
//		previous_active_sum = active_sum;


		for(int i = idx_start; i <= idx_end; i++)// pour tous les beta entre 0 et beta_max
		{
			if(change_pre)
			{
				// save previous things
				previous_clip_idx   = clip_idx;
				previous_zero_idx   = zero_idx;
				previous_active_sum = active_sum;
			}
			change_pre = false;

			beta = T_out[i];
			if(order_out[i] <= r)
			{
				clip_idx--;
				active_sum += llr[order_out[i]];
			}
			else
			{
				zero_idx++;
				active_sum -= llr[order_out[i]];
			}


			if (i < length - 1)
			{
				if (beta != T_out[i+1])
				{
					total_sum = (clip_idx+1) + active_sum - beta * (zero_idx - clip_idx - 1);
					change_pre = true;
					if(total_sum < r)
						break;
				}

			}
			else if (i == length - 1)
			{
				total_sum = (clip_idx + 1)  + active_sum - beta * (zero_idx - clip_idx - 1);
				change_pre = true;
			}
		}

		if (total_sum > r)
		{
			beta = -(r - clip_idx - 1 - active_sum)/(zero_idx - clip_idx - 1);
		}
		else
		{
			beta = -(r - previous_clip_idx - 1 - previous_active_sum)/(previous_zero_idx - previous_clip_idx - 1);
		}

#if 0
			if (i <= r)
			{
				results[zSort[i].index] = min(max(zSort[i].value - beta, 0.0f),1.0f);
			}
			else
			{
				results[zSort[i].index] = min(max(zSort[i].value + beta, 0.0f),1.0f);
			}
#else
		#pragma unroll
		for(int i = 0; i < length; i++)
		{
			const auto vA = llr[i];
			const auto vB = vA - beta;
			const auto vC = vA + beta;
			const auto vD = (i <= r) ? vB : vC;
			const auto vMin = std::max(vD,   0.0f);
			const auto vMax = std::min(vMin, 1.0f);
			results[zSorti[i]] = vMax;
		}
#endif
		#ifdef PERFORMANCE_ANALYSIS
            counter_4 += 1;
            temps_4   += (timer() - start);
		#endif
		return results;
	}


    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////

		#ifdef __AVX2__
    inline __m256 projection_deg6(const __m256 llrS)
    {
        const auto length = 6;
        const auto zero     = _mm256_setzero_ps();//_mm256_set1_ps( 0.0f );
        const auto un       = _mm256_set1_ps( 1.0f );
        const auto allZeros = _mm256_cmp_ps( llrS, zero, _CMP_GT_OS );
        const auto allOnes  = _mm256_cmp_ps(   un, llrS, _CMP_GT_OS );
        int _zeros = (_mm256_movemask_ps( allZeros ) & 0x3F) == 0x00; // degree 6
        if( _zeros == true )
        {
            return zero;
        }
        int _ones  = (_mm256_movemask_ps( allOnes  ) & 0x3F) == 0x00; // degree 6
        if( _ones  )
        {
            return un;
        }
        float t_llrS[8]      __attribute__((aligned(64)));
        int   Indices4Tri[8] __attribute__((aligned(64)));
        int   Indices4Avx[8] __attribute__((aligned(64)));
        _mm256_store_ps   (t_llrS, llrS);
        Sort4Deg6(t_llrS, Indices4Tri, Indices4Avx);
        const auto imask_6  = _mm256_set_epi32(0x00000000, 0x00000000, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF);
        const auto  mask_6  = _mm256_castsi256_ps(imask_6);
        const auto sL       = _mm256_and_ps( llrS, mask_6 );
        const auto t1       = _mm256_max_ps( sL,   zero   );
        const auto llrClip  = _mm256_min_ps( t1,   un     );
        T constituent = 0.0f;
        const auto hsum = _mm256_hadd_ps(llrClip, llrClip);
        const auto psum = _mm256_add_ps(hsum, _mm256_permute2f128_ps(hsum, hsum, 0x1));
        _mm_store_ss(&constituent, _mm_hadd_ps( _mm256_castps256_ps128(psum), _mm256_castps256_ps128(psum) ) );
        const int  p = (int)constituent;
        const auto r = p - (p & 0x01);
        const auto perm  = _mm256_loadu_si256      ((const __m256i*)Indices4Tri);
        const auto cllr  = _mm256_permutevar8x32_ps(llrS,    perm); // non scale
        const auto llrs  = _mm256_permutevar8x32_ps(llrClip, perm); // scale
        _mm256_storeu_ps(t_llrS, cllr);
        float sum_Clip = 0;
        const auto fmults  = iClipTab[r];
        const auto s_clip  = _mm256_xor_ps  ( fmults, llrs  );
        const auto hclip   = _mm256_hadd_ps (s_clip, s_clip);
        const auto pclip   = _mm256_add_ps  (hclip, _mm256_permute2f128_ps(hclip, hclip, 0x1));
        _mm_store_ss(&sum_Clip, _mm_hadd_ps ( _mm256_castps256_ps128(pclip), _mm256_castps256_ps128(pclip) ) );
        if (sum_Clip <= r) // then done, return projection, beta = 0
        {
            return llrClip;
        }
        int clip_idx = 0, zero_idx= 0;
        T beta_max = 0;
        if (r + 2 <= length)
            beta_max = (t_llrS[r] - t_llrS[r+1])/2; // assign beta_max
        else
            beta_max = t_llrS[r];
        float T_in[8]      __attribute__((aligned(64)));
        float T_out[8]     __attribute__((aligned(64)));
        int   order_out[8] __attribute__((aligned(64)));
        const auto pmults  = iClipTab[r];
        const auto pdiffs  = minusOne[r];
        const auto temp_1  = _mm256_add_ps( cllr,   pdiffs  );
        const auto temp_2  = _mm256_xor_ps( pmults, temp_1  );
        _mm256_storeu_ps(T_in, temp_2);
        SecondTriDesDonnesDeg6(T_in, order_out);
        const auto permu = _mm256_loadu_si256      ((const __m256i*)order_out);
        const auto dperm = _mm256_permutevar8x32_ps(temp_2, permu);
        _mm256_storeu_ps(T_out, dperm);
        const auto small      = _mm256_set1_ps( -1e-10f );
        const auto t_clip_idx = _mm256_cmp_ps( llrS, un,    _CMP_GT_OS );
        const auto t_zero_idx = _mm256_cmp_ps( llrS, small, _CMP_GE_OS );
        const int _clip_idx   = (_mm256_movemask_ps( t_clip_idx ) & 0x3F);
        const int _zero_idx   = (_mm256_movemask_ps( t_zero_idx ) & 0x3F);
        clip_idx              = _mm_popcnt_u32( _clip_idx ) - 1;
        zero_idx              = _mm_popcnt_u32( _zero_idx );
        const auto t_out     = _mm256_loadu_ps( T_in );
        const auto vbeta     = _mm256_set1_ps ( beta_max );
        const auto t_start   = _mm256_cmp_ps  ( t_out, small, _CMP_LT_OS );
        const auto t_end     = _mm256_cmp_ps  ( t_out, vbeta, _CMP_LT_OS );
        const int _start     = (_mm256_movemask_ps( t_start ) & 0x3F);
        const int _end       = (_mm256_movemask_ps( t_end   ) & 0x3F);
        const int idx_start  = _mm_popcnt_u32( _start );
        const int idx_end    = _mm_popcnt_u32( _end   ) - 1;
        T active_sum = 0;
        #pragma unroll
        for(int i = 0;i < length; i++)
        {
            if(i > clip_idx && i <= r)
                active_sum += t_llrS[i];
            if(i > r && i < zero_idx)
                active_sum -= t_llrS[i];
        }
        T total_sum = 0;
        total_sum   = active_sum + clip_idx + 1;
        int previous_clip_idx = 0;
        int previous_zero_idx = 0;
        T previous_active_sum = 0.0f;
        bool change_pre     = true;
        T beta     = 0;
        for(int i = idx_start; i <= idx_end; i++)// pour tous les beta entre 0 et beta_max
        {
            if(change_pre)
            {
                previous_clip_idx   = clip_idx;
                previous_zero_idx   = zero_idx;
                previous_active_sum = active_sum;
            }
            change_pre = false;
            beta = T_out[i];
            if(order_out[i] <= r)
            {
                clip_idx--;
                active_sum += t_llrS[order_out[i]];
            }
            else
            {
                zero_idx++;
                active_sum -= t_llrS[order_out[i]];
            }
            if (i < length - 1)
            {
                if (beta != T_out[i+1])
                {
                    total_sum = (clip_idx + 1) + active_sum - beta * (zero_idx - clip_idx - 1);
                    change_pre = true;
                    if(total_sum < r)
                        break;
                }

            }
            else if (i == length - 1)
            {
                total_sum = (clip_idx + 1)  + active_sum - beta * (zero_idx - clip_idx - 1);
                change_pre = true;
            }
        }

        if (total_sum > r)
        {
            beta = -(r - clip_idx - 1 - active_sum)/(zero_idx - clip_idx - 1);
        }
        else
        {
            beta = -(r - previous_clip_idx - 1 - previous_active_sum) / (previous_zero_idx - previous_clip_idx - 1);
        }
        const auto facteur = invTab[r];
        const auto f_llrs  = _mm256_loadu_ps(    t_llrS        );
        const auto f_beta  = _mm256_set1_ps (   beta           );
        const auto s_beta  = _mm256_xor_ps  ( f_beta, facteur  );
        const auto r_llrs  = _mm256_add_ps  ( f_llrs, s_beta   );
        const auto s_llrs  = _mm256_max_ps  ( r_llrs,   zero   );
        const auto t_llrs  = _mm256_min_ps  ( s_llrs,   un     );
        const auto permut  = _mm256_loadu_si256((const __m256i*)Indices4Avx);
        const auto rFinal  = _mm256_permutevar8x32_ps(t_llrs, permut);
        return rFinal;
    }
		#endif
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
		#ifdef __AVX2__
    inline __m256 projection_deg7(const __m256 llrS)
    {
        const auto length = 7;
        const auto mask7  = 0x7F;
        const auto zero     = _mm256_setzero_ps();
        const auto un       = _mm256_set1_ps( 1.0f );
        const auto allZeros = _mm256_cmp_ps( llrS, zero, _CMP_GT_OS );
        const auto allOnes  = _mm256_cmp_ps(   un, llrS, _CMP_GT_OS );

        /////////////////////////////////////////////////////////////////////////////////////////////////////////

        int _zeros = (_mm256_movemask_ps( allZeros ) & mask7) == 0x00; // degree 6
        if( _zeros == true )
        {
            return zero;
        }

        /////////////////////////////////////////////////////////////////////////////////////////////////////////

        int _ones  = (_mm256_movemask_ps( allOnes  ) & mask7) == 0x00; // degree 6
        if( _ones && ( (length&0x01) == 0) )
        {
            return un;
        }

        /////////////////////////////////////////////////////////////////////////////////////////////////////////
        //
        // ON STOCKE ET ON TRIE
        //
        float t_llrS     [8] __attribute__((aligned(64)));
        int   Indices4Tri[8] __attribute__((aligned(64)));
        int   Indices4Avx[8] __attribute__((aligned(64)));
        _mm256_store_ps   (t_llrS, llrS);
        Sort4Deg7         (t_llrS, Indices4Tri, Indices4Avx);                                // <== DEGREEE RELATED CODE

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //
        // ON CLIP LES DONNEES
        //
        const auto imask_7  = _mm256_set_epi32(0x00000000, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF);
        const auto  mask_7  = _mm256_castsi256_ps(imask_7);
        const auto sL       = _mm256_and_ps( llrS, mask_7 );
        const auto t1       = _mm256_max_ps( sL,   zero   );
        const auto llrClip  = _mm256_min_ps( t1,   un     );
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //
        // ON LES ACCUMULE
        //

        T constituent = 0.0f;
        const auto hsum = _mm256_hadd_ps(llrClip, llrClip);
        const auto psum = _mm256_add_ps(hsum, _mm256_permute2f128_ps(hsum, hsum, 0x1));
        _mm_store_ss(&constituent, _mm_hadd_ps( _mm256_castps256_ps128(psum), _mm256_castps256_ps128(psum) ) );
        const int  p = (int)constituent;
        const auto r = p - (p & 0x01);

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //
        // ON LES PERMUTE
        //

        const auto perm  = _mm256_loadu_si256      ((const __m256i*)Indices4Tri);
        const auto cllr  = _mm256_permutevar8x32_ps(llrS,    perm); // non scale
        const auto llrs  = _mm256_permutevar8x32_ps(llrClip, perm); // scale
        _mm256_storeu_ps(t_llrS, cllr);


        ////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //
        // ON CALCULE LA VALEUR DE sum_Clip
        //

        float sum_Clip = 0;
        const auto fmults  = iClipTab[r];
        const auto s_clip  = _mm256_xor_ps  ( fmults, llrs  );
        const auto hclip   = _mm256_hadd_ps (s_clip, s_clip);
        const auto pclip   = _mm256_add_ps  (hclip, _mm256_permute2f128_ps(hclip, hclip, 0x1));
        _mm_store_ss(&sum_Clip, _mm_hadd_ps ( _mm256_castps256_ps128(pclip), _mm256_castps256_ps128(pclip) ) );

        if (sum_Clip <= r) // then done, return projection, beta = 0
        {
            return llrClip;
        }
        int clip_idx = 0, zero_idx= 0;

        /////////////////////////////////////////////////////////////////////////////////////////////////////////

        T beta_max = 0;
        if (r + 2 <= length)
            beta_max = (t_llrS[r] - t_llrS[r+1])/2; // assign beta_max
        else
            beta_max = t_llrS[r];

        float T_in[8]      __attribute__((aligned(64)));
        float T_out[8]     __attribute__((aligned(64)));
        int   order_out[8] __attribute__((aligned(64)));

        const auto pmults  = iClipTab[r];
        const auto pdiffs  = minusOne[r];
        const auto temp_1  = _mm256_add_ps( cllr,   pdiffs  );
        const auto temp_2  = _mm256_xor_ps( pmults, temp_1  );
        _mm256_storeu_ps(T_in, temp_2);

        SecondTriDesDonnesDeg7(T_in, order_out);                                // <== DEGREEE RELATED CODE
        const auto permu = _mm256_loadu_si256      ((const __m256i*)order_out);
        const auto dperm = _mm256_permutevar8x32_ps(temp_2, permu);
        _mm256_storeu_ps(T_out, dperm);

        const auto small      = _mm256_set1_ps( -1e-10f );
        const auto t_clip_idx = _mm256_cmp_ps( llrS, un,    _CMP_GT_OS );
        const auto t_zero_idx = _mm256_cmp_ps( llrS, small, _CMP_GE_OS );
        const int _clip_idx   = (_mm256_movemask_ps( t_clip_idx ) & mask7);
        const int _zero_idx   = (_mm256_movemask_ps( t_zero_idx ) & mask7);
        clip_idx              = _mm_popcnt_u32( _clip_idx ) - 1;
        zero_idx              = _mm_popcnt_u32( _zero_idx );

        const auto t_out     = _mm256_loadu_ps( T_in );
        const auto vbeta     = _mm256_set1_ps ( beta_max );
        const auto t_start   = _mm256_cmp_ps  ( t_out, small, _CMP_LT_OS );
        const auto t_end     = _mm256_cmp_ps  ( t_out, vbeta, _CMP_LT_OS );
        const int _start     = (_mm256_movemask_ps( t_start ) & mask7);
        const int _end       = (_mm256_movemask_ps( t_end   ) & mask7);
        const int idx_start  = _mm_popcnt_u32( _start );
        const int idx_end    = _mm_popcnt_u32( _end   ) - 1;
        //#endif

        T active_sum = 0;
        #pragma unroll
        for(int i = 0;i < length; i++)
        {
            if(i > clip_idx && i <=        r) active_sum += t_llrS[i];
            if(i > r        && i <  zero_idx) active_sum -= t_llrS[i];
        }

        T total_sum = 0;
        total_sum   = active_sum + clip_idx + 1;

        int previous_clip_idx = 0;
        int previous_zero_idx = 0;
        T previous_active_sum = 0.0f;
        bool change_pre     = true;
        T beta     = 0;
        for(int i = idx_start; i <= idx_end; i++)// pour tous les beta entre 0 et beta_max
        {
            if(change_pre)
            {
                // save previous things
                previous_clip_idx   = clip_idx;
                previous_zero_idx   = zero_idx;
                previous_active_sum = active_sum;
            }
            change_pre = false;

            beta = T_out[i];
            if(order_out[i] <= r)
            {
                clip_idx--;
                active_sum += t_llrS[order_out[i]];
            }
            else
            {
                zero_idx++;
                active_sum -= t_llrS[order_out[i]];
            }


            if (i < length - 1)
            {
                if (beta != T_out[i+1])
                {
                    total_sum = (clip_idx + 1) + active_sum - beta * (zero_idx - clip_idx - 1);
                    change_pre = true;
                    if(total_sum < r)
                        break;
                }

            }
            else if (i == length - 1)
            {
                total_sum = (clip_idx + 1)  + active_sum - beta * (zero_idx - clip_idx - 1);
                change_pre = true;
            }
        }

        if (total_sum > r)
        {
            beta = -(r - clip_idx - 1 - active_sum)/(zero_idx - clip_idx - 1);
        }
        else
        {
            beta = -(r - previous_clip_idx - 1 - previous_active_sum) / (previous_zero_idx - previous_clip_idx - 1);
        }

        const auto facteur = invTab[r];
        const auto f_llrs  = _mm256_loadu_ps(    t_llrS        );
        const auto f_beta  = _mm256_set1_ps (   beta           );
        const auto s_beta  = _mm256_xor_ps  ( f_beta, facteur  );
        const auto r_llrs  = _mm256_add_ps  ( f_llrs, s_beta   );
        const auto s_llrs  = _mm256_max_ps  ( r_llrs,   zero   );
        const auto t_llrs  = _mm256_min_ps  ( s_llrs,   un     );
        const auto permut  = _mm256_loadu_si256((const __m256i*)Indices4Avx);
        const auto rFinal  = _mm256_permutevar8x32_ps(t_llrs, permut);
        return rFinal;
    }
		#endif

    /////////////////////////////////////////////////////////////////////////////////////////////////////////

		#ifdef __AVX2__
    inline __m256 projection_deg8(const __m256 llrS)
    {
        const auto length = 8;

        ///////////// VERSION ACCELEREEE DU CALCUL DES ZEROS ET UN //////////////

        const auto zero     = _mm256_setzero_ps();
        const auto un       = _mm256_set1_ps( 1.0f );
        const auto allZeros = _mm256_cmp_ps( llrS, zero, _CMP_GT_OS );
        const auto allOnes  = _mm256_cmp_ps(   un, llrS, _CMP_GT_OS );

        /////////////////////////////////////////////////////////////////////////////////////////////////////////

        int _zeros = (_mm256_movemask_ps( allZeros )) == 0x00; // degree 6
        if( _zeros == true )
        {
            return zero;
        }

        /////////////////////////////////////////////////////////////////////////////////////////////////////////

        int _ones  = (_mm256_movemask_ps( allOnes  )) == 0x00; // degree 6
        if( _ones && ( (length&0x01) == 0) )
        {
            return un;
        }

        /////////////////////////////////////////////////////////////////////////////////////////////////////////
        //
        // ON STOCKE ET ON TRIE
        //
        float t_llrS     [8] __attribute__((aligned(64)));
        int   Indices4Tri[8] __attribute__((aligned(64)));
        int   Indices4Avx[8] __attribute__((aligned(64)));
        _mm256_store_ps   (t_llrS, llrS);
        Sort4Deg8         (t_llrS, Indices4Tri, Indices4Avx);                                // <== DEGREEE RELATED CODE

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //
        // ON CLIP LES DONNEES
        //
        //const auto sL       = _mm256_and_ps( llrS, mask_7 );
        const auto t1       = _mm256_max_ps( llrS,   zero   );
        const auto llrClip  = _mm256_min_ps( t1,   un     );
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //
        // ON LES ACCUMULE
        //

        T constituent = 0.0f;
        const auto hsum = _mm256_hadd_ps(llrClip, llrClip);
        const auto psum = _mm256_add_ps(hsum, _mm256_permute2f128_ps(hsum, hsum, 0x1));
        _mm_store_ss(&constituent, _mm_hadd_ps( _mm256_castps256_ps128(psum), _mm256_castps256_ps128(psum) ) );
        const int  p = (int)constituent;
        const auto r = p - (p & 0x01);

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //
        // ON LES PERMUTE
        //

        const auto perm  = _mm256_loadu_si256      ((const __m256i*)Indices4Tri);
        const auto cllr  = _mm256_permutevar8x32_ps(llrS,    perm); // non scale
        const auto llrs  = _mm256_permutevar8x32_ps(llrClip, perm); // scale
        _mm256_storeu_ps(t_llrS, cllr);


        ////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //
        // ON CALCULE LA VALEUR DE sum_Clip
        //

        float sum_Clip = 0;
        const auto fmults  = iClipTab[r];
        const auto s_clip  = _mm256_xor_ps  ( fmults, llrs  );
        const auto hclip   = _mm256_hadd_ps (s_clip, s_clip);
        const auto pclip   = _mm256_add_ps  (hclip, _mm256_permute2f128_ps(hclip, hclip, 0x1));
        _mm_store_ss(&sum_Clip, _mm_hadd_ps ( _mm256_castps256_ps128(pclip), _mm256_castps256_ps128(pclip) ) );

        if (sum_Clip <= r) // then done, return projection, beta = 0
        {
            return llrClip;
        }
        int clip_idx = 0, zero_idx= 0;

        /////////////////////////////////////////////////////////////////////////////////////////////////////////

        T beta_max = 0;
        if (r + 2 <= length)
            beta_max = (t_llrS[r] - t_llrS[r+1])/2; // assign beta_max
        else
            beta_max = t_llrS[r];

        float T_in[8]      __attribute__((aligned(64)));
        float T_out[8]     __attribute__((aligned(64)));
        int   order_out[8] __attribute__((aligned(64)));

        const auto pmults  = iClipTab[r];
        const auto pdiffs  = minusOne[r];
        const auto temp_1  = _mm256_add_ps( cllr,   pdiffs  );
        const auto temp_2  = _mm256_xor_ps( pmults, temp_1  );
        _mm256_storeu_ps(T_in, temp_2);

        SecondTriDesDonnesDeg8(T_in, order_out);                                // <== DEGREEE RELATED CODE
        const auto permu = _mm256_loadu_si256      ((const __m256i*)order_out);
        const auto dperm = _mm256_permutevar8x32_ps(temp_2, permu);
        _mm256_storeu_ps(T_out, dperm);

        const auto small      = _mm256_set1_ps( -1e-10f );
        const auto t_clip_idx = _mm256_cmp_ps( llrS, un,    _CMP_GT_OS );
        const auto t_zero_idx = _mm256_cmp_ps( llrS, small, _CMP_GE_OS );
        const int _clip_idx   = (_mm256_movemask_ps( t_clip_idx ) );
        const int _zero_idx   = (_mm256_movemask_ps( t_zero_idx ) );
        clip_idx              = _mm_popcnt_u32( _clip_idx ) - 1;
        zero_idx              = _mm_popcnt_u32( _zero_idx );

        const auto t_out     = _mm256_loadu_ps( T_in );
        const auto vbeta     = _mm256_set1_ps ( beta_max );
        const auto t_start   = _mm256_cmp_ps  ( t_out, small, _CMP_LT_OS );
        const auto t_end     = _mm256_cmp_ps  ( t_out, vbeta, _CMP_LT_OS );
        const int _start     = (_mm256_movemask_ps( t_start ) );
        const int _end       = (_mm256_movemask_ps( t_end   ) );
        const int idx_start  = _mm_popcnt_u32( _start );
        const int idx_end    = _mm_popcnt_u32( _end   ) - 1;
        //#endif

        T active_sum = 0;
#pragma unroll
        for(int i = 0;i < length; i++)
        {
            if(i > clip_idx && i <=        r) active_sum += t_llrS[i];
            if(i > r        && i <  zero_idx) active_sum -= t_llrS[i];
        }

        T total_sum = 0;
        total_sum   = active_sum + clip_idx + 1;

        int previous_clip_idx, previous_zero_idx;
        T previous_active_sum;
        bool change_pre     = true;
        T beta     = 0;
        for(int i = idx_start; i <= idx_end; i++)// pour tous les beta entre 0 et beta_max
        {
            if(change_pre)
            {
                // save previous things
                previous_clip_idx   = clip_idx;
                previous_zero_idx   = zero_idx;
                previous_active_sum = active_sum;
            }
            change_pre = false;

            beta = T_out[i];
            if(order_out[i] <= r)
            {
                clip_idx--;
                active_sum += t_llrS[order_out[i]];
            }
            else
            {
                zero_idx++;
                active_sum -= t_llrS[order_out[i]];
            }


            if (i < length - 1)
            {
                if (beta != T_out[i+1])
                {
                    total_sum = (clip_idx + 1) + active_sum - beta * (zero_idx - clip_idx - 1);
                    change_pre = true;
                    if(total_sum < r)
                        break;
                }

            }
            else if (i == length - 1)
            {
                total_sum = (clip_idx + 1)  + active_sum - beta * (zero_idx - clip_idx - 1);
                change_pre = true;
            }
        }

        if (total_sum > r)
        {
            beta = -(r - clip_idx - 1 - active_sum)/(zero_idx - clip_idx - 1);
        }
        else
        {
            beta = -(r - previous_clip_idx - 1 - previous_active_sum) / (previous_zero_idx - previous_clip_idx - 1);
        }

        const auto facteur = invTab[r];
        const auto f_llrs  = _mm256_loadu_ps(    t_llrS        );
        const auto f_beta  = _mm256_set1_ps (   beta           );
        const auto s_beta  = _mm256_xor_ps  ( f_beta, facteur  );
        const auto r_llrs  = _mm256_add_ps  ( f_llrs, s_beta   );
        const auto s_llrs  = _mm256_max_ps  ( r_llrs,   zero   );
        const auto t_llrs  = _mm256_min_ps  ( s_llrs,   un     );
        const auto permut  = _mm256_loadu_si256((const __m256i*)Indices4Avx);
        const auto rFinal  = _mm256_permutevar8x32_ps(t_llrs, permut);
        return rFinal;
    }
		#endif


};

#endif /* DECODER_HPP_ */
