/*
 * Decoder.hpp
 *
 *  Created on: 17 avr. 2015
 *      Author: legal
 *  Modified: june 2015
 */

inline void sort6_rank_order_reg(float llr[], int pos[])
{
		//register float x0,x1,x2,x3,x4,x5;
		const auto x0 = llr[0];
    const auto x1 = llr[1];
    const auto x2 = llr[2];
    const auto x3 = llr[3];
    const auto x4 = llr[4];
    const auto x5 = llr[5];
    const int o0 = (x0<x1) +(x0<x2)+(x0<x3)+(x0<x4)+(x0<x5);
    const int o1 = (x1<=x0)+(x1<x2)+(x1<x3)+(x1<x4)+(x1<x5);
    const int o2 = (x2<=x0)+(x2<=x1)+(x2<x3)+(x2<x4)+(x2<x5);
    const int o3 = (x3<=x0)+(x3<=x1)+(x3<=x2)+(x3<x4)+(x3<x5);
    const int o4 = (x4<=x0)+(x4<=x1)+(x4<=x2)+(x4<=x3)+(x4<x5);
    const int o5 = 15-(o0+o1+o2+o3+o4);
    llr[o0]=x0; llr[o1]=x1; llr[o2]=x2; llr[o3]=x3; llr[o4]=x4; llr[o5]=x5;
    pos[o0]= 0; pos[o1]= 1; pos[o2]= 2; pos[o3]= 3; pos[o4]= 4; pos[o5]= 5;
}

void sort6_rank_order_reg_modif(float illr[], float rllr[], int ipos[], int rpos[])
{
		const auto x0 = illr[0];
		const auto x1 = illr[1];
		const auto x2 = illr[2];
		const auto x3 = illr[3];
		const auto x4 = illr[4];
		const auto x5 = illr[5];
		const int  o0 = (x0>x1) +(x0>x2) +(x0>x3) +(x0>x4) +(x0>x5);
		const int  o1 = (x1>=x0)+(x1>x2) +(x1>x3) +(x1>x4) +(x1>x5);
		const int  o2 = (x2>=x0)+(x2>=x1)+(x2>x3) +(x2>x4) +(x2>x5);
		const int  o3 = (x3>=x0)+(x3>=x1)+(x3>=x2)+(x3>x4) +(x3>x5);
		const int  o4 = (x4>=x0)+(x4>=x1)+(x4>=x2)+(x4>=x3)+(x4>x5);
		const int  o5 = 15-(o0+o1+o2+o3+o4);
		rllr[o0]=x0;      rllr[o1]=x1;       rllr[o2]=x2;       rllr[o3]=x3;       rllr[o4]=x4;       rllr[o5]=x5;
		rpos[o0]=ipos[0]; rpos[o1]=ipos[1];  rpos[o2]=ipos[2];  rpos[o3]=ipos[3];  rpos[o4]=ipos[4];  rpos[o5]=ipos[5];
}

void sort6_rank_order_reg_modif(double illr[], double rllr[], int ipos[], int rpos[])
{
		const double x0 = illr[0];
		const double x1 = illr[1];
		const double x2 = illr[2];
		const double x3 = illr[3];
		const double x4 = illr[4];
		const double x5 = illr[5];
		const int o0 = (x0>x1) +(x0>x2) +(x0>x3) +(x0>x4) +(x0>x5);
		const int o1 = (x1>=x0)+(x1>x2) +(x1>x3) +(x1>x4) +(x1>x5);
		const int o2 = (x2>=x0)+(x2>=x1)+(x2>x3) +(x2>x4) +(x2>x5);
		const int o3 = (x3>=x0)+(x3>=x1)+(x3>=x2)+(x3>x4) +(x3>x5);
		const int o4 = (x4>=x0)+(x4>=x1)+(x4>=x2)+(x4>=x3)+(x4>x5);
		const int o5 = 15-(o0+o1+o2+o3+o4);

		rllr[o0]=x0;      rllr[o1]=x1;       rllr[o2]=x2;       rllr[o3]=x3;       rllr[o4]=x4;       rllr[o5]=x5;
		rpos[o0]=ipos[0]; rpos[o1]=ipos[1];  rpos[o2]=ipos[2];  rpos[o3]=ipos[3];  rpos[o4]=ipos[4];  rpos[o5]=ipos[5];
}

inline void Sort4Deg6(float llr[], int ipos[])
{
	//register float x0,x1,x2,x3,x4,x5;
	const auto x0 = llr[0];
  const auto x1 = llr[1];
  const auto x2 = llr[2];
  const auto x3 = llr[3];
  const auto x4 = llr[4];
  const auto x5 = llr[5];
  const int  o0 = (x0<x1) +(x0<x2)+(x0<x3)+(x0<x4)+(x0<x5);
  const int  o1 = (x1<=x0)+(x1<x2)+(x1<x3)+(x1<x4)+(x1<x5);
  const int  o2 = (x2<=x0)+(x2<=x1)+(x2<x3)+(x2<x4)+(x2<x5);
  const int  o3 = (x3<=x0)+(x3<=x1)+(x3<=x2)+(x3<x4)+(x3<x5);
  const int  o4 = (x4<=x0)+(x4<=x1)+(x4<=x2)+(x4<=x3)+(x4<x5);
  const int  o5 = 15-(o0+o1+o2+o3+o4);
  llr[o0] = x0; llr[o1]=x1; llr[o2]=x2; llr[o3]=x3; llr[o4]=x4; llr[o5]=x5;
  ipos[0] = o0; ipos[1]=o1; ipos[2]=o2; ipos[3]=o3; ipos[4]=o4; ipos[5]=o5; ipos[6]=6; ipos[7]=7;
}


    //
    //
    // ON NE FAIT QUE PRECACLULER LA VALEUR DES POSITIONS POST-TRI
    //
    //
    inline void Sort4Deg6(float llr[], int pos[], int ipos[])
    {
        //register float x0,x1,x2,x3,x4,x5;
        const auto x0 = llr[0];
        const auto x1 = llr[1];
        const auto x2 = llr[2];
        const auto x3 = llr[3];
        const auto x4 = llr[4];
        const auto x5 = llr[5];
        const int  vv1 = (x0<x1);
        const int  vv2 = (x0<x2);
        const int  vv3 = (x0<x3);
        const int  vv4 = (x0<x4);
        const int  vv5 = (x1<x2);
        const int  vv6 = (x1<x3);
        const int  vv7 = (x1<x4);
        const int  vv8 = (x2<x3);
        const int  vv9 = (x2<x4);
        const int  vvA = (x3<x4);
        const int   o0 =  (vv1) +  (vv2) +  (vv3) +  (vv4) + (x0<x5);
        const int   o1 = (!vv1) +  (vv5) +  (vv6) +  (vv7) + (x1<x5);
        const int   o2 = (!vv2) + (!vv5) +  (vv8) +  (vv9) + (x2<x5);
        const int   o3 = (!vv3) + (!vv6) + (!vv8) +  (vvA) + (x3<x5);
        const int   o4 = (!vv4) + (!vv7) + (!vv9) + (!vvA) + (x4<x5);
        const int   o5 = 15 - (o0 + o1 + o2 + o3 + o4);
         pos[o0] =  0;  pos[o1] =  1;  pos[o2] =  2;  pos[o3] =  3;  pos[o4] =  4;  pos[o5] =  5;  pos[6] = 6;  pos[7] =7;
        ipos[ 0] = o0; ipos[ 1] = o1; ipos[ 2] = o2; ipos[ 3] = o3; ipos[ 4] = o4; ipos[ 5] = o5; ipos[6] = 6; ipos[7] = 7;
    }

		#ifdef __AVX2__
    inline void Sort4Deg6(__m256 llrI, int pos[], int ipos[])
    {
        int llr[8] __attribute__((aligned(64)));
        const auto v1 = _mm256_set1_ps( 67108864.0f );
        const auto v2 = _mm256_mul_ps( v1, llrI );
        _mm256_store_si256((__m256i *)llr,  _mm256_cvttps_epi32(v2));

        //register float x0,x1,x2,x3,x4,x5;
        const auto x0 = llr[0];
        const auto x1 = llr[1];
        const auto x2 = llr[2];
        const auto x3 = llr[3];
        const auto x4 = llr[4];
        const auto x5 = llr[5];
        int o0 = (x0<x1) +(x0<x2)+(x0<x3)+(x0<x4)+(x0<x5);
        int o1 = (x1<=x0)+(x1<x2)+(x1<x3)+(x1<x4)+(x1<x5);
        int o2 = (x2<=x0)+(x2<=x1)+(x2<x3)+(x2<x4)+(x2<x5);
        int o3 = (x3<=x0)+(x3<=x1)+(x3<=x2)+(x3<x4)+(x3<x5);
        int o4 = (x4<=x0)+(x4<=x1)+(x4<=x2)+(x4<=x3)+(x4<x5);
        int o5 = 15-(o0+o1+o2+o3+o4);
        pos[o0] =  0;  pos[o1]= 1;  pos[o2]= 2;  pos[o3]= 3;  pos[o4]= 4;  pos[o5]= 5;  pos[6]=6;  pos[7]=7;
        ipos[ 0] = o0; ipos[ 1]=o1; ipos[ 2]=o2; ipos[ 3]=o3; ipos[ 4]=o4; ipos[ 5]=o5; ipos[6]=6; ipos[7]=7;
    }
		#endif

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    inline void SecondTriDesDonnesDeg6(float illr[], float rllr[], int rpos[])
    {
        const auto x0 = illr[0];
        const auto x1 = illr[1];
        const auto x2 = illr[2];
        const auto x3 = illr[3];
        const auto x4 = illr[4];
        const auto x5 = illr[5];

        const int o0 = ((x0>x1) +(x0>x2) )+((x0>x3) +(x0>x4) )+(x0>x5);
        const int o1 = ((x1>=x0)+(x1>x2) )+((x1>x3) +(x1>x4) )+(x1>x5);
        const int o2 = ((x2>=x0)+(x2>=x1))+((x2>x3) +(x2>x4) )+(x2>x5);
        const int o3 = ((x3>=x0)+(x3>=x1))+((x3>=x2)+(x3>x4) )+(x3>x5);
        const int o4 = ((x4>=x0)+(x4>=x1))+((x4>=x2)+(x4>=x3))+(x4>x5);
        const int o5 = 15-(o0+o1+o2+o3+o4);

        rllr[o0]=x0; rllr[o1]=x1; rllr[o2]=x2; rllr[o3]=x3; rllr[o4]=x4; rllr[o5]=x5;
        rpos[o0]=0;  rpos[o1]=1;  rpos[o2]=2;  rpos[o3]=3;  rpos[o4]=4;  rpos[o5]=5;   rpos[6]=6;   rpos[7]=7;
    }

    inline void SecondTriDesDonnesDeg6(float llr[], int pos[])
    {
        const auto x0 = llr[0];
        const auto x1 = llr[1];
        const auto x2 = llr[2];
        const auto x3 = llr[3];
        const auto x4 = llr[4];
        const auto x5 = llr[5];
#if 1
        const int  vv1 = (x0>x1);
        const int  vv2 = (x0>x2);
        const int  vv3 = (x0>x3);
        const int  vv4 = (x0>x4);
        const int  vv5 = (x1>x2);
        const int  vv6 = (x1>x3);
        const int  vv7 = (x1>x4);
        const int  vv8 = (x2>x3);
        const int  vv9 = (x2>x4);
        const int  vvA = (x3>x4);
        const int   o0 =  (vv1) +  (vv2) +  (vv3) +  (vv4) + (x0>x5);
        const int   o1 = (!vv1) +  (vv5) +  (vv6) +  (vv7) + (x1>x5);
        const int   o2 = (!vv2) + (!vv5) +  (vv8) +  (vv9) + (x2>x5);
        const int   o3 = (!vv3) + (!vv6) + (!vv8) +  (vvA) + (x3>x5);
        const int   o4 = (!vv4) + (!vv7) + (!vv9) + (!vvA) + (x4>x5);
#else
        const int   o0 = ((x0>x1) +(x0>x2) )+((x0>x3) +(x0>x4) )+(x0>x5);
        const int   o1 = ((x1>=x0)+(x1>x2) )+((x1>x3) +(x1>x4) )+(x1>x5);
        const int   o2 = ((x2>=x0)+(x2>=x1))+((x2>x3) +(x2>x4) )+(x2>x5);
        const int   o3 = ((x3>=x0)+(x3>=x1))+((x3>=x2)+(x3>x4) )+(x3>x5);
        const int   o4 = ((x4>=x0)+(x4>=x1))+((x4>=x2)+(x4>=x3))+(x4>x5);
#endif
        const int   o5 = 15-(o0+o1+o2+o3+o4);
         pos[o0] =  0;  pos[o1]= 1;  pos[o2]= 2;  pos[o3]= 3;  pos[o4]= 4;  pos[o5]= 5;  pos[6]=6;  pos[7]=7;
    }

		void function_sort_6VNs(NODE* v){
			const int nE = 6;
			int i, j;
			for (i = 1; i < nE; i++)
			{
				NODE tmp = v[i]; // ON RECUPERE LE NOEUD
				float value = tmp.value;
				for (j = i; (j >= 1) && (value > v[j - 1].value); j--)
				{
					v[j] = v[j - 1];
				}
				v[j] = tmp;
			}
		}
