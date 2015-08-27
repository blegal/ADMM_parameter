
//
// NOT TESTED !
//
//void sort7_rank_order_reg(float llr[], int pos[])
//{
//    const float x0 = llr[0];
//    const float x1 = llr[1];
//    const float x2 = llr[2];
//    const float x3 = llr[3];
//    const float x4 = llr[4];
//    const float x5 = llr[5];
//    const float x6 = llr[6];
//    const int o0 = (x0< x1)+(x0 <x2)+(x0 <x3)+(x0 <x4)+(x0< x5)+(x0<x6);
//    const int o1 = (x1<=x0)+(x1 <x2)+(x1 <x3)+(x1 <x4)+(x1< x5)+(x1<x6);
//    const int o2 = (x2<=x0)+(x2<=x1)+(x2 <x3)+(x2 <x4)+(x2< x5)+(x2<x6);
//    const int o3 = (x3<=x0)+(x3<=x1)+(x3<=x2)+(x3 <x4)+(x3< x5)+(x3<x6);
//    const int o4 = (x4<=x0)+(x4<=x1)+(x4<=x2)+(x4<=x3)+(x4< x5)+(x4<x6);
//    const int o5 = (x5<=x0)+(x5<=x1)+(x5<=x2)+(x5<=x3)+(x5<=x4)+(x5<x6);
//    const int o6 = 21 - (o0 + o1 + o2 + o3 + o4 + o5);
//    llr[o0]=x0; llr[o1]=x1; llr[o2]=x2; llr[o3]=x3;
//    llr[o4]=x4; llr[o5]=x5; llr[o6]=x6;
//    pos[o0]= 0; pos[o1]= 1; pos[o2]= 2; pos[o3]= 3;
//    pos[o4]= 4; pos[o5]= 5; pos[o6]= 6;
//}
void sort8_rank_order_reg(float llr[], int pos[])
{
    const float x0 = llr[0]; const float x1 = llr[1];
    const float x2 = llr[2]; const float x3 = llr[3];
    const float x4 = llr[4]; const float x5 = llr[5];
    const float x6 = llr[6]; const float x7 = llr[7];
    const int o0 = (x0< x1)+(x0 <x2)+(x0 <x3)+(x0 <x4)+(x0 <x5)+(x0 <x6)+(x0<x7);
    const int o1 = (x1<=x0)+(x1 <x2)+(x1 <x3)+(x1 <x4)+(x1 <x5)+(x1 <x6)+(x1<x7);
    const int o2 = (x2<=x0)+(x2<=x1)+(x2 <x3)+(x2 <x4)+(x2 <x5)+(x2 <x6)+(x2<x7);
    const int o3 = (x3<=x0)+(x3<=x1)+(x3<=x2)+(x3 <x4)+(x3 <x5)+(x3 <x6)+(x3<x7);
    const int o4 = (x4<=x0)+(x4<=x1)+(x4<=x2)+(x4<=x3)+(x4 <x5)+(x4 <x6)+(x4<x7);
    const int o5 = (x5<=x0)+(x5<=x1)+(x5<=x2)+(x5<=x3)+(x5<=x4)+(x5 <x6)+(x5<x7);
    const int o6 = (x6<=x0)+(x6<=x1)+(x6<=x2)+(x6<=x3)+(x6<=x4)+(x6<=x5)+(x6<x7);
    const int o7 = 28 - (o0 + o1 + o2 + o3 + o4 + o5 + o6);
    llr[o0]=x0; llr[o1]=x1; llr[o2]=x2; llr[o3]=x3;
    llr[o4]=x4; llr[o5]=x5; llr[o6]=x6; llr[o7]=x7;
    pos[o0]= 0; pos[o1]= 1; pos[o2]= 2; pos[o3]= 3;
    pos[o4]= 4; pos[o5]= 5; pos[o6]= 6; pos[o7]= 7;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//void sort7_rank_order_reg_modif(double illr[], double rllr[], int ipos[], int rpos[])
//{
//    const double x0 = illr[0];
//    const double x1 = illr[1];
//    const double x2 = illr[2];
//    const double x3 = illr[3];
//    const double x4 = illr[4];
//    const double x5 = illr[5];
//    const double x6 = illr[6];
//    const int o0 = (x0>x1) +(x0>x2) +(x0>x3) +(x0>x4) +(x0>x5) + (x0>x6);
//    const int o1 = (x1>=x0)+(x1>x2) +(x1>x3) +(x1>x4) +(x1>x5) + (x1>x6);
//    const int o2 = (x2>=x0)+(x2>=x1)+(x2>x3) +(x2>x4) +(x2>x5) + (x2>x6);
//    const int o3 = (x3>=x0)+(x3>=x1)+(x3>=x2)+(x3>x4) +(x3>x5) + (x3>x6);
//    const int o4 = (x4>=x0)+(x4>=x1)+(x4>=x2)+(x4>=x3)+(x4>x5) + (x4>x6);
//    const int o5 = (x5>=x0)+(x5>=x1)+(x5>=x2)+(x5>=x3)+(x5>=x4)+ (x5>x6);
//    const int o6 = 21 - (o0 + o1 + o2 + o3 + o4 + o5);
//    rllr[o0]=x0; rllr[o1]=x1; rllr[o2]=x2; rllr[o3]=x3;
//    rllr[o4]=x4; rllr[o5]=x5; rllr[o6]=x6;
//    rpos[o0]=ipos[0]; rpos[o1]=ipos[1]; rpos[o2]=ipos[2]; rpos[o3]=ipos[3];
//    rpos[o4]=ipos[4]; rpos[o5]=ipos[5]; rpos[o6]=ipos[6];
//}
void sort8_rank_order_reg_modif(double llr[], double rllr[], int ipos[], int rpos[])
{
    const float x0 = llr[0]; const float x1 = llr[1];
    const float x2 = llr[2]; const float x3 = llr[3];
    const float x4 = llr[4]; const float x5 = llr[5];
    const float x6 = llr[6]; const float x7 = llr[7];
    const int o0 = (x0> x1)+(x0 >x2)+(x0 >x3)+(x0 >x4)+(x0 >x5)+(x0 >x6)+(x0>x7);
    const int o1 = (x1>=x0)+(x1 >x2)+(x1 >x3)+(x1 >x4)+(x1 >x5)+(x1 >x6)+(x1>x7);
    const int o2 = (x2>=x0)+(x2>=x1)+(x2 >x3)+(x2 >x4)+(x2 >x5)+(x2 >x6)+(x2>x7);
    const int o3 = (x3>=x0)+(x3>=x1)+(x3>=x2)+(x3 >x4)+(x3 >x5)+(x3 >x6)+(x3>x7);
    const int o4 = (x4>=x0)+(x4>=x1)+(x4>=x2)+(x4>=x3)+(x4 >x5)+(x4 >x6)+(x4>x7);
    const int o5 = (x5>=x0)+(x5>=x1)+(x5>=x2)+(x5>=x3)+(x5>=x4)+(x5 >x6)+(x5>x7);
    const int o6 = (x6>=x0)+(x6>=x1)+(x6>=x2)+(x6>=x3)+(x6>=x4)+(x6>=x5)+(x6>x7);
    const int o7 = 28 - (o0 + o1 + o2 + o3 + o4 + o5 + o6);
    rllr[o0]=x0; rllr[o1]=x1; rllr[o2]=x2; rllr[o3]=x3;
    rllr[o4]=x4; rllr[o5]=x5; rllr[o6]=x6; rllr[o7]=x7;
    rpos[o0]=ipos[0]; rpos[o1]=ipos[1]; rpos[o2]=ipos[2]; rpos[o3]=ipos[3];
    rpos[o4]=ipos[4]; rpos[o5]=ipos[5]; rpos[o6]=ipos[6]; rpos[o7]=ipos[7];
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//void sort7_rank_order_reg_modif(float illr[], float rllr[], int ipos[], int rpos[])
//{
//    const auto x0 = illr[0];
//    const auto x1 = illr[1];
//    const auto x2 = illr[2];
//    const auto x3 = illr[3];
//    const auto x4 = illr[4];
//    const auto x5 = illr[5];
//    const auto x6 = illr[6];
//    const int o0 = (x0>x1) +(x0>x2) +(x0>x3) +(x0>x4) +(x0>x5) + (x0>x6);
//    const int o1 = (x1>=x0)+(x1>x2) +(x1>x3) +(x1>x4) +(x1>x5) + (x1>x6);
//    const int o2 = (x2>=x0)+(x2>=x1)+(x2>x3) +(x2>x4) +(x2>x5) + (x2>x6);
//    const int o3 = (x3>=x0)+(x3>=x1)+(x3>=x2)+(x3>x4) +(x3>x5) + (x3>x6);
//    const int o4 = (x4>=x0)+(x4>=x1)+(x4>=x2)+(x4>=x3)+(x4>x5) + (x4>x6);
//    const int o5 = (x5>=x0)+(x5>=x1)+(x5>=x2)+(x5>=x3)+(x5>=x4)+ (x5>x6);
//    const int o6 = 21 - (o0 + o1 + o2 + o3 + o4 + o5);
//    rllr[o0]=x0; rllr[o1]=x1; rllr[o2]=x2; rllr[o3]=x3;
//    rllr[o4]=x4; rllr[o5]=x5; rllr[o6]=x6;
//    rpos[o0]=ipos[0]; rpos[o1]=ipos[1]; rpos[o2]=ipos[2]; rpos[o3]=ipos[3];
//    rpos[o4]=ipos[4]; rpos[o5]=ipos[5]; rpos[o6]=ipos[6];
//}
void sort8_rank_order_reg_modif(float llr[], float rllr[], int ipos[], int rpos[])
{
    const float x0 = llr[0]; const float x1 = llr[1];
    const float x2 = llr[2]; const float x3 = llr[3];
    const float x4 = llr[4]; const float x5 = llr[5];
    const float x6 = llr[6]; const float x7 = llr[7];
    const int o0 = (x0> x1)+(x0 >x2)+(x0 >x3)+(x0 >x4)+(x0 >x5)+(x0 >x6)+(x0>x7);
    const int o1 = (x1>=x0)+(x1 >x2)+(x1 >x3)+(x1 >x4)+(x1 >x5)+(x1 >x6)+(x1>x7);
    const int o2 = (x2>=x0)+(x2>=x1)+(x2 >x3)+(x2 >x4)+(x2 >x5)+(x2 >x6)+(x2>x7);
    const int o3 = (x3>=x0)+(x3>=x1)+(x3>=x2)+(x3 >x4)+(x3 >x5)+(x3 >x6)+(x3>x7);
    const int o4 = (x4>=x0)+(x4>=x1)+(x4>=x2)+(x4>=x3)+(x4 >x5)+(x4 >x6)+(x4>x7);
    const int o5 = (x5>=x0)+(x5>=x1)+(x5>=x2)+(x5>=x3)+(x5>=x4)+(x5 >x6)+(x5>x7);
    const int o6 = (x6>=x0)+(x6>=x1)+(x6>=x2)+(x6>=x3)+(x6>=x4)+(x6>=x5)+(x6>x7);
    const int o7 = 28 - (o0 + o1 + o2 + o3 + o4 + o5 + o6);
    rllr[o0]=x0; rllr[o1]=x1; rllr[o2]=x2; rllr[o3]=x3;
    rllr[o4]=x4; rllr[o5]=x5; rllr[o6]=x6; rllr[o7]=x7;
    rpos[o0]=ipos[0]; rpos[o1]=ipos[1]; rpos[o2]=ipos[2]; rpos[o3]=ipos[3];
    rpos[o4]=ipos[4]; rpos[o5]=ipos[5]; rpos[o6]=ipos[6]; rpos[o7]=ipos[7];
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//void Sort4Deg7(float llr[], int ipos[])
//{
//    const auto x0 = llr[0];
//    const auto x1 = llr[1];
//    const auto x2 = llr[2];
//    const auto x3 = llr[3];
//    const auto x4 = llr[4];
//    const auto x5 = llr[5];
//    const auto x6 = llr[6];
//    const int o0 = (x0<x1) +(x0<x2) +(x0<x3) +(x0<x4) +(x0<x5) + (x0<x6);
//    const int o1 = (x1<=x0)+(x1<x2) +(x1<x3) +(x1<x4) +(x1<x5) + (x1<x6);
//    const int o2 = (x2<=x0)+(x2<=x1)+(x2<x3) +(x2<x4) +(x2<x5) + (x2<x6);
//    const int o3 = (x3<=x0)+(x3<=x1)+(x3<=x2)+(x3<x4) +(x3<x5) + (x3<x6);
//    const int o4 = (x4<=x0)+(x4<=x1)+(x4<=x2)+(x4<=x3)+(x4<x5) + (x4<x6);
//    const int o5 = (x5<=x0)+(x5<=x1)+(x5<=x2)+(x5<=x3)+(x5<=x4)+ (x5<x6);
//    const int o6 = 21 - (o0 + o1 + o2 + o3 + o4 + o5);
//    llr[o0]=x0; llr[o1]=x1; llr[o2]=x2; llr[o3]=x3;
//    llr[o4]=x4; llr[o5]=x5; llr[o6]=x6;
//    ipos[0] = o0; ipos[1] = o1; ipos[2] = o2; ipos[3] = o3;
//    ipos[4] = o4; ipos[5] = o5; ipos[6] = o6; ipos[7] =  7;
//}

void Sort4Deg8(float llr[], int ipos[])
{
    const float x0 = llr[0]; const float x1 = llr[1];
    const float x2 = llr[2]; const float x3 = llr[3];
    const float x4 = llr[4]; const float x5 = llr[5];
    const float x6 = llr[6]; const float x7 = llr[7];
    const int o0 = (x0< x1)+(x0 <x2)+(x0 <x3)+(x0 <x4)+(x0 <x5)+(x0 <x6)+(x0<x7);
    const int o1 = (x1<=x0)+(x1 <x2)+(x1 <x3)+(x1 <x4)+(x1 <x5)+(x1 <x6)+(x1<x7);
    const int o2 = (x2<=x0)+(x2<=x1)+(x2 <x3)+(x2 <x4)+(x2 <x5)+(x2 <x6)+(x2<x7);
    const int o3 = (x3<=x0)+(x3<=x1)+(x3<=x2)+(x3 <x4)+(x3 <x5)+(x3 <x6)+(x3<x7);
    const int o4 = (x4<=x0)+(x4<=x1)+(x4<=x2)+(x4<=x3)+(x4 <x5)+(x4 <x6)+(x4<x7);
    const int o5 = (x5<=x0)+(x5<=x1)+(x5<=x2)+(x5<=x3)+(x5<=x4)+(x5 <x6)+(x5<x7);
    const int o6 = (x6<=x0)+(x6<=x1)+(x6<=x2)+(x6<=x3)+(x6<=x4)+(x6<=x5)+(x6<x7);
    const int o7 = 28 - (o0 + o1 + o2 + o3 + o4 + o5 + o6);
    llr[o0] = x0; llr[o1]=x1; llr[o2]=x2; llr[o3]=x3; llr[o4]=x4; llr[o5]=x5; llr[o6]=x6; llr[o7]=x7;
    ipos[0] = o0; ipos[1]=o1; ipos[2]=o2; ipos[3]=o3; ipos[4]=o4; ipos[5]=o5; ipos[6]=o6; ipos[7]=o7;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//inline void Sort4Deg7(float llr[], int pos[], int ipos[])
//{
//    const auto x0 = llr[0];
//    const auto x1 = llr[1];
//    const auto x2 = llr[2];
//    const auto x3 = llr[3];
//    const auto x4 = llr[4];
//    const auto x5 = llr[5];
//    const auto x6 = llr[6];
//    const int o0 = (x0<x1) +(x0<x2) +(x0<x3) +(x0<x4) +(x0<x5) + (x0<x6);
//    const int o1 = (x1<=x0)+(x1<x2) +(x1<x3) +(x1<x4) +(x1<x5) + (x1<x6);
//    const int o2 = (x2<=x0)+(x2<=x1)+(x2<x3) +(x2<x4) +(x2<x5) + (x2<x6);
//    const int o3 = (x3<=x0)+(x3<=x1)+(x3<=x2)+(x3<x4) +(x3<x5) + (x3<x6);
//    const int o4 = (x4<=x0)+(x4<=x1)+(x4<=x2)+(x4<=x3)+(x4<x5) + (x4<x6);
//    const int o5 = (x5<=x0)+(x5<=x1)+(x5<=x2)+(x5<=x3)+(x5<=x4)+ (x5<x6);
//    const int o6 = 21 - (o0 + o1 + o2 + o3 + o4 + o5);
//    pos[o0] =  0;  pos[o1] =  1;  pos[o2] = 2;  pos[o3] = 3;
//    pos[o4] =  4;  pos[o5] =  5;  pos[o6] = 6;  pos[ 7] = 7;
//    ipos[0] = o0; ipos[1] = o1; ipos[2] = o2; ipos[3] = o3;
//    ipos[4] = o4; ipos[5] = o5; ipos[6] = o6; ipos[7] = 7;
//}

inline void Sort4Deg8(float llr[], int pos[], int ipos[])
{
    const float x0 = llr[0]; const float x1 = llr[1];
    const float x2 = llr[2]; const float x3 = llr[3];
    const float x4 = llr[4]; const float x5 = llr[5];
    const float x6 = llr[6]; const float x7 = llr[7];
    const int o0 = (x0< x1)+(x0 <x2)+(x0 <x3)+(x0 <x4)+(x0 <x5)+(x0 <x6)+(x0<x7);
    const int o1 = (x1<=x0)+(x1 <x2)+(x1 <x3)+(x1 <x4)+(x1 <x5)+(x1 <x6)+(x1<x7);
    const int o2 = (x2<=x0)+(x2<=x1)+(x2 <x3)+(x2 <x4)+(x2 <x5)+(x2 <x6)+(x2<x7);
    const int o3 = (x3<=x0)+(x3<=x1)+(x3<=x2)+(x3 <x4)+(x3 <x5)+(x3 <x6)+(x3<x7);
    const int o4 = (x4<=x0)+(x4<=x1)+(x4<=x2)+(x4<=x3)+(x4 <x5)+(x4 <x6)+(x4<x7);
    const int o5 = (x5<=x0)+(x5<=x1)+(x5<=x2)+(x5<=x3)+(x5<=x4)+(x5 <x6)+(x5<x7);
    const int o6 = (x6<=x0)+(x6<=x1)+(x6<=x2)+(x6<=x3)+(x6<=x4)+(x6<=x5)+(x6<x7);
    const int o7 = 28 - (o0 + o1 + o2 + o3 + o4 + o5 + o6);
    pos[o0] =  0;  pos[o1] =  1;  pos[o2] = 2;  pos[o3] = 3;
    pos[o4] =  4;  pos[o5] =  5;  pos[o6] = 6;  pos[o7] = 7;
    ipos[0] = o0;  ipos[1] = o1;  ipos[2] = o2; ipos[3] = o3;
    ipos[4] = o4;  ipos[5] = o5;  ipos[6] = o6; ipos[7] = o7;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//inline void SecondTriDesDonnesDeg7(float llr[], float rllr[], int rpos[])
//{
//    const auto x0 = llr[0];
//    const auto x1 = llr[1];
//    const auto x2 = llr[2];
//    const auto x3 = llr[3];
//    const auto x4 = llr[4];
//    const auto x5 = llr[5];
//    const auto x6 = llr[6];
//    const int o0 = (x0>x1) +(x0>x2) +(x0>x3) +(x0>x4) +(x0>x5) + (x0>x6);
//    const int o1 = (x1>=x0)+(x1>x2) +(x1>x3) +(x1>x4) +(x1>x5) + (x1>x6);
//    const int o2 = (x2>=x0)+(x2>=x1)+(x2>x3) +(x2>x4) +(x2>x5) + (x2>x6);
//    const int o3 = (x3>=x0)+(x3>=x1)+(x3>=x2)+(x3>x4) +(x3>x5) + (x3>x6);
//    const int o4 = (x4>=x0)+(x4>=x1)+(x4>=x2)+(x4>=x3)+(x4>x5) + (x4>x6);
//    const int o5 = (x5>=x0)+(x5>=x1)+(x5>=x2)+(x5>=x3)+(x5>=x4)+ (x5>x6);
//    const int o6 = 21 - (o0 + o1 + o2 + o3 + o4 + o5);
//    rllr[o0]=x0; rllr[o1]=x1; rllr[o2]=x2; rllr[o3]=x3;
//    rllr[o4]=x4; rllr[o5]=x5; rllr[o6]=x6;
//    rpos[o0]=0;  rpos[o1]=1;  rpos[o2]=2;  rpos[o3]=3;
//    rpos[o4]=4;  rpos[o5]=5;  rpos[o6]=6;  rpos[ 7]=7;
//}

inline void SecondTriDesDonnesDeg8(float llr[], float rllr[], int rpos[])
{
    const float x0 = llr[0]; const float x1 = llr[1];
    const float x2 = llr[2]; const float x3 = llr[3];
    const float x4 = llr[4]; const float x5 = llr[5];
    const float x6 = llr[6]; const float x7 = llr[7];
    const int o0 = (x0> x1)+(x0 >x2)+(x0 >x3)+(x0 >x4)+(x0 >x5)+(x0 >x6)+(x0>x7);
    const int o1 = (x1>=x0)+(x1 >x2)+(x1 >x3)+(x1 >x4)+(x1 >x5)+(x1 >x6)+(x1>x7);
    const int o2 = (x2>=x0)+(x2>=x1)+(x2 >x3)+(x2 >x4)+(x2 >x5)+(x2 >x6)+(x2>x7);
    const int o3 = (x3>=x0)+(x3>=x1)+(x3>=x2)+(x3 >x4)+(x3 >x5)+(x3 >x6)+(x3>x7);
    const int o4 = (x4>=x0)+(x4>=x1)+(x4>=x2)+(x4>=x3)+(x4 >x5)+(x4 >x6)+(x4>x7);
    const int o5 = (x5>=x0)+(x5>=x1)+(x5>=x2)+(x5>=x3)+(x5>=x4)+(x5 >x6)+(x5>x7);
    const int o6 = (x6>=x0)+(x6>=x1)+(x6>=x2)+(x6>=x3)+(x6>=x4)+(x6>=x5)+(x6>x7);
    const int o7 = 28 - (o0 + o1 + o2 + o3 + o4 + o5 + o6);
	rllr[o0]=x0; rllr[o1]=x1; rllr[o2]=x2; rllr[o3]=x3;
	rllr[o4]=x4; rllr[o5]=x5; rllr[o6]=x6; rllr[o7]=x7;
	rpos[o0]=0;  rpos[o1]=1;  rpos[o2]=2;  rpos[o3]=3;
	rpos[o4]=4;  rpos[o5]=5;  rpos[o6]=6;  rpos[o7]=7;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//inline void SecondTriDesDonnesDeg7(float llr[], int pos[])
//{
//    const auto x0 = llr[0];
//    const auto x1 = llr[1];
//    const auto x2 = llr[2];
//    const auto x3 = llr[3];
//    const auto x4 = llr[4];
//    const auto x5 = llr[5];
//    const auto x6 = llr[6];
//    const int o0 = (x0>x1) +(x0>x2) +(x0>x3) +(x0>x4) +(x0>x5) + (x0>x6);
//    const int o1 = (x1>=x0)+(x1>x2) +(x1>x3) +(x1>x4) +(x1>x5) + (x1>x6);
//    const int o2 = (x2>=x0)+(x2>=x1)+(x2>x3) +(x2>x4) +(x2>x5) + (x2>x6);
//    const int o3 = (x3>=x0)+(x3>=x1)+(x3>=x2)+(x3>x4) +(x3>x5) + (x3>x6);
//    const int o4 = (x4>=x0)+(x4>=x1)+(x4>=x2)+(x4>=x3)+(x4>x5) + (x4>x6);
//    const int o5 = (x5>=x0)+(x5>=x1)+(x5>=x2)+(x5>=x3)+(x5>=x4)+ (x5>x6);
//    const int o6 = 21 - (o0 + o1 + o2 + o3 + o4 + o5);
//    pos[o0] = 0;  pos[o1]= 1;  pos[o2]= 2;  pos[o3]= 3;
//    pos[o4] = 4;  pos[o5]= 5;  pos[o6]= 6;  pos[7]=7;
//}

inline void SecondTriDesDonnesDeg8(float llr[], int pos[])
{
    const float x0 = llr[0]; const float x1 = llr[1];
    const float x2 = llr[2]; const float x3 = llr[3];
    const float x4 = llr[4]; const float x5 = llr[5];
    const float x6 = llr[6]; const float x7 = llr[7];
    const int o0 = (x0> x1)+(x0 >x2)+(x0 >x3)+(x0 >x4)+(x0 >x5)+(x0 >x6)+(x0>x7);
    const int o1 = (x1>=x0)+(x1 >x2)+(x1 >x3)+(x1 >x4)+(x1 >x5)+(x1 >x6)+(x1>x7);
    const int o2 = (x2>=x0)+(x2>=x1)+(x2 >x3)+(x2 >x4)+(x2 >x5)+(x2 >x6)+(x2>x7);
    const int o3 = (x3>=x0)+(x3>=x1)+(x3>=x2)+(x3 >x4)+(x3 >x5)+(x3 >x6)+(x3>x7);
    const int o4 = (x4>=x0)+(x4>=x1)+(x4>=x2)+(x4>=x3)+(x4 >x5)+(x4 >x6)+(x4>x7);
    const int o5 = (x5>=x0)+(x5>=x1)+(x5>=x2)+(x5>=x3)+(x5>=x4)+(x5 >x6)+(x5>x7);
    const int o6 = (x6>=x0)+(x6>=x1)+(x6>=x2)+(x6>=x3)+(x6>=x4)+(x6>=x5)+(x6>x7);
    const int o7 = 28 - (o0 + o1 + o2 + o3 + o4 + o5 + o6);
    pos[o0] = 0;  pos[o1] = 1;  pos[o2] = 2;  pos[o3] = 3;
    pos[o4] = 4;  pos[o5] = 5;  pos[o6] = 6;  pos[o7] = 7;
}
